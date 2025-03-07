# Callaghan parametrization
import numpy as np
from statsmodels.sandbox.regression.tools import tstd_dlldy

Parameters = {
    'Phi_0': 'variable',
    'Phi_1': 'variable',
    'conc': {'other': {'dynamic': 'model_cylinder', 'params': ['c_amp', 'c_Dpara', 'c_R']},
             'Mac': {'dynamic': 'model_lin', 'params': ['c_amp', 'c_slope']}},
    'eps': 'variable',
    'gamma': 'variable',
    'baseline': {'dynamic': 'model_exp_offset', 'params': ['b_amp', 'b_adc', 'b_off']}
}

# Optionally define bounds on the parameters
Bounds = {
    'c_amp': (1e-12, None),
    'c_Dpara': (0.005, 3),
    'c_R': (0.001, 15),
    'gamma': (0, None),
    'b_amp': (None, None),
    'b_adc': (1E-5, 3),
    'b_off': (None, None)
}

# Dynamic models here
# These define how the parameters of the dynamic models change as a function
# of the time variable (in dwMRS, that is the bvalue)
from numpy import exp, sum
from numpy import sqrt
from numpy import pi
from numpy import asarray
from numpy import ones_like
from scipy.special import erf, j1, jnp_zeros, jvp
from mpmath import nsum, inf, besselj, besseljzero, exp

# Mono-exponential model with offset
def model_exp_offset(p,t):
    # p = [amp,adc,off]
    return p[2]+p[0]*exp(-p[1]*t)

# Exponential model
def model_exp(p, t):
    # p = [amp,adc]
    return p[0] * exp(-p[1] * t)

# Linear model
def model_lin(p, t):
    # p = [amp,adc]
    return p[0] - p[1] * t

# Cylinder model
## model parameters
t_Delta = 43.2 # ms diffusion time capital Delta
t_delta = 3 # ms diffusion gradient duration
## helper functions
def q(b):
    return np.sqrt(b / (t_Delta - t_delta/3))

def b_nm(n, m, bessel_arg, D, r):
    # inner summand in general
    alpha_nm2 = jnp_zeros(n, m)[-1] ** 2
    return alpha_nm2 / (alpha_nm2 - n ** 2) / (alpha_nm2 - bessel_arg ** 2) ** 2 * np.exp(
        -D * alpha_nm2 * t_Delta / r ** 2)

def E_theta(q, D, r, theta):
    # q value projections
    q_perp = q * np.sin(theta)
    q_para = q * np.cos(theta)
    # helpers
    bessel_arg = q_perp * r

    # calculate inner sum
    def a_n(n):
        # define summand dependent on n
        def this_b_m(m):
            #print(n,m,bessel_arg,D,r)
            return b_nm(n, m, bessel_arg, D, r)

        B = nsum(this_b_m, [1, inf])

        # outer summand
        return jvp(int(n), bessel_arg) ** 2 * B / (1 + (int(n) == 0))

    # caculate outer sum
    return ((2 * j1(bessel_arg)) ** 2 / bessel_arg ** 2 + 8 * bessel_arg ** 2 * nsum(a_n, [0, inf])) * np.exp(
        -D * q_para ** 2 * t_Delta)

def E(q, D, r, n_theta=1e0):
    delta_theta = np.pi/n_theta
    thetas = np.arange(1e-15,np.pi,delta_theta)
    return np.sum([np.sin(theta)*E_theta(q, D, r, theta)*delta_theta for theta in thetas])/2

def model_cylinder(p, t):
    # p = [amp, Dpara, R]
    return p[0]*np.array([E(qi, p[1], p[2]) for qi in q(t)])

# ------------------------------------------------------------------------
# Gradients
# For each of the models defined above, specify the gradient
# And call these functions using the same names as above with
# '_grad' appended in the end

def model_exp_offset_grad(p,t):
    e1 = exp(-p[1]*t)
    g0 = e1
    g1 = -t*p[0]*e1
    g2 = ones_like(t)
    return asarray([g0,g1,g2], dtype=object)

def model_exp_grad(p, t):
    e1 = exp(-p[1] * t)
    g0 = e1
    g1 = -t * p[0] * e1
    return asarray([g0, g1], dtype=object)

def model_lin_grad(p, t):
    g0 = ones_like(t)
    g1 = -t
    return asarray([g0, g1], dtype=object)

def model_cylinder_grad(p, t): # start with numerical gradient
    this_q = q(t)
    this_E = np.array([E(q, p[1], p[2]) for q in this_q])
    g0 = this_E
    delta_p1 = 1e-10
    this_E_dDplus = np.array([E(q, p[1]+delta_p1/2, p[2]) for q in this_q])
    this_E_dDminus = np.array([E(q, p[1] - delta_p1 / 2, p[2]) for q in this_q])
    g1 = p[0]*(this_E_dDplus - this_E_dDminus)/delta_p1
    delta_p2 = 1e-10
    this_E_dRplus = np.array([E(q, p[1], p[2]+delta_p2/2) for q in this_q])
    this_E_dRminus = np.array([E(q, p[1], p[2]-delta_p2/2) for q in this_q])
    g2 = p[0] * (this_E_dRplus - this_E_dRminus) / delta_p2

    return asarray([g0, g1, g2], dtype=object)