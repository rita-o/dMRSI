# Dynamic models here
# These define how the parameters of the dynamic models change as a function
# of the time variable (in dwMRS, that is the bvalue)
from numpy import exp
from numpy import sqrt
from numpy import pi
from numpy import asarray
from numpy import ones_like
from scipy.special import erf

# Mono-exponential model with offset
def model_exp_offset(p,t):
    # p = [amp,adc,off]
    return p[2]+p[0]*exp(-p[1]*t)

# Bi-exponential model
def model_biexp(p,t):
    # p = [amp,adc1,adc2,frac]
    return p[0]*(p[3]*exp(-p[1]*t)+(1-p[3])*exp(-p[2]*t))

# Exponential model
def model_exp(p, t):
    # p = [amp,adc]
    return p[0] * exp(-p[1] * t)

# Linear model
def model_lin(p, t):
    # p = [amp,adc]
    return p[0] - p[1] * t

# Callaghan model
def model_callaghan(p, t):
    # p = [amp, adc]
    sqrtbD = sqrt(t*p[1])
    return p[0]*sqrt(pi)/(2*sqrtbD)*erf(sqrtbD)

def model_dki(p, t):
    # p = [amp, adc, kurtosis]
    mbD = -t*p[1]
    return p[0]*exp(mbD + 1/6*mbD**2*p[2])

# ------------------------------------------------------------------------
# Gradients
# For each of the models defined above, specify the gradient
# And call these functions using the same names as above with
# '_grad' appended in the end
def model_biexp_grad(p,t):
    e1 = exp(-p[1]*t)
    e2 = exp(-p[2]*t)
    g0 = p[3]*e1+(1-p[3])*e2
    g1 = p[0]*(-p[3]*t*e1)
    g2 = p[0]*(-(1-p[3])*t*e2)
    g3 = p[0]*(e1-e2)
    return asarray([g0,g1,g2,g3])

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

def model_callaghan_grad(p, t): # validated by comparison to numerical derivative (scratch)
    bD = t*p[1]
    sqrtbD = sqrt(bD)
    erfsqrtbD = erf(sqrtbD)
    g0 = sqrt(pi)/(2*sqrtbD)*erfsqrtbD
    g1 = p[0]*(exp(-bD)  / (2 * p[1]) - sqrt(pi)*t*erfsqrtbD / (4*(t*p[1])**(3/2)))
    return asarray([g0, g1], dtype=object)

def model_dki_grad(p, t):
    mbD = -t * p[1]
    g0 = exp(mbD + 1/6*mbD**2*p[2])
    g1 = p[0]*exp(mbD + 1/6*mbD**2*p[2])*(-t + 1/3*t**2 *p[1]*p[2])
    g2 = p[0]*exp(mbD + 1/6*mbD**2*p[2])*(1/6*mbD**2)
    return asarray([g0, g1, g2], dtype=object)