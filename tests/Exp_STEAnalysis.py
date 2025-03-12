#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 10:14:38 2025

@author: localadmin
"""

import numpy as np
from scipy.optimize import curve_fit, least_squares
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, least_squares
from types import SimpleNamespace
from scipy.special import erf  # Import the function

def dtd_gamma_1d_data2fit(signal, xps, opt):
    """
    Converts 1D diffusion tensor distribution (DTD) data to fit parameters.
    Equivalent to the MATLAB function `dtd_gamma_1d_data2fit`.
    
    Parameters:
    - signal: np.array of shape (N,) -> The signal intensities.
    - xps: SimpleNamespace -> Contains `b` (b-values) and `n` (number of points).
    - opt: SimpleNamespace -> Contains fitting options.
    
    Returns:
    - m: np.array of fitted parameters [s0, MD, Vi, Va]
    """
    
    if opt.dtd_gamma.do_multiple_s0 and hasattr(xps, 's_ind'):
        ns = len(np.unique(xps.s_ind)) - 1
    else:
        ns = 0

    unit_to_SI = np.array([max(signal) + np.finfo(float).eps, 1e-9, (1e-9)**2, 1] + [1] * ns)

    def t2m(t):
        """Convert internal parameters to output format."""
        s0, d_iso, mu2_iso, mu2_aniso = t[:4]
        sw = t[4:] if len(t) > 4 else []
        return np.array([s0, d_iso, mu2_iso, mu2_aniso] + list(sw)) * unit_to_SI

    def my_1d_fit2data(t):
        """Convert model parameters to signal predictions."""
        m = t2m(t)
        s = dtd_gamma_1d_fit2data(m, xps)
        return s * weight  # Element-wise multiplication

    def weightfun(sthresh, mdthresh, wthresh):
        """Soft heaviside weighting function."""
        bthresh = -np.log(sthresh) / mdthresh
        return 0.5 * (1 - erf(wthresh * (xps.b - bthresh) / bthresh))

    def calc_weight_from_signal_samples():
        """Handle heteroscedasticity weighting."""
        if not hasattr(xps, 'pa_w'):
            return np.ones_like(xps.b)
        return np.sqrt(xps.pa_w / np.max(xps.pa_w))

    # Set fitting bounds
    m_lb = np.array(opt.dtd_gamma.fit_lb + [0.5] * ns)
    m_ub = np.array(opt.dtd_gamma.fit_ub + [2.0] * ns)

    m_lb[0] *= max(signal) + np.finfo(float).eps
    m_ub[0] *= max(signal) + np.finfo(float).eps

    t_lb = m_lb / unit_to_SI
    t_ub = m_ub / unit_to_SI

    r_thr = np.inf

    for _ in range(opt.dtd_gamma.fit_iters):
        weight = np.ones_like(xps.b)

        if opt.dtd_gamma.do_weight:
            weight *= weightfun(opt.dtd_gamma.weight_sthresh, opt.dtd_gamma.weight_mdthresh, opt.dtd_gamma.weight_wthresh)

        if opt.dtd_gamma.do_pa_weight:
            weight *= calc_weight_from_signal_samples()

        # Generate a random guess
        m_guess = np.random.uniform(m_lb, m_ub)
        t_guess = m_guess / unit_to_SI

        # Fit the model
        res = least_squares(my_1d_fit2data, t_guess, bounds=(t_lb, t_ub))
        t = res.x
        m = t2m(t)

        # Redo fit with updated weighting
        if opt.dtd_gamma.do_weight:
            weight *= weightfun(opt.dtd_gamma.weight_sthresh, m[1], opt.dtd_gamma.weight_wthresh)
            if opt.dtd_gamma.do_pa_weight:
                weight *= calc_weight_from_signal_samples()
            res = least_squares(my_1d_fit2data, t, bounds=(t_lb, t_ub))
            t = res.x
            m = t2m(t)

        # Compute residual
        s_fit = dtd_gamma_1d_fit2data(m, xps)
        res_value = np.sum(((signal - s_fit) * weight) ** 2)

        if res_value < r_thr:
            r_thr = res_value
            m_keep = m

    if opt.dtd_gamma.do_plot:
        plt.semilogy(xps.b, signal, '.', label='Measured')
        plt.semilogy(xps.b, s_fit, 'o', label='Fitted')
        plt.xlabel('b-value')
        plt.ylabel('Signal')
        plt.legend()
        plt.show()

    return m_keep

def dtd_gamma_1d_fit2data(m, xps):
    """
    Compute the signal based on the gamma distribution model.

    Parameters:
    - m: array-like, model parameters [s0, d_iso, mu2_iso, mu2_aniso, (optional sw)]
    - xps: SimpleNamespace, containing `b`, `b_delta`, `b_eta`, and `s_ind` (if needed)

    Returns:
    - s: Computed signal
    """
    # Convert to readable parameters
    s0, d_iso, mu2_iso, mu2_aniso = m[:4]
    rs = np.array([1] + list(m[4:]))  # Relative signal scaling factors

    # Signal baseline coefficient adjustment for multiple series
    if len(rs) > 1:
        sw = s0 * np.sum(
            np.outer(np.ones_like(xps.b), rs) * (xps.s_ind[:, None] == np.arange(1, len(rs) + 1)),
            axis=1
        )
    else:
        sw = s0

    # Compute total diffusional variance (mu2)
    mu2 = mu2_iso + (mu2_aniso * xps.b_delta**2 * (xps.b_eta**2 + 3) / 3)

    # Ensure b is a 1D array and perform broadcasting
    xps_b_reshaped = np.reshape(xps.b, (-1))  # Flatten to a 1D array

    # Now perform the operation
    s = sw * (1 + xps_b_reshaped * mu2 / d_iso) ** (-1 * (d_iso**2 / mu2))
    
    # Ensure the signal is real (important for negative variances)
    return np.real(s)

