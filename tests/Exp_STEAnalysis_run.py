#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 10:16:29 2025

@author: localadmin
"""

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.path.expanduser('~'), 'Documents', 'Rita','Codes_GitHub','dMRSI'))
sys.path.append(os.path.join(os.path.expanduser('~'), 'Documents', 'Rita','Codes_GitHub','dMRSI','analysis_dwi'))


from Exp_STEAnalysis import dtd_gamma_1d_data2fit

# Example: Simulated signal data
signal_ste = [-0.41398666, -0.83259541, -1.41036642]  # Simulated signal with noise
signal_lte = [-0.71039153, -1.2283481 , -1.73936573, -2.06857871, -2.29722797,
       -0.72125336, -1.25311263, -1.77195249, -2.08851685, -2.33610801,
       -0.71043807, -1.250951  , -1.77059149, -2.08862232, -2.33959565]

b_vals_lte = [1., 2. , 3.5, 5., 7., 1., 2. , 3.5, 5., 7., 1., 2. , 3.5, 5., 7.]

b_vals_ste = [0.5, 1. ,1.7]

from types import SimpleNamespace  # Allows xps.b instead of xps['b']

bvals = b_vals_lte
signal = signal_lte

# Define xps using SimpleNamespace (so you can do xps.b, xps.n, etc.)
xps = SimpleNamespace(
    b= bvals,  # Example b-values
    n= len(bvals),  # Number of data points
    b_delta = np.ones(len(bvals)),
    b_eta = np.ones(len(bvals)),
   # b_eta=np.array([0, 0.3, 0.6, 0.9]),  # Ensure this exists
   # s_ind=np.array([1, 1, 2, 2])  # Series indices
)

# Example: Optimization options
opt = SimpleNamespace(
    dtd_gamma=SimpleNamespace(
        do_multiple_s0=False,
        fit_lb=[0, 0, 0, 0],
        fit_ub=[2, 2e-9, (1e-9) ** 2, 1],
        fit_iters=5,
        do_weight=True,
        weight_sthresh=0.2,
        weight_mdthresh=1e-9,
        weight_wthresh=5,
        do_pa_weight=False,
        lsq_opts={'max_nfev': 1000},
        do_plot=True,  # Enable plotting
    )
)

# Call the function
m_lte = dtd_gamma_1d_data2fit(signal, xps, opt)

print("Fitted Parameters:", m)
