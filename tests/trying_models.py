#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 16:55:37 2025

@author: localadmin
"""

############################## ADD PATHS ##############################

import os
import sys
import importlib

sys.path.append(os.path.join(os.path.expanduser('~'), 'Documents', 'Rita','Codes_GitHub','dMRSI'))
sys.path.append(os.path.join(os.path.expanduser('~'), 'Documents', 'Rita','Codes_GitHub','dMRSI','simulations_dwi'))
from utils_sim import *
importlib.reload(sys.modules['utils_sim'])

############################## MAKE SCHEME FILE ##############################

from dmipy.core.acquisition_scheme import acquisition_scheme_from_bvalues
simulator_folder = os.path.join(os.path.expanduser('~'), 'Documents', 'Rita','Codes_GitHub','Simulator_WM','Simulator_GMWM')
scheme_name      = os.path.join(simulator_folder,'instructions/scheme/PGSE_70_dir_high_b.scheme')
b_values    = get_bvals(scheme_name)
b_values_SI = b_values * 1e6  # now given in SI units as s/m^2
b_vecs      = get_bvecs(scheme_name)
delta       = 0.0165  
Delta       = 0.05 

acq_scheme = acquisition_scheme_from_bvalues(b_values_SI, b_vecs, delta, Delta)
acq_scheme.print_acquisition_info

############################## PREPARE MODEL ##############################

from dmipy.signal_models import cylinder_models, gaussian_models
ball = gaussian_models.G1Ball()
stick = cylinder_models.C1Stick()
zeppelin = gaussian_models.G2Zeppelin()

from dmipy.distributions.distribute_models import SD1WatsonDistributed
watson_dispersed_bundle = SD1WatsonDistributed(models=[stick, zeppelin])

watson_dispersed_bundle.set_tortuous_parameter('G2Zeppelin_1_lambda_perp','C1Stick_1_lambda_par','partial_volume_0')
watson_dispersed_bundle.set_equal_parameter('G2Zeppelin_1_lambda_par', 'C1Stick_1_lambda_par')
watson_dispersed_bundle.set_fixed_parameter('G2Zeppelin_1_lambda_par', 1.7e-9)

from dmipy.core.modeling_framework import MultiCompartmentModel
NODDI_mod = MultiCompartmentModel(models=[ball, watson_dispersed_bundle])

NODDI_mod.set_fixed_parameter('G1Ball_1_lambda_iso', 3e-9)


############################## FIT MODEL ##############################

# Get data
outpath     = '/home/localadmin/Documents/Rita/Data/Simulations_GMWM/metab_in_axons/'
dwi_filename    = os.path.join(outpath,'_DWI_img.bfloat') 
dwi      = get_dwi(dwi_filename)

# Create an array
affine = np.eye(4)
img_nii = nib.Nifti1Image(dwi, affine)
img = img_nii.get_fdata()

# Fit
NODDI_fit_hcp = NODDI_mod.fit(acq_scheme, img, solver='mix')
fitted_parameters = NODDI_fit_hcp.fitted_parameters

fig, axs = plt.subplots(2, 2, figsize=[15, 15])
axs = axs.ravel()

counter = 0
for name, values in fitted_parameters.items():
    if values.squeeze().ndim != 2:
        continue
    cf = axs[counter].imshow(values.squeeze().T, origin=True, interpolation='nearest')
    axs[counter].set_axis_off()
    axs[counter].set_title(name)
    fig.colorbar(cf, ax=axs[counter], shrink=0.7)
    counter += 1


