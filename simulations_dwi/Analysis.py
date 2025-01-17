#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 11:30:50 2025

@author: localadmin
"""
import os
import sys
import matplotlib.pyplot as plt


plt.close('all');
os.system('clear')
os.system('cls')

############################## ADD CODE PATH ##############################
sys.path.append(os.path.join(os.path.expanduser('~'), 'Documents', 'Rita','Codes_GitHub','dMRSI'))
sys.path.append(os.path.join(os.path.expanduser('~'), 'Documents', 'Rita','Codes_GitHub','dMRSI','simulations_dwi'))

simulator_folder = os.path.join(os.path.expanduser('~'), 'Documents', 'Rita','Codes_GitHub','Simulator_WM','Simulator_GMWM')
sys.path.append(simulator_folder)

import numpy as np
import pandas as pd
import nibabel as nib
import importlib
from utils_sim import *
from custom_functions import *
importlib.reload(sys.modules['utils_sim'])
importlib.reload(sys.modules['custom_functions'])

scheme_folder = os.path.join(simulator_folder,'instructions/scheme')

############################## RUN SIMULATIONS ##############################
scheme_name = 'DTI_multi_shell.scheme'
outpath     = '/home/localadmin/Bureau/Rita/Data/Simulations_GMWM/20250116/'
create_directory(outpath) 

create_conf_MCSim(10, 115500, 0.077, 2.5e-9, 1.5e-9, scheme_name, outpath, 'voxel', simulator_folder)
run_sim(simulator_folder)
    
############################## USER INPUT ##############################
dwi_filename    = os.path.join(outpath,'_DWI_img.bfloat') 
scheme_filename = os.path.join(scheme_folder,'DTI_multi_shell.scheme')

############################## ANALYSIS ##############################

b_values = get_bvals(scheme_filename)
b_vecs   = get_bvecs(scheme_filename)
dwi      = get_dwi(dwi_filename)

(FA, MD, AD, RD, MK, AK, RK) = calculate_DKI(scheme_filename, dwi)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(b_vecs[:, 0], b_vecs[:, 1], b_vecs[:, 2])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()



