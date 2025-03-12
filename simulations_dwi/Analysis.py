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

common_folder        = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','common','substracts_sims')

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
scheme_name = 'SMI.scheme'
outpath     = '/home/localadmin/Documents/Rita/Data/Simulations_GMWM/metab_in_astrocytes_SMI/'
substract   = os.path.join(common_folder,'astrocytes_0.03.swc')
create_directory(outpath) 

Di      = 0.3e-9
steps   = round(1/(np.pow(1e-7,2)/(6*Di*0.077)))
with open(substract.replace('.swc','_info.txt'), 'r') as file:
    content = file.read()
match = re.search(r'Voxel (\d+)', content)
vx_size = int(match.group(1))

create_conf_MCSim(int(10e4), steps, 0.077,  Di, 0.3e-9, scheme_name, outpath, substract, simulator_folder, vx_size)
run_sim(simulator_folder)
  
plot_traj(outpath, substract)
  
############################## USER INPUT ##############################
dwi_filename    = os.path.join(outpath,'_DWI_img.bfloat') 
scheme_filename = os.path.join(scheme_folder, scheme_name)

############################## ANALYSIS ##############################

b_values = get_bvals(scheme_filename)
b_vecs   = get_bvecs(scheme_filename)
dwi      = get_dwi(dwi_filename)

#(FA1, MD1, AD1, RD1, MK1, AK1, RK1) = calculate_DKI(scheme_filename, dwi)
(FA2, MD2, AD2, RD2, MK2, AK2, RK2) = calculate_DKI(scheme_filename, dwi)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(b_vecs[:, 0], b_vecs[:, 1], b_vecs[:, 2])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()



