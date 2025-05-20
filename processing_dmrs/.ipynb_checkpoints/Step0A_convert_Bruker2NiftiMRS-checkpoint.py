#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert mat lab file to fsl mrs 
Not to be inside a pipeline, just a helper function to convert one dataset from EPFL to nifti files
Needs a lot of manual input

Last changed Jan 2025
@author: Rita O
"""
from fsl_mrs.utils import mrs_io
from fsl_mrs.utils import plotting as splot
from fsl_mrs.utils.misc import parse_metab_groups
from fsl_mrs.utils import fitting
from fsl_mrs.utils import preproc as proc
from fsl_mrs.core import nifti_mrs as ntools
import glob
import os
import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.io import loadmat
import copy
from custom_functions import *
from mrs_plots import *
import importlib
from bids_structure import *

importlib.reload(sys.modules['custom_functions'])
importlib.reload(sys.modules['mrs_plots'])


# Load data
dmrs_list = []
for b in list([43, 44, 47, 48, 52, 45]): # MANUAL INPUT
    mat_data   = loadmat(f'/home/malte/Documents/Projects/dMRS_starting_data_cristina/CristinasTestData/TheirFolder/Data/processed/sum/SUM_Bruker_2022-11-18_{b}_ser_processed.mat')  # MANUAL INPUT: Replace 'your_file.mat' with your actual file name
    real_field = mat_data['study'][0]['data'][0]['real'][0][0]
    real_field = np.squeeze(real_field)
    imag_field = mat_data['study'][0]['data'][0]['imag'][0][0]
    imag_field = np.squeeze(imag_field)
    nt         = mat_data['study'][0]['params'][0]['nt'][0][0][0][0]
    fid_data      = (real_field+1j*imag_field)
    fid_data        = fid_data.conj()
    cf = mat_data['study'][0]['params'][0]['sfrq'][0][0][0][0]
    bw = mat_data['study'][0]['params'][0]['sw'][0][0][0][0]
    dwell_time=1/bw
    dmrs_list.append(fid_data)

# Plot
# fig, (ax1, ax2) = plt.subplots(1, 2)
# ax1.plot(np.fft.fft(fid_data.real))
# ax2.plot(np.fft.fft(fid_data.imag))
# plt.show()
 
# Build nifti file
dmrs_list = np.stack(dmrs_list).T
dmrs_list = dmrs_list.reshape(((1, 1, 1,) + dmrs_list.shape))
bvals = list([50, 1000, 3000, 5000, 10000]) # MANUAL INPUT
bvals = [x*0.001 for x in bvals]

nifti_obj = ntools.create_nmrs.gen_nifti_mrs(
       dmrs_list,
       dwell_time,
       cf,
       dim_tags=['DIM_USER_0', None, None])

nifti_obj.set_dim_tag(
    'DIM_USER_0',
    'DIM_USER_0',
    'b-value increment',
    {'b_value': {'Value': bvals, 'Description': 'b-value in ms.μm^-2'}})


# Save nifti - # MANUAL INPUT
data_path = os.path.join(os.path.expanduser('~'), 'Documents','Projects','dMRS_starting_data_cristina','CristinasTestData')
bids_strc = create_bids_structure(subj='sub-02', sess=2, datatype='dmrs', root=data_path, 
                                            folderlevel='derivatives', workingdir='preprocessed')     
create_directory(bids_strc.get_path())
      
nifti_obj.save(bids_strc.get_path('dmrs.nii.gz'))

## Plot QA
data = mrs_io.read_FID(bids_strc.get_path('dmrs.nii.gz'))
dmrs_list = data.mrs()

QA_dir = os.path.join(bids_strc.get_path(),'QA')
plot_spectrum(mrs_list=dmrs_list, time_var=[x*1000 for x in bvals], ppmlim=(0, 5), proj='real', output_folder=QA_dir)
splot.plotly_dynMRS(dmrs_list, time_var=bvals, ppmlim=(0, 15))
   