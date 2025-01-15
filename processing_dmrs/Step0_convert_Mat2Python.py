#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert mat lab file to fsl mrs
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

dwell_time=1/7.14e3
bw = 1/dwell_time
cf = 14.1*42.576
# dwell_time = 0.0002
# bw = 1/dwell_time
# cf = 500.3


## Load data
dmrs_list = []
for b in list(range(22, 27)):
    mat_data   = loadmat(f'/home/localadmin/Documents/Rita/Data/CristinasTestData/TheirFolder/Data/processed/sum/SUM_Bruker_2022-10-31_{b}_1_ser_processed.mat')  # Replace 'your_file.mat' with your actual file name
    real_field = mat_data['study'][0]['data'][0]['real'][0][0]
    real_field = np.squeeze(real_field)
    imag_field = mat_data['study'][0]['data'][0]['imag'][0][0]
    imag_field = np.squeeze(imag_field)
    nt         = mat_data['study'][0]['params'][0]['nt'][0][0][0][0]
    fid_data      = (real_field+1j*imag_field)
    fid_data        = fid_data
    #fid_data = proc.hlsvd(fid_data, dwell_time, cf, [4, 5])
    #fid_data = proc.phaseCorrect(fid_data, bw, cf)[0]
    # fid_data[0] *= 2.0
    # fid_data = proc.applyLinPhase(
    #         fid_data,
    #         np.linspace(-bw/2, bw, len(fid_data)),
    #         -0.061E-3)
    dmrs_list.append(fid_data)

# dmrs_list, phs, freq = proc.phase_freq_align(
#         dmrs_list,
#         bw,
#         cf,
#         ppmlim=(0.2, 4.0),
#         niter=4,
#         apodize=0)



fixed_shift = proc.shiftToRef(
    dmrs_list[0],
    4.65,
    bw,
    cf,
    ppmlim=(0.2, 4.0))[1]

dmrs_list = [proc.freqshift(ff, dwell_time, -fixed_shift * cf) for ff in dmrs_list]


a =  np.fft.fft(dmrs_list[0].real)
plt.figure()
plt.plot(a)
plt.show()


dmrs_list = np.stack(dmrs_list).T
dmrs_list = dmrs_list.reshape(((1, 1, 1,) + dmrs_list.shape))
bvals = list([50, 1000, 300, 5000, 10000])
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


data_path = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','CristinasTestData')
ss_ctr=0
bids_strc = create_bids_structure(subj='sub-01', sess=1, datatype='dmrs', root=data_path, 
                                            folderlevel='derivatives', workingdir='preprocessed')
         
   
nifti_obj.save(bids_strc.get_path('dmrs.nii.gz'))

## Plot QA
data = mrs_io.read_FID(bids_strc.get_path('dmrs.nii.gz'))
dmrs_list = data.mrs()

QA_dir = os.path.join(bids_strc.get_path(),'QA')
#plot_spectrum(mrs_list=dmrs_list, time_var=bvals, ppmlim=(-5, 15), proj='real', output_folder=QA_dir)
splot.plotly_dynMRS(dmrs_list, time_var=bvals, ppmlim=(0, 15))
   