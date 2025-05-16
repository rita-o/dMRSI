#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main script to analyse dMRS data

Last changed Jan 2025
@author: Rita, Malte
"""
import os
import sys
import matplotlib.pyplot as plt

plt.close('all');
os.system('clear')

############################## ADD CODE PATH ##############################
dmrsi_path = os.path.join(os.path.expanduser('~'), 'Documents','dMRSI')
sys.path.append(dmrsi_path)
sys.path.append(os.path.join(dmrsi_path,'processing_dmrs'))


import importlib, sys
from custom_functions import *
from bids_structure import *
importlib.reload(sys.modules['custom_functions'])
importlib.reload(sys.modules['bids_structure'])
from Step0_convert_brukerraw2niimrs import *
from Step1_Fitting import *

########################## DATA PATH AND SUBJECTS ##########################
subj_list = ['sub-01']#['sub-01','sub-02','sub-03']#

cfg                         = {}
cfg['data_path']            = os.path.join('/media','localadmin','DATA','data','20250424')
cfg['water_ref_seq_number'] = 41
cfg['min_seq_number']       = 42
cfg['max_seq_number']       = 62
cfg['fixed_phase_shift']    = 245
cfg['prep_foldername']      = 'preprocessed_bruker_cc'
cfg['coil_combination_method'] = 'Bruker header info' # 'FSL MRS'
cfg['analysis_foldername']  = 'analysis'
cfg['common_folder']        = os.path.join(dmrsi_path,'common')
cfg['basis_filename']       = os.path.join(cfg['common_folder'], 'mrs_basis')
cfg['diffusion_times']      = 'all'
cfg['atlas']                = 'Atlas_WHS_v4'
cfg['diffusion_models']     =  ['exp']#['biexp','callaghan','dki','exp'] # #[]
cfg['ppm_lim']              = [.2, 4.3]
cfg['baseline']             = 'poly, 0' # 'spline, moderate' #
cfg['save_fit']             = True
cfg['model']                = 'free_shift'

#### STEP 1. Fitting of data >>> Use fsl_mrs env
Step0_convert_bruker(subj_list, cfg)
#Step0_combine_diffusion_times(subj_list, cfg)

#### STEP 1. Fitting of data >>> Use fsl_mrs env
#Step1_Fitting(subj_list, cfg)
