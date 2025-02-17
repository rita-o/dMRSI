#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main script to preprocess and analyse dMRI data

Last changed Jan 2025
@author: Rita O
"""

import os
import sys
import matplotlib.pyplot as plt
import json

plt.close('all');
os.system('clear')
os.system('cls')

from custom_functions import antsreg_simple
importlib.reload(sys.modules['custom_functions'])
from custom_functions import antsreg_simple

########################## SCRIPT CONFIGURATION ##########################
################### STEP 1 DATA PATH AND SUBJECTS ###################
subj_list = ['sub-01']

cfg                         = {}
cfg['subj_list']            = subj_list
cfg['data_path']            = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','dMRI_Pilot_20250207')
cfg['code_path']            = os.path.join(os.path.expanduser('~'),  'Documents','Rita','Codes_GitHub','dMRSI')
cfg['prep_foldername']      = 'preprocessed_designer'
cfg['analysis_foldername']  = 'analysis'
cfg['common_folder']        = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','common')
cfg['scan_list_name']       = 'ScanList.xlsx'
cfg['atlas']                = 'Atlas_WHS_v4'

################### ADD CODE PATH ###################
sys.path.append(cfg['code_path'] )
sys.path.append(os.path.join(cfg['code_path'], 'processing_dwi'))

import importlib
from bids_structure import *
from custom_functions import *

importlib.reload(sys.modules['custom_functions'])
importlib.reload(sys.modules['bids_structure'])

################### STEP 3 DWI PREPROCESSING CONFIG ###################
cfg['do_topup']             = 1
cfg['redo_all']             = 0
cfg['redo_bet_anat']        = 0
cfg['redo_b0_extract']      = 0
cfg['redo_merge_dwi']       = 0
cfg['redo_denoise']         = 0
cfg['redo_gibbs']           = 0
cfg['redo_topup']           = 0
cfg['redo_eddy']            = 0
cfg['redo_final_mask']      = 0
cfg['preproc_type']         = 'combined' #  'individual' or'combined'

################### STEP 4 DWI MODELING CONFIG ###################
cfg['model_list_GM']        =  ['Nexi','Sandi']
cfg['model_list_WM']        =  ['SMI']

################### STEP 5 BRAIN REGION ESTIMATES CONFIG ###################
cfg['ROIs_GM']       = ['hippocampus','M1','M2','S1','S2', 'V1', 'PL','CG', 'Thal', 'WB']
cfg['ROIs_WM']       = ['CC']

# store config file for subprocesses calling python scripts in other environments
cfg = update_cfg(cfg)

#### STEP 1. COHORT DEFINITION
from Step1_fill_study_excel import *
Step1_fill_study_excel(cfg)   

#### STEP 2. NIFTI CONVERT SUBJECT
run_script_in_conda_environment(os.path.join(cfg['code_path'], 'processing_dwi','Step2_run.py') + ' ' + cfg['data_path'],'Dicomifier')

#### STEP 3. PREPROCESS SUBJECT
from Step3_preproc import *
Step3_preproc(subj_list,cfg) 

#### STEP 4. MODELLING SUBJECT
run_script_in_conda_environment(os.path.join(cfg['code_path'], 'processing_dwi','Step4_run.py') + ' ' + cfg['data_path'],'SwissKnife')

#### STEP 5. GET VALUES 
from Step5_GetEstimates import *
Step5_GetEstimates(subj_list,cfg) 