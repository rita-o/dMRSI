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

########################## SCRIPT CONFIGURATION ##########################

#### DATA PATH AND SUBJECTS ####
subj_list = ['sub-01']
cfg                         = {}
cfg['subj_list']            = subj_list
cfg['data_path']            = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','dMRI_dMRSI_Pilot_20250428')
cfg['code_path']            = os.path.join(os.path.expanduser('~'),  'Documents','Rita','Codes_GitHub','dMRSI')
cfg['code_path2']           = os.path.join(os.path.expanduser('~'),  'Documents','Rita','Codes_GitHub','dMRSI','processing_dwi')
cfg['toolboxes']            = os.path.join(os.path.expanduser('~'),  'Documents','Rita','Toolboxes')
cfg['prep_foldername']      = 'preprocessed_tMPPCA_5D_old'
cfg['analysis_foldername']  = 'analysis_tMPPCA_5D_old'
cfg['common_folder']        = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','common')
cfg['scan_list_name']       = 'ScanList.xlsx'
cfg['atlas']                = 'Atlas_WHS_v4'
cfg['atlas_TPM']            = 'TPM_C57Bl6'

#### ADD CODE PATH ####     
sys.path.append(cfg['code_path'] )
sys.path.append(os.path.join(cfg['code_path'], 'processing_dwi'))

import importlib
from bids_structure import *
from custom_functions import *

importlib.reload(sys.modules['custom_functions'])
importlib.reload(sys.modules['bids_structure'])

#### DWI PREPROCESSING CONFIG #### 
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
cfg['algo_denoising']       = 'tMPPCA_5D'     # 'MPPCA', or ''tMPPCA''
cfg['anat_thr']             = '4000'     # 2100, 4000 depending on your data

#### DWI MODEL CONFIG ####
cfg['model_list_GM']        =  ['Nexi','Sandi']
cfg['model_list_WM']        =  ['SMI']
cfg['LTEDelta_for_microFA'] =  38 
cfg['redo_modelling']       =  0

#### ROIS CONFIG ####
cfg['ROIs_GM']       = ['hippocampus','M1','M2','S1','S2', 'V1', 'PL','CG', 'Thal', 'WB']
cfg['ROIs_WM']       = ['CC']

#### SAVE CONFIG FILE ####
cfg = update_cfg(cfg)


########################## DATA PROCESSING ##########################

#### STEP 1. COHORT DEFINITION ####
from Step1_fill_study_excel import *
Step1_fill_study_excel(cfg)   

#### STEP 2. NIFTI CONVERT SUBJECT  ####

# 2.1 Use Dicomifier to convert to nifi (pass as argument the datapath to load the cfg file)
subprocess.run( ["conda", "run", "-n", "Dicomifier", "python", 
                 os.path.join(cfg['code_path'], 'processing_dwi','Step2_raw2nii2bids.py')] 
                + [cfg['data_path']] , check=True)

# 2.2 Correct orientation from bruker system to be consisten with normal atlas and everything else
from Step2_correct_orientation import *
Step2_correct_orientation(subj_list, cfg)  

#### STEP 3. PREPROCESS SUBJECT ####

# 3.1 Process normal linear encoding data with PSGE (LTE), along with anatomical image
from Step3_preproc import *
Step3_preproc(subj_list,cfg) 

# 3.2 Process spherical encoding (STE) data, assumes an anatomical image has already been processed before
from Step3_preproc_STE import *
Step3_preproc_STE(subj_list,cfg) 

# 3.3 Register anatomical T2w to dwi for individual or combined diffusion times. 
# If exists also fits STE to LTE. Does not register across sessions for now
from Step3_registrations import *
Step3_registrations(subj_list, cfg)

#### STEP 4. MODELLING SUBJECT ####

# 4.1 Fit the dwi signal with models like Nexi, Sandi, SMI, ....
from Step4_modelling import *
Step4_modelling(subj_list,cfg)

#### STEP 5. GET VALUES ####

# 5.1 Retreives parameter estimates from the model fits, making summary figures and excel with data in certain ROIs
from Step5_get_estimates import *
Step5_get_estimates(subj_list,cfg) 
