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

plt.close('all');
os.system('clear')
os.system('cls')

# from custom_functions import QA_gc
# importlib.reload(sys.modules['custom_functions'])
# from custom_functions import QA_gc

############################## ADD CODE PATH ##############################
sys.path.append(os.path.join(os.path.expanduser('~'),  'Documents', 'Rita','Codes_GitHub','dMRSI','processing_dwi'))
sys.path.append(os.path.join(os.path.expanduser('~'),  'Documents', 'Rita','Codes_GitHub','dMRSI'))

from Step1_fill_study_excel import *
from Step2_raw2nii2bids import *
from Step2_correct_orientation import *
from Step3_preproc import *

from bids_structure import *
from custom_functions import *

importlib.reload(sys.modules['custom_functions'])
importlib.reload(sys.modules['bids_structure'])


########################## DATA PATH AND SUBJECTS ##########################
subj_list = ['sub-01','sub-02','sub-03','sub-04']
subj_list = ['sub-01']

cfg                         = {}
cfg['data_path']            = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','dMRI_Pilot_20220116')
cfg['data_path']            = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','dMRI_Pilot_20220121')

cfg['prep_foldername']      = 'preprocessed'
cfg['analysis_foldername']  = 'analysis'
cfg['common_folder']        = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','common')
cfg['scan_list_name']       = 'ScanList.xlsx'
cfg['atlas']                = 'Atlas_WHS_v4'

    
#### STEP 1. COHORT DEFINITION >>> Use Base env
#Step1_fill_study_excel(cfg)   ## Do once or if new data is added to the excel study file


#### STEP 2. NIFTI CONVERT SUBJECT >>> Use Dicomifier env
#Step2_raw2nii2bids(subj_list,cfg) # Do once for subject
#Step2_correct_orientation(subj_list,cfg) # Do once for subject

#### STEP 3. PREPROCESS SUBJECT >>> Use Base env
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
Step3_preproc(subj_list,cfg) ### Do more than once if needed


#### STEP 4. MODELLING SUBJECT >>> Use Swissknife env
from Step4_modelling_GM import *
cfg['model_list_GM'] =  ['Nexi','Sandi']
Step4_modelling_GM(subj_list,cfg) ### Do more than once if needed

from Step4_modelling_WM import *
cfg['model_list_WM'] =  ['SMI']
Step4_modelling_WM(subj_list,cfg) ### Do more than once if needed


#### STEP 5. GET VALUES - not finished yet
from Step5_GetEstimates import *
cfg['model_list'] = ['Nexi']
cfg['ROIs']       = ['Primary visual area']
Step5_GetEstimates(subj_list,cfg) ### Do more than once if needed
