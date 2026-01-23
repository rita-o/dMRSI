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
import json

plt.close('all');
os.system('clear')
os.system('cls')

########################## SCRIPT CONFIGURATION (EDIT AS APPPROPRIATE) ##########################

#### DATA PATH AND SUBJECTS ####
subj_list = ['sub-15']    # list of subjects to analyse

cfg                         = {}
cfg['subj_list']            = subj_list
cfg['data_path']            = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','data_CTD')          # path to where the data from the cohort is
cfg['code_path']            = os.path.join(os.path.expanduser('~'),  'Documents','Rita','Codes_GitHub','dMRSI')                     # path to code folder
cfg['code_path2']           = os.path.join(os.path.expanduser('~'),  'Documents','Rita','Codes_GitHub','dMRSI','processing_dmrs')    # path to code subfolder
cfg['toolboxes']            = os.path.join(os.path.expanduser('~'),  'Documents','Rita','Toolboxes')                                # path to where some toolboxes from matlab are (including MPPCA and tMPPCA)
cfg['prep_foldername']      = 'preprocessed'    # name of the preprocessed folder (keep 'preprocessed' as default)
cfg['analysis_foldername']  = 'analysis'        # name of the analysis folder (keep 'analysis' as default)
cfg['common_folder']        = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','common')  # path to the common folder with files needed throught the pipeline
cfg['scan_list_name']       = 'ScanList.xlsx'   # name of the excel file containing the metadata of the cohort 

#### ADD CODE PATH ####     
sys.path.append(cfg['code_path'])
sys.path.append(cfg['code_path2'])

import importlib
from bids_structure import *
from custom_functions import *

importlib.reload(sys.modules['custom_functions'])
importlib.reload(sys.modules['bids_structure'])

#### DMRS PREPROCESSING CONFIG #### 

cfg['LC_model']             = os.path.join(cfg['toolboxes'], 'LCModel')
cfg['basis_set']            = os.path.join(cfg['common_folder'], 'mrs_basis_sets','Basis_Set_dSPECIAL_differentTM')
cfg['coil_type']            = 'rat' # Options are: 'rat' for P30 rats scanned with rat cryo prob; or 'mouse' for P10, rat pups scanned with moise cryo probe

#### STEP 0. Process bruker data
from Step1_preproc import *
Step1_preproc(subj_list, cfg)

#### STEP 1. Fitting of data >>> Use fsl_mrs env
Step1_Fitting(subj_list, cfg)
