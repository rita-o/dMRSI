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
from datetime import datetime

plt.close('all');
os.system('clear')
os.system('cls')

########################## SCRIPT CONFIGURATION (EDIT AS APPPROPRIATE) ##########################

#### DATA PATH AND SUBJECTS ####
subj_list = ['sub-08','sub-09','sub-10','sub-11','sub-12','sub-13','sub-14','sub-15']    # list of subjects to analyse
subj_list = ['sub-08']    # list of subjects to analyse

cfg                         = {}
cfg['subj_list']            = subj_list
cfg['data_path']            = os.path.join(os.path.expanduser('~'),'Documents','Rita','Sinergia_dMRI_dMRS','data_CTD')          # path to where the data from the cohort is
cfg['code_path']            = os.path.join(os.path.expanduser('~'),'Documents','Rita','Codes_GitHub','dMRI-dMRS-Processing-Toolbox')                     # path to code folder
cfg['code_path2']           = os.path.join(os.path.expanduser('~'),'Documents','Rita','Codes_GitHub','dMRI-dMRS-Processing-Toolbox','processing_dmrs')    # path to code subfolder
cfg['toolboxes']            = os.path.join(os.path.expanduser('~'),'Documents','Rita','Toolboxes')                                # path to where some toolboxes from matlab are (including MPPCA and tMPPCA)
cfg['prep_foldername']      = 'preprocessed'    # name of the preprocessed folder (keep 'preprocessed' as default)
cfg['analysis_foldername']  = 'analysis'        # name of the analysis folder (keep 'analysis' as default)
cfg['common_folder']        = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','common')  # path to the common folder with files needed throught the pipeline
cfg['scan_list_name']       = 'ScanList_CTD.xlsx'   # name of the excel file containing the metadata of the cohort 

#### ADD CODE PATH ####     
sys.path.append(cfg['code_path'])
sys.path.append(cfg['code_path2'])

import importlib
from bids_structure import *
from custom_functions import *

importlib.reload(sys.modules['custom_functions'])
importlib.reload(sys.modules['bids_structure'])

#### DMRS PREPROCESSING CONFIG ####  

cfg['LC_model']             = os.path.join(cfg['toolboxes'], 'LCModel','binaries','linux')                              # path to LC model executable
cfg['basis_set']            = os.path.join(cfg['common_folder'], 'mrs_basis_sets','Basis_Set_dSPECIAL_differentTM')     # path to where the basis set are
cfg['coil_type']            = 'rat' # Options are: 'rat' for P30 rats scanned with rat cryo prob; or 'mouse' for P10, rat pups scanned with moise cryo probe
cfg['models']               = ["dti","stick", "dki","cylinder","cylinder_sphere","stick_sphere"]    # models used for fitting
cfg['metabolites']          = ['NAA+NAAG','Glu','Ins','GPC+PCho','Cr+PCr','Tau','Gln']              # metabolites for analysis
cfg['redo_processing']      = 1  # 1 to remove previous file and redo all processing; 0 to process only missing steps

#### SAVE CONFIG FILE ####
with open(cfg['data_path'] + '/.config_mrs.json', 'w') as f:
    json.dump(cfg, f)

#### STEP 1. Process and quantify bruker data
from Step1_preproc import *
Step1_preproc(cfg)

#### STEP 2. Fitting of data (needs SwissKnife environment)
env = os.environ.copy()
env["QT_QPA_PLATFORM"] = "offscreen"
env["XDG_RUNTIME_DIR"] = "/tmp"
timestamp = datetime.now().strftime("%Y%m%d_%H%M")
path_log = os.path.join(cfg['data_path'],'derivatives',cfg['analysis_foldername'],f"Step2_fitting_{timestamp}.log")
logfile  = open(path_log, "w")
subprocess.run( ["conda", "run", "-n", "SwissKnife", "python", 
                 os.path.join(cfg['code_path'], 'processing_dmrs','Step2_fitting.py')] 
                + [cfg['data_path']] , env=env,stdout=logfile,stderr=logfile, check=True)
logfile.close()