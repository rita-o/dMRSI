#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is the main script used to preprocess and analyze diffusion MRI (dMRI) data.  
This pipeline is designed to process multi-shell diffusion data with multiple diffusion times, 
supporting both Linear Tensor Encoding (LTE) and Spherical Tensor Encoding (STE) 
for processing and analysis, along with an anatomical reference image (T1- or T2-weighted).

Current parameters are set up for HUMAN data by default.

Please open each processing step script (StepX.py) to understand better what 
is being done at each step.
I don't advise just clicking run on this script, but rather running each step 
individually and checking each step outputs.

Last changed June 2025
@author: Rita O
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
subj_list = ['sub-01']   # list of subjects to analyse

cfg                         = {}
cfg['subj_list']            = subj_list
cfg['data_path']            = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','Human_NEXI_aim1_pilot_short')        # path to where the data from the cohort is
cfg['code_path']            = os.path.join(os.path.expanduser('~'),  'Documents','Rita','Codes_GitHub','dMRSI')                     # path to code folder
cfg['code_path2']           = os.path.join(os.path.expanduser('~'),  'Documents','Rita','Codes_GitHub','dMRSI','processing_dwi')    # path to code subfolder
cfg['toolboxes']            = os.path.join(os.path.expanduser('~'),  'Documents','Rita','Toolboxes')                                # path to where some toolboxes from matlab are (including MPPCA and tMPPCA)
cfg['prep_foldername']      = 'preprocessed'    # name of the preprocessed folder (keep 'preprocessed' as default)
cfg['analysis_foldername']  = 'analysis'        # name of the analysis folder (keep 'analysis' as default)
cfg['common_folder']        = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','common')  # path to the common folder with files needed throught the pipeline
cfg['scan_list_name']       = 'ScanList_secondRec.xlsx'   # name of the excel file containing the metadata of the cohort 
cfg['atlas']                = 'Atlas_DKT'                 # name of the brain atlas to be used in the analysis. Possible options are 'Atlas_Juelich', 'Atlas_DKT', Atlas_Neuromorphometrics. This atlas needs to exists in the common folder
cfg['atlas_TPM']            = 'TPM_human_spm'             # name of the tissue probability map (tpm) to be used to threshold GM and WM to define more precisly the ROIs. This atlas needs to exists in the common folder

#### ADD CODE PATH ####     
sys.path.append(cfg['code_path'] )
sys.path.append(os.path.join(cfg['code_path2']))

import importlib
from bids_structure import *
from custom_functions import *

importlib.reload(sys.modules['custom_functions'])
importlib.reload(sys.modules['bids_structure'])

#### DWI PREPROCESSING CONFIG #### 
cfg['do_topup']             = 1   # If there is data to do topup correction, set to 1, otherwise to 0
cfg['redo_all']             = 0   # If you want to redo everything set to 1 (will delete previous results)
cfg['redo_bet_anat']        = 0   # If you want to redo the processing from the brain extraction set to 1
cfg['redo_b0_extract']      = 0   # If you want to redo the processing from the dwi b0 extraction set to 1
cfg['redo_merge_dwi']       = 0   # If you want to redo the processing from the merging of dwi files from different diffusion times set to 1
cfg['redo_denoise']         = 0   # If you want to redo the processing from the denoising set to 1
cfg['redo_gibbs']           = 0   # If you want to redo the processing from the gibbs unringing set to 1
cfg['redo_topup']           = 0   # If you want to redo the processing from the topup correction set to 1
cfg['redo_eddy']            = 0   # If you want to redo the processing from the eddy correction set to 1
cfg['redo_final_mask']      = 0   # If you want to redo the processing from the creation of the final brain masks set to 1

cfg['algo_denoising']       = 'matlab_tMPPCA_4D'  # Options are: 'matlab_MPPCA', or 'matlab_tMPPCA_4D' or 'matlab_tMPPCA_5D' or 'mrtrix_MPPCA' or 'designer_tMPPCA'. Note that designer sigma output map is not caculated the same as for the other methods
cfg['algo_brainextract']    = 'BET'             # Options are: 'BET' or 'RATS'
cfg['anat_thr']             = '4000'            # 2100, 4000 depending on your data
cfg['anat_format']          = 'T1w'             # Depends on you anatomical image. Common options are: 'T1w' or 'T2w'
cfg['subject_type']         = 'human'           # Options are: 'human' or 'rat'
cfg['is_alive']             = 'in_vivo'         # Options are: 'in_vivo' or 'ex_vivo'
cfg['individual_rev']       = 0                 # If there is one rev direction acquired for each diffusion time write 1, otherwise 0
cfg['topup_cfg_name']       = 'mycnf_fmri.cnf'  # name of the file with parameter details for topup (should be in the common folder)

#### DWI MODEL CONFIG ####
cfg['model_list_GM']        =  ['Nexi','Smex']  # List of model names to use for fitting the GM signal
cfg['model_list_WM']        =  []               # List of model names to use for fitting the WM signal
cfg['LTEDelta_for_microFA'] =  38               # Diffusion time (in ms) from the LTE (linear tensor encoding) acquisition that needs to be used together with STE (spherical tensor encoding) data to compute microFA
cfg['redo_modelling']       =  0                # If you want to redo the modelling set to 1, otherwise will redo just the models that didnt do before (because it crashed or something)

#### ROIS CONFIG ####
cfg['ROIs_GM']       = ['hippocampus','V1','V2','M1','premotor','parietal', 'S1', 'S2','Broca'] # List of ROIs to analyse (in GM) for Atlas_Juelich. Defined previously for each atlas in custom_functions. Please read instructions of Step3_registrations
cfg['ROIs_GM']       = ['frontal','precentral','postcentral','occipital','parietal', 'temporal'] # List of ROIs to analyse (in GM) for Atlas_Neuromorphometrics. Defined previously for each atlas in custom_functions. Please read instructions of Step3_registrations
cfg['ROIs_GM']       = ['frontal','precentral','postcentral','occipital','parietal', 'temporal'] # List of ROIs to analyse (in GM) for Atlas_DKT. Defined previously for each atlas in custom_functions. Please read instructions of Step3_registrations

cfg['ROIs_WM']       = []    # List of ROIs to analyse (in WM). Defined previously for each atlas in custom_functions. Please read instructions of Step3_registrations
cfg['tpm_thr']       = 0.2   # Threshold to be used for the tissue probability map (tpm) to define the different tissues

cfg['mrs_vx']        = 0                        # Does the dataset include mrs. 1 if yes, 0 if no
cfg['lat_ROIS']      = 0                        # Do you want to have ROIs in left and right hemispheres separately? 1 if yes, 0 if no. It requires adding a column VoxMidHem in the excel with the voxel of the middle plane that separates the hemisphere for each subject. It assumes a given orientation in the data order so it might not work for human and organoid data.

#### SAVE CONFIG FILE ####
cfg = update_cfg(cfg)


########################## DATA PROCESSING ##########################

#### STEP 2. NIFTI CONVERT SUBJECT  ####

subprocess.run( ["conda", "run", "-n", "niix2bids", "python", 
                 os.path.join(cfg['code_path'], 'processing_dwi','Step2_raw2nii2bids_human.py')] 
                + [cfg['data_path']] ,  check=True)


#### STEP 3. PREPROCESS SUBJECT ####

# 3.1 Process normal linear encoding data with PSGE (LTE), along with anatomical image
from Step3_preproc import *
Step3_preproc(subj_list,cfg) 

# 3.3 Register atlas to anatomical image and then to dwi for individual or combined diffusion times. 
# If exists also fits STE to LTE. Does not register across sessions for now
# Takes a long time. comment out if no regional estimations are needed
from Step3_registrations import *
Step3_registrations(subj_list, cfg)

# #### STEP 4. MODELLING SUBJECT ####

# 4.1 Fit the dwi signal with models like Nexi, Sandi, SMI, ....
# Diffusion Tensor (DTI) and Kurtusis (DKI) is always done by default. 
# To be faster and perform only DTI and DKI fitting just leave cfg['model_list_GM'] 
# and cfg['model_list_WM'] empty.
from Step4_modelling import *
Step4_modelling(subj_list,cfg)

# #### STEP 5. GET VALUES ####

# 5.1 Retreives parameter estimates from the model fits, making summary figures and excel with data in certain ROIs
# Needs the registration step 3 to be done before.
from Step5_get_estimates import *
Step5_get_estimates(subj_list,cfg) 
