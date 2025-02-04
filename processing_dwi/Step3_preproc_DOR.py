"""
Script to convert preprocess dMRI data.
Includes: denoising, gibbs correction, top up, eddy
It can be run for the combined dataset (all diffusion times) or 
for one diffusion time only (the longest, for WM fitting).
It does not use a particular python environment (base).

Last changed Jan 2025
@author: Rita O
"""

import os
import sys
import pandas as pd
import platform
import math
import importlib, sys
import numpy as np
from custom_functions import *
from bids_structure import *
import copy
import glob
import matplotlib.pyplot as plt
import shutil
import keyboard
plt.close('all')

def Step3_preproc_DOR(subj_list, cfg):
    
    data_path   = cfg['data_path']     
    scan_list   = pd.read_excel(os.path.join(data_path, cfg['scan_list_name']))
    topupon     = cfg['do_topup'] 
    
    # update cfg in case one of the steps was 1 and the next ones 0, it changes the subsequents to 1 as well
    update_cfg(cfg)
    
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in subj_list:
    
        print('Preprocessing ' + subj + '...')
        
        # Extract data for subject
        subj_data      = scan_list[(scan_list['newstudyName'] == subj)].reset_index(drop=True)
        
        # List of acquisition sessions
        sess_list    = [x for x in list(subj_data['blockNo'].unique()) if not math.isnan(x)] # clean NaNs
    
        ######## SESSION-WISE OPERATIONS ########
        for sess in sess_list:
          
            # Define anat bids structure
            bids_strc_anat = create_bids_structure(subj=subj, sess=sess, datatype="anat", root=data_path, 
                                        folderlevel='derivatives', workingdir=cfg['prep_foldername'])
            
            ###### DWI SCAN-WISE OPERATIONS ######
            bids_strc = create_bids_structure(subj=subj, sess=sess, datatype="dwi_DOR", root=data_path, 
                                         folderlevel='derivatives', workingdir=cfg['prep_foldername'])
            # Index of diff scans for this session 
            dwi_indices = np.where(
                (np.array(subj_data['acqType']) == 'DOR') &
                (np.array(subj_data['scanQA']) == 'ok') &
                (np.array(subj_data['blockNo']) == sess))[0]

            # Generate paths for fwd and rev acquisition types
            masks_paths = []; paths_to_process = []; paths_b0_fwd =[];  paths_dwi_fwd = []; paths_b0_rev =[]; paths_dwi_rev =[]; 
            for scn_ctr in dwi_indices:    
                bids_strc.set_param(description=subj_data['phaseDir'][scn_ctr])
                if subj_data['phaseDir'][scn_ctr] == 'fwd' :
                    paths_dwi_fwd = bids_strc.get_path('dwi.nii.gz')
                    paths_b0_fwd = bids_strc.get_path('b0.nii.gz')
                elif subj_data['phaseDir'][scn_ctr] == 'rev':
                    paths_dwi_rev = bids_strc.get_path('dwi.nii.gz')
                    paths_b0_rev = bids_strc.get_path('b0.nii.gz')
                                  
            # B0 extract for each scan and REGISTRATION of brain mask with avg B0
            if paths_dwi_fwd:
                paths_to_process = [paths_dwi_fwd]
            else:
                paths_to_process = []
            if paths_dwi_rev:
                paths_to_process.append(paths_dwi_rev)
            for paths_dwi in paths_to_process:
                if paths_dwi:
                    
                    # Extract b0 volume
                    extract_vols(paths_dwi, paths_dwi.replace('dwi.nii.gz', 'b0.nii.gz'), 0, 1)
                    # Bias field correct b0
                    N4_unbias(paths_dwi.replace('dwi.nii.gz', 'b0.nii.gz'), paths_dwi.replace('dwi.nii.gz', 'b0_bc.nii.gz'))
            
                    # Register dwi --> T2w
                    antsreg(bids_strc_anat.get_path('T2w.nii.gz'),  # fixed
                            paths_dwi.replace('dwi.nii.gz', 'b0_bc.nii.gz'),  # moving
                            paths_dwi.replace('dwi.nii.gz', 'dwi2T2w'))
            
                    # Apply inverse transform to put T2w in dwi space
                    ants_apply_transforms([bids_strc_anat.get_path('T2w.nii.gz'),
                                           bids_strc_anat.get_path('T2w_brain.nii.gz')],  # input 
                                        paths_dwi.replace('dwi.nii.gz', 'b0_bc.nii.gz'),  # moving
                                        [paths_dwi.replace('dwi.nii.gz', 'T2w_in_dwi.nii.gz'),
                                         paths_dwi.replace('dwi.nii.gz', 'T2w_brain_in_dwi.nii.gz')],  # output
                                        [paths_dwi.replace('dwi.nii.gz', 'dwi2T2w0GenericAffine.mat'), 1],  # transform 1
                                        paths_dwi.replace('dwi.nii.gz', 'dwi2T2w1InverseWarp.nii.gz'))  # transform 2
            
                    # QA
                    QA_reg(paths_dwi.replace('dwi.nii.gz', 'T2w_brain_in_dwi.nii.gz'),
                           paths_dwi.replace('dwi.nii.gz', 'b0_bc.nii.gz'),
                           os.path.join(os.path.dirname(paths_dwi), 'QA_reg'))
            
                    # Create mask
                    make_mask(paths_dwi.replace('dwi.nii.gz', 'T2w_brain_in_dwi.nii.gz'),
                              paths_dwi.replace('dwi.nii.gz', 'b0_mask.nii.gz'), 100)
                    
                    # save masks path to be used later
                    masks_paths.append(paths_dwi_fwd.replace('dwi.nii.gz', 'b0_mask.nii.gz'))

            ###### DWI COMBINED OPERATIONS ######
            if paths_to_process: # only analyse if there is that type of data
                
                # Set output path
                bids_strc.set_param(description='fwd')
                create_directory(bids_strc.get_path())
                    
                # Create deformed mask
                union_niftis(masks_paths, bids_strc.get_path('mask.nii.gz'))
                filter_clusters_by_size(bids_strc.get_path('mask.nii.gz'), bids_strc.get_path('mask.nii.gz'), 200)
                dilate_im(bids_strc.get_path('mask.nii.gz'), bids_strc.get_path('mask_dil.nii.gz'), '1.5')
    
                # DENOISE
                if not os.path.exists(bids_strc.get_path('dwi_dn.nii.gz')) or cfg['redo_denoise']:
                    denoise_vols_default_kernel(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('dwi_dn_sigma.nii.gz'))
                    calc_noise_floor(bids_strc.get_path('dwi_dn_sigma.nii.gz'), bids_strc.get_path('dwi_nf.nii.gz') ) # added by rita
                    calc_snr(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn_sigma.nii.gz'),bids_strc.get_path('dwi_snr.nii.gz'))
    
                # GIBBS UNRINGING
                if not os.path.exists(bids_strc.get_path('dwi_dn_gc.nii.gz')) or cfg['redo_gibbs']:
                    gibbs_corr(bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('dwi_dn_gc.nii.gz'))
        
                # TOPUP
                if (not os.path.exists(bids_strc.get_path('dwi_dn_gc_topup.nii.gz')) or cfg['redo_topup']) and topupon and paths_dwi_rev:
                    topup_routine(paths_b0_fwd, paths_b0_rev, bids_strc.get_path('dwi_dn_gc.nii.gz'), bids_strc,  os.path.join(cfg['common_folder'],'mycnf_fmri.cnf'))
                    QA_topup(bids_strc, 'dwi_dn_gc.nii.gz', 'dwi_dn_gc_topup.nii.gz', os.path.join(output_path, 'QA_topup'))
    
                # Quality analysis
                output_path = bids_strc.get_path();
                QA_plotSNR(bids_strc, 'dwi_snr.nii.gz', 'dwi_nf.nii.gz', 'mask.nii.gz', 'bvalsNom.txt',os.path.join(output_path, 'QA_acquisition'))
                QA_denoise(bids_strc, 'dwi_dn_res.nii.gz','dwi_dn_sigma.nii.gz',os.path.join(output_path, 'QA_denoise'))
                QA_gc(bids_strc, 'dwi_dn.nii.gz', 'dwi_dn_gc.nii.gz', os.path.join(output_path, 'QA_gc'))
                plt.close('all')   
                    
