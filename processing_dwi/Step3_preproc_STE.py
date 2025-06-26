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

def Step3_preproc_STE(subj_list, cfg):
    
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
        
        # Check that data exists
        if not np.any(np.array(subj_data['acqType']) == 'STE'):
                print("No dwi scans with STE found — exiting.")
                return  
        
        # Copy nifti data to preprocessing folder
        nifti_path      = os.path.join(data_path, 'nifti_data', 'sorted', subj)
        preproc_path    = os.path.join(data_path, 'derivatives', cfg['prep_foldername'], subj)
        for sess in sess_list:
            preproc_path_sess    = os.path.join(data_path, 'derivatives', cfg['prep_foldername'], subj, f"ses-{sess:02}",'dwi_STE')
            nifti_path_sess      = os.path.join(data_path, 'nifti_data', 'sorted', subj, f"ses-{sess:02}",'dwi_STE')

            if not os.path.exists(preproc_path_sess) or cfg['redo_all']:
                if os.path.exists(preproc_path_sess):
                    print("Your previous results will be deleted, are you sure? Press any key to continue or abort.")
                    input()
                    shutil.rmtree(preproc_path_sess)
                
                shutil.copytree(nifti_path_sess, preproc_path_sess)
            
      
        ######## SESSION-WISE OPERATIONS ########
        for sess in sess_list:
          
            print('Working on session ' + str(sess) + '...')

            # Define anat bids structure
            bids_strc_anat = create_bids_structure(subj=subj, sess=sess, datatype="anat", root=data_path, 
                                        folderlevel='derivatives', workingdir=cfg['prep_foldername'])
            
            ########################## DWI PROCESSING PREPARATION ##########################
            bids_strc = create_bids_structure(subj=subj, sess=sess, datatype="dwi_STE", root=data_path, 
                                         folderlevel='derivatives', workingdir=cfg['prep_foldername'])
            # Index of diff scans for this session 
            dwi_indices = np.where(
                (np.array(subj_data['acqType']) == 'STE') &
                (np.array(subj_data['scanQA']) == 'ok') &
                (np.array(subj_data['blockNo']) == sess))[0]

            # Generate paths for fwd and rev acquisition types
            masks_paths = []; paths_to_process = []; paths_b0_fwd =[];  paths_dwi_fwd = []; paths_b0_rev =[]; paths_dwi_rev =[]; 
            for scn_ctr in dwi_indices:    
                bids_strc.set_param(description='STE_' + subj_data['phaseDir'][scn_ctr])
                if subj_data['phaseDir'][scn_ctr] == 'fwd' :
                    paths_dwi_fwd.append(bids_strc.get_path('dwi.nii.gz'))
                    paths_b0_fwd.append(bids_strc.get_path('b0.nii.gz'))
                elif subj_data['phaseDir'][scn_ctr] == 'rev':
                    paths_dwi_rev.append(bids_strc.get_path('dwi.nii.gz'))
                    paths_b0_rev.append(bids_strc.get_path('b0.nii.gz'))
                                  
            # B0 extract for each scan and REGISTRATION of brain mask with avg B0 - fwd scans
            for paths_dwi in paths_dwi_fwd:
                if paths_dwi and not os.path.exists(paths_dwi.replace('dwi.nii.gz', 'b0_mask.nii.gz')):
                    
                    # Extract b0 volume
                    extract_vols(paths_dwi, paths_dwi.replace('dwi.nii.gz', 'b0.nii.gz'), 0, 1)
                    
                    # Bias field correct b0
                    N4_unbias(paths_dwi.replace('dwi.nii.gz', 'b0.nii.gz'), paths_dwi.replace('dwi.nii.gz', 'b0_bc.nii.gz'))
            
                    # Register dwi --> T2w
                    antsreg_full(bids_strc_anat.get_path('T2w_bc.nii.gz'),  # fixed
                            paths_dwi.replace('dwi.nii.gz', 'b0_bc.nii.gz'),  # moving
                            paths_dwi.replace('dwi.nii.gz', 'dwi2T2w'))
            
                    # Apply inverse transform to put T2w in dwi space
                    ants_apply_transforms([bids_strc_anat.get_path('T2w_bc.nii.gz'),
                                           bids_strc_anat.get_path('T2w_bc_brain.nii.gz')],  # input 
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
                masks_paths.append(paths_dwi.replace('dwi.nii.gz', 'b0_mask.nii.gz'))
           
            # B0 extract for each scan - rev scans
            for paths_dwi in paths_dwi_rev:
                if paths_dwi and not os.path.exists(paths_dwi.replace('dwi.nii.gz', 'b0_mask.nii.gz')):
                    
                    # Extract b0 volume
                    extract_vols(paths_dwi, paths_dwi.replace('dwi.nii.gz', 'b0.nii.gz'), 0, 1)
            
            # Copy files to working folder
            bids_strc.set_param(description='STE_fwd')
            if paths_dwi_fwd:
                concat_niftis(paths_b0_fwd, bids_strc.get_path('b0_fwd.nii.gz'), 1)
            if paths_dwi_rev:
                 concat_niftis(paths_b0_rev, bids_strc.get_path('b0_rev.nii.gz'), 1) # assumes only one B0 value was collected in rev direction 
                    
            ########################## DWI PROCESSING ##########################       
            if paths_dwi_fwd: # only analyse if there is that type of data
                
                # Set output path
                bids_strc.set_param(description='STE_fwd')
                create_directory(bids_strc.get_path())
                output_path = bids_strc.get_path();
         
                # Create deformed mask
                union_niftis(masks_paths, bids_strc.get_path('mask_before_preproc.nii.gz'))
                #filter_clusters_by_size(bids_strc.get_path('mask.nii.gz'), bids_strc.get_path('mask.nii.gz'), 200)
    
                # DENOISE
                if not os.path.exists(bids_strc.get_path('dwi_dn.nii.gz')) or cfg['redo_denoise']:
                    if cfg['algo_denoising']=='matlab_MPPCA':
                        denoise_matlab(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('DiffTime.txt'), cfg['code_path2'], cfg['toolboxes'],'MPPCA')
                    elif cfg['algo_denoising']=='mrtrix_MPPCA':
                        denoise_vols_default_kernel(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('dwi_dn_sigma.nii.gz'))
                    elif cfg['algo_denoising']=='matlab_tMPPCA_4D':
                        denoise_matlab(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('DiffTime.txt'), cfg['code_path2'], cfg['toolboxes'],'tMPPCA-4D')
                    elif cfg['algo_denoising']=='designer_tMPPCA' or cfg['algo_denoising']=='matlab_tMPPCA_5D':
                         denoise_designer(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('bvecs_fake.txt'), bids_strc.get_path('bvalsNom.txt'), bids_strc.get_path('dwi_dn.nii.gz'), data_path, 'jespersen')

                    calc_snr(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn_sigma.nii.gz'),bids_strc.get_path('dwi_snr.nii.gz'))
                    QA_denoise(bids_strc, 'dwi_dn_res.nii.gz','dwi_dn_sigma.nii.gz',os.path.join(output_path, 'QA_denoise'))

                # GIBBS UNRINGING
                if not os.path.exists(bids_strc.get_path('dwi_dn_gc.nii.gz')) or cfg['redo_gibbs']:
                    gibbs_corr(bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('dwi_dn_gc.nii.gz'))
                    QA_gc(bids_strc, 'dwi_dn.nii.gz', 'dwi_dn_gc.nii.gz', os.path.join(output_path, 'QA_gc'))

                # TOPUP
                if (not os.path.exists(bids_strc.get_path('dwi_dn_gc_topup.nii.gz')) or cfg['redo_topup']) and topupon and paths_dwi_rev:
                    topup_routine(bids_strc.get_path('dwi_dn_gc.nii.gz'), bids_strc,  os.path.join(cfg['common_folder'],'mycnf_fmri.cnf'))
                    QA_topup(bids_strc, 'dwi_dn_gc.nii.gz', 'dwi_dn_gc_topup.nii.gz', os.path.join(output_path, 'QA_topup'))
                    
                    extract_vols(bids_strc.get_path('dwi_dn_gc_topup.nii.gz'), bids_strc.get_path('b0_dn_gc_topup.nii.gz'), 0, 1)

                    # Generate non-deformed masks
                    if not os.path.exists(bids_strc.get_path('mask.nii.gz')) or cfg['redo_final_mask']:
                       
                        # average b0
                        make_avg(3, [bids_strc.get_path('b0_dn_gc_topup.nii.gz')], [bids_strc.get_path('b0_dn_gc_topup_avg.nii.gz')])
                        
                        # bias field correct b0
                        N4_unbias(bids_strc.get_path('b0_dn_gc_topup_avg.nii.gz'),bids_strc.get_path('b0_dn_gc_topup_avg_bc.nii.gz'))
 
                        # get brain
                        binary_op(bids_strc.get_path('b0_dn_gc_topup_avg_bc.nii.gz'),bids_strc.get_path('mask_before_preproc.nii.gz'), '-mul', bids_strc.get_path('b0_dn_gc_topup_avg_bc_brain_before_preproc.nii.gz'))

                        # register dwi --> T2w
                        antsreg_simple(bids_strc_anat.get_path('T2w_bc_brain.nii.gz'), # fixed
                                bids_strc.get_path('b0_dn_gc_topup_avg_bc_brain_before_preproc.nii.gz'),  # moving
                                bids_strc.get_path('dwiafterpreproc2T2w'))
                        
                        # apply inverse transform to put T2w in dwi space
                        ants_apply_transforms_simple([bids_strc_anat.get_path('T2w_bc_brain.nii.gz')],  # input 
                                              bids_strc.get_path('b0_dn_gc_topup_avg_bc_brain_before_preproc.nii.gz'), # moving
                                              [bids_strc.get_path('T2w_brain_in_dwiafterpreproc.nii.gz')], # output
                                              [bids_strc.get_path('dwiafterpreproc2T2w0GenericAffine.mat'), 1]) # transform 1
         
                        
                        make_mask(bids_strc.get_path('T2w_brain_in_dwiafterpreproc.nii.gz'), bids_strc.get_path('mask.nii.gz'), 100)                
                        #filter_clusters_by_size(bids_strc.get_path('mask.nii.gz'), bids_strc.get_path('mask.nii.gz'), 200)
                        dilate_im(bids_strc.get_path('mask.nii.gz'), bids_strc.get_path('mask_dil.nii.gz'), '1')
        
                        # update mask brain
                        binary_op(bids_strc.get_path('b0_dn_gc_topup_avg_bc.nii.gz'),bids_strc.get_path('mask.nii.gz'), '-mul', bids_strc.get_path('b0_dn_gc_topup_avg_bc_brain.nii.gz'))

                    # Convert to mif in case
                    nifti_to_mif(bids_strc.get_path('dwi_dn_gc_topup.nii.gz'), bids_strc.get_path('bvecs_fake.txt'), bids_strc.get_path('bvalsNom.txt'), bids_strc.get_path('dwi_dn_gc_topup.mif'))

                # Quality analysis
                output_path = bids_strc.get_path();
                QA_plotSNR(bids_strc,'dwi.nii.gz', 'dwi_snr.nii.gz', 'dwi_dn_sigma.nii.gz', 'mask.nii.gz', 'bvalsNom.txt',os.path.join(output_path, 'QA_acquisition'))
                plt.close('all')
                make_summary_pdf(bids_strc.get_path(),bids_strc.get_path('summary.pdf'))

                
                    
