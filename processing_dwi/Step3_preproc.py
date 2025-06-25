"""
Script to preprocess dMRI data.

Includes: 
- Denoising. several options exist
    Options are: 'matlab_MPPCA' - uses matlab and performs normal MPPCA on matrix with this format [x, y, z, bvals/bvecs x diffTime] (4D)
                 'matlab_tMPPCA_4D' - uses matlab and performs tMPPCA on matrix with this format [x, y, z, bvals/bvecs x diffTime] (4D)
                 'matlab_tMPPCA_5D'- uses matlab and performs tMPPCA on matrix with this format [x, y, z, bvals/bvecs, diffTime] (5D)
                 'mrtrix_MPPCA' - uses mrtrix to perform MPPCA on matrix with this format [x, y, z, bvals/bvecs x diffTime] (4D)
                 'designer_tMPPCA' - uses mrtrix to perform tMPPCA on matrix with this format [x, y, z, bvals/bvecs x diffTime] (4D)
    Note that designer sigma output map is not caculated on the software, so it's calculated here by hand but it's the same formula as for the other methods.
    Note that the implementation of MPPCA on mrtrix and on matlab are slightly different (matlab version denoises more the data and creates smoother maps)
    Note that for the 5D version you have to have a complete matrix (this is with the same number of diffusion times for each bvals/bvecs)
- Gibbs correction with mrtrix 
- Top up with FSL
- Eddy with FSL

It runs for the combined dataset (all diffusion times) or 
for each diffusion time alone.

It does not use a particular python environment (base).
It plots some quality analysis (QA) results for checking.

Last changed June 2025
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

def Step3_preproc(subj_list, cfg):
    
    data_path   = cfg['data_path']     
    scan_list   = pd.read_excel(os.path.join(data_path, cfg['scan_list_name']))
    topupon     = cfg['do_topup']
    anat_format = cfg['anat_format']

    ######## SUBJECT-WISE OPERATIONS ########
    for subj in subj_list:
    
        print('Preprocessing ' + subj + '...')
        
        # Copy nifti data to preprocessing folder
        nifti_path      = os.path.join(data_path, 'nifti_data', 'sorted', subj)
        preproc_path    = os.path.join(data_path, 'derivatives', cfg['prep_foldername'], subj)
        if not os.path.exists(preproc_path) or cfg['redo_all']:
            if os.path.exists(preproc_path):
                print("Your previous results will be deleted, are you sure? Press any key to continue or abort.")
                input()
                shutil.rmtree(preproc_path)
            
            shutil.copytree(nifti_path, preproc_path)
    
        # Extract data for subject
        subj_data      = scan_list[(scan_list['newstudyName'] == subj)].reset_index(drop=True)
        
        # List of acquisition sessions
        sess_list    = [x for x in list(subj_data['blockNo'].unique()) if not math.isnan(x)] # clean NaNs
    
        ######## SESSION-WISE OPERATIONS ########
        for sess in sess_list:
            
            print('Working on session ' + str(sess) + '...')

            ########################## ANATATOMICAL PROCESSING ##########################
            bids_strc_anat = create_bids_structure(subj=subj, sess=sess, datatype="anat", root=data_path, 
                                        folderlevel='derivatives', workingdir=cfg['prep_foldername'])
            # BET
            if not os.path.exists(bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz')) or cfg['redo_bet_anat']:
                print('Processing anatomical data')
                
                # Check that data exists
                if not os.path.exists(bids_strc_anat.get_path(f'{anat_format}.nii.gz')):
                    print("No anat scans found — exiting.")
                    return
                
                # If human data ensure it's oriented in a standard way
                if cfg['subject_type']=='human':
                    fsl_reorient(bids_strc_anat.get_path(f'{anat_format}.nii.gz'))
                    
                # Correct for bias field   
                N4_unbias(bids_strc_anat.get_path(f'{anat_format}.nii.gz'),bids_strc_anat.get_path(f'{anat_format}_bc.nii.gz'))
                
                # Brain extraction
                if cfg['subject_type']=='human':
                    brain_extract_BET(bids_strc_anat.get_path(f'{anat_format}_bc.nii.gz'))
                elif cfg['subject_type']=='rat':
                    brain_extract_RATS(bids_strc_anat.get_path(f'{anat_format}_bc.nii.gz'),cfg['anat_thr'])
                elif cfg['subject_type']=='organoid':
                    #brain_extract_organoids(bids_strc_anat.get_path(f'{anat_format}_bc.nii.gz'),cfg['anat_thr'])
                    
                    # make mask of the organoids themselves
                    make_mask(bids_strc_anat.get_path(f'{anat_format}_bc.nii.gz'), bids_strc_anat.get_path('organoids_mask.nii.gz'), 3000)
                    copy_files([bids_strc_anat.get_path('organoids_mask.nii.gz')],[bids_strc_anat.get_path(f'{anat_format}_bc_brain_mask.nii.gz')])    
                    call = [f'fslmaths',
                            f'{bids_strc_anat.get_path(f'{anat_format}_bc.nii.gz')}',
                            f'-mul {bids_strc_anat.get_path(f'{anat_format}_bc_brain_mask.nii.gz')}',
                            f'{bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz')}']
                    os.system(' '.join(call))
                    create_inverse_mask(bids_strc_anat.get_path('organoids_mask.nii.gz'),bids_strc_anat.get_path(f'{anat_format}_bc_brain_mask.nii.gz'),bids_strc_anat.get_path())
                    make_atlas_label_organoid(bids_strc_anat.get_path('organoids_mask.nii.gz'),bids_strc_anat.get_path('organoids_inv_mask.nii.gz'),bids_strc_anat.get_path(f"{cfg['atlas']}.label"))
                    

                # QA
                QA_brain_extract(bids_strc_anat.get_path(f'{anat_format}_bc.nii.gz'),os.path.join(bids_strc_anat.get_path(),'QA_brain_extract'),anat_format)

            ########################## DWI PROCESSING PREPARATION ##########################
            bids_strc = create_bids_structure(subj=subj, sess=sess, datatype="dwi", root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['prep_foldername'])
            
            print('Processing diffusion data')

            # Index of diff scans for this session 
            dwi_indices = np.where(
                (np.array(subj_data['acqType']) == 'PGSE') &
                (np.array(subj_data['scanQA']) == 'ok') &
                (np.array(subj_data['blockNo']) == sess))[0]

            # Check that data exists
            if dwi_indices.size == 0:
                   print("No dwi scans found — exiting.")
                   return
     
            # Generate paths for fwd and rev acquisition types
            paths_dwi_fwd=[]; paths_b0_fwd=[]; paths_dwi_rev=[]; paths_b0_rev=[]
            paths_bvals_fwd =[];paths_bvecs_fwd =[]; paths_bvals_rev =[]; paths_bvecs_rev =[]; 
            diffTimes = []; diffDurations =[];
            for scn_ctr in dwi_indices:    
                bids_strc.set_param(description='Delta_'+str(int(subj_data['diffTime'][scn_ctr]))+'_'+subj_data['phaseDir'][scn_ctr])
                if subj_data['phaseDir'][scn_ctr] == 'fwd' :
                    target_lists = (paths_dwi_fwd, paths_b0_fwd, paths_bvals_fwd, paths_bvecs_fwd)     
                    diffDurations.append(subj_data['diffDuration'][scn_ctr])
                    diffTimes.append(subj_data['diffTime'][scn_ctr])
                elif  subj_data['phaseDir'][scn_ctr] == 'rev':
                    target_lists = (paths_dwi_rev, paths_b0_rev, paths_bvals_rev, paths_bvecs_rev)
                                  
                for lst, suffix in zip(target_lists, ['dwi.nii.gz', 'b0.nii.gz', 'bvalsNom.txt', 'bvecs.txt']):
                    lst.append(bids_strc.get_path(suffix))
           
            # B0 extract for each scan and REGISTRATION of brain mask with avg B0 - fwd scans
            print('B0 extract for each scan and REGISTRATION of brain mask with avg B0 - fwd scans ...')
            masks_paths = []
            for kk in range(len(paths_dwi_fwd)):
                
                if not os.path.exists(paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg_mask.nii.gz')) or cfg['redo_b0_extract']:
                    # extract b0
                    extract_b0(    [paths_dwi_fwd[kk]], \
                                    [paths_dwi_fwd[kk].replace('dwi.nii.gz', 'bvecs.txt')], \
                                    [paths_dwi_fwd[kk].replace('dwi.nii.gz', 'bvalsNom.txt')], \
                                    [paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0.nii.gz')]    )
                    # average b0
                    make_avg(3, [paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0.nii.gz')], [paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg.nii.gz')])
                    
                    # bias field correct b0
                    N4_unbias(paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg.nii.gz'),paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg_bc.nii.gz'))
                    
                    if cfg['subject_type']=='organoid':
                         # do simple mask of the flask
                         #brain_extract_organoids(paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg_bc.nii.gz'),cfg['anat_thr'])
                         #make_mask(paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg_bc_brain.nii.gz'), paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg_mask.nii.gz'), 1e4)                

                         # pad image temporarily for registration
                         pad_image(bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'), bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'))
                         pad_image(bids_strc_anat.get_path(f'{anat_format}_bc.nii.gz'), bids_strc_anat.get_path(f'{anat_format}_bc.nii.gz'))
                         pad_image(paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg_bc.nii.gz'), paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg_bc.nii.gz'))
                         
                    # register dwi --> T2w
                    antsreg_full(bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'), # fixed
                            paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg_bc.nii.gz'),  # moving
                            paths_dwi_fwd[kk].replace('dwi.nii.gz', f'dwi2{anat_format}'))
                    
                    # apply inverse transform to put T2w in dwi space
                    ants_apply_transforms([bids_strc_anat.get_path(f'{anat_format}_bc.nii.gz'),bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz')],  # input 
                                          paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg_bc.nii.gz'), # moving
                                          [paths_dwi_fwd[kk].replace('dwi.nii.gz', f'{anat_format}_in_dwi.nii.gz'),paths_dwi_fwd[kk].replace('dwi.nii.gz', f'{anat_format}_brain_in_dwi.nii.gz')], # output
                                          [ paths_dwi_fwd[kk].replace('dwi.nii.gz', f'dwi2{anat_format}0GenericAffine.mat'), 1], # transform 1
                                          paths_dwi_fwd[kk].replace('dwi.nii.gz', f'dwi2{anat_format}1InverseWarp.nii.gz'))   # transform 2

                    if cfg['subject_type']=='organoid':
                        # unpad the images previousy padded
                        unpad_image(bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'), bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'))
                        unpad_image(bids_strc_anat.get_path(f'{anat_format}_bc.nii.gz'), bids_strc_anat.get_path(f'{anat_format}_bc.nii.gz'))
                        unpad_image(paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg_bc.nii.gz'), paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg_bc.nii.gz'))
                        unpad_image(paths_dwi_fwd[kk].replace('dwi.nii.gz', f'{anat_format}_in_dwi.nii.gz'), paths_dwi_fwd[kk].replace('dwi.nii.gz', f'{anat_format}_in_dwi.nii.gz'))
                        unpad_image(paths_dwi_fwd[kk].replace('dwi.nii.gz', f'{anat_format}_brain_in_dwi.nii.gz'), paths_dwi_fwd[kk].replace('dwi.nii.gz', f'{anat_format}_brain_in_dwi.nii.gz'))
                        unpad_image(paths_dwi_fwd[kk].replace('dwi.nii.gz', f'dwi2{anat_format}.nii.gz'), paths_dwi_fwd[kk].replace('dwi.nii.gz', f'dwi2{anat_format}.nii.gz'))

                    #else:
                    # create mask
                    make_mask(paths_dwi_fwd[kk].replace('dwi.nii.gz', f'{anat_format}_brain_in_dwi.nii.gz'), paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg_mask.nii.gz'), 100)                
                    
                    # QA
                    QA_reg(paths_dwi_fwd[kk].replace('dwi.nii.gz', f'{anat_format}_brain_in_dwi.nii.gz'),paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg_bc.nii.gz'),os.path.join(os.path.dirname(paths_dwi_fwd[kk]), 'QA_reg'))
                 
                 
                # save masks path to be used later
                masks_paths.append(paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg_mask.nii.gz'))

            # B0 extract for each scan and REGISTRATION of brain mask with avg B0 - rev scans
            print('B0 extract for each scan and REGISTRATION of brain mask with avg B0 - rev scans ...')
            for kk in range(len(paths_dwi_rev)):
     
                if not os.path.exists(paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg.nii.gz')) or cfg['redo_b0_extract']:
                  
                    # since rev is only b0, just copy the file and rename it
                    copy_files([paths_dwi_rev[kk]],[paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0.nii.gz')])
                   
                    # average b0 - if its just one volume it will be the same, but it's needed if more than one volume are acquired
                    make_avg(3, [paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0.nii.gz')], [paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg.nii.gz')])
                    
                    # bias field correct b0
                    N4_unbias(paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg.nii.gz'),paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg_bc.nii.gz'))
                    
                    if cfg['subject_type']=='organoid':
                         # do simple mask of the flask
                        # brain_extract_organoids(paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg_bc.nii.gz'),cfg['anat_thr'])
                         #make_mask(paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg_bc_brain.nii.gz'), paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg_mask.nii.gz'), 1e4)                

                         # pad image temporarily for registration
                         pad_image(bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'), bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'))
                         pad_image(bids_strc_anat.get_path(f'{anat_format}_bc.nii.gz'), bids_strc_anat.get_path(f'{anat_format}_bc.nii.gz'))
                         pad_image(paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg_bc.nii.gz'), paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg_bc.nii.gz'))
                         

                    # register dwi --> T2w
                    antsreg_full(bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'), # fixed
                            paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg_bc.nii.gz'),  # moving
                            paths_dwi_rev[kk].replace('dwi.nii.gz', f'dwi2{anat_format}'))
                    
                    # apply inverse transform to put T2w in dwi space
                    ants_apply_transforms([bids_strc_anat.get_path(f'{anat_format}_bc.nii.gz'),bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz')],  # input 
                                          paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg_bc.nii.gz'), # moving
                                          [paths_dwi_rev[kk].replace('dwi.nii.gz', f'{anat_format}_in_dwi.nii.gz'),paths_dwi_rev[kk].replace('dwi.nii.gz', f'{anat_format}_brain_in_dwi.nii.gz')], # output
                                          [paths_dwi_rev[kk].replace('dwi.nii.gz', f'dwi2{anat_format}0GenericAffine.mat'), 1], # transform 1
                                          paths_dwi_rev[kk].replace('dwi.nii.gz', f'dwi2{anat_format}1InverseWarp.nii.gz'))   # transform 2
               
                    if cfg['subject_type']=='organoid':
                        # unpad the images previousy padded
                        unpad_image(bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'), bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'))
                        unpad_image(bids_strc_anat.get_path(f'{anat_format}_bc.nii.gz'), bids_strc_anat.get_path(f'{anat_format}_bc.nii.gz'))
                        unpad_image(paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg_bc.nii.gz'), paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg_bc.nii.gz'))
                        unpad_image(paths_dwi_rev[kk].replace('dwi.nii.gz', f'{anat_format}_in_dwi.nii.gz'), paths_dwi_rev[kk].replace('dwi.nii.gz', f'{anat_format}_in_dwi.nii.gz'))
                        unpad_image(paths_dwi_rev[kk].replace('dwi.nii.gz', f'{anat_format}_brain_in_dwi.nii.gz'), paths_dwi_rev[kk].replace('dwi.nii.gz', f'{anat_format}_brain_in_dwi.nii.gz'))
                        unpad_image(paths_dwi_rev[kk].replace('dwi.nii.gz', f'dwi2{anat_format}.nii.gz'),paths_dwi_rev[kk].replace('dwi.nii.gz', f'dwi2{anat_format}.nii.gz'))

                    # create mask
                    make_mask(paths_dwi_rev[kk].replace('dwi.nii.gz', f'{anat_format}_brain_in_dwi.nii.gz'), paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg_mask.nii.gz'), 100)                
        
                    # QA
                    QA_reg(paths_dwi_rev[kk].replace('dwi.nii.gz', f'{anat_format}_brain_in_dwi.nii.gz'),paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg_bc.nii.gz'),os.path.join(os.path.dirname(paths_dwi_rev[kk]), 'QA_reg'))
                            
                       
                # save masks path to be used later
                masks_paths.append(paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg_mask.nii.gz'))

            # Define which data to process (all diffusion times together and individual scans)
            data_to_process = []
            for preproc_type in ['combined']: # if desidered add 'individual' to process each Delta individually
                              
                # Combine data from multiple diffusion times
                if preproc_type == 'combined':
                  
                    # Generate combined output path
                    bids_strc.set_param(description='allDelta-allb')
                    create_directory(bids_strc.get_path())
                    
                    if not os.path.exists(bids_strc.get_path('dwi.nii.gz')) or cfg['redo_merge_dwi']:
                        
                        # Combine niftis, bvals, bvecs
                        concat_niftis(paths_dwi_fwd,  bids_strc.get_path('dwi.nii.gz'), 'all')
                        concat_files(paths_bvals_fwd, bids_strc.get_path('bvalsNom.txt'))
                        paths_bvalsEff_fwd = [path.replace('Nom', 'Eff') for path in paths_bvals_fwd]
                        concat_files(paths_bvalsEff_fwd, bids_strc.get_path('bvalsEff.txt'))
                        concat_files(paths_bvecs_fwd, bids_strc.get_path('bvecs.txt')) 
                        
                        # Create the diffusion times and diffusion durations file
                        concat_param(np.array(diffTimes),paths_bvals_fwd,bids_strc.get_path('DiffTime.txt'))
                        concat_param(np.array(diffDurations),paths_bvals_fwd,bids_strc.get_path('DiffDuration.txt'))
                        
                        # Concatenate b0 files from both fwd and rev directions
                        if paths_b0_fwd:
                            concat_niftis(paths_b0_fwd, bids_strc.get_path('b0_fwd.nii.gz'), 1)
                        if paths_b0_rev:
                            concat_niftis(paths_b0_rev, bids_strc.get_path('b0_rev.nii.gz'), 1) 
                        
                        # QA
                        output_path = bids_strc.get_path();
                        QA_plotbvecs(bids_strc.get_path('bvecs.txt'), bids_strc.get_path('bvalsNom.txt'),os.path.join(output_path, 'QA_acquisition'))
    
                        # CREATE low b vals dataset (containing only low b values data for Nexi models eg) 
                        print('Create low b vals dataset...')
                        o_bids_strc = copy.deepcopy(bids_strc)
                        bids_strc.set_param(description='allDelta-lowb')
                        create_directory(bids_strc.get_path())
                        desiredbvals = [0, 1000]
                        old_dataset = {"dwi":   o_bids_strc.get_path('dwi.nii.gz'), 
                                        "bvals": o_bids_strc.get_path('bvalsNom.txt'),
                                        "bvecs": o_bids_strc.get_path('bvecs.txt')}
                        new_dataset = {"dwi":   bids_strc.get_path('dwi.nii.gz'), 
                                        "bvals": bids_strc.get_path('bvalsNom.txt'),
                                        "bvecs": bids_strc.get_path('bvecs.txt')}              
                        dwi_extract(old_dataset, new_dataset,','.join(map(str, desiredbvals)))
                        param_extract(o_bids_strc.get_path('DiffTime.txt'), o_bids_strc.get_path('bvalsNom.txt'), desiredbvals, bids_strc.get_path('DiffTime.txt'))
                        param_extract(o_bids_strc.get_path('DiffDuration.txt'), o_bids_strc.get_path('bvalsNom.txt'), desiredbvals, bids_strc.get_path('DiffDuration.txt'))
    
                    # DENOISE low b values dataset
                    if not os.path.exists(bids_strc.get_path('dwi_dn.nii.gz')) or cfg['redo_denoise']:
                        print('Denoise low b vals dataset...')
                        if cfg['algo_denoising']=='matlab_MPPCA':
                            denoise_matlab(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('DiffTime.txt'), cfg['code_path2'], cfg['toolboxes'],'MPPCA')
                        elif cfg['algo_denoising']=='mrtrix_MPPCA':
                             denoise_vols_default_kernel(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('dwi_dn_sigma.nii.gz'))                        
                        elif cfg['algo_denoising']=='matlab_tMPPCA_4D':
                            denoise_matlab(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('DiffTime.txt'), cfg['code_path2'], cfg['toolboxes'],'tMPPCA-4D')
                        elif cfg['algo_denoising']=='designer_tMPPCA':
                             denoise_designer(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('bvecs.txt'), bids_strc.get_path('bvalsNom.txt'), bids_strc.get_path('dwi_dn.nii.gz'), data_path, 'jespersen')                         
                        elif cfg['algo_denoising']=='matlab_tMPPCA_5D':
                            denoise_matlab(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('DiffTime.txt'), cfg['code_path2'], cfg['toolboxes'],'tMPPCA-5D')

                        # Calculates SNR
                        calc_snr(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn_sigma.nii.gz'),bids_strc.get_path('dwi_snr.nii.gz'))
                        
                        # QA
                        output_path = bids_strc.get_path();
                        QA_plotbvecs(bids_strc.get_path('bvecs.txt'), bids_strc.get_path('bvalsNom.txt'),os.path.join(output_path, 'QA_acquisition'))
                        QA_denoise(bids_strc, 'dwi_dn_res.nii.gz','dwi_dn_sigma.nii.gz',os.path.join(output_path, 'QA_denoise'))

                    # Generate combined output path
                    bids_strc.set_param(description='allDelta-allb')
                    data_to_process.append(bids_strc)
    
                # Choose individual datasets of single diffusion times
                elif preproc_type == 'individual':
                    print('Prepare each diffusion time dataset...')

                    # if individual, get the diffusion times (Deltas) of the individual datasets
                    filtered_data = subj_data[(subj_data['phaseDir'] == 'fwd') & (subj_data['blockNo'] == sess) & (subj_data['noBval'] > 1)]
                    Delta_list = filtered_data["diffTime"].dropna().astype(int).tolist()
                    delta_list = filtered_data["diffDuration"].dropna().astype(int).tolist()

                    for Delta, delta in zip(Delta_list,delta_list):
                        
                        # Create a file for the fwd direction (basically a copy)
                        fwd_bids_strc = copy.deepcopy(bids_strc)
                        fwd_bids_strc.set_param(description=f'Delta_{Delta}_fwd')
                        copy_file([fwd_bids_strc.get_path('b0.nii.gz')],[fwd_bids_strc.get_path('b0_fwd.nii.gz')])
                        
                        # Create a file for the rev direction depending on the acquired data:
                        rev_bids_strc = copy.deepcopy(bids_strc)
                            # If there is one reverse encoding for each diffusion time, use the specific one for each diffusion time
                        if cfg['individual_rev']==1:
                            rev_bids_strc.set_param(description=f'Delta_{Delta}_rev')
                            # If there is only one reverse encoding acquired for all diffusion times, use that for each diffusion time
                        else:
                            delta_rev = int((subj_data[(subj_data['phaseDir'] == 'rev') & (subj_data['blockNo'] == sess)]['diffTime']).iloc[0])
                            rev_bids_strc.set_param(description='Delta_'+str(delta_rev)+'_rev')
                        if os.path.exists(rev_bids_strc.get_path('b0.nii.gz')):
                            copy_file([rev_bids_strc.get_path('b0.nii.gz')],[fwd_bids_strc.get_path('b0_rev.nii.gz')])
                      
                        # Create the diffusion times and diffusion durations file
                        concat_param(np.array([Delta]),[fwd_bids_strc.get_path('bvalsNom.txt')],fwd_bids_strc.get_path('DiffTime.txt'))
                        concat_param(np.array([delta]),[fwd_bids_strc.get_path('bvalsNom.txt')],fwd_bids_strc.get_path('DiffDuration.txt'))

                        # QA
                        QA_plotbvecs(fwd_bids_strc.get_path('bvecs.txt'), fwd_bids_strc.get_path('bvalsNom.txt'), os.path.join( fwd_bids_strc.get_path(), 'QA_acquisition'))
                        data_to_process.append(fwd_bids_strc)
             
            ########################## DWI PROCESSING ##########################       
            for data in data_to_process:
                bids_strc = data
                output_path = bids_strc.get_path();
                print(f'Processing {os.path.basename(output_path)}....')

                # Create deformed mask
                if 'allDelta-allb' in output_path:
                    # If processing the combined dataset, merge the masks from the different diffusion times
                    union_niftis(masks_paths, bids_strc.get_path('mask_before_preproc.nii.gz'))
                else:
                    # If processing the individual datasets, just rename the mask file for consistency
                    copy_file([find_files_with_pattern(bids_strc,'b0_avg_mask')[0]],[ bids_strc.get_path('mask_before_preproc.nii.gz')])
                filter_clusters_by_size(bids_strc.get_path('mask_before_preproc.nii.gz'), bids_strc.get_path('mask_before_preproc.nii.gz'), 200)
                dilate_im(bids_strc.get_path('mask_before_preproc.nii.gz'), bids_strc.get_path('mask_before_preproc.nii.gz'), '1.5')
               
                # DENOISE
                if not os.path.exists(bids_strc.get_path('dwi_dn.nii.gz')) or cfg['redo_denoise']:
                    if cfg['algo_denoising']=='matlab_MPPCA':
                        denoise_matlab(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('DiffTime.txt'), cfg['code_path2'], cfg['toolboxes'],'MPPCA')
                    elif cfg['algo_denoising']=='mrtrix_MPPCA':
                        denoise_vols_default_kernel(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('dwi_dn_sigma.nii.gz'))
                    elif cfg['algo_denoising']=='matlab_tMPPCA_4D':
                        denoise_matlab(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('DiffTime.txt'), cfg['code_path2'], cfg['toolboxes'],'tMPPCA-4D')
                    elif cfg['algo_denoising']=='designer_tMPPCA':
                         denoise_designer(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('bvecs.txt'), bids_strc.get_path('bvalsNom.txt'), bids_strc.get_path('dwi_dn.nii.gz'), data_path, 'jespersen')
                    elif cfg['algo_denoising']=='matlab_tMPPCA_5D' and os.path.basename(output_path)=='allDelta-allb':
                        denoise_matlab(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('DiffTime.txt'), cfg['code_path2'], cfg['toolboxes'],'tMPPCA-5D')
                    elif cfg['algo_denoising']=='matlab_tMPPCA_5D' and not os.path.basename(output_path)=='allDelta-allb':
                        denoise_designer(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('bvecs.txt'), bids_strc.get_path('bvalsNom.txt'), bids_strc.get_path('dwi_dn.nii.gz'), data_path, 'jespersen')

                    calc_snr(bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('dwi_dn_sigma.nii.gz'),bids_strc.get_path('dwi_snr.nii.gz'))
                    QA_denoise(bids_strc, 'dwi_dn_res.nii.gz','dwi_dn_sigma.nii.gz',os.path.join(output_path, 'QA_denoise'))

                # GIBBS UNRINGING
                if not os.path.exists(bids_strc.get_path('dwi_dn_gc.nii.gz')) or cfg['redo_gibbs']:
                    gibbs_corr(bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('dwi_dn_gc.nii.gz'))
                    QA_gc(bids_strc, 'dwi_dn.nii.gz', 'dwi_dn_gc.nii.gz', os.path.join(output_path, 'QA_gc'))

                # TOPUP
                if (not os.path.exists(bids_strc.get_path('dwi_dn_gc_topup.nii.gz')) or cfg['redo_topup']) and topupon:
                    topup_routine(bids_strc.get_path('dwi_dn_gc.nii.gz'), bids_strc,  os.path.join(cfg['common_folder'],cfg['topup_cfg_name']))
                    QA_topup(bids_strc, 'dwi_dn_gc.nii.gz', 'dwi_dn_gc_topup.nii.gz', os.path.join(output_path, 'QA_topup'))

                # EDDY
                if not os.path.exists(bids_strc.get_path('dwi_dn_gc_ec.nii.gz')) or cfg['redo_eddy']:
                    eddy_routine(bids_strc.get_path('dwi_dn_gc.nii.gz'), bids_strc.get_path('dwi_dn_gc_ec.nii.gz'), bids_strc.get_path('mask_before_preproc.nii.gz'), bids_strc.get_path('bvalsNom.txt'), bids_strc.get_path('bvecs.txt'), topupon) # added the corr
                    extract_b0([bids_strc.get_path('dwi_dn_gc_ec.nii.gz')], [bids_strc.get_path('dwi_dn_gc_ec.eddy_rotated_bvecs')], \
                                [bids_strc.get_path('bvalsNom.txt')], [bids_strc.get_path('b0_dn_gc_ec.nii.gz')]) # check this image for motion
                   
                # Generate non-deformed masks
                if (not os.path.exists(bids_strc.get_path('mask.nii.gz')) or cfg['redo_final_mask']) and os.path.exists(bids_strc.get_path('dwi_dn_gc_ec.nii.gz')):
    
                    # average b0
                    make_avg(3, [bids_strc.get_path('b0_dn_gc_ec.nii.gz')], [bids_strc.get_path('b0_dn_gc_ec_avg.nii.gz')])
                    #threshold_image(bids_strc.get_path('b0_dn_gc_ec_avg.nii.gz'), bids_strc.get_path('mask.nii.gz'), 1e4, 9e4)
                    
                    # bias field correct b0
                    N4_unbias(bids_strc.get_path('b0_dn_gc_ec_avg.nii.gz'),bids_strc.get_path('b0_dn_gc_ec_avg_bc.nii.gz'))
                       
                    #if cfg['subject_type']=='organoid':
                         # do simple mask of the flask
                         #brain_extract_organoids(bids_strc.get_path('b0_dn_gc_ec_avg_bc.nii.gz'),cfg['anat_thr'])
                         #make_mask(bids_strc.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'), bids_strc.get_path('mask.nii.gz'), 0.8e4)                

                    # get b0 avg bias field correct brain
                    binary_op(bids_strc.get_path('b0_dn_gc_ec_avg_bc.nii.gz'),bids_strc.get_path('mask_before_preproc.nii.gz'), '-mul', bids_strc.get_path('b0_dn_gc_ec_avg_bc_brain_before_preproc.nii.gz'))

                    if cfg['subject_type']=='organoid':
                        
                         # pad image temporarily for registration
                         pad_image(bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'), bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'))
                         pad_image(bids_strc.get_path('b0_dn_gc_ec_avg_bc_brain_before_preproc.nii.gz'), bids_strc.get_path('b0_dn_gc_ec_avg_bc_brain_before_preproc.nii.gz'))
  
                         # register dwi --> T2w
                         antsreg_full(bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'), # fixed
                                bids_strc.get_path('b0_dn_gc_ec_avg_bc_brain_before_preproc.nii.gz'),  # moving
                                bids_strc.get_path(f'dwiafterpreproc2{anat_format}'))
                        
                         # apply inverse transform to put T2w in dwi space
                         ants_apply_transforms([bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz')],  # input 
                                              bids_strc.get_path('b0_dn_gc_ec_avg_bc_brain_before_preproc.nii.gz'), # moving
                                              [bids_strc.get_path(f'{anat_format}_brain_in_dwiafterpreproc.nii.gz')], # output
                                              [bids_strc.get_path(f'dwiafterpreproc2{anat_format}0GenericAffine.mat'), 1], # transform 1
                                              bids_strc.get_path(f'dwiafterpreproc2{anat_format}1InverseWarp.nii.gz'))   # transform 2
                                                            
                         # unpad the images previousy padded
                         unpad_image(bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'), bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'))
                         unpad_image(bids_strc.get_path('b0_dn_gc_ec_avg_bc_brain_before_preproc.nii.gz'), bids_strc.get_path('b0_dn_gc_ec_avg_bc_brain_before_preproc.nii.gz'))
                         unpad_image(bids_strc.get_path(f'{anat_format}_brain_in_dwiafterpreproc.nii.gz'), bids_strc.get_path(f'{anat_format}_brain_in_dwiafterpreproc.nii.gz'))
                         unpad_image(bids_strc.get_path(f'dwiafterpreproc2{anat_format}.nii.gz'), bids_strc.get_path(f'dwiafterpreproc2{anat_format}.nii.gz'))

                    else:
                        # register dwi --> T2w
                        antsreg_simple(bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'), # fixed
                                bids_strc.get_path('b0_dn_gc_ec_avg_bc_brain_before_preproc.nii.gz'),  # moving
                                bids_strc.get_path(f'dwiafterpreproc2{anat_format}'))
                        
                        # apply inverse transform to put T2w in dwi space
                        ants_apply_transforms_simple([bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz')],  # input 
                                              bids_strc.get_path('b0_dn_gc_ec_avg_bc_brain_before_preproc.nii.gz'), # moving
                                              [bids_strc.get_path(f'{anat_format}_brain_in_dwiafterpreproc.nii.gz')], # output
                                              [bids_strc.get_path(f'dwiafterpreproc2{anat_format}0GenericAffine.mat'), 1]) # transform 1
                                  
                    # make mask
                    make_mask(bids_strc.get_path(f'{anat_format}_brain_in_dwiafterpreproc.nii.gz'), bids_strc.get_path('mask.nii.gz'), 0)                
                        #filter_clusters_by_size(bids_strc.get_path('mask.nii.gz'), bids_strc.get_path('mask.nii.gz'), 200)
                    
                    dilate_im(bids_strc.get_path('mask.nii.gz'), bids_strc.get_path('mask_dil.nii.gz'), '1')
    
                    # get brain with good non-deformed mask
                    binary_op(bids_strc.get_path('b0_dn_gc_ec_avg_bc.nii.gz'),bids_strc.get_path('mask.nii.gz'), '-mul', bids_strc.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'))

                    # QA eddy
                    QA_eddy(bids_strc.get_path('mask.nii.gz'),bids_strc.get_path('mask_dil.nii.gz'), bids_strc.get_path('dwi_dn_gc.nii.gz'), bids_strc.get_path('dwi_dn_gc_ec.nii.gz'), os.path.join(output_path, 'QA_eddy'),bids_strc.get_path('bvalsNom.txt'),bids_strc)
                    QA_DTI_fit(bids_strc.get_path('dwi_dn_gc.nii.gz'), bids_strc.get_path('bvalsNom.txt'), bids_strc.get_path('bvecs.txt'), bids_strc.get_path('mask.nii.gz'), os.path.join(output_path, 'QA_dti_before_eddy'))
                    QA_DTI_fit(bids_strc.get_path('dwi_dn_gc_ec.nii.gz'), bids_strc.get_path('bvalsNom.txt'), bids_strc.get_path('dwi_dn_gc_ec.eddy_rotated_bvecs'), bids_strc.get_path('mask.nii.gz'),os.path.join(output_path, 'QA_dti_after_eddy'))
                    QA_plotSNR(bids_strc,'dwi.nii.gz', 'dwi_snr.nii.gz', 'dwi_dn_sigma.nii.gz', 'mask.nii.gz', 'bvalsNom.txt',os.path.join(output_path, 'QA_acquisition'))
                    QA_mask(bids_strc.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'), bids_strc.get_path('mask.nii.gz'), bids_strc.get_path('mask_dil.nii.gz'),bids_strc.get_path('mask_before_preproc.nii.gz'),os.path.join(output_path, 'QA_mask'))

                    # Convert to mif in case
                    nifti_to_mif(bids_strc.get_path('dwi_dn_gc_ec.nii.gz'), bids_strc.get_path('dwi_dn_gc_ec.eddy_rotated_bvecs'), bids_strc.get_path('bvalsNom.txt'), bids_strc.get_path('dwi_dn_gc_ec.mif'))

                # Uncombine data for each diffusion time
                if not os.path.exists(os.path.join(os.path.dirname(bids_strc.get_path('bvalsNom.txt')), f'Delta_{int(diffTimes[0])}')) or cfg['redo_final_mask']:
                    unconcat_files(bids_strc.get_path('bvalsNom.txt'),bids_strc.get_path('DiffTime.txt'))
                    unconcat_files(bids_strc.get_path('bvalsEff.txt'),bids_strc.get_path('DiffTime.txt'))
                    unconcat_files(bids_strc.get_path('dwi_dn_gc_ec.eddy_rotated_bvecs'),bids_strc.get_path('DiffTime.txt'))
                    unconcat_files(bids_strc.get_path('DiffTime.txt'),bids_strc.get_path('DiffTime.txt'))
                    unconcat_files(bids_strc.get_path('DiffDuration.txt'),bids_strc.get_path('DiffTime.txt'))
                    unconcat_niftis(bids_strc.get_path('dwi_dn_gc_ec.nii.gz'),bids_strc.get_path('DiffTime.txt'))



                plt.close('all') 
               
                make_summary_pdf(bids_strc.get_path(),bids_strc.get_path('summary.pdf'))
                    
              