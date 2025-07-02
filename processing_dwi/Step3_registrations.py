#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to do registrations.
A. Registers atlas and tissue probability map (TPM) to anatomical space. Then it registers that to diffusion space.
B. Register sperical tensor encoding (STE) to one of the linear tensor encoding (LTE), defined in cfg which Diffusion time of LTE is used.

If new atlas arrives and something needs to be done to them to adjust it, please edit prepare_atlas in custom_functions.
Otherwise it assumes that:
    - In the atlas folder there are two nifiti files containing the strings: 
            '*atlas*', '*template_brain*'
    - In the TPM folder there are two nifiti files containing the strings: 
            '*TPM*' and '*template_brain*'. 
 
TPM is a 5D volume containing the tissue probability maps for the following tissues in order:
    GM (1), WM (2), CSF (3). The 4th and 5th dimension usually correspond to skull and outside brain but migth be different depending on the atlas of TPM used.
    
You don't need to have both an atlas and a TPM to do this registration. 
    It will register any of these if the corresponding fields in config file 
    (cfg['atlas'] and cfg['atlas_TPM']) are not empty.

Last changed April 2025
@author: Rita O
"""


import os
import sys
import pandas as pd
import platform
import math
import importlib, sys
from custom_functions import *
from bids_structure import *
import glob
import numpy as np
import SimpleITK as sitk
import numpy.ma as ma
import matplotlib.pyplot as plt
import nibabel as nib
import imutils
import nibabel
import nibabel.processing
import csv
import copy

plt.close('all')

def Step3_registrations(subj_list, cfg):
    
    
    data_path   = cfg['data_path']     
    scan_list   = pd.read_excel(os.path.join(data_path, cfg['scan_list_name'] ))
    atlas_path  = cfg['common_folder']
    anat_format = cfg['anat_format']
        
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in subj_list:
        
        print('Doing registration on ' + subj + '...')
    
        # Extract data for subject
        subj_data    = scan_list[scan_list['newstudyName'] == subj].reset_index(drop=True)
        
        # List of acquisition sessions
        sess_list    = [x for x in list(subj_data['blockNo'].unique()) if not math.isnan(x)] # clean NaNs
        
        ######## SESSION-WISE OPERATIONS ########
        for sess in sess_list:

           ########################## A. REGISTRATION OF ATLAS TO DWI ##########################
           for dossier in  [cfg['atlas_TPM'], cfg['atlas']]:
         
               if dossier:
                    print('Registering ' + dossier + ' ...')
                   
                    ########################## 1. PREPARATION OF ATLAS/TPM ##########################
                    ## Please edit function prepare_atlas to add as much conditions as atlas used!
                    if 'TPM' in dossier:
                        atlas, template = prepare_atlas(dossier, cfg['common_folder'],'TPM')
                    elif 'Atlas' in dossier:
                        atlas, template = prepare_atlas(dossier, cfg['common_folder'],'atlas')
                    elif 'anat_space_organoids' in dossier:
                        bids = create_bids_structure(subj=subj, sess=sess, datatype='anat', root=data_path, 
                                                   folderlevel='derivatives', workingdir=cfg['prep_foldername'])
                        template = bids.get_path(f'{anat_format}_bc_brain.nii.gz')
                        atlas = bids.get_path('organoids_atlas.nii.gz')

                    ########################## 2. REGISTRATION ATLAS TO ANAT ##########################
    
                    # Create BIDS structure
                    bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=dossier+f'-To-{anat_format}', root=data_path, 
                                               folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                    bids_strc_reg.set_param(base_name='')
                    bids_strc_anat = create_bids_structure(subj=subj, sess=sess, datatype='anat', root=data_path, 
                                               folderlevel='derivatives', workingdir=cfg['prep_foldername'])
        
                    if not os.path.exists(bids_strc_reg.get_path(f'{anat_format}2atlas.nii.gz')):

                        create_directory(bids_strc_reg.get_path())
                        # Copy ref file
                        shutil.copyfile(bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'),bids_strc_reg.get_path(f'ref_anat.nii.gz'))

                        if 'anat_space_organoids' not in dossier:
                           # Register anat --> template
                           # if os.path.exists(bids_strc_anat.get_path('lesion_inv_mask.nii.gz')):
                           #           lesion_mask = bids_strc_anat.get_path('lesion_inv_mask.nii.gz')
                           # else:
                           #           lesion_mask=None
                           antsreg_full(template, # fixed
                                   bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'),  # moving
                                   bids_strc_reg.get_path(f'{anat_format}2atlas'))
                     
                           # Apply inverse transform to put template in anat
                           ants_apply_transforms([template],  # input 
                                           bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'), # reference
                                           [bids_strc_reg.get_path(f'template_in_{anat_format}.nii.gz')], # output
                                           [bids_strc_reg.get_path(f'{anat_format}2atlas0GenericAffine.mat'), 1], # transform 1
                                           bids_strc_reg.get_path(f'{anat_format}2atlas1InverseWarp.nii.gz'))   # transform 2
            
                           # Apply inverse transform to put atlas in anat, make sure if it's a label atlas that the labels are still integers
                           if 'TPM' in atlas:
                               ants_apply_transforms([atlas.replace('.nii', '_vol_0000.nii.gz'),atlas.replace('.nii', '_vol_0001.nii.gz'),atlas.replace('.nii', '_vol_0002.nii.gz')],  # input 
                                               bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'), # reference
                                              [bids_strc_reg.get_path(f'atlas_in_{anat_format}_GM.nii.gz'),bids_strc_reg.get_path(f'atlas_in_{anat_format}_WM.nii.gz'),bids_strc_reg.get_path(f'atlas_in_{anat_format}_CSF.nii.gz')], # output
                                              [bids_strc_reg.get_path(f'{anat_format}2atlas0GenericAffine.mat'), 1], # transform 1
                                              bids_strc_reg.get_path(f'{anat_format}2atlas1InverseWarp.nii.gz'))  # transform 2  
                           
                           else:
                               ants_apply_transforms([atlas],  # input 
                                              bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'), # reference
                                             [bids_strc_reg.get_path(f'atlas_in_{anat_format}.nii.gz')], # output
                                             [bids_strc_reg.get_path(f'{anat_format}2atlas0GenericAffine.mat'), 1], # transform 1
                                             bids_strc_reg.get_path(f'{anat_format}2atlas1InverseWarp.nii.gz'),  # transform 2
                                             '--interpolation NearestNeighbor -u int')  
                               # remove lesion mask
                               if os.path.exists(bids_strc_anat.get_path('lesion_inv_mask.nii.gz')):
                                   fsl_mult(bids_strc_reg.get_path(f'atlas_in_{anat_format}.nii.gz'),
                                            bids_strc_anat.get_path('lesion_inv_mask.nii.gz'),
                                            bids_strc_reg.get_path(f'atlas_in_{anat_format}.nii.gz'))
                                   
               
                        else: # if the "atlas" was derived on the anatomical space of each anat_space_organoids, just copy those files
                             shutil.copyfile(template,bids_strc_reg.get_path(f'template_in_{anat_format}.nii.gz'))
                             shutil.copyfile(atlas,bids_strc_reg.get_path(f'atlas_in_{anat_format}.nii.gz'))

                    ########################## 3. REGISTRATION (ATLAS TO ANAT) TO DWI ##########################
        
                    # Define dwi data to be used for registration
                    filtered_data = subj_data[(subj_data['acqType'] == 'PGSE') & (subj_data['phaseDir'] == 'fwd') & (subj_data['blockNo'] == sess) & (subj_data['noBval'] > 1)]
                    Delta_list = [f'Delta_{int(x)}_fwd' for x in filtered_data["diffTime"].dropna()]
    
                    # Loop through the different dwi data
                    for data_type in ['allDelta-allb']: #  + Delta_list
                      
                        bids_strc_prep = create_bids_structure(subj=subj, sess=sess, datatype='dwi', description=data_type, root=data_path, 
                                                    folderlevel='derivatives', workingdir=cfg['prep_foldername'])
                        bids_strc_reg_dwi  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=dossier+'-To-'+data_type, root=data_path, 
                                                    folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                        bids_strc_reg_dwi.set_param(base_name='')
                       
                        if not os.path.exists(bids_strc_reg.get_path('template_in_dwi.nii.gz')):
                        
                            create_directory(bids_strc_reg_dwi.get_path())
                            # Copy ref file
                            shutil.copyfile(bids_strc_prep.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'),bids_strc_reg_dwi.get_path(f'ref_dwi.nii.gz'))

                            if cfg['subject_type']=='organoid' :
                                
                                # pad image temporarily for registration
                                pad_image(bids_strc_reg.get_path(f'template_in_{anat_format}.nii.gz'), bids_strc_reg.get_path(f'template_in_{anat_format}.nii.gz'))
                                pad_image(bids_strc_prep.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'),  bids_strc_prep.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'))
                                pad_image(bids_strc_reg.get_path(f'atlas_in_{anat_format}.nii.gz'), bids_strc_reg.get_path(f'atlas_in_{anat_format}.nii.gz'))

                                # apply full transform like in preprocessing because it works better than just affine
                                ants_apply_transforms([bids_strc_reg.get_path(f'template_in_{anat_format}.nii.gz')],  # input 
                                                      bids_strc_prep.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'), # moving
                                                      [bids_strc_reg_dwi.get_path('template_in_dwi.nii.gz')], # output
                                                      [bids_strc_prep.get_path(f'dwiafterpreproc2{anat_format}0GenericAffine.mat'), 1], # transform 1
                                                      bids_strc_prep.get_path(f'dwiafterpreproc2{anat_format}1InverseWarp.nii.gz'))   # transform 2
                                if 'TPM' not in atlas: # ensure atlas is binary, assumes no TPM
                                    ants_apply_transforms([bids_strc_reg.get_path(f'atlas_in_{anat_format}.nii.gz')],  # input 
                                                          bids_strc_prep.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'), # moving
                                                          [bids_strc_reg_dwi.get_path('atlas_in_dwi.nii.gz')], # output
                                                          [bids_strc_prep.get_path(f'dwiafterpreproc2{anat_format}0GenericAffine.mat'), 1], # transform 1
                                                          bids_strc_prep.get_path(f'dwiafterpreproc2{anat_format}1InverseWarp.nii.gz'),  # transform 2
                                                          '--interpolation NearestNeighbor -u int')  
                               
                                # unpad the images previousy padded
                                unpad_image(bids_strc_reg.get_path(f'template_in_{anat_format}.nii.gz'), bids_strc_reg.get_path(f'template_in_{anat_format}.nii.gz'))
                                unpad_image(bids_strc_prep.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'),  bids_strc_prep.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'))
                                unpad_image(bids_strc_reg.get_path(f'atlas_in_{anat_format}.nii.gz'), bids_strc_reg.get_path(f'atlas_in_{anat_format}.nii.gz'))
                                unpad_image(bids_strc_reg_dwi.get_path('atlas_in_dwi.nii.gz'), bids_strc_reg_dwi.get_path('atlas_in_dwi.nii.gz'))
                                unpad_image(bids_strc_reg_dwi.get_path('template_in_dwi.nii.gz'), bids_strc_reg_dwi.get_path('template_in_dwi.nii.gz'))
                                   
                                if 'TPM' not in atlas: # ensure atlas is binary
                                    # ensure mask is binary
                                    call = [
                                        "fslmaths",
                                        bids_strc_reg_dwi.get_path("atlas_in_dwi.nii.gz"),
                                        "-mul", "1",  # multiply by 1 → forces rounding to nearest integer label
                                        bids_strc_reg_dwi.get_path("atlas_in_dwi.nii.gz")
                                    ]
                                    os.system(' '.join(call))
                                
                            else:
         
                                # Apply inverse transform to put template anat in dwi
                                ants_apply_transforms_simple([bids_strc_reg.get_path(f'template_in_{anat_format}.nii.gz')],  # input 
                                                     bids_strc_prep.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'), # moving
                                                     [bids_strc_reg_dwi.get_path('template_in_dwi.nii.gz')], # output
                                                     [bids_strc_prep.get_path(f'dwiafterpreproc2{anat_format}0GenericAffine.mat'), 1]) # # transform 1
                           
                                # Apply inverse transform to put atlas anat in dwi
                                if 'TPM' in atlas: # ensure atlas is binary
                                    ants_apply_transforms_simple([bids_strc_reg.get_path(f'atlas_in_{anat_format}_GM.nii.gz'),bids_strc_reg.get_path(f'atlas_in_{anat_format}_WM.nii.gz'),bids_strc_reg.get_path(f'atlas_in_{anat_format}_CSF.nii.gz')],  # input 
                                                     bids_strc_prep.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'), # moving
                                                     [bids_strc_reg_dwi.get_path('atlas_TPM_GM_in_dwi.nii.gz'),bids_strc_reg_dwi.get_path('atlas_TPM_WM_in_dwi.nii.gz'),bids_strc_reg_dwi.get_path('atlas_TPM_CSF_in_dwi.nii.gz')], # output
                                                     [bids_strc_prep.get_path(f'dwiafterpreproc2{anat_format}0GenericAffine.mat'), 1]) # # transform 1
                            
                                else: # but not TPM
                                    # Apply inverse transform to put anat in dwi
                                    ants_apply_transforms_simple([bids_strc_reg.get_path(f'atlas_in_{anat_format}.nii.gz')],  # input 
                                                     bids_strc_prep.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'), # moving
                                                     [bids_strc_reg_dwi.get_path('atlas_in_dwi.nii.gz')], # output
                                                     [bids_strc_prep.get_path(f'dwiafterpreproc2{anat_format}0GenericAffine.mat'), 1],
                                                     '--interpolation NearestNeighbor -u int') # # transform 1
                              
                                
           ########################## B. REGISTRATION STE TO LTE ##########################
             
           #data_type =f"Delta_{cfg['LTEDelta_for_microFA']}" # Diffusion time of LTE we will compare the STE to
           
           # Create BIDS structures
           bids_LTE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=cfg['data_path'] , 
                         folderlevel='derivatives', workingdir=cfg['prep_foldername'],description='allDelta-allb')
           #extract_vols(find_files_with_pattern(bids_LTE,'pwd_avg_norm.nii.gz')[0], bids_LTE.get_path('b0.nii.gz'), 0, 1)
           bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                         folderlevel='derivatives', workingdir=cfg['prep_foldername'],description='STE_fwd')
           #extract_vols(find_files_with_pattern(bids_STE,'pwd_avg_norm.nii.gz')[0], bids_STE.get_path('b0.nii.gz'), 0, 1)
           bids_strc_reg_ste  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description='STE-To-LTE_allDelta-allb', root=data_path, 
                                          folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
           bids_strc_reg_ste.set_param(base_name='')
        

           # Register STE to LTE
           if os.path.exists(bids_STE.get_path('b0_bc.nii.gz')):                

               binary_op(bids_STE.get_path('b0_bc.nii.gz'),bids_STE.get_path('b0_mask.nii.gz'), '-mul', bids_STE.get_path('b0_bc_brain.nii.gz'))
    
               create_directory(bids_strc_reg_ste.get_path())
               antsreg_simple(bids_LTE.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'),  # fixed
                        bids_STE.get_path('b0_dn_gc_topup_avg_bc_brain.nii.gz'),# moving
                        bids_strc_reg_ste.get_path('STE2dwi'))
               
         
               # Apply inverse transform to put anat in dwi space
               ants_apply_transforms_simple([bids_STE.get_path('b0_dn_gc_topup_avg_bc_brain.nii.gz')],  # input
                                   bids_LTE.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'),# reference
                                   [bids_strc_reg_ste.get_path('STE_in_LTE_b0_brain.nii.gz')],  # output
                                   [bids_strc_reg_ste.get_path('STE2dwi0GenericAffine.mat'), 0])  # transform 1

               ants_apply_transforms_simple_4D([bids_STE.get_path('dwi_dn_gc_topup.nii.gz')],  # input
                                   bids_LTE.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'),# reference
                                   [bids_strc_reg_ste.get_path('STE_in_LTE_dn_gc_topup.nii.gz')],  # output
                                   [bids_strc_reg_ste.get_path('STE2dwi0GenericAffine.mat'), 0])  # transform 1

