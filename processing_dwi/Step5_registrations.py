#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to do registrations

* not finished yet * 
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

def Step5_registrations(subj_list, cfg):
    
    
    data_path   = cfg['data_path']     
    scan_list   = pd.read_excel(os.path.join(data_path, 'ScanList.xlsx'))
    atlas_path  = cfg['common_folder']
                 
    ########################## SUBJECT-WISE OPERATIONS ##########################
    for subj in subj_list:
        
        print('Doing registration on ' + subj + '...')
    
        # Extract data for subject
        subj_data    = scan_list[scan_list['newstudyName'] == subj].reset_index(drop=True)
        
        # List of acquisition sessions
        sess_list    = [x for x in list(subj_data['blockNo'].unique()) if not math.isnan(x)] # clean NaNs
        
        ###### SESSION-WISE OPERATIONS 
        for sess in sess_list:

          
           ### REGISTRATION ATLAS TO T2W ###
          
           # Create BIDS structure
           bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas'] +'-To-T2w', root=data_path, 
                                       folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
           bids_strc_reg.set_param(base_name='')
           bids_strc_anat = create_bids_structure(subj=subj, sess=sess, datatype='anat', root=data_path, 
                                       folderlevel='derivatives', workingdir=cfg['prep_foldername'])
       
           # Define atlas and make small improvements for betetr registration
           atlas      = glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*atlas.nii.gz'))[0]
           template   = glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*T2s_brain.nii.gz'))[0]
           
           if not os.path.exists(template.replace('.nii.gz', '_crop_lowres.nii.gz')):
               # remove regions to make it more like brain
               img  = nib.load(atlas)
               data_a = img.get_fdata()
               data_a[data_a == 42] = 0
               data_a[data_a == 41] = 0
               data_a[data_a == 45] = 0
               data_a[data_a == 76] = 0
               img  = nib.load(template)
               data_t = img.get_fdata()           
               mask = (data_a != 0).astype(np.uint8)
               data_t = data_t*mask
               nib.save(nib.Nifti1Image(data_a, img.affine), atlas.replace('.nii.gz', '_crop.nii.gz'))
               nib.save(nib.Nifti1Image(data_t, img.affine), template.replace('.nii.gz', '_crop.nii.gz'))
    
               for image in (atlas,template):
                    
                    # Crop template/atlas - otherwise too much to register
                    img  = nib.load(image.replace('.nii.gz', '_crop.nii.gz'))
                    data = img.get_fdata()
                    masked_data = np.zeros_like(data)
                    masked_data[:, 270:840, :] = data[:, 270:840, :]
                    nib.save(nib.Nifti1Image(masked_data, img.affine), image.replace('.nii.gz', '_crop.nii.gz'))
                    
    
                    # Downsample template/atlas 
                    input_img = nibabel.load(image.replace('.nii.gz', '_crop.nii.gz'))
                    if image==atlas:
                        resampled_img = nibabel.processing.resample_to_output(input_img, [0.05, 0.5, 0.05],order=0  )
                    elif image==template:
                        resampled_img = nibabel.processing.resample_to_output(input_img, [0.05, 0.5, 0.05])
                    nibabel.save(resampled_img,  image.replace('.nii.gz', '_crop_lowres.nii.gz')) 
                
           # Register T2w --> template
           if not os.path.exists(bids_strc_reg.get_path('T2w2atlas.nii.gz')):
               create_directory(bids_strc_reg.get_path())
               antsreg_full(template.replace('.nii.gz', '_crop_lowres.nii.gz'), # fixed
                       bids_strc_anat.get_path('T2w_bc_brain.nii.gz'),  # moving
                       bids_strc_reg.get_path('T2w2atlas'))
         
               # Apply inverse transform to put template and atlas in T2w
               ants_apply_transforms([template.replace('.nii.gz', '_crop_lowres.nii.gz')],  # input 
                               bids_strc_anat.get_path('T2w_bc_brain.nii.gz'), # reference
                               [bids_strc_reg.get_path('template_in_T2w.nii.gz')], # output
                               [bids_strc_reg.get_path('T2w2atlas0GenericAffine.mat'), 1], # transform 1
                               bids_strc_reg.get_path('T2w2atlas1InverseWarp.nii.gz'))   # transform 2

               ants_apply_transforms([atlas.replace('.nii.gz', '_crop_lowres.nii.gz')],  # input 
                               bids_strc_anat.get_path('T2w_bc_brain.nii.gz'), # reference
                              [bids_strc_reg.get_path('atlas_in_T2w.nii.gz')], # output
                              [bids_strc_reg.get_path('T2w2atlas0GenericAffine.mat'), 1], # transform 1
                              bids_strc_reg.get_path('T2w2atlas1InverseWarp.nii.gz'),  # transform 2
                              '--interpolation NearestNeighbor -u int')  
               
           ### REGISTRATION ATLAS TO T2W TO DWI ###
           
           # Define dwi data to be used for registration
           filtered_data = subj_data[(subj_data['acqType'] == 'PGSE') & (subj_data['phaseDir'] == 'fwd') & (subj_data['blockNo'] == sess) & (subj_data['noBval'] > 1)]
           ind_folder = getattr(filtered_data["diffTime"], 'idxmin')()
           data_Deltamin = 'Delta_'+str(int(filtered_data['diffTime'][ind_folder]))+'_fwd'  
                
           filtered_data = subj_data[(subj_data['acqType'] == 'PGSE') & (subj_data['phaseDir'] == 'fwd') & (subj_data['blockNo'] == sess) & (subj_data['noBval'] > 1)]
           ind_folder = getattr(filtered_data["diffTime"], 'idxmax')()
           data_Deltamax = 'Delta_'+str(int(filtered_data['diffTime'][ind_folder]))+'_'+filtered_data['phaseDir'][ind_folder]  

           # Loop through the different dwi data
           for data_type in ['allDelta-allb', data_Deltamin, data_Deltamax]:
              
                bids_strc_prep = create_bids_structure(subj=subj, sess=sess, datatype='dwi', description=data_type, root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['prep_foldername'])
                bids_strc_reg_dwi  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas']+'_To_'+data_type, root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                bids_strc_reg_dwi.set_param(base_name='')

                if not os.path.exists(bids_strc_reg.get_path('template_in_dwi.nii.gz')):
                
                    create_directory(bids_strc_reg_dwi.get_path())
 
                    # Apply inverse transform to put T2w in dwi
                    ants_apply_transforms_simple([bids_strc_reg.get_path('template_in_T2w.nii.gz'),bids_strc_reg.get_path('atlas_in_T2w.nii.gz')],  # input 
                                     bids_strc_prep.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'), # moving
                                     [bids_strc_reg_dwi.get_path('template_in_dwi.nii.gz'),bids_strc_reg_dwi.get_path('atlas_in_dwi.nii.gz')], # output
                                     [bids_strc_prep.get_path('dwiafterpreproc2T2w0GenericAffine.mat'), 1],
                                     '--interpolation NearestNeighbor -u int') # # transform 1
                                    # bids_strc_prep.get_path('dwi2T2w1InverseWarp.nii.gz'),'NearestNeighbor')# transform 2
                
           ### REGISTRATION STE to LTE 
             
           data_type ='Delta_38_fwd' # Diffusion time of LTE we will compare the STE to
           
           # Create BIDS structures
           bids_LTE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=cfg['data_path'] , 
                         folderlevel='derivatives', workingdir=cfg['prep_foldername'],description=data_type)
           #extract_vols(find_files_with_pattern(bids_LTE,'pwd_avg_norm.nii.gz')[0], bids_LTE.get_path('b0.nii.gz'), 0, 1)
           bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                         folderlevel='derivatives', workingdir=cfg['prep_foldername'],description='STE_fwd')
           #extract_vols(find_files_with_pattern(bids_STE,'pwd_avg_norm.nii.gz')[0], bids_STE.get_path('b0.nii.gz'), 0, 1)
           bids_strc_reg_ste  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description='STE_To_LTE_'+data_type, root=data_path, 
                                          folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
           bids_strc_reg_ste.set_param(base_name='')
        

           # Register STE to LTE
           if os.path.exists(bids_STE.get_path('b0_bc.nii.gz')):                

               binary_op(bids_STE.get_path('b0_bc.nii.gz'),bids_STE.get_path('b0_mask.nii.gz'), '-mul', bids_STE.get_path('b0_bc_brain.nii.gz'))
    
               create_directory(bids_strc_reg_ste.get_path())
               antsreg_simple(bids_LTE.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'),  # fixed
                        bids_STE.get_path('b0_dn_gc_topup_avg_bc_brain.nii.gz'),# moving
                        bids_strc_reg_ste.get_path('STE2dwi'))
               
         
               # Apply inverse transform to put T2w in dwi space
               ants_apply_transforms_simple([bids_STE.get_path('b0_dn_gc_topup_avg_bc_brain.nii.gz')],  # input
                                   bids_LTE.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'),# reference
                                   [bids_strc_reg_ste.get_path('STE_in_LTE.nii.gz')],  # output
                                   [bids_strc_reg_ste.get_path('STE2dwi0GenericAffine.mat'), 0])  # transform 1

             
