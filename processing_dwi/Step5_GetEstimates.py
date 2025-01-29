"""
Script to retreive model estimates within regions of interest.
It does not use a particular python environment.

* not finished yet * 
Last changed Jan 2025
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

plt.close('all')

def Step5_GetEstimates(subj_list, cfg):
    
    
    data_path   = cfg['data_path']     
    scan_list   = pd.read_excel(os.path.join(data_path, 'ScanList.xlsx'))
    atlas_path  = cfg['common_folder']
                 
    ########################## SUBJECT-WISE OPERATIONS ##########################
    for subj in subj_list:
        
        print('Getting model estimates ' + subj + '...')
    
        # Extract data for subject
        subj_data    = scan_list[scan_list['newstudyName'] == subj].reset_index(drop=True)
        
        # List of acquisition sessions
        sess_list    = [x for x in list(subj_data['blockNo'].unique()) if not math.isnan(x)] # clean NaNs
        
        ###### SESSION-WISE OPERATIONS 
        for sess in list(subj_data['blockNo'].unique()):

           for model in cfg['model_list_GM'] + cfg['model_list_WM']:
               
                # REGISTRATION atlas
                if model in cfg['model_list_GM']:
                    data_used = 'allDelta-allb'
                    
                elif model in cfg['model_list_WM']:
                    filtered_data = subj_data[(subj_data['phaseDir'] == 'fwd') & (subj_data['blockNo'] == sess) & (subj_data['noBval'] > 1)]
                    ind_folder =filtered_data["diffTime"].idxmax()
                    data_used = 'Delta_'+str(int(filtered_data['diffTime'][ind_folder]))+'_'+filtered_data['phaseDir'][ind_folder]
                                  
                # Define BIDS structure for the analysis data
                bids_strc_prep = create_bids_structure(subj=subj, sess=sess, datatype='dwi', description=data_used, root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                bids_strc_anat = create_bids_structure(subj=subj, sess=sess, datatype='anat', root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['prep_foldername'])
            
                if not os.path.exists(bids_strc_reg.get_path('template_in_dwi.nii.gz')):
                   
                    # define atlas and make small improvements for betetr registration
                    atlas      = glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*atlas.nii.gz'))[0]
                    template   = glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*T2s_brain.nii.gz'))[0]
                    for image in (atlas,template):
                        
                        # crop template/atlas - otherwise too much to register
                        img = nib.load(image)
                        cropped_img = img.slicer[:, 270:800, :]
                        nib.save(cropped_img, image.replace('.nii.gz', '_crop.nii.gz'))
                        
                        # downsample template/atlas 
                        input_img = nibabel.load(image.replace('.nii.gz', '_crop.nii.gz'))
                        resampled_img = nibabel.processing.resample_to_output(input_img, [0.07, 0.07, 0.07])
                        nibabel.save(resampled_img,  image.replace('.nii.gz', '_crop_lowres.nii.gz')) 
                        
                    # register T2w --> template
                    create_directory(bids_strc_reg.get_path())
                    antsreg(template.replace('.nii.gz', '_crop_lowres.nii.gz'), # fixed
                            bids_strc_anat.get_path('T2w_brain.nii.gz'),  # moving
                            bids_strc_reg.get_path('T2w2atlas'))
              
                    # apply inverse transform to put template in T2w
                    ants_apply_transforms([template.replace('.nii.gz', '_crop_lowres.nii.gz'),atlas.replace('.nii.gz', '_crop_lowres.nii.gz')],  # input 
                                    bids_strc_anat.get_path('T2w.nii.gz'), # moving
                                    [bids_strc_reg.get_path('template_in_T2w.nii.gz'),bids_strc_reg.get_path('atlas_in_T2w.nii.gz')], # output
                                    [bids_strc_reg.get_path('T2w2atlas0GenericAffine.mat'), 1], # transform 1
                                    bids_strc_reg.get_path('T2w2atlas1InverseWarp.nii.gz'))   # transform 2

                    # apply inverse transform to put T2w in dwi
                    ants_apply_transforms([bids_strc_reg.get_path('template_in_T2w.nii.gz'),bids_strc_reg.get_path('atlas_in_T2w.nii.gz')],  # input 
                                     bids_strc_prep.get_path('b0_dn_gc_ec_avg_bc.nii.gz'), # moving
                                     [bids_strc_reg.get_path('template_in_dwi.nii.gz'),bids_strc_reg.get_path('atlas_in_dwi.nii.gz')], # output
                                     [bids_strc_prep.get_path('dwi2T2w0GenericAffine.mat'), 1], # transform 1
                                     bids_strc_prep.get_path('dwi2T2w1InverseWarp.nii.gz'))   # transform 2

                    # get b0 image, bias correct, brain skull extract
                    # b0       = find_files_with_pattern(bids_strc_analysis,'b0_dn_gc_ec_avg.nii.gz')[0]
                    # mask       = find_files_with_pattern(bids_strc_analysis,'mask_dil.nii.gz')[0]
                    # b0_bc    = b0.replace('.nii.gz','_bc.nii.gz')
                    # N4_unbias(b0,b0_bc)
                    # multiply_by_mask(b0_bc,bids_strc_analysis.get_path(),mask)
                    # dwi = find_files_with_pattern(bids_strc_analysis,'b0_dn_gc_ec_avg_bc_masked.nii.gz')[0]
                    
                    # # find T2 weighted
                    # T2w       = find_files_with_pattern(bids_strc_analysis,'T2w_brain.nii.gz')[0]

                    # output_path = os.path.dirname(b0_bc)
                    # antsreg_atlas(template_cropped_rs,T2w, os.path.join(output_path,'T2w2template'))

               
                # atlas = nib.load(os.path.join(output_path,'atlas_in_dwi.nii.gz')).get_fdata()

                # atlas_labels = pd.read_csv(
                #     glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*atlas.label'))[0],
                #     delim_whitespace=True,  # Use whitespace as the delimiter
                #     skiprows=14,  # Skip the header lines (modify as per actual file structure)
                #     header=None,  # No column headers in the file
                #     names=['IDX', 'R', 'G', 'B', 'A', 'VIS', 'MSH', 'LABEL'],  # Assign column names
                #     quotechar='"',  # Handle quoted strings for the LABEL column
                # )
   
                # SUMMARY PLOT
                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path, 
                             folderlevel='derivatives', workingdir=cfg['analysis_foldername'], description=model)
    
                # Define output and parameters names and values
                output_path = os.path.join(bids_strc_analysis.get_path(),'Masked')
                patterns, lims = get_param_names_model(model)

                # Start figure
                if len(patterns)<=4:
                    fig, axs = plt.subplots(1, len(patterns), figsize=(12, 4))  
                else:
                    fig, axs = plt.subplots(2, math.ceil(len(patterns)/2), figsize=(12, 4))  

                axs = axs.flatten()
                fig.subplots_adjust(wspace=0.05,hspace=0.02, top=0.95, bottom=0.05, left=0.05, right=0.95)  
                k = 0
                for pattern, lim in zip(patterns, lims):
                    
                    # Load the first file matching the pattern
                    parameter_data = nib.load(glob.glob(os.path.join(output_path, pattern))[0]).get_fdata()
                
                    # Extract the relevant slice and rotate it
                    img = parameter_data[:, 10, :]
                    fact = int((max(img.shape) - min(img.shape)) / 2)
                    img = imutils.rotate(img, angle=90)[fact:-fact, :]
                    img[np.isnan(img)]=0
                                 
                    # Display the image
                    im = axs[k].imshow(img, cmap='jet',vmin=lim[0], vmax=lim[1])
                    axs[k].axis('off') 
                    axs[k].set_title(pattern[1:-1])
                    cbar = plt.colorbar(im, ax=axs[k], orientation='horizontal', pad=0.02)  
                    labels = [lim[0],lim[1]]
                    cbar.set_ticks(labels) 
                    cbar.ax.set_xticklabels(labels, rotation=-45) 
                    cbar.ax.tick_params(labelsize=10)
                    k += 1
                    
                for ax in axs[len(patterns):]:
                    ax.set_visible(False)
                plt.tight_layout(rect=[0, 0, 1, 1])
                plt.show()
                plt.savefig(os.path.join(bids_strc_analysis.get_path(),'output_summary.png'))

                ## EXTRACT ESTIMATES
                # outfile = os.path.join(output_path,'output_ROIs.txt')
                # Data = np.zeros((len(cfg['ROIs']), len(patterns)))
                # ROI_ctr=0; pattern_ctr=0
                # for ROI in cfg['ROIs']:
                    
                #     ind_list=atlas_labels["LABEL"].str.find(ROI)
                #     match_idx = ind_list[ind_list != -1].index
                #     match_idx= atlas_labels["IDX"][match_idx].to_numpy()
                #     mask_indexes = (atlas==match_idx)
                    
                #     # Loop through each specified pattern and process the first matching file
                #     for pattern in patterns:
                        
                #         parameter_data = nib.load(glob.glob(os.path.join(bids_strc_analysis.get_path(),'Masked', pattern))[0]).get_fdata()  
                #         param_masked = parameter_data * mask_indexes

                #         indices = np.nonzero(param_masked)  
                #         data = param_masked[indices]
                #         data = data[~np.isnan(data)]
                #         Data[ROI_ctr,pattern_ctr] = np.mean(data)
                #         pattern_ctr=pattern_ctr+1
                # ROI_ctr=ROI_ctr+1
 
                        
                       