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
import csv
import copy

plt.close('all')

def Step5_get_estimates(subj_list, cfg):
    
    
    data_path   = cfg['data_path']     
    scan_list   = pd.read_excel(os.path.join(data_path, 'ScanList.xlsx'))
    atlas_path  = cfg['common_folder']
    cfg['model_list'] = cfg['model_list_GM'] + cfg['model_list_WM']
       
    ########################## SUBJECT-WISE OPERATIONS ##########################
    for subj in subj_list:
        
        print('Getting model estimates ' + subj + '...')
    
        # Extract data for subject
        subj_data    = scan_list[scan_list['newstudyName'] == subj].reset_index(drop=True)
        
        # List of acquisition sessions
        sess_list    = [x for x in list(subj_data['blockNo'].unique()) if not math.isnan(x)] # clean NaNs
        
        ###### SESSION-WISE OPERATIONS 
        for sess in sess_list:
           
           print('Working on session ' + str(sess) + '...')

           for model in cfg['model_list']:
               
                print('Getting model estimates from ' + model + '...')

                # Define BIDS structure for the analysis data depending on the input
                if model=='Nexi':
                    data_used = 'allDelta-allb'
                elif model=='Sandi': # lowest diff time
                    filtered_data = subj_data[(subj_data['acqType'] == 'PGSE') & (subj_data['phaseDir'] == 'fwd') & (subj_data['blockNo'] == sess) & (subj_data['noBval'] > 1)]
                    ind_folder = getattr(filtered_data["diffTime"], 'idxmin')()
                    data_used = 'Delta_'+str(int(filtered_data['diffTime'][ind_folder]))+'_fwd'  
                elif model in 'SMI' or model=='SMI_wSTE': # largest diff time
                    filtered_data = subj_data[(subj_data['acqType'] == 'PGSE') & (subj_data['phaseDir'] == 'fwd') & (subj_data['blockNo'] == sess) & (subj_data['noBval'] > 1)]
                    ind_folder = getattr(filtered_data["diffTime"], 'idxmax')()
                    data_used = 'Delta_'+str(int(filtered_data['diffTime'][ind_folder]))+'_'+filtered_data['phaseDir'][ind_folder]  
                
                bids_strc_prep = create_bids_structure(subj=subj, sess=sess, datatype='dwi', description=data_used, root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['prep_foldername'])
                bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas']+'_To_'+data_used, root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                bids_strc_reg.set_param(base_name='')
                bids_strc_anat = create_bids_structure(subj=subj, sess=sess, datatype='anat', root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['prep_foldername'])
            
                # Get atlas in dwi space and atlas labels
                atlas = bids_strc_reg.get_path('atlas_in_dwi.nii.gz')
                atlas_labels = pd.read_csv(
                    glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*atlas.label'))[0],
                    sep=r'\s+',
                    skiprows=14,  
                    header=None,  
                    names=['IDX', 'R', 'G', 'B', 'A', 'VIS', 'MSH', 'LABEL'], 
                    quotechar='"',)  
                
                ### SUMMARY PLOT ###
                
                # Define bids structure, output and parameters names and values
                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path, 
                             folderlevel='derivatives', workingdir=cfg['analysis_foldername'], description=model)
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
                    
                    # Load data
                    parameter_data = nib.load(glob.glob(os.path.join(output_path, pattern))[0]).get_fdata()
                
                    # Extract the relevant slice and rotate it
                    img = parameter_data[:, 10, :]
                    fact = int((max(img.shape) - min(img.shape)) / 2)
                    img = imutils.rotate(img, angle=90)[fact:-fact, :]
                    img[np.isnan(img)]=0
                                 
                    # Display 
                    im = axs[k].imshow(img, cmap='jet',vmin=lim[0], vmax=lim[1])
                    axs[k].axis('off') 
                    axs[k].set_title(pattern[1:-1])
                    cbar = plt.colorbar(im, ax=axs[k], orientation='horizontal', pad=0.02)  
                    labels = [lim[0],lim[1]]
                    cbar.set_ticks(labels) 
                    cbar.ax.set_xticklabels(labels, rotation=-45) 
                    cbar.ax.tick_params(labelsize=10)
                    k += 1
                
                # finish figure and save
                for ax in axs[len(patterns):]:
                    ax.set_visible(False)
                plt.tight_layout(rect=[0, 0, 1, 1])
                plt.savefig(os.path.join(bids_strc_analysis.get_path(),'output_summary.png'))
                plt.show()
                plt.close(fig)
                
                ### EXTRACT MODEL ESTIMATES ###
                
                if model in cfg['model_list_GM']:
                    outfile   = os.path.join(os.path.dirname(os.path.dirname(output_path)),f'output_ROIs_GM_{model}.xlsx')
                    ROI_list  = cfg['ROIs_GM']

                elif model in cfg['model_list_WM']:
                    outfile   = os.path.join(os.path.dirname(os.path.dirname(output_path)),f'output_ROIs_WM_{model}.xlsx')
                    ROI_list  = cfg['ROIs_WM']
                
                # Loop through ROIs
                ROI_ctr = 0
                Data = np.zeros((len(ROI_list), len(patterns)))  
                
                for ROI in ROI_list:
                    
                    mask_indexes = create_ROI_mask(atlas, atlas_labels, ROI, bids_strc_reg)
                    
                    # Loop through each parameter outputed in the model
                    pattern_ctr = 0
                    for pattern in patterns:
                        
                        # Load parameter data
                        parameter_data = nib.load(glob.glob(os.path.join(bids_strc_analysis.get_path(), 'Masked', pattern))[0]).get_fdata()
                        
                        # Mask the parameter with the ROI
                        param_masked = parameter_data * 1 * mask_indexes
                        data = param_masked.reshape(param_masked.shape[0]*param_masked.shape[1]*param_masked.shape[2], 1)
                        data = data[~(np.isnan(data).any(axis=1) | (data == 0).any(axis=1))]
                        
                        # Store the mean of the masked data
                        Data[ROI_ctr, pattern_ctr] = np.nanmean(data) if len(data) > 0 else np.nan
           
                        pattern_ctr += 1
                    
                    ROI_ctr += 1
                
                # Save excel      
                df_data = pd.DataFrame(Data, columns=patterns)
                df_data.insert(0, 'ROI Name', ROI_list)
                df_data.to_excel(outfile, index=False)

            
                ### PLOT SNR FOR EACH ROI ###
                
                if model in cfg['model_list_GM']:
                    
                    # load data
                    bvals = read_numeric_txt(os.path.join(bids_strc_analysis.get_path(),'powderaverage.bval'))
                    S_S0  = nib.load(os.path.join(bids_strc_analysis.get_path(),'powderaverage_dwi.nii.gz')).get_fdata()

                    nf = nib.load(os.path.join(bids_strc_analysis.get_path(),'normalized_sigma.nii.gz')).get_fdata()*np.sqrt(np.pi/2)
 
                    # Loop through ROIs                    
                    fig, axs = plt.subplots(1, len(ROI_list), figsize=(12, 4))  
                    k=0
                    for ROI in ROI_list:
     
                        mask_indexes = create_ROI_mask(atlas, atlas_labels, ROI, bids_strc_reg)

                        S_S0_masked = copy.deepcopy(S_S0)
                        for v in range(S_S0_masked.shape[-1]):
                            S_S0_masked[:, :, :, v] = np.multiply(S_S0_masked[:, :, :, v], mask_indexes)
    
                        data = S_S0_masked.reshape(S_S0_masked.shape[0]*S_S0_masked.shape[1]*S_S0_masked.shape[2], S_S0_masked.shape[3])
                        data = data[~(np.isnan(data).any(axis=1) | (data == 0).any(axis=1))]
                                   
                        nf_masked = nf * mask_indexes
                        nf_masked = nf_masked.reshape(nf_masked.shape[0]*nf_masked.shape[1]*nf_masked.shape[2], 1)
                        nf_masked = nf_masked[~(np.isnan(nf_masked).any(axis=1) | (nf_masked == 0).any(axis=1))]

                        axs[k].plot(np.transpose(bvals), np.nanmean(data, axis=0), 'bo', markersize=3)

                        axs[k].plot(np.transpose(bvals), np.repeat(np.nanmean(nf_masked), np.transpose(bvals).shape[0]))

                        # Set axes
                        if k==0:
                            axs[k].set_ylabel(r'$S / S_0$', fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})  # Correct LaTeX formatting
                        axs[k].set_xlabel('b-val',
                                      fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})
                        axs[k].grid(True)
                        axs[k].set_title(ROI)
                        k += 1
                    plt.savefig(os.path.join(bids_strc_analysis.get_path(),'SignalDecay_summary.png'))
                    plt.close(fig)


                    
                 