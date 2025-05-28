"""
Script to retrieve model estimates within regions of interest.
Compatible with any Python environment.

Last updated: Jan 2025
@author: Rita O
"""

import os
import sys
import glob
import math
import copy
import csv
import platform
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import nibabel as nib
import imutils
import numpy.ma as ma
import SimpleITK as sitk

from custom_functions import *
from bids_structure import *

plt.close('all')

def Step5_get_estimates(subj_list, cfg):
    data_path = cfg['data_path']
    scan_list = pd.read_excel(os.path.join(data_path, 'ScanList.xlsx'))
    cfg['model_list'] = cfg['model_list_GM'] + cfg['model_list_WM']

    ######## SUBJECT-WISE OPERATIONS ########
    for subj in subj_list:
        print(f'Getting model estimates for {subj}...')
        
        subj_data = scan_list[scan_list['newstudyName'] == subj].reset_index(drop=True)
        sess_list = [x for x in subj_data['blockNo'].unique() if not math.isnan(x)]

        ######## SESSION-WISE OPERATIONS ########
        for sess in sess_list:
            print(f'Session {sess}...')
            
            ######## MODEL-WISE OPERATIONS ########
            for model in cfg['model_list']:
                print(f'Getting model estimates from {model}...')

                filtered_data = subj_data[(subj_data['acqType'] == 'PGSE') &
                                          (subj_data['phaseDir'] == 'fwd') &
                                          (subj_data['blockNo'] == sess) &
                                          (subj_data['noBval'] > 1)]
                
                # Define BIDS structure for the analysis data depending on the input
                if model == 'Nexi':
                    data_used = 'allDelta-allb'
                elif model == 'Sandi':
                    idx = filtered_data['diffTime'].idxmin()
                    data_used = f"Delta_{int(filtered_data['diffTime'][idx])}_fwd"
                elif model in ('SMI', 'SMI_wSTE'):
                    idx = filtered_data['diffTime'].idxmax()
                    data_used = f"Delta_{int(filtered_data['diffTime'][idx])}_{filtered_data['phaseDir'][idx]}"

                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path, 
                             folderlevel='derivatives', workingdir=cfg['analysis_foldername'], description=model)
                output_path = os.path.join(bids_strc_analysis.get_path(), 'Masked')

                bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas']+'_To_'+data_used, root=data_path, 
                                             folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                bids_strc_reg.set_param(base_name='')
                
                # Define atlas
                atlas = bids_strc_reg.get_path('atlas_in_dwi.nii.gz')

                ######## OPERATIONS INVOLVING THE NEED OF AN ATLAS  ########
                if os.path.exists(atlas):
                    
                    # Define atlas labels 
                    atlas_labels = prepare_atlas_labels(cfg['atlas'], cfg['common_folder'])

                    # Define TPMs
                    bids_strc_reg_TPM  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas_TPM']+'_To_'+data_used, root=cfg['data_path'] , 
                                                 folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                    bids_strc_reg_TPM.set_param(base_name='')
                    TPMs = []
                    for tissue in ['GM', 'WM', 'CSF']:
                        path = bids_strc_reg_TPM.get_path(f'atlas_TPM_{tissue}_in_dwi.nii.gz')
                        TPMs.append(path if os.path.exists(path) else '')

                    # Determine ROI list and output file
                    if model in cfg['model_list_GM']:
                        ROI_list = cfg['ROIs_GM']
                        outfile = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_GM_{model}.xlsx")
                    else:
                        ROI_list = cfg['ROIs_WM']
                        outfile = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_WM_{model}.xlsx")

                    ######## EXTRACT MODEL ESTIMATES ########
                    # Extract estimates
                    patterns, lims, maximums = get_param_names_model(model)
                    Data = np.zeros((len(ROI_list), len(patterns)))
                    for i, ROI in enumerate(ROI_list):
                        mask = create_ROI_mask(atlas, atlas_labels, TPMs, ROI, cfg['tpm_thr'], bids_strc_reg)
                        for j, (pattern, maximum) in enumerate(zip(patterns, maximums)):
                            param_img = nib.load(glob.glob(os.path.join(output_path, pattern))[0]).get_fdata()
                            masked = param_img[mask > 0]  # Select only voxels inside the ROI
                            masked_clean = masked[~np.isnan(masked) & (masked > maximum[0]) & (masked < maximum[1])]
                            Data[i, j] = np.nanmean(masked_clean) if len(masked_clean) > 0 else np.nan
    
                    # Create table structure with results
                    df_data = pd.DataFrame(Data, columns=patterns)
                    df_data.insert(0, 'ROI Name', ROI_list)
                    
                    # Extract estimates from all ROIs (dont save nifiti files of mask)
                    Data2 = np.zeros((len(atlas_labels['IDX'].to_numpy()), len(patterns)))
                    for i, indx in enumerate(atlas_labels['IDX'].to_numpy()):
                        mask = create_ROI_mask_fromindx(atlas, atlas_labels, indx, bids_strc_reg)

                        for j, (pattern, maximum) in enumerate(zip(patterns, maximums)):
                            param_img = nib.load(glob.glob(os.path.join(output_path, pattern))[0]).get_fdata()
                            masked = param_img[mask > 0]  # Select only voxels inside the ROI
                            #masked_clean = masked[~np.isnan(masked) & (masked > maximum[0]) & (masked < maximum[1])]
                            masked_clean = masked[~np.isnan(masked)]
                            Data2[i, j] = np.median(masked_clean) if len(masked_clean) > 0 else np.nan

                    # Add mean of medians
                    df_data.loc[len(df_data)] = ['mean of medians'] + np.nanmean(Data2, axis=0).tolist()
                    
                    # Save in excel
                    df_data.to_excel(outfile, index=False)


                    ######## PLOT SNR ########
                    if model == 'Nexi':
                        
                        # Load data
                        bvals = read_numeric_txt(os.path.join(bids_strc_analysis.get_path(),'powderaverage.bval'))
                        S_S0  = nib.load(os.path.join(bids_strc_analysis.get_path(),'powderaverage_dwi.nii.gz')).get_fdata()
                        nf = nib.load(os.path.join(bids_strc_analysis.get_path(),'normalized_sigma.nii.gz')).get_fdata()*np.sqrt(np.pi/2)
     
                        # Loop through ROIs                    
                        fig, axs = plt.subplots(1, len(ROI_list), figsize=(12, 4))  
                        k=0
                        for ROI in ROI_list:
         
                            mask_indexes = create_ROI_mask(atlas, atlas_labels, TPMs, ROI, cfg['tpm_thr'], bids_strc_reg)
    
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
                        plt.savefig(os.path.join(os.path.dirname(os.path.dirname(output_path)), 'SignalDecay_summary.png'))
                        plt.close(fig)
                    
                ######## PLOT SUMMARY MAP PLOT ########
     
                # Make colorbar
                import matplotlib.colors as mcolors
                import matplotlib.cm as cm
                jet = cm.get_cmap('jet', 256)
                jet_colors = jet(np.linspace(0, 1, 256))
                fade_len = 20
                fade = np.linspace(0, 1, fade_len).reshape(-1, 1)
                jet_colors[:fade_len, :3] *= fade  # Keep alpha (4th channel) unchanged
                jet_colors[-fade_len:, :3] *= fade[::-1]  # Reverse fade for the end
                custom_jet_black = mcolors.ListedColormap(jet_colors)
 
                # Define BIDS structure and output path
                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype='dwi',  root=data_path,
                                            folderlevel='derivatives', workingdir=cfg['analysis_foldername'], description=model )
                output_path = os.path.join(bids_strc_analysis.get_path(), 'Masked')
            
                # Get model parameter names and display ranges
                patterns, lims, maximums = get_param_names_model(model)
                
                # Create subplot grid
                n_params = len(patterns)
                n_rows = 1 if n_params <= 4 else 2
                n_cols = math.ceil(n_params / n_rows)
                fig, axs = plt.subplots(n_rows, n_cols, figsize=(12, 4))
                axs = axs.flatten()
                fig.subplots_adjust(wspace=0.05, hspace=0.02, top=0.95, bottom=0.05, left=0.05, right=0.95)
                
                # Display each parameter slice
                for ax, pattern, lim in zip(axs, patterns, lims):
                    # Load parameter image
                    param_path = glob.glob(os.path.join(output_path, pattern))[0]
                    param_data = nib.load(param_path).get_fdata()
                
                    # Extract and process middle slice
                    if cfg['subject_type'] =='rat':
                        slicee = int(np.ceil(nib.load(param_path).shape[1]/2))
                        img = imutils.rotate(param_data[:,slicee, :], angle=90)
                        fact = int((max(img.shape) - min(img.shape)) / 2)
                        img = img[fact:-fact, :]
                        img[np.isnan(img)] = 0
                    elif cfg['subject_type'] =='human':
                        slicee = int(np.ceil(nib.load(param_path).shape[2]/2))
                        img = imutils.rotate(param_data[:,:, slicee], angle=90)
                        img[np.isnan(img)] = 0
                
                    # Show slice
                    im = ax.imshow(img, cmap=custom_jet_black, vmin=lim[0], vmax=lim[1])
                    ax.set_title(pattern[1:-1])
                    ax.axis('off')
                
                    # Add colorbar
                    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.02)
                    cbar.set_ticks([lim[0], lim[1]])
                    cbar.ax.set_xticklabels([lim[0], lim[1]], rotation=-45)
                    cbar.ax.tick_params(labelsize=10)
                
                # Hide unused axes
                for ax in axs[n_params:]:
                    ax.set_visible(False)
                
                # Save and close figure
                plt.tight_layout(rect=[0, 0, 1, 1])
                plt.savefig(os.path.join(bids_strc_analysis.get_path(), 'output_summary.png'))
                plt.close(fig)
                
            ######## EXTRACT DTI,DKI ESTIMATES ########
            ######## OPERATIONS INVOLVING THE NEED OF AN ATLAS  ########

            # Shortened and improved DTI/DKI Delta loop
            ROI_list = cfg['ROIs_GM'] + cfg['ROIs_WM']
            
            filtered_data = subj_data[
                (subj_data['phaseDir'] == 'fwd') &
                (subj_data['blockNo'] == sess) &
                (subj_data['noBval'] > 1) &
                (subj_data['acqType'] == 'PGSE') &
                (subj_data['scanQA'] == 'ok')
            ]
            Delta_list = sorted(filtered_data['diffTime'].dropna().astype(int).unique())
            
            patterns, lims, maximums = get_param_names_model('DTI_DKI_short')
            Data_DTIDKI = np.zeros((len(Delta_list), len(ROI_list), len(patterns)))
            
            for d_idx, Delta in enumerate(Delta_list):
                data_used = f'Delta_{Delta}'
                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path,
                    folderlevel='derivatives', workingdir=cfg['analysis_foldername'], description=f'DTI_DKI_{data_used}')
                output_path = os.path.join(bids_strc_analysis.get_path(), 'Masked')
            
                bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas']+'_To_'+data_used+'_fwd', root=data_path, 
                                         folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                bids_strc_reg.set_param(base_name='')
            
                # Define atlas  
                atlas = bids_strc_reg.get_path('atlas_in_dwi.nii.gz')
            
                # Define atlas labels 
                atlas_labels = prepare_atlas_labels(cfg['atlas'], cfg['common_folder'])
     
                # Define TPMs
                bids_strc_reg_TPM  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas_TPM']+'_To_'+data_used+'_fwd', root=cfg['data_path'] , 
                                         folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                bids_strc_reg_TPM.set_param(base_name='')
                TPMs = []
                for tissue in ['GM', 'WM', 'CSF']:
                    path = bids_strc_reg_TPM.get_path(f'atlas_TPM_{tissue}_in_dwi.nii.gz')
                    TPMs.append(path if os.path.exists(path) else '')
            
                # Get data
                if os.path.exists(atlas) and os.path.exists(output_path):
                    for r_idx, ROI in enumerate(ROI_list):
                        mask = create_ROI_mask(atlas, atlas_labels, TPMs, ROI, cfg['tpm_thr'], bids_strc_reg)
                        for p_idx, (pattern, maximum) in enumerate(zip(patterns, maximums)):
                            param_data = nib.load(glob.glob(os.path.join(output_path, pattern))[0]).get_fdata()
                            masked = param_img[mask > 0]  # Select only voxels inside the ROI
                            masked_clean = masked[~np.isnan(masked) & (masked > maximum[0]) & (masked < maximum[1])]
                            Data_DTIDKI[d_idx, r_idx, p_idx] = np.nanmean(masked_clean) if len(masked_clean) > 0 else np.nan
            
            # Plot results
            if os.path.exists(atlas) and os.path.exists(output_path):
                import distinctipy
                color_list = distinctipy.get_colors(len(ROI_list), pastel_factor=0.5)
                fig, ax = plt.subplots(3, 1, figsize=(10, 6))
                fig.subplots_adjust(wspace=0.01, hspace=0.10, top=0.91, bottom=0.14, left=0.15, right=0.95)
                
                line_handles = []
                for i, ROI in enumerate(ROI_list):
                    line, = ax[0].plot(Delta_list, Data_DTIDKI[:, i, 0], linestyle='--', marker='o', c=color_list[i])
                    ax[1].plot(Delta_list, Data_DTIDKI[:, i, 1], linestyle='--', marker='o', c=color_list[i])
                    ax[2].plot(Delta_list, Data_DTIDKI[:, i, 2], linestyle='--', marker='o', c=color_list[i])
                    line_handles.append(line)
                
                ax[0].set_ylabel('$MD$'); ax[0].set_ylim([0, 1]); ax[0].tick_params(labelbottom=False)
                ax[1].set_ylabel('$MK$'); ax[1].set_ylim([0, 1]); ax[1].tick_params(labelbottom=False)
                ax[2].set_ylabel('$FA$'); ax[2].set_ylim([0, 1]); ax[2].set_xlabel('Diffusion time [ms]')
                
                plt.legend(handles=line_handles, labels=ROI_list, loc='upper right', fontsize=7)
                plt.suptitle('DKI')
                plt.rc('font', size=9)
                plt.savefig(os.path.join(os.path.dirname(os.path.dirname(output_path)), 'DKI.png'))
                plt.close(fig)
                
