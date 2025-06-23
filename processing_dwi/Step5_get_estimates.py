"""
Script to retrieve model estimates within regions of interest.
Needs registration of an atlas to work.
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
from sklearn.linear_model import LinearRegression
from custom_functions import *
from bids_structure import *

plt.close('all')

def Step5_get_estimates(subj_list, cfg):
    data_path = cfg['data_path']
    scan_list = pd.read_excel(os.path.join(data_path, cfg['scan_list_name'] ))
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
                if model == 'Nexi' or model == 'Smex':
                    data_used = 'allDelta-allb'
                elif model == 'Sandi': # lowest diff time
                    idx = filtered_data['diffTime'].idxmin()
                    data_used = f"Delta_{int(filtered_data['diffTime'][idx])}_fwd"
                elif model in ('SMI', 'SMI_wSTE'):  # largest diff time
                    idx = filtered_data['diffTime'].idxmax()
                    data_used = f"Delta_{int(filtered_data['diffTime'][idx])}_{filtered_data['phaseDir'][idx]}"
               
                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path, 
                             folderlevel='derivatives', workingdir=cfg['analysis_foldername'], description=model)
                output_path = os.path.join(bids_strc_analysis.get_path(), 'Output_masked')

                bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas']+'-To-'+data_used, root=data_path, 
                                             folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                bids_strc_reg.set_param(base_name='')
                
                # Define atlas
                atlas = bids_strc_reg.get_path('atlas_in_dwi.nii.gz')

                ######## OPERATIONS INVOLVING THE NEED OF AN ATLAS  ########
                if os.path.exists(atlas):
                    
                    # Define atlas labels 
                    if 'anat_space_organoids' not in  cfg['atlas'] :
                        atlas_labels = prepare_atlas_labels(cfg['atlas'], glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*label*'))[0])
                    else:
                        bids_strc_anat = create_bids_structure(subj=subj, sess=sess, datatype='anat', root=data_path, 
                                                   folderlevel='derivatives', workingdir=cfg['prep_foldername'])   
                        atlas_labels = prepare_atlas_labels(cfg['atlas'], glob.glob(os.path.join(bids_strc_anat.get_path(), '*label*'))[0])

                    # Define TPMs
                    if cfg['atlas_TPM']:
                        bids_strc_reg_TPM  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas_TPM']+'-'+data_used, root=cfg['data_path'] , 
                                                     folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                        bids_strc_reg_TPM.set_param(base_name='')
                        TPMs = []
                        for tissue in ['GM', 'WM', 'CSF']:
                            path = bids_strc_reg_TPM.get_path(f'atlas_TPM_{tissue}_in_dwi.nii.gz')
                            TPMs.append(path if os.path.exists(path) else '')
                    else:
                        TPMs = []

                    # Determine ROI list and output file
                    if model in cfg['model_list_GM']:
                        ROI_list = cfg['ROIs_GM']
                        outfile = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_GM_{model}.xlsx")
                    else:
                        ROI_list = cfg['ROIs_WM']
                        outfile = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_WM_{model}.xlsx")

                    ######## EXTRACT MODEL ESTIMATES ########
                    # Extract estimates
                    patterns, lims, maximums = get_param_names_model(model,cfg['subject_type'])
                    cleaned_patterns = [re.sub(r'^\*\*(.*?)\*\*$', r'*\1*', re.sub(model, '', p, flags=re.IGNORECASE)) for p in patterns] 
                    Data = np.zeros((len(ROI_list), len(patterns)))
                    for i, ROI in enumerate(ROI_list):
                        mask = create_ROI_mask(atlas, atlas_labels, TPMs, ROI, cfg['tpm_thr'], bids_strc_reg)
                        for j, (pattern, maximum) in enumerate(zip(patterns, maximums)): 
                            param_img = nib.load(glob.glob(os.path.join(output_path, pattern))[0]).get_fdata()
                            masked = param_img[mask > 0]  # Select only voxels inside the ROI
                            masked_clean = masked[~np.isnan(masked) & (masked > maximum[0]) & (masked < maximum[1])]
                            Data[i, j] = np.nanmean(masked_clean) if len(masked_clean) > 0 else np.nan
    
                    # Create table structure with results
                    df_data = pd.DataFrame(Data, columns=cleaned_patterns)
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
                        if len(ROI_list) == 1:
                                axs = [axs]  # ensure axs is always iterable
                        k=0
                        for ROI in ROI_list:
         
                            # Plot data
                            mask_indexes = create_ROI_mask(atlas, atlas_labels, TPMs, ROI, cfg['tpm_thr'], bids_strc_reg)
    
                            S_S0_masked = copy.deepcopy(S_S0)
                            for v in range(S_S0_masked.shape[-1]):
                                S_S0_masked[:, :, :, v] = np.multiply(S_S0_masked[:, :, :, v], mask_indexes)
        
                            data = S_S0_masked.reshape(S_S0_masked.shape[0]*S_S0_masked.shape[1]*S_S0_masked.shape[2], S_S0_masked.shape[3])
                            data = data[~(np.isnan(data).any(axis=1) | (data == 0).any(axis=1))]
                                       
                            nf_masked = nf * mask_indexes
                            nf_masked = nf_masked.reshape(nf_masked.shape[0]*nf_masked.shape[1]*nf_masked.shape[2], 1)
                            nf_masked = nf_masked[~(np.isnan(nf_masked).any(axis=1) | (nf_masked == 0).any(axis=1))]
    
                            # Plot data
                            axs[k].plot(np.transpose(bvals), np.nanmean(data, axis=0), 'bo', markersize=3)
    
                            # Plot noise floor
                            axs[k].plot(np.transpose(bvals), np.repeat(np.nanmean(nf_masked), np.transpose(bvals).shape[0]))
    
                            # # Make linear fit
                            # bval_mask = bvals < 4
                            # bvals_subset = bvals[bval_mask]
                            # log_signal = np.log(np.nanmean(data[:,bval_mask[0]], axis=0))
                            # bvals_reshaped = np.transpose(bvals_subset).reshape(-1, 1)
                            
                            # # Fit linear model: log(S/S0) = -b * ADC
                            # lin_reg = LinearRegression()
                            # lin_reg.fit(bvals_reshaped, log_signal)
                            
                            # # Predict the fitted line in linear space
                            # fitted_log = lin_reg.predict(bvals_reshaped)
                            # fitted_signal = np.exp(fitted_log)
                            
                            # # Plot fitted line
                            # axs[k].plot(np.transpose(bvals_subset), fitted_signal, 'r--', label='Linear fit')
                            
                            # Set axes
                            if k==0:
                                axs[k].set_ylabel(r'$S / S_0$', fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})  # Correct LaTeX formatting
                            axs[k].set_xlabel('b-val',
                                          fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})
                            axs[k].grid(True)
                            axs[k].set_title(ROI)
                            axs[k].set_yscale('log')
                            axs[k].set_yticks([0.05, 0.1, 1])

                            k += 1
                        plt.savefig(os.path.join(os.path.dirname(os.path.dirname(output_path)), 'SignalDecay_summary.png'))
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
            
            patterns, lims, maximums = get_param_names_model('DTI_DKI',cfg['subject_type'])
            Data_DTIDKI = np.zeros((len(Delta_list), len(ROI_list), len(patterns)))

            for d_idx, Delta in enumerate(Delta_list):
                data_used = f'Delta_{Delta}'
                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path,
                    folderlevel='derivatives', workingdir=cfg['analysis_foldername'], description=f'DTI_DKI_{data_used}')
                output_path = os.path.join(bids_strc_analysis.get_path(), 'Output_masked')
            
                bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas']+'-To-'+data_used+'_fwd', root=data_path, 
                                         folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                bids_strc_reg.set_param(base_name='')
            
                # Define atlas  
                atlas = bids_strc_reg.get_path('atlas_in_dwi.nii.gz')
            
                # Define atlas labels 
                if 'anat_space_organoids' not in  cfg['atlas'] :
                    atlas_labels = prepare_atlas_labels(cfg['atlas'], glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*label*'))[0])
                else:
                    bids_strc_anat = create_bids_structure(subj=subj, sess=sess, datatype='anat', root=data_path, 
                                               folderlevel='derivatives', workingdir=cfg['prep_foldername'])   
                    atlas_labels = prepare_atlas_labels(cfg['atlas'], glob.glob(os.path.join(bids_strc_anat.get_path(), '*label*'))[0])
     
                # Define TPMs
                if cfg['atlas_TPM']:
                    bids_strc_reg_TPM  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas_TPM']+'-To-'+data_used+'_fwd', root=cfg['data_path'] , 
                                             folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                    bids_strc_reg_TPM.set_param(base_name='')
                    TPMs = []
                    for tissue in ['GM', 'WM', 'CSF']:
                        path = bids_strc_reg_TPM.get_path(f'atlas_TPM_{tissue}_in_dwi.nii.gz')
                        TPMs.append(path if os.path.exists(path) else '')
                else:
                    TPMs = []
                
                # Get data
                if os.path.exists(atlas) and os.path.exists(output_path):
                    for r_idx, ROI in enumerate(ROI_list):
                        mask = create_ROI_mask(atlas, atlas_labels, TPMs, ROI, cfg['tpm_thr'], bids_strc_reg)
                        for p_idx, (pattern, maximum) in enumerate(zip(patterns, maximums)):
                            param_img = nib.load(glob.glob(os.path.join(output_path, pattern))[0]).get_fdata()
                            masked = param_img[mask > 0]  # Select only voxels inside the ROI
                            masked_clean = masked[~np.isnan(masked) & (masked > maximum[0]) & (masked < maximum[1])]
                            Data_DTIDKI[d_idx, r_idx, p_idx] = np.nanmean(masked_clean) if len(masked_clean) > 0 else np.nan
            
                # Create table structure with results
                df_data = []
                df_data = pd.DataFrame(Data_DTIDKI[d_idx,:,:], columns=patterns)
                df_data.insert(0, 'ROI Name', ROI_list)
 
                # Save in excel
                outfile = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_DTI_DKI_Delta_{Delta}.xlsx")
                df_data.to_excel(outfile, index=False)
                
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
                
                ax[0].set_ylabel('$MD$'); ax[0].set_ylim([0, 1.5]); ax[0].tick_params(labelbottom=False)
                ax[1].set_ylabel('$MK$'); ax[1].set_ylim([0, 1]); ax[1].tick_params(labelbottom=False)
                ax[2].set_ylabel('$FA$'); ax[2].set_ylim([0, 1]); ax[2].set_xlabel('Diffusion time [ms]')
                
                plt.legend(handles=line_handles, labels=ROI_list, loc='upper right', fontsize=7)
                plt.suptitle('DKI')
                plt.rc('font', size=9)
                plt.savefig(os.path.join(os.path.dirname(os.path.dirname(output_path)), 'DKI.png'))
                plt.close(fig)
                
