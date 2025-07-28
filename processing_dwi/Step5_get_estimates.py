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
    import distinctipy
    color_list = distinctipy.get_colors(len(cfg['ROIs_GM'] + cfg['ROIs_WM']), pastel_factor=0.5)

    ######## SUBJECT-WISE OPERATIONS ########
    for subj in subj_list:
        print(f'Getting model estimates for {subj}...')
        
        subj_data = scan_list[scan_list['newstudyName'] == subj].reset_index(drop=True)
        sess_list = [x for x in subj_data['blockNo'].unique() if not math.isnan(x)]

        ######## SESSION-WISE OPERATIONS ########
        for sess in sess_list:
            print(f'Session {sess}...')
            
            filtered_data = subj_data[
                (subj_data['phaseDir'] == 'fwd') &
                (subj_data['blockNo'] == sess) &
                (subj_data['noBval'] > 1) &
                (subj_data['acqType'] == 'PGSE') &
                (subj_data['scanQA'] == 'ok')
            ]
            Delta_list = sorted(filtered_data['diffTime'].dropna().astype(int).unique())
            
            ######## MODEL-WISE OPERATIONS ########
            for model in cfg['model_list']:
                print(f'Getting model estimates from {model}...')

                filtered_data = subj_data[(subj_data['acqType'] == 'PGSE') &
                                          (subj_data['phaseDir'] == 'fwd') &
                                          (subj_data['blockNo'] == sess) &
                                          (subj_data['noBval'] > 1)]
                
               
                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path, 
                             folderlevel='derivatives', workingdir=cfg['analysis_foldername'], description=model)
                output_path = os.path.join(bids_strc_analysis.get_path(), 'Output_masked')

                bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas']+'-To-'+'allDelta-allb', root=data_path, 
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
                        bids_strc_reg_TPM  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas_TPM']+'-To-'+'allDelta-allb', root=cfg['data_path'] , 
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
                        outfile2 = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_GM_{model}.npy")

                    else:
                        ROI_list = cfg['ROIs_WM']
                        outfile = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_WM_{model}.xlsx")
                        outfile2 = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_WM_{model}.npy")

                    ######## EXTRACT MODEL ESTIMATES ########
                    # Option 1. Extract estimates with user defined ROIs
                    patterns, lims, maximums = get_param_names_model(model,cfg['is_alive'])
                    cleaned_patterns = [re.sub(r'\*{2,}', '*', re.sub(model, '', p, flags=re.IGNORECASE).replace('[^s]', '')).strip('*') for p in patterns]          
                    Data = np.zeros((len(ROI_list), len(patterns)))
                    Data_all = np.empty((len(ROI_list), len(patterns)), dtype=object)
                    for i, ROI in enumerate(ROI_list):
                        mask = create_ROI_mask(atlas, atlas_labels, TPMs, ROI, cfg['tpm_thr'], bids_strc_reg)
                        for j, (pattern, maximum) in enumerate(zip(patterns, maximums)): 
                            matched_file = glob.glob(os.path.join(output_path, pattern))

                            # Filter out files where 'fs' appears in the filename for sandi when looking for f
                            if pattern=='*sandi*f*':
                                matched_file = [
                                    f for f in matched_file
                                    if 'fs' not in os.path.basename(f).lower()  
                                ]
                            
                            param_img = nib.load(matched_file[0]).get_fdata()
                            masked = param_img[mask > 0]  # Select only voxels inside the ROI
                            masked_clean = masked[~np.isnan(masked) & (masked > maximum[0]) & (masked < maximum[1])]
                            Data[i, j]     = np.nanmean(masked_clean) if len(masked_clean) > 0 else np.nan
                            Data_all[i, j] = masked_clean if len(masked_clean) > 0 else np.nan

                    # Create table structure with results
                    df_data = pd.DataFrame(Data, columns=cleaned_patterns)
                    df_data.insert(0, 'ROI Name', ROI_list)
                    
                    # Option 2. Extract estimates from all labels in the atlas (dont save nifiti files of mask)
                    # Personally think it's not very relevant
                    # Data2 = np.zeros((len(atlas_labels['IDX'].to_numpy()), len(patterns)))
                    # for i, indx in enumerate(atlas_labels['IDX'].to_numpy()):
                    #     mask = create_ROI_mask_fromindx(atlas, atlas_labels, indx, bids_strc_reg)

                    #     for j, (pattern, maximum) in enumerate(zip(patterns, maximums)):
                    #         param_img = nib.load(glob.glob(os.path.join(output_path, pattern))[0]).get_fdata()
                    #         masked = param_img[mask > 0]  # Select only voxels inside the ROI
                    #         #masked_clean = masked[~np.isnan(masked) & (masked > maximum[0]) & (masked < maximum[1])]
                    #         masked_clean = masked[~np.isnan(masked)]
                    #         Data2[i, j] = np.median(masked_clean) if len(masked_clean) > 0 else np.nan

                    # Add mean of medians
                    #df_data.loc[len(df_data)] = ['mean of medians'] + np.nanmean(Data2, axis=0).tolist()
                    
                    # Save in excel summary
                    df_data.to_excel(outfile, index=False)

                    # Save all data
                    df_data_all = pd.DataFrame(Data_all, columns=cleaned_patterns)
                    df_data_all.insert(0, 'ROI Name', ROI_list)
                    np.save(outfile2, df_data_all)

                    ######## PLOT SNR ########
                    color_list_snr  =  [(0.3462032775519162, 0.3531070236388303, 0.8723545410491512), (0.4324776646215159, 0.9594894437563749, 0.33498906309524773), (0.9406821354408894, 0.3893640567923122, 0.3692878637990247), (0.4072131417883373, 0.8783467338575792, 0.9732166499618664), (0.883784919700095, 0.5084984818886036, 0.9831871536574889), (0.9725395306324506, 0.954894866579424, 0.37771765906993093)]

                    if model == 'Nexi':
                        print(f'Plotting powder average signal within ROIs...')

                        # Load data
                        bvals = read_numeric_txt(os.path.join(bids_strc_analysis.get_path(),'powderaverage.bval'))
                        S_S0  = nib.load(os.path.join(bids_strc_analysis.get_path(),'powderaverage_dwi.nii.gz')).get_fdata()
                        nf = nib.load(os.path.join(bids_strc_analysis.get_path(),'normalized_sigma.nii.gz')).get_fdata()*np.sqrt(np.pi/2)
     
                        # organize
                        bvals_split = np.split(bvals[0], len(Delta_list))

                        # Loop through ROIs     
                        n_params = len(ROI_list)
                        n_rows = 1 if n_params <= 4 else 2
                        n_cols = math.ceil(n_params / n_rows)
                        
                        fig, axs = plt.subplots(n_rows, n_cols, figsize=(12, 4))
                        if len(ROI_list) != 1:
                            axs = axs.flatten()
                        fig.subplots_adjust(wspace=0.05, hspace=0.2, top=0.92, bottom=0.15, left=0.2, right=0.95)

                        if len(ROI_list) == 1:
                                axs = [axs]  # ensure axs is always iterable
                        k=0
                            
                        for ROI in ROI_list:
     
                            mask_indexes = create_ROI_mask(atlas, atlas_labels, TPMs, ROI, cfg['tpm_thr'], bids_strc_reg)
    
                            S_S0_masked = copy.deepcopy(S_S0)
                            for v in range(S_S0_masked.shape[-1]):
                                S_S0_masked[:, :, :, v] = np.multiply(S_S0_masked[:, :, :, v], mask_indexes)
                                
                            data = S_S0_masked.reshape(S_S0_masked.shape[0]*S_S0_masked.shape[1]*S_S0_masked.shape[2], S_S0_masked.shape[3])
                            data = data[~(np.isnan(data).any(axis=1) | (data == 0).any(axis=1))]
                            data_split = np.split(data, len(Delta_list), axis=1)
                            
                            nf_masked = nf * mask_indexes
                            nf_masked = nf_masked.reshape(nf_masked.shape[0]*nf_masked.shape[1]*nf_masked.shape[2], 1)
                            nf_masked = nf_masked[~(np.isnan(nf_masked).any(axis=1) | (nf_masked == 0).any(axis=1))]
    
                            # Plot data
                            for i in range(len(Delta_list)):
                                axs[k].plot(
                                    np.transpose(bvals_split[i]),
                                    np.nanmean(data_split[i], axis=0),
                                    marker='o',
                                    linestyle='--',
                                    color=color_list_snr[i],
                                    label=f'$\\Delta$= {Delta_list[i]}',
                                    markersize=4
                                )
                            axs[k].plot(np.transpose(bvals_split[0]), np.repeat(np.nanmean(nf_masked), np.transpose(bvals_split[0]).shape[0]),color='black',label='nf')

                            # Settings
                            axs[k].set_yscale('log')
                            row = k // n_cols
                            col = k % n_cols
                            axs[k].set_ylim([0.02, 1])
                            axs[k].set_xticks(bvals_split[i])
                            axs[k].set_yticks([0.02, 0.1, 1])
                            if col == 0:
                               axs[k].set_ylabel(r'$S / S_0$', fontdict={'size': 10, 'weight': 'bold', 'style': 'italic'})
                               if row==0:
                                   axs[k].legend()
                            else:
                                axs[k].set_yticklabels([])
                            if row == n_rows -1:
                               axs[k].set_xlabel(r'$b$ $[ms/µm^2]$', fontdict={'size': 10, 'weight': 'bold', 'style': 'italic'})
                               axs[k].set_xticklabels(np.round(bvals_split[i]).astype(int))
                            else:
                               axs[k].set_xticklabels([])
                            axs[k].grid(True)
                            axs[k].set_title(ROI)
                            k += 1
                        n_used = n_params 
                        if len(axs) > n_used:
                             for ax in axs[n_used:]:
                                 ax.set_visible(False)  
                       
                        plt.savefig(os.path.join(os.path.dirname(os.path.dirname(output_path)), 'SignalDecay_summary.png'))
                        plt.close(fig)
                        
                    ######## PLOT Signal vs b^-1/2 ########

                    if model == 'Nexi':
                        print(f'Plotting powder average signal within ROIs against 1/sqrt(b)...')

                        # Load data
                        bvals = read_numeric_txt(os.path.join(bids_strc_analysis.get_path(),'powderaverage.bval'))
                        S_S0  = nib.load(os.path.join(bids_strc_analysis.get_path(),'powderaverage_dwi.nii.gz')).get_fdata()
     
                        # organize
                        bvals_split = np.split(bvals[0], len(Delta_list))

                        # Loop through ROIs     
                        n_params = len(ROI_list)
                        n_rows = 1 if n_params <= 4 else 2
                        n_cols = math.ceil(n_params / n_rows)
                        
                        fig, axs = plt.subplots(n_rows, n_cols, figsize=(12, 4))
                        if len(ROI_list) != 1:
                           axs = axs.flatten()
                        fig.subplots_adjust(wspace=0.05, hspace=0.45, top=0.90, bottom=0.19, left=0.2, right=0.95)

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
                            data_split = np.split(data, len(Delta_list), axis=1)
                            
                            # Plot data
                            for i in range(len(Delta_list)):
                                axs[k].plot(
                                    1 / np.sqrt(bvals_split[i]),
                                    np.nanmean(data_split[i], axis=0),
                                    marker='o',
                                    linestyle='--',
                                    color=color_list_snr[i],
                                    label=f'$\\Delta$= {Delta_list[i]}',
                                    markersize=4
                                )
                                
                            axs[k].set_yscale('log')
                            row = k // n_cols
                            col = k % n_cols
                            axs[k].set_ylim([0.02, 1])
                            axs[k].set_xticks(1 / np.sqrt(bvals_split[i]))
                            axs[k].set_yticks([0.02, 0.1, 1])
                            if col == 0:
                               axs[k].set_ylabel(r'$S / S_0$', fontdict={'size': 10, 'weight': 'bold', 'style': 'italic'})
                               if row==0:
                                   axs[k].legend()
                            else:
                                axs[k].set_yticklabels([])
                            if row == n_rows -1:
                               axs[k].set_xlabel(r'$1/b^{-1/2}$ $[µm/ms^{1/2}]$', fontdict={'size': 10, 'weight': 'bold', 'style': 'italic'})
                               axs[k].set_xticks( 1 / np.sqrt(bvals_split[i]))
                               axs[k].set_xticklabels(np.round(1 / np.sqrt(bvals_split[i]), 2))
                               axs[k].set_xticklabels(axs[k].get_xticklabels(), rotation=45)
                               ax_top = axs[k].twiny()
                               ax_top.set_xlim(axs[k].get_xlim())
                               tick_positions = 1 / np.sqrt(bvals_split[i])
                               ax_top.set_xticks(tick_positions)
                               ax_top.set_xticklabels(np.round(bvals_split[i], 0).astype(int))


                            else:
                               axs[k].set_xticklabels([])
                            axs[k].grid(True)
                            axs[k].set_title(ROI)
                            k += 1
                        n_used = n_params 
                        if len(axs) > n_used:
                             for ax in axs[n_used:]:
                                 ax.set_visible(False)
                                 
                        plt.savefig(os.path.join(os.path.dirname(os.path.dirname(output_path)), 'SignalDecay_summary2.png'))
                        plt.close(fig)
                    
                
            ######## EXTRACT DTI,DKI ESTIMATES ########
            ######## OPERATIONS INVOLVING THE NEED OF AN ATLAS  ########

            ROI_list = cfg['ROIs_GM'] + cfg['ROIs_WM']
            
            patterns, lims, maximums = get_param_names_model('DTI_DKI',cfg['is_alive'])
            cleaned_patterns = [re.sub(r'\*{2,}', '*', re.sub(model, '', p, flags=re.IGNORECASE).replace('[^s]', '')).strip('*') for p in patterns]          
            Data_DTIDKI      = np.zeros((len(Delta_list), len(ROI_list), len(patterns)))
            Data_DTIDKI_all  = np.empty((len(Delta_list), len(ROI_list), len(patterns)), dtype=object)

            for d_idx, Delta in enumerate(Delta_list):
                print(f'Getting model estimates from DTI/DKI for Delta={Delta}...')

                data_used = f'Delta_{Delta}'
                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path,
                    folderlevel='derivatives', workingdir=cfg['analysis_foldername'], description=f'DTI_DKI_{data_used}')
                output_path = os.path.join(bids_strc_analysis.get_path(), 'Output_masked')
            
                bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas']+'-To-'+'allDelta-allb', root=data_path, 
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
                    bids_strc_reg_TPM  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas_TPM']+'-To-'+'allDelta-allb', root=cfg['data_path'] , 
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
                            Data_DTIDKI_all[d_idx, r_idx, p_idx] = masked_clean if len(masked_clean) > 0 else np.nan

                # Create table structure with summary results
                df_data = []
                df_data = pd.DataFrame(Data_DTIDKI[d_idx,:,:], columns=patterns)
                df_data.insert(0, 'ROI Name', ROI_list)
 
                # Save in excel
                outfile = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_DTI_DKI_Delta_{Delta}.xlsx")
                df_data.to_excel(outfile, index=False)
                
                # Save all data
                outfile2 = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_DTI_DKI_Delta_{Delta}.npy")
                df_data_all = pd.DataFrame(Data_DTIDKI_all[d_idx,:,:], columns=cleaned_patterns)
                df_data_all.insert(0, 'ROI Name', ROI_list)
                np.save(outfile2, df_data_all)
                
            # Plot results
            if os.path.exists(atlas) and os.path.exists(output_path):
                
                fig, ax = plt.subplots(3, 1, figsize=(3, 6))
                fig.subplots_adjust(wspace=0.01, hspace=0.10, top=0.91, bottom=0.14, left=0.2, right=0.9)
                
                line_handles = []
                for i, ROI in enumerate(ROI_list):
                    line, = ax[0].plot(Delta_list, Data_DTIDKI[:, i, 0], linestyle='--', marker='o', c=color_list[i])
                    ax[1].plot(Delta_list, Data_DTIDKI[:, i, 1], linestyle='--', marker='o', c=color_list[i])
                    ax[2].plot(Delta_list, Data_DTIDKI[:, i, 2], linestyle='--', marker='o', c=color_list[i])
                    line_handles.append(line)
                
                ax[0].set_ylabel('$MD$'); ax[0].tick_params(labelbottom=False)
                ax[1].set_ylabel('$MK$'); ax[1].tick_params(labelbottom=False)
                ax[2].set_ylabel('$FA$'); ax[2].set_xlabel('Diffusion time [ms]')
                
                for i, a in enumerate(ax):
                   
                    lower = round(np.floor(np.nanmin(Data_DTIDKI[:, :, i])/0.05) * 0.05,2)
                    upper = round(np.ceil(np.nanmax(Data_DTIDKI[:, :, i])/0.05) * 0.05,2)

                    a.set_ylim([lower, upper])
                    a.set_yticks([lower, upper])
                    
                plt.legend(handles=line_handles, labels=ROI_list, loc='upper right', fontsize=7)
                plt.suptitle('DKI')
                plt.rc('font', size=9)
                plt.savefig(os.path.join(os.path.dirname(os.path.dirname(output_path)), 'DKI.png'))
                plt.close(fig)
              
                
            ######## EXTRACT MicroFA ESTIMATES ########
            ######## OPERATIONS INVOLVING THE NEED OF AN ATLAS  ########

            bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                          folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='microFA')               
            bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas']+'-To-'+'allDelta-allb', root=data_path, 
                                          folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
            bids_strc_reg.set_param(base_name='')
             
            # Define atlas
            atlas = bids_strc_reg.get_path('atlas_in_dwi.nii.gz')

            if os.path.exists(atlas) and os.path.exists(bids_STE.get_path()):
                print(f'Getting model estimates from Micro FA...')

                # Define atlas labels 
                if 'anat_space_organoids' not in  cfg['atlas'] :
                    atlas_labels = prepare_atlas_labels(cfg['atlas'], glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*label*'))[0])
                else:
                    bids_strc_anat = create_bids_structure(subj=subj, sess=sess, datatype='anat', root=data_path, 
                                               folderlevel='derivatives', workingdir=cfg['prep_foldername'])   
                    atlas_labels = prepare_atlas_labels(cfg['atlas'], glob.glob(os.path.join(bids_strc_anat.get_path(), '*label*'))[0])

                # Define TPMs
                if cfg['atlas_TPM']:
                    bids_strc_reg_TPM  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas_TPM']+'-To-'+'allDelta-allb', root=cfg['data_path'] , 
                                                 folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                    bids_strc_reg_TPM.set_param(base_name='')
                    TPMs = []
                    for tissue in ['GM', 'WM', 'CSF']:
                        path = bids_strc_reg_TPM.get_path(f'atlas_TPM_{tissue}_in_dwi.nii.gz')
                        TPMs.append(path if os.path.exists(path) else '')
                else:
                    TPMs = []

                # Determine ROI list and output file
                outfile = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_Micro_FA.xlsx")
                outfile2 = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_Micro_FA.npy")

                # Option 1. Extract estimates with user defined ROIs
                patterns, lims, maximums = get_param_names_model('Micro_FA',cfg['is_alive'])
                cleaned_patterns = [re.sub(r'\*{2,}', '*', re.sub('Micro_FA', '', p, flags=re.IGNORECASE).replace('[^s]', '')).strip('*') for p in patterns]          
                Data = np.zeros((len(ROI_list), len(patterns)))
                Data_all = np.empty((len(ROI_list), len(patterns)), dtype=object)
                for i, ROI in enumerate(ROI_list):
                    if 'CSF' in ROI:
                        bids_STE.set_param(description='microFA_lowb')
                    else:
                        bids_STE.set_param(description='microFA')
                    output_path = bids_STE.get_path()
                    mask = create_ROI_mask(atlas, atlas_labels, TPMs, ROI, cfg['tpm_thr'], bids_strc_reg)
                    for j, (pattern, maximum) in enumerate(zip(patterns, maximums)): 
                        param_img = nib.load(glob.glob(os.path.join(output_path, pattern))[0]).get_fdata()
                        masked = param_img[mask > 0]  # Select only voxels inside the ROI
                        masked_clean = masked[~np.isnan(masked) & (masked > maximum[0]) & (masked < maximum[1])]
                        Data[i, j]     = np.nanmean(masked_clean) if len(masked_clean) > 0 else np.nan
                        Data_all[i, j] = masked_clean if len(masked_clean) > 0 else np.nan

                # Create table structure with results
                df_data = pd.DataFrame(Data, columns=cleaned_patterns)
                df_data.insert(0, 'ROI Name', ROI_list)
                
                # Save in excel summary
                df_data.to_excel(outfile, index=False)

                # Save all data
                df_data_all = pd.DataFrame(Data_all, columns=cleaned_patterns)
                df_data_all.insert(0, 'ROI Name', ROI_list)
                np.save(outfile2, df_data_all)
                
            ######## MAKE ROI MAP ########
            bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas']+'-To-'+'allDelta-allb', root=data_path, 
                                          folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
            bids_strc_reg.set_param(base_name='')
            
            roi_paths =  []
            for ROI in cfg['ROIs_GM'] + cfg['ROIs_WM']:
                roi_paths.append(bids_strc_reg.get_path(f'mask_{ROI}.nii.gz'))
                
            QA_ROIs(roi_paths, bids_strc_reg.get_path('ref_dwi.nii.gz'), bids_strc_reg.get_path())

