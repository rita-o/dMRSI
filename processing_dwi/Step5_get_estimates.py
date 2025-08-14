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
    color_list = distinctipy.get_colors(len(cfg['ROIs_GM'] + cfg['ROIs_WM'])+1, pastel_factor=0.5)

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
            if cfg['lat_ROIS']==1:
                vx_middle = filtered_data['VoxMidHem'].dropna().astype(int).unique()[0]
            else:
                vx_middle = 0
                print('You are not going to have lateralized ROIs. Please ignore the files _left and _right because they are not valid in this case')
            
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
                        ROI_list = cfg['ROIs_GM'].copy() 
                        outfile = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_GM_{model}.xlsx")
                        outfile2 = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_GM_{model}.npy")

                    else:
                        ROI_list = cfg['ROIs_WM'].copy() 
                        outfile = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_WM_{model}.xlsx")
                        outfile2 = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_WM_{model}.npy")

                    # If there is mrs voxel, add it as ROI
                    bids_mrs  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=f"dmrs-to-allDelta-allb", root=data_path, 
                                                   folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                    bids_mrs.set_param(base_name='')
                    if cfg['mrs_vx'] == 1 and os.path.exists(bids_mrs.get_path()):
                        ROI_list.append('voxel_mrs')
                                  
                    ######## EXTRACT MODEL ESTIMATES ########
                    # Option 1. Extract estimates with user defined ROIs
                    
                    # Initialize variables
                    patterns, lims, maximums = get_param_names_model(model,cfg['is_alive'])
                    cleaned_patterns = [re.sub(r'\*{2,}', '*', re.sub(model, '', p, flags=re.IGNORECASE).replace('[^s]', '')).strip('*') for p in patterns]          
                   
                    # Get values of parameters inside the ROIs        
                    Data, Data_all, Data_l, Data_r = get_values_within_ROI(
                        ROI_list, atlas, atlas_labels, TPMs, cfg['tpm_thr'], 
                        vx_middle, patterns, maximums, bids_strc_reg, bids_mrs, output_path)
                    
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
                    
                    # Save summary of means in excel format
                    df_data.to_excel(outfile, index=False)

                    # Save all data in npy format
                    df_data_all = pd.DataFrame(Data_all, columns=cleaned_patterns)
                    df_data_all.insert(0, 'ROI Name', ROI_list)
                    np.save(outfile2, df_data_all)
                    
                    if cfg['lat_ROIS']==1:
                        df_data_l = pd.DataFrame(Data_l, columns=cleaned_patterns)
                        df_data_l.insert(0, 'ROI Name', ROI_list)
                        np.save(outfile2.replace('.npy','_left.npy'), df_data_l)
                                
                        df_data_r = pd.DataFrame(Data_r, columns=cleaned_patterns)
                        df_data_r.insert(0, 'ROI Name', ROI_list)
                        np.save(outfile2.replace('.npy','_right.npy'), df_data_r)

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
                        
                        fig, axs = plt.subplots(n_rows, n_cols, figsize=(8, 4))
                        if len(ROI_list) != 1:
                            axs = axs.flatten()
                        fig.subplots_adjust(wspace=0.05, hspace=0.18, top=0.92, bottom=0.15, left=0.1, right=0.95)

                        if len(ROI_list) == 1:
                                axs = [axs]  # ensure axs is always iterable
                        k=0
                            
                        for ROI in ROI_list:
     
                            if ROI == 'voxel_mrs':
                                mask_indexes = nib.load(bids_mrs.get_path('voxel_mrs.nii.gz')).get_fdata()
                            else:
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
                                    markersize=3
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
                        
                        fig, axs = plt.subplots(n_rows, n_cols, figsize=(8, 4))
                        if len(ROI_list) != 1:
                           axs = axs.flatten()
                        fig.subplots_adjust(wspace=0.05, hspace=0.48, top=0.90, bottom=0.19, left=0.1, right=0.95)

                        if len(ROI_list) == 1:
                                axs = [axs]  # ensure axs is always iterable
                        k=0
                            
                        for ROI in ROI_list:
     
                            # Plot data
                            if ROI == 'voxel_mrs':
                                mask_indexes = nib.load(bids_mrs.get_path('voxel_mrs.nii.gz')).get_fdata()
                            else:
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
                                    markersize=3
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
            
            # If there is mrs voxel, add it as ROI
            bids_mrs  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=f"dmrs-to-allDelta-allb", root=data_path, 
                                           folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
            bids_mrs.set_param(base_name='')
            if cfg['mrs_vx'] == 1 and os.path.exists(bids_mrs.get_path()):
                ROI_list.append('voxel_mrs')

            Data_DTIDKI    = np.empty((len(Delta_list), len(ROI_list), len(patterns)), dtype=object)
            Data_DTIDKI_l  = np.empty((len(Delta_list), len(ROI_list), len(patterns)), dtype=object)
            Data_DTIDKI_r  = np.empty((len(Delta_list), len(ROI_list), len(patterns)), dtype=object)

            # Loop through deltas
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
                    
                    # Get values of parameters inside the ROIs        
                    Data, Data_all, Data_all_l, Data_all_r = get_values_within_ROI(
                        ROI_list, atlas, atlas_labels, TPMs, cfg['tpm_thr'], 
                        vx_middle, patterns, maximums, bids_strc_reg, bids_mrs, output_path)
                    
                    Data_DTIDKI[d_idx, :, :]     = Data_all
                    Data_DTIDKI_l[d_idx, :, :]   = Data_all_l
                    Data_DTIDKI_r[d_idx, :, :]   = Data_all_r
                                    
                # Save summary of means in excel format
                df_data = []
                df_data = pd.DataFrame( np.vectorize(np.nanmean, otypes=[float])(Data_all), columns=cleaned_patterns)
                df_data.insert(0, 'ROI Name', ROI_list)
                outfile = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_DTI_DKI_Delta_{Delta}.xlsx")
                df_data.to_excel(outfile, index=False)
                
                # Save all data in npy format
                outfile2 = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_DTI_DKI_Delta_{Delta}.npy")
                df_data_all = pd.DataFrame(Data_all, columns=cleaned_patterns)
                df_data_all.insert(0, 'ROI Name', ROI_list)
                np.save(outfile2, df_data_all)
                
                if cfg['lat_ROIS']==1:
                    df_data_l = pd.DataFrame(Data_all_l, columns=cleaned_patterns)
                    df_data_l.insert(0, 'ROI Name', ROI_list)
                    np.save(outfile2.replace('.npy','_left.npy'), df_data_l)
                            
                    df_data_r = pd.DataFrame(Data_all_r, columns=cleaned_patterns)
                    df_data_r.insert(0, 'ROI Name', ROI_list)
                    np.save(outfile2.replace('.npy','_right.npy'), df_data_r)
                
            # Plot results
            if os.path.exists(atlas) and os.path.exists(output_path):
                
                # Start figure
                fig, ax = plt.subplots(3, 3, figsize=(9, 6))
                fig.subplots_adjust(wspace=0.02, hspace=0.10, top=0.91, bottom=0.14, left=0.2, right=0.9)
                
                # Main figure
                line_handles = []
                for i, ROI in enumerate(ROI_list):
                    y0 =  np.vectorize(np.nanmean, otypes=[float])(Data_DTIDKI[:, i, 0])
                    s0 =  np.vectorize(np.nanstd, otypes=[float])(Data_DTIDKI[:, i, 0])
                    y1 =  np.vectorize(np.nanmean, otypes=[float])(Data_DTIDKI[:, i, 1])
                    s1 =  np.vectorize(np.nanstd, otypes=[float])(Data_DTIDKI[:, i, 1])
                    y2 =  np.vectorize(np.nanmean, otypes=[float])(Data_DTIDKI[:, i, 2])
                    s2 =  np.vectorize(np.nanstd, otypes=[float])(Data_DTIDKI[:, i, 2])

                    line = plot_with_shade(ax[0,0], Delta_list, y0, s0, color_list[i],
                           linestyle='--', marker='o')
                    plot_with_shade(ax[1,0], Delta_list, y1, s1, color_list[i],
                                    linestyle='--', marker='o')
                    plot_with_shade(ax[2,0], Delta_list, y2, s2, color_list[i],
                                    linestyle='--', marker='o')
                    ax[0, 0].set_title('All', fontsize=12, fontweight='bold')
                    
                    
                    y0 =  np.vectorize(np.nanmean, otypes=[float])(Data_DTIDKI_l[:, i, 0])
                    s0 =  np.vectorize(np.nanstd, otypes=[float])(Data_DTIDKI_l[:, i, 0])
                    y1 =  np.vectorize(np.nanmean, otypes=[float])(Data_DTIDKI_l[:, i, 1])
                    s1 =  np.vectorize(np.nanstd, otypes=[float])(Data_DTIDKI_l[:, i, 1])
                    y2 =  np.vectorize(np.nanmean, otypes=[float])(Data_DTIDKI_l[:, i, 2])
                    s2 =  np.vectorize(np.nanstd, otypes=[float])(Data_DTIDKI_l[:, i, 2])
                    
                    plot_with_shade(ax[0,1], Delta_list, y0, s0, color_list[i],
                    linestyle='--', marker='o')
                    plot_with_shade(ax[1,1], Delta_list, y1, s1, color_list[i],
                                    linestyle='--', marker='o')
                    plot_with_shade(ax[2,1], Delta_list, y2, s2, color_list[i],
                                    linestyle='--', marker='o')
                    ax[0, 1].set_title('Left', fontsize=12, fontweight='bold')
                    
                    y0 =  np.vectorize(np.nanmean, otypes=[float])(Data_DTIDKI_r[:, i, 0])
                    s0 =  np.vectorize(np.nanstd, otypes=[float])(Data_DTIDKI_r[:, i, 0])
                    y1 =  np.vectorize(np.nanmean, otypes=[float])(Data_DTIDKI_r[:, i, 1])
                    s1 =  np.vectorize(np.nanstd, otypes=[float])(Data_DTIDKI_r[:, i, 1])
                    y2 =  np.vectorize(np.nanmean, otypes=[float])(Data_DTIDKI_r[:, i, 2])
                    s2 =  np.vectorize(np.nanstd, otypes=[float])(Data_DTIDKI_r[:, i, 2])
                    
                    plot_with_shade(ax[0,2], Delta_list, y0, s0, color_list[i],
                    linestyle='--', marker='o')
                    plot_with_shade(ax[1,2], Delta_list, y1, s1, color_list[i],
                                    linestyle='--', marker='o')
                    plot_with_shade(ax[2,2], Delta_list, y2, s2, color_list[i],
                                    linestyle='--', marker='o')
                    ax[0, 2].set_title('Right', fontsize=12, fontweight='bold')
    
                   
                    line_handles.append(line)
                
                # Put everything with the same ylim
                for (i, j), a in np.ndenumerate(ax):
                    
                    if cfg['lat_ROIS']==1:
                        m1=np.vectorize(np.nanmean, otypes=[float])(Data_DTIDKI[:, :, i])
                        m2=np.vectorize(np.nanmean, otypes=[float])(Data_DTIDKI_l[:, :, i])
                        m3=np.vectorize(np.nanmean, otypes=[float])(Data_DTIDKI_r[:, :, i])

                        temp = np.min([np.nanmin(m1), np.nanmin(m2), np.nanmin(m3)])
                        lower = round(np.floor(temp/0.05) * 0.05,2)
                        temp = np.min([np.nanmax(m1), np.nanmax(m2), np.nanmax(m3)])
                        upper = round(np.ceil(temp/0.05) * 0.05,2)
                    else:
                        lower = round(np.floor(np.nanmin(np.vectorize(np.nanmean, otypes=[float])(Data_DTIDKI[:, :, i]))/0.05) * 0.05,2)
                        upper = round(np.ceil(np.nanmax(np.vectorize(np.nanmean, otypes=[float])(Data_DTIDKI[:, :, i]))/0.05) * 0.05,2)
                        
                    a.set_ylim([lower, upper])
                    a.set_yticks([lower, upper])
                    
                # Set axis and labels for first column
                ax[0,0].set_ylabel('$MD$'); ax[0,0].tick_params(labelbottom=False)
                ax[1,0].set_ylabel('$MK$'); ax[1,0].tick_params(labelbottom=False)
                ax[2,0].set_ylabel('$FA$'); ax[2,0].set_xlabel('Diffusion time [ms]')
                
                # Remove y-axis labels and ticks for 2nd and 3rd columns
                for col in (1, 2):
                    for r in range(3):
                        ax[r, col].set_yticklabels([])   
                        ax[r, col].set_ylabel('')        
                        ax[r, col].tick_params(labelleft=False, labelbottom=False)  
                ax[2,1].set_xlabel('Diffusion time [ms]'); ax[2,2].set_xlabel('Diffusion time [ms]')
                ax[2,1].tick_params(labelbottom=True); ax[2,2].tick_params(labelbottom=True)  
               
                # Remove axis completly if there is no left/right information
                if cfg['lat_ROIS'] == 0:
                    for col in (1, 2):  
                        for r in (0,1,2):
                            ax[r, col].remove()
            
                # Plot legends, title and save
                plt.legend(handles=line_handles, labels=ROI_list, loc='center left', fontsize=7, bbox_to_anchor=(0.98, 0.5))
                plt.suptitle('DKI')
                plt.rc('font', size=9)
                plt.savefig(os.path.join(os.path.dirname(os.path.dirname(output_path)), 'DKI.png'))
                plt.close(fig)
              
                
            ######## EXTRACT MicroFA ESTIMATES ########
            ######## OPERATIONS INVOLVING THE NEED OF AN ATLAS  ########
            ROI_list = cfg['ROIs_GM'] + cfg['ROIs_WM']

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
                    
                # If there is mrs voxel, add it as ROI
                bids_mrs  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=f"dmrs-to-allDelta-allb", root=data_path, 
                                               folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                bids_mrs.set_param(base_name='')
                if cfg['mrs_vx'] == 1 and os.path.exists(bids_mrs.get_path()):
                    ROI_list.append('voxel_mrs')
                    
                # Determine output file
                outfile = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_Micro_FA.xlsx")
                outfile2 = os.path.join(os.path.dirname(os.path.dirname(output_path)), f"output_ROIs_{cfg['atlas']}_Micro_FA.npy")

                # Option 1. Extract estimates with user defined ROIs
                patterns, lims, maximums = get_param_names_model('Micro_FA',cfg['is_alive'])
                cleaned_patterns = [re.sub(r'\*{2,}', '*', re.sub('Micro_FA', '', p, flags=re.IGNORECASE).replace('[^s]', '')).strip('*') for p in patterns]          
                
                # Get values of parameters inside the ROIs        
                Data, Data_all, Data_l, Data_r = get_values_within_ROI(
                    ROI_list, atlas, atlas_labels, TPMs, cfg['tpm_thr'], 
                    vx_middle, patterns, maximums, bids_strc_reg, bids_mrs, output_path, bids_STE)
                
                # Save summary of means in excel format
                df_data = []
                df_data = pd.DataFrame(Data, columns=cleaned_patterns)
                df_data.insert(0, 'ROI Name', ROI_list)
                df_data.to_excel(outfile, index=False)
                
                # Save all data in npy format
                df_data_all = pd.DataFrame(Data_all, columns=cleaned_patterns)
                df_data_all.insert(0, 'ROI Name', ROI_list)
                np.save(outfile2, df_data_all)
                
                if cfg['lat_ROIS']==1:
                    df_data_l = pd.DataFrame(Data_l, columns=cleaned_patterns)
                    df_data_l.insert(0, 'ROI Name', ROI_list)
                    np.save(outfile2.replace('.npy','_left.npy'), df_data_l)
                            
                    df_data_r = pd.DataFrame(Data_r, columns=cleaned_patterns)
                    df_data_r.insert(0, 'ROI Name', ROI_list)
                    np.save(outfile2.replace('.npy','_right.npy'), df_data_r)
                
            ######## MAKE ROI MAP ########
            bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas']+'-To-'+'allDelta-allb', root=data_path, 
                                          folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
            bids_strc_reg.set_param(base_name='')
            
            roi_paths =  []
            for ROI in cfg['ROIs_GM'] + cfg['ROIs_WM']:
                roi_paths.append(bids_strc_reg.get_path(f'mask_{ROI}.nii.gz'))
           
            QA_ROIs(roi_paths, bids_strc_reg.get_path('ref_dwi.nii.gz'), bids_strc_reg.get_path())