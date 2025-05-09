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

                bids_strc_analysis = create_bids_structure(subj, sess, 'dwi', model, data_path, 'derivatives', cfg['analysis_foldername'])
                output_path = os.path.join(bids_strc_analysis.get_path(), 'Masked')

                bids_strc_reg = create_bids_structure(subj, sess, 'registration', f"{cfg['atlas']}_To_{data_used}", data_path, 'derivatives', cfg['analysis_foldername'])
                bids_strc_reg.set_param(base_name='')
                atlas = bids_strc_reg.get_path('atlas_in_dwi.nii.gz')

                atlas_label_path = glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*atlas.label'))[0]
                atlas_labels = pd.read_csv(atlas_label_path, sep=r'\s+', skiprows=14, header=None,
                                           names=['IDX', 'R', 'G', 'B', 'A', 'VIS', 'MSH', 'LABEL'], quotechar='"')

                bids_strc_reg_TPM = create_bids_structure(subj, 1, 'registration', f"{cfg['atlas_TPM']}_To_{data_used}_fwd", cfg['data_path'], 'derivatives', cfg['analysis_foldername'])
                bids_strc_reg_TPM.set_param(base_name='')
                TPMs = [bids_strc_reg_TPM.get_path(f'atlas_TPM_{tissue}_in_dwi.nii.gz') for tissue in ['GM', 'WM', 'CSF']]

                # Determine ROI list and output file
                if model in cfg['model_list_GM']:
                    ROI_list = cfg['ROIs_GM']
                    outfile = os.path.join(os.path.dirname(os.path.dirname(output_path)), f'output_ROIs_GM_{model}.xlsx')
                else:
                    ROI_list = cfg['ROIs_WM']
                    outfile = os.path.join(os.path.dirname(os.path.dirname(output_path)), f'output_ROIs_WM_{model}.xlsx')

                ######## EXTRACT MODEL ESTIMATES ########

                patterns, lims = get_param_names_model(model)
                Data = np.zeros((len(ROI_list), len(patterns)))
                for i, ROI in enumerate(ROI_list):
                    mask = create_ROI_mask(atlas, atlas_labels, TPMs, ROI, bids_strc_reg)
                    for j, pattern in enumerate(patterns):
                        param_img = nib.load(glob.glob(os.path.join(output_path, pattern))[0]).get_fdata()
                        masked = param_img * mask
                        flat = masked.reshape(-1, 1)
                        valid = flat[~(np.isnan(flat) | (flat == 0)).any(axis=1)]
                        Data[i, j] = np.nanmean(valid) if len(valid) > 0 else np.nan

                df_data = pd.DataFrame(Data, columns=patterns)
                df_data.insert(0, 'ROI Name', ROI_list)
                df_data.to_excel(outfile, index=False)


                ######## PLOT SNR ########
                if model == 'Nexi':
                    bvals = read_numeric_txt(os.path.join(bids_strc_analysis.get_path(), 'powderaverage.bval'))
                    S_S0 = nib.load(os.path.join(bids_strc_analysis.get_path(), 'powderaverage_dwi.nii.gz')).get_fdata()
                    nf = nib.load(os.path.join(bids_strc_analysis.get_path(), 'normalized_sigma.nii.gz')).get_fdata() * np.sqrt(np.pi / 2)

                    fig, axs = plt.subplots(1, len(ROI_list), figsize=(12, 4))
                    for k, ROI in enumerate(ROI_list):
                        mask = create_ROI_mask(atlas, atlas_labels, TPMs, ROI, bids_strc_reg)

                        masked_S_S0 = np.array([S_S0[..., v] * mask for v in range(S_S0.shape[-1])]).transpose(1,2,3,0)
                        data = masked_S_S0.reshape(-1, S_S0.shape[-1])
                        data = data[~(np.isnan(data).any(axis=1) | (data == 0).any(axis=1))]

                        nf_masked = (nf * mask).reshape(-1, 1)
                        nf_masked = nf_masked[~(np.isnan(nf_masked) | (nf_masked == 0)).any(axis=1)]

                        axs[k].plot(bvals, np.nanmean(data, axis=0), 'bo', markersize=3)
                        axs[k].plot(bvals, [np.nanmean(nf_masked)] * len(bvals))
                        axs[k].set_title(ROI)
                        axs[k].set_xlabel('b-val', fontsize=12, fontstyle='italic', weight='bold')
                        if k == 0:
                            axs[k].set_ylabel(r'$S / S_0$', fontsize=12, fontstyle='italic', weight='bold')
                        axs[k].grid(True)

                    plt.tight_layout()
                    outfile = os.path.join(os.path.dirname(os.path.dirname(output_path)), 'SignalDecay_summary.png')
                    plt.savefig(outfile)
                    plt.close(fig)
                    
                ######## PLOT SUMMARY MAP PLOT ########
     
                # Define BIDS structure and output path
                bids_strc_analysis = create_bids_structure(
                    subj=subj,
                    sess=sess,
                    datatype='dwi',
                    root=data_path,
                    folderlevel='derivatives',
                    workingdir=cfg['analysis_foldername'],
                    description=model
                )
                output_path = os.path.join(bids_strc_analysis.get_path(), 'Masked')
            
                # Get model parameter names and display ranges
                patterns, lims = get_param_names_model(model)
                
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
                    img = imutils.rotate(param_data[:, 10, :], angle=90)
                    fact = int((max(img.shape) - min(img.shape)) / 2)
                    img = img[fact:-fact, :]
                    img[np.isnan(img)] = 0
                
                    # Show slice
                    im = ax.imshow(img, cmap='jet', vmin=lim[0], vmax=lim[1])
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
                fig_path = os.path.join(bids_strc_analysis.get_path(), 'output_summary.png')
                plt.savefig(fig_path)
                plt.close(fig)
                
            ######## EXTRACT DTI,DKI ESTIMATES ########

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
            
            patterns, lims = get_param_names_model('DTI_DKI_short')
            Data_DTIDKI = np.zeros((len(Delta_list), len(ROI_list), len(patterns)))
            
            for d_idx, Delta in enumerate(Delta_list):
                desc = f'Delta_{Delta}'
                bids_strc_analysis = create_bids_structure(subj, sess, 'dwi', f'DTI_DKI_{desc}', data_path, 'derivatives', cfg['analysis_foldername'])
                output_path = os.path.join(bids_strc_analysis.get_path(), 'Masked')
            
                bids_strc_reg = create_bids_structure(subj, sess, 'registration', f"{cfg['atlas']}_To_{desc}_fwd", data_path, 'derivatives', cfg['analysis_foldername'])
                bids_strc_reg.set_param(base_name='')
                atlas = bids_strc_reg.get_path('atlas_in_dwi.nii.gz')
            
                atlas_labels = pd.read_csv(
                    glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*atlas.label'))[0],
                    sep=r'\s+', skiprows=14, header=None,
                    names=['IDX', 'R', 'G', 'B', 'A', 'VIS', 'MSH', 'LABEL'], quotechar='"')
            
                bids_strc_reg_TPM = create_bids_structure(subj, 1, 'registration', f"{cfg['atlas_TPM']}_To_{desc}_fwd", cfg['data_path'], 'derivatives', cfg['analysis_foldername'])
                bids_strc_reg_TPM.set_param(base_name='')
                TPMs = [bids_strc_reg_TPM.get_path(f'atlas_TPM_{t}_in_dwi.nii.gz') for t in ['GM', 'WM', 'CSF']]
            
                for r_idx, ROI in enumerate(ROI_list):
                    mask = create_ROI_mask(atlas, atlas_labels, TPMs, ROI, bids_strc_reg)
                    for p_idx, pattern in enumerate(patterns):
                        param_data = nib.load(glob.glob(os.path.join(output_path, pattern))[0]).get_fdata()
                        masked = param_data * mask
                        flat = masked.reshape(-1, 1)
                        valid = flat[~((np.isnan(flat)) | (flat == 0) | (flat > 100) | (flat < 0)).any(axis=1)]
                        Data_DTIDKI[d_idx, r_idx, p_idx] = np.nanmean(valid) if len(valid) > 0 else np.nan
            
            # Plot results
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
                
