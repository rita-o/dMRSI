#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 13 16:57:50 2025

@author: localadmin
"""


#### ADD LIBRARIES AND CODE PATHS ####  

import os
import sys

code_path           = os.path.join(os.path.expanduser('~'),  'Documents','Rita','Codes_GitHub','dMRSI')   
sys.path.append(code_path )
sys.path.append(os.path.join(code_path, 'processing_dwi'))

import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
import imutils
from scipy.stats import shapiro
import statsmodels.api as sm 
import pandas as pd
import glob
from bids_structure import *
from custom_functions import *
import itertools

plt.close('all')

########################## CONFIG ##########################       

cfg = {
    'atlas': 'Atlas_WHS_v4',
    'atlas_TPM': 'TPM_C57Bl6',
    'ROIs': ['M1','M2','S1','S2', 'V1', 'CC','CG', 'Thal', 'Cereb GM','Cereb WM'],
    'common_folder': os.path.join(os.path.expanduser('~'), 'Documents', 'Rita', 'Data', 'common'),
    'bvals_pershell': [4, 12, 16, 24, 30, 40],
    'methods': ['MPPCA', 'tMPPCA','tMPPCA_5D'],
    'tpm_thr': 0.8
}


color_list = distinctipy.get_colors(len(cfg['methods']),pastel_factor=0)

subj = 'sub-01'
sess = 1
data_path = os.path.join(os.path.expanduser('~'), 'Documents', 'Rita', 'Data', 'dMRI_dMRSI_Pilot_20250428')
save_path = os.path.join(data_path, 'results', 'denoising')
create_directory(save_path)

########################## PREPARE DATA ##########################       

bids_structs = {}
volumes = {}
residuals = {}
#SNR = {}
SNR_gain = {}

# Create BIDS structures and load data
for method in cfg['methods']:
    bids_structs[method] = create_bids_structure(
        subj=subj,
        sess=sess,
        datatype='dwi',
        root=data_path,
        description='allDelta-allb',
        folderlevel='derivatives',
        workingdir=f'preprocessed_{method}'
    )

# Load original and denoised data
vol_dwi = nib.load(bids_structs['MPPCA'].get_path('dwi.nii.gz')).get_fdata()
mask = nib.load(bids_structs['MPPCA'].get_path('mask.nii.gz')).get_fdata()

for method in cfg['methods']:
    volumes[method] = nib.load(bids_structs[method].get_path('dwi_dn.nii.gz')).get_fdata()
    SNR_gain[method] = nib.load(bids_structs[method].get_path('dwi_dn_SNR_gain.nii')).get_fdata()

# Order data by bvals
bvals_path = bids_structs['MPPCA'].get_path('bvalsNom.txt')
bvals = np.loadtxt(bvals_path)
bvals_ordered, bvals_indx = order_bvals(bvals)
bvals_unique = np.unique(bvals_ordered)
vol_dwi = vol_dwi[:, :, :, bvals_indx]

for method in cfg['methods']:
    volumes[method]  = volumes[method][:, :, :, bvals_indx]

# Mask data
for method in cfg['methods']:
    volumes[method] = volumes[method] * mask[..., np.newaxis]
    SNR_gain[method] = SNR_gain[method] * mask
vol_dwi = vol_dwi * mask[..., np.newaxis]
   
# Compute residuals and SNR
for method in cfg['methods']:
    res = volumes[method] - vol_dwi
    residuals[method] = res
    #sigma = np.std(res, axis=3)
    #SNR[method] = np.stack([(volumes[method][:,:,:,i] / sigma) for i in range(volumes[method].shape[-1])], axis=3)

# Create BIDS structures and load data
for method in cfg['methods']:
    bids_structs[method] = create_bids_structure(
        subj=subj,
        sess=sess,
        datatype='dwi',
        root=data_path,
        description='Nexi',
        folderlevel='derivatives',
        workingdir=f'analysis_{method}'
    )
    bids_structs[method].set_param(base_name='')
    
volumes_pwd = {}
for method in cfg['methods']:
    volumes_pwd[method] = nib.load(bids_structs[method].get_path('powderaverage_dwi.nii.gz')).get_fdata()
bvals_pwd = np.loadtxt(bids_structs['MPPCA'].get_path('powderaverage.bval'))
   
########################## PLOT IMAGES DENOISING ##########################       

# Define methods and visualization configuration
vols_to_plot = [
    sum(cfg['bvals_pershell'][0:2]) * 3,
    sum(cfg['bvals_pershell'][0:3]) * 3,
    sum(cfg['bvals_pershell'][0:4]) * 3,
    sum(cfg['bvals_pershell'][0:5]) * 3
]
slice_to_plot = vol_dwi.shape[1] // 2

fig, axes = plt.subplots(len(vols_to_plot), 1 + 2 * len(cfg['methods']), figsize=(15, 10))
fig.subplots_adjust(wspace=0.05, hspace=0.02, top=0.95, bottom=0.05, left=0.05, right=0.95)

for i, vol_idx in enumerate(vols_to_plot):
    title = str(int(bvals_ordered[vol_idx]))
    col = 0

    # Original
    axes[i, col].imshow(imutils.rotate(vol_dwi[:, slice_to_plot, :, vol_idx], angle=90), cmap='gray', vmin=0, vmax=1e4)
    axes[i, col].set_ylabel(f'b-val {title}, vol {vol_idx}')
    if i == 0:
        axes[i, col].set_title('Original')
    axes[i, col].set_xticks([])
    axes[i, col].set_yticks([])
    col += 1

    for method in cfg['methods']:
        vol_s = imutils.rotate(volumes[method][:, slice_to_plot, :, vol_idx], angle=90)
        res_s = imutils.rotate(residuals[method][:, slice_to_plot, :, vol_idx], angle=90)
        data_SNR_gain = SNR_gain[method]
        data_SNR_gain = round(np.mean(data_SNR_gain[mask.astype(bool)]),2)

        axes[i, col].imshow(vol_s, cmap='gray', vmin=0, vmax=1e4)
        if i == 0:
            axes[i, col].set_title(f'{method}, \n SNR gain = {data_SNR_gain}')
        axes[i, col].axis('off')
        col += 1

        axes[i, col].imshow(res_s, cmap='gray', vmin=-3e3, vmax=3e3)
        if i == 0:
            axes[i, col].set_title(f'{method} Res')
        axes[i, col].axis('off')
        col += 1

plt.savefig(os.path.join(save_path, 'Denoising_MPPCAvstMPPCA.png'), bbox_inches='tight', dpi=300)


########################## PLOT RESIDUALS ##########################       

fig, axes = plt.subplots(2, len(cfg['methods']), figsize=(10, 6))
fig.subplots_adjust(top=0.88)
for i, method in enumerate(cfg['methods']):
    
    data = residuals[method]
    data = data[mask.astype(bool), :].reshape(-1)
    data = data[~np.isnan(data)]  

    stat, p = shapiro(data)
    hist_vals, bins, _ = axes[0, i].hist(data, bins=100, color='lightgray', edgecolor='black')

    # Calculate and plot mode
    mode_val = bins[np.argmax(hist_vals)]
    axes[0, i].axvline(mode_val, color='r', linestyle='--', label=f'Mode = {mode_val:.1f}')

    # Plot vertical line at 0
    axes[0, i].axvline(0, color='k', linestyle='-', label='Res = 0')

    axes[0, i].set_title(f'{method} Residual Histogram')
    axes[0, i].legend(fontsize=8)
    axes[0, i].set_xlim(-2000, 2000)

    # Optional: QQ plot
    sm.qqplot(data, ax=axes[1, i], line='45', fit=True)
    axes[1, i].set_title(f'{method} Residual QQ Plot')

fig.suptitle('Residuals Distribution', fontsize=14, weight='bold')
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(os.path.join(save_path, 'Denoising_MPPCAvstMPPCA_residuals.png'), bbox_inches='tight', dpi=300)

########################## PLOT SIGNAL DECAY ##########################     

marker_list = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'h', 'H', 'X']
marker_cycle = itertools.cycle(marker_list)

fig, axes = plt.subplots(2, int(np.ceil(len(cfg['ROIs']) / 2)), figsize=(10, 6))
fig.subplots_adjust(wspace=0.05, hspace=0.2, top=0.9, bottom=0.1, left=0.09, right=0.95)  
axes = axes.flatten()

n_rows = 2
n_cols = int(np.ceil(len(cfg['ROIs']) / 2))

# Loop through ROIs
for i, ROI in enumerate(cfg['ROIs']):
    ax = axes[i]
    mctr = 0

    for mctr, method in enumerate(cfg['methods']):
        
        bids_strc = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path, description= 'Nexi',
                                        folderlevel='derivatives', workingdir='analysis_'+method)
        output_path = os.path.join(bids_strc.get_path(), 'Masked')

        
        # Define atlas
        bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas']+'_To_'+'allDelta-allb', root=data_path, 
                                                folderlevel='derivatives', workingdir='analysis_'+method)
        bids_strc_reg.set_param(base_name='')
        atlas = bids_strc_reg.get_path('atlas_in_dwi.nii.gz')
        
        # Define atlas labels 
        atlas_labels = prepare_atlas_labels(cfg['atlas'], cfg['common_folder'])
        
        # Define TPMs
        bids_strc_reg_TPM  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas_TPM']+'_To_'+'allDelta-allb', root=data_path , 
                                      folderlevel='derivatives', workingdir='analysis_'+method)
        bids_strc_reg_TPM.set_param(base_name='')
        TPMs = []
        for tissue in ['GM', 'WM', 'CSF']:
            path = bids_strc_reg_TPM.get_path(f'atlas_TPM_{tissue}_in_dwi.nii.gz')
            TPMs.append(path if os.path.exists(path) else '')
            
        # Define ROI
        mask_indexes = create_ROI_mask(atlas, atlas_labels, TPMs, ROI, cfg['tpm_thr'], bids_strc_reg)
        
        # Get data in ROI and plot
        marker = next(marker_cycle)
        data = volumes_pwd[method][mask_indexes.astype(bool), :]
        ax.plot(
            bvals_pwd,
            np.nanmean(data, axis=0),
            marker,
            color=color_list[mctr],
            markersize=4,
            alpha=0.6,
            label=method
        )
        mctr += 1

    # Add legend only to the last subplot
    if i == len(cfg['ROIs']) - 1:
        ax.legend()

    ax.set_title(f'{ROI}')

    # Y-axis label only for leftmost plots
    if i % n_cols == 0:
        ax.set_ylabel('S', fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})
    else:
        ax.set_yticklabels([])

    # X-axis label only for bottom row
    if i >= n_cols * (n_rows - 1):
        ax.set_xlabel('Nominal b-val', fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})
    else:
        ax.set_xticklabels([])

    ax.grid(True)

plt.savefig(os.path.join(save_path,'Denoising_MPPCAvstMPPCA_SignalDecay.png'))


########################## MODEL ESTIMATES ##########################   


# Get values
patterns, lims, maximums = get_param_names_model('Nexi')
Data = np.zeros((len(cfg['methods']),len(cfg['ROIs']), len(patterns)))  
k = 0

# Loop through methods
for method in cfg['methods']:
    
    bids_strc = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path, description= 'Nexi',
                                    folderlevel='derivatives', workingdir='analysis_'+method)
    output_path = os.path.join(bids_strc.get_path(), 'Masked')

    bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas']+'_To_'+'allDelta-allb', root=data_path, 
                                            folderlevel='derivatives', workingdir='analysis_'+method)
    bids_strc_reg.set_param(base_name='')
               
    # Define atlas
    atlas = bids_strc_reg.get_path('atlas_in_dwi.nii.gz')
    
    # Define atlas labels 
    atlas_labels = prepare_atlas_labels(cfg['atlas'], cfg['common_folder'])
    
    # Define TPMs
    bids_strc_reg_TPM  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas_TPM']+'_To_'+'allDelta-allb', root=data_path , 
                                 folderlevel='derivatives', workingdir='analysis_'+method)
    bids_strc_reg_TPM.set_param(base_name='')
    bids_strc_reg_TPM.set_param(base_name='')
    TPMs = []
    for tissue in ['GM', 'WM', 'CSF']:
        path = bids_strc_reg_TPM.get_path(f'atlas_TPM_{tissue}_in_dwi.nii.gz')
        TPMs.append(path if os.path.exists(path) else '')

    # Loop through ROIs
    ROI_ctr = 0
    for ROI in cfg['ROIs']:
        
        # Define ROI
        mask_indexes = create_ROI_mask(atlas, atlas_labels, TPMs, ROI, cfg['tpm_thr'], bids_strc_reg)
        
        # Loop through each parameter outputed in the model
        pattern_ctr = 0
        for (pattern, maximum) in zip(patterns, maximums):
            
            # Load parameter data
            parameter_data = nib.load(glob.glob(os.path.join(output_path, pattern))[0]).get_fdata()
            
            # Mask the parameter with the ROI
            data = parameter_data[mask_indexes.astype(bool)]
            data = data[~np.isnan(data) & (data > maximum[0]) & (data < maximum[1])]

            # Store the mean of the masked data
            Data[k,ROI_ctr, pattern_ctr] = np.nanmean(data) if len(data) > 0 else np.nan
    
            pattern_ctr += 1
        
        ROI_ctr += 1
   
    k += 1
    
num_patterns = Data.shape[2]  
num_ROIs = Data.shape[1]  
num_methods = len(cfg['methods'])
bar_width = 0.8 / num_methods  # dynamic width

fig, axes = plt.subplots(num_patterns, 1, figsize=(10, 4 * num_patterns), sharex=True)

for param_idx in range(num_patterns):  
    ax = axes[param_idx] if num_patterns > 1 else axes

    x = np.arange(num_ROIs)  # ROIs on x-axis

    for mctr, method in enumerate(cfg['methods']):
        # Offset x for each method
        offset = (mctr - num_methods / 2) * bar_width + bar_width / 2

        data_to_plot = Data[mctr, :, param_idx]
        ax.bar(x + offset, data_to_plot, width=bar_width, color=color_list[mctr],label=method if param_idx == 0 else "")
        
        
    # Label and formatting
    paramname = patterns[param_idx][1:-1]
    ax.set_ylabel(f'{paramname}', fontsize=12, fontweight='bold')
    
    if param_idx == 1:
        ax.set_ylim(1, 4)
    elif param_idx == 2:
        ax.set_ylim(0.5, 1.5)
    elif param_idx == 3:
        ax.set_ylim(0, 1)

    ax.grid(True)

# X-axis tick labels
axes[-1].set_xticks(x)
axes[-1].set_xticklabels(cfg['ROIs'], rotation=45, ha='right')

# Add legend
axes[0].legend()

plt.tight_layout()

plt.subplots_adjust(wspace=0.05,hspace=0.05, top=0.95, bottom=0.1, left=0.1, right=0.90)  
plt.suptitle('Parameter Values Across ROIs', fontsize=14, fontweight='bold')
plt.savefig(os.path.join(save_path,'Denoising_MPPCAvstMPPCA_NexiEstimates.png'))