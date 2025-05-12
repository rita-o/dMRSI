#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  5 16:45:15 2025

@author: localadmin
"""

"""
Compare two images 
Last changed Feb 2025
@author: Rita O
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

from bids_structure import *
from custom_functions import *


########################## PREPARE DATA ##########################       

# Define BIDs and paths
plt.close('all')
cfg                         = {}
cfg['atlas']                = 'Atlas_WHS_v4'
cfg['atlas_TPM']            = 'TPM_C57Bl6'
cfg['ROIs']                 = ['M1','M2','S1','S2', 'V1', 'CC','CG', 'Thal', 'cerebellum']
cfg['common_folder']        = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','common')
subj  ='sub-01'
sess  = 1
data_path   = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','dMRI_dMRSI_Pilot_20250428')
save_path = os.path.join(data_path,'results','denoising')
create_directory(save_path)
bvals_pershell =  [4, 12, 16, 24, 30, 40]

bids_strc_MPPCA = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path, description= 'allDelta-allb',
                                folderlevel='derivatives', workingdir='preprocessed_MPPCA')
bids_strc_tMPPCA = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path, description= 'allDelta-allb',
                                folderlevel='derivatives', workingdir='preprocessed_tMPPCA')

# Load images
vol_dwi     = nib.load(bids_strc_MPPCA.get_path('dwi.nii.gz')).get_fdata()
vol_MPPCA   = nib.load(bids_strc_MPPCA.get_path('dwi_dn.nii.gz')).get_fdata()
vol_tMPPCA  = nib.load(bids_strc_tMPPCA.get_path('dwi_dn.nii.gz')).get_fdata()
mask        = nib.load(bids_strc_MPPCA.get_path('mask.nii.gz')).get_fdata()

# Order data by bvals as they were acquired randomly
bvals_path = bids_strc_MPPCA.get_path('bvalsNom.txt')
with open(bvals_path, 'r') as file:
    data = file.read()
    bvals = [float(num) for num in data.split()]

bvals_ordered, bvals_indx = order_bvals(bvals)
bvals_unique  = np.unique(bvals_ordered)
vol_dwi       = vol_dwi[:,:,:,bvals_indx]
vol_MPPCA     = vol_MPPCA[:,:,:,bvals_indx]
vol_tMPPCA    = vol_tMPPCA[:,:,:,bvals_indx]
   

# calculate residuals and SNR map
res_MPPCA   = vol_MPPCA-vol_dwi
res_tMPPCA  = vol_tMPPCA-vol_dwi

sigma_MPPCA      = np.std(res_MPPCA,3)
SNR_MPPCA        = np.zeros(vol_MPPCA.shape)
for dim in range(vol_MPPCA.shape[-1]):
    SNR_MPPCA[:,:,:,dim]  = vol_MPPCA[:,:,:,dim] / sigma_MPPCA

sigma_tMPPCA      = np.std(res_tMPPCA,3)
SNR_tMPPCA        = np.zeros(res_tMPPCA.shape)
for dim in range(vol_tMPPCA.shape[-1]):
    SNR_tMPPCA[:,:,:,dim]  = vol_tMPPCA[:,:,:,dim] / sigma_tMPPCA
    
 
########################## PLOT SNR ##########################       

# multiply SNR by mask
SNR_MPPCA_mask  = np.zeros_like(SNR_MPPCA)
SNR_tMPPCA_mask = np.zeros_like(SNR_tMPPCA)
for v in range(vol_dwi.shape[-1]):
    SNR_MPPCA_mask[:, :, :, v]  = SNR_MPPCA[:, :, :, v] * mask
    SNR_tMPPCA_mask[:, :, :, v] = SNR_tMPPCA[:, :, :, v] * mask
    
# Reshape and compute mean SNR for MPPCA
data_SNR_MPPCA = SNR_MPPCA_mask.reshape(-1, SNR_MPPCA.shape[3])
data_SNR_MPPCA[data_SNR_MPPCA == 0] = np.nan
SNR_MPPCA_mean = np.nanmean(data_SNR_MPPCA[:, :bvals_pershell[0]*3])

# Reshape and compute mean SNR for tMPPCA
data_SNR_tMPPCA = SNR_tMPPCA_mask.reshape(-1, SNR_tMPPCA.shape[3])
data_SNR_tMPPCA[data_SNR_tMPPCA == 0] = np.nan
SNR_tMPPCA_mean = np.nanmean(data_SNR_tMPPCA[:, :bvals_pershell[0]*3])

# Plot
fig, ax = plt.subplots()
ax.plot(bvals_ordered, np.nanmean(data_SNR_MPPCA, axis=0), 'bo', markersize=3,label='MPPCA')
ax.plot(bvals_ordered, np.nanmean(data_SNR_tMPPCA, axis=0), 'ro', markersize=3,label='t-MPPCA')
ax.legend()
ax.set_xlabel('Nominal b-val',
              fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})
ax.set_ylabel('Mean SNR', fontdict={
              'size': 12, 'weight': 'bold', 'style': 'italic'})
ax.grid(True)
plt.savefig(os.path.join(save_path, 'Denoising_MPPCAvstMPPCA_SNR.png'),
            bbox_inches='tight', dpi=300)
    

########################## PLOT IMAGES DENOISING ##########################       

vols_to_plot  = [(bvals_pershell[0]+bvals_pershell[1])*3,
                 (bvals_pershell[0]+bvals_pershell[1]+bvals_pershell[2])*3,
                 (bvals_pershell[0]+bvals_pershell[1]+bvals_pershell[2]+bvals_pershell[3])*3,
                 (bvals_pershell[0]+bvals_pershell[1]+bvals_pershell[2]+bvals_pershell[3]+bvals_pershell[4])*3]

slice_to_plot = vol_MPPCA.shape[1] // 2  

fig, axes = plt.subplots(len(vols_to_plot), 5, figsize=(10, 10))
fig.subplots_adjust(wspace=0.05,hspace=0.02, top=0.95, bottom=0.05, left=0.05, right=0.95)  
k = 0
for vol_idx in vols_to_plot:
    title = str(int(bvals_ordered[vol_idx]))
    
    vol_dwi_s = imutils.rotate(vol_dwi[:, slice_to_plot,:, vol_idx], angle=90)
    vol_MPPCA_s = imutils.rotate(vol_MPPCA[:, slice_to_plot,:, vol_idx], angle=90)
    vol_tMPPCA_s = imutils.rotate(vol_tMPPCA[:, slice_to_plot,:, vol_idx], angle=90)
    res_MPPCA_s = imutils.rotate(res_MPPCA[:, slice_to_plot,:, vol_idx], angle=90)
    res_tMPPCA_s = imutils.rotate(res_tMPPCA[:, slice_to_plot,:, vol_idx], angle=90)

    axes[k,0].imshow(vol_dwi_s, cmap='gray',vmin=0, vmax=1e4)
    axes[k,0].set_ylabel(f'b-val {title}, vol {vol_idx}')
    axes[k,0].set_xticks([])
    axes[k,0].set_yticks([]) 
    
    axes[k,1].imshow(vol_MPPCA_s, cmap='gray',vmin=0, vmax=1e4)
    axes[k,1].axis('off') 
    
    axes[k,2].imshow(res_MPPCA_s, cmap='gray',vmin=-3e3, vmax=3e3)
    axes[k,2].axis('off') 
    
    axes[k,3].imshow(vol_tMPPCA_s, cmap='gray',vmin=0, vmax=1e4)
    axes[k,3].axis('off') 
    
    axes[k,4].imshow(res_tMPPCA_s, cmap='gray',vmin=-3e3, vmax=3e3)
    axes[k,4].axis('off') 
    
    if k ==0:
        axes[k,0].set_title('Original')
        axes[k,1].set_title('MP-PCA')
        axes[k,2].set_title('MP-PCA Res')
        axes[k,3].set_title('t-MP-PCA')
        axes[k,4].set_title('t-MP-PCA Res')

    k += 1
    
plt.savefig(os.path.join(save_path,'Denoising_MPPCAvstMPPCA.png'))

########################## PLOT RESIDUALS ##########################       


fig, axes = plt.subplots(2, 2, figsize=(10, 10))

# MP-PCA residuals
data = res_MPPCA.reshape(-1, 1)
data = data[~(np.isnan(data).any(axis=1) | (data == 0).any(axis=1))]
data = np.squeeze(data)
stat, p1 = shapiro(data)
axes[0, 0].hist(data, bins=200)
axes[0, 0].set_title(f'MP-PCA Residual Histogram\nShapiro-Wilk p={p1:.3e}')
sm.qqplot(data, ax=axes[1, 0], line='45', fit=True)
axes[1, 0].set_title('MP-PCA Residual QQ Plot')

# tMP-PCA residuals
data = res_tMPPCA.reshape(-1, 1)
data = data[~(np.isnan(data).any(axis=1) | (data == 0).any(axis=1))]
data = np.squeeze(data)
stat, p2 = shapiro(data)
axes[0, 1].hist(data, bins=200)
axes[0, 1].set_title(f'tMP-PCA Residual Histogram\nShapiro-Wilk p={p2:.3e}')
sm.qqplot(data, ax=axes[1, 1], line='45', fit=True)
axes[1, 1].set_title('tMP-PCA Residual QQ Plot')

plt.tight_layout()
fig.suptitle('Residuals distribution') # or plt.suptitle('Main title')
plt.savefig(os.path.join(save_path,'Denoising_MPPCAvstMPPCA_residuals.png'))


########################## COMPUTE RATIOS ##########################   

# # Load images
# vol_dwi     = nib.load(bids_strc_MPPCA.get_path('dwi.nii.gz')).get_fdata()
# vol_MPPCA   = nib.load(bids_strc_MPPCA.get_path('dwi_dn.nii.gz')).get_fdata()
# vol_tMPPCA  = nib.load(bids_strc_tMPPCA.get_path('dwi_dn.nii.gz')).get_fdata()
# mask_dil    = nib.load(bids_strc_MPPCA.get_path('mask_dil.nii.gz')).get_fdata()

# Multiply by mask
inverted_mask = 1 - mask
for i in range(0,vol_dwi.shape[3]):
    vol_dwi[:,:,:,i] =  vol_dwi[:,:,:,i]*mask
    vol_MPPCA[:,:,:,i] =  vol_MPPCA[:,:,:,i]*mask
    vol_tMPPCA[:,:,:,i] =  vol_tMPPCA[:,:,:,i]*mask


## ratio
bvals_ordered = np.array(bvals_ordered)  # Ensure it's a NumPy array
indices_min = np.where(bvals_ordered == min(bvals_ordered[bvals_ordered>0]))[0]
indices_max = np.where(bvals_ordered == max(bvals_ordered[bvals_ordered>0]))[0]

sigma_dwi = np.mean(np.std(vol_dwi[:,:,:,indices_min].reshape(vol_dwi.shape[0]*vol_dwi.shape[1]*vol_dwi.shape[2],len(indices_min)),axis=0)) / np.mean(np.std(vol_dwi[:,:,:,indices_max].reshape(vol_dwi.shape[0]*vol_dwi.shape[1]*vol_dwi.shape[2],len(indices_max)),axis=0))
sigma_MPPCA = np.mean(np.std(vol_MPPCA[:,:,:,indices_min].reshape(vol_MPPCA.shape[0]*vol_MPPCA.shape[1]*vol_MPPCA.shape[2],len(indices_min)),axis=0)) / np.mean(np.std(vol_MPPCA[:,:,:,indices_max].reshape(vol_MPPCA.shape[0]*vol_MPPCA.shape[1]*vol_MPPCA.shape[2],len(indices_max)),axis=0))
sigma_tMPPCA = np.mean(np.std(vol_tMPPCA[:,:,:,indices_min].reshape(vol_tMPPCA.shape[0]*vol_tMPPCA.shape[1]*vol_tMPPCA.shape[2],len(indices_min)),axis=0)) / np.mean(np.std(vol_tMPPCA[:,:,:,indices_max].reshape(vol_tMPPCA.shape[0]*vol_tMPPCA.shape[1]*vol_tMPPCA.shape[2],len(indices_max)),axis=0))




########################## MODEL ESTIMATES ##########################   

# Get values
k = 0
patterns, lims = get_param_names_model('Nexi')
Data = np.zeros((2,len(cfg['ROIs']), len(patterns)))  

# Loop through standard vs random sessions
for analysis_folder in ['analysis_MPPCA','analysis_tMPPCA']:
    
    bids_strc = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path, description= 'Nexi',
                                    folderlevel='derivatives', workingdir=analysis_folder)
    output_path = os.path.join(bids_strc.get_path(), 'Masked')

    bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas']+'_To_'+'allDelta-allb', root=data_path, 
                                            folderlevel='derivatives', workingdir=analysis_folder)
    bids_strc_reg.set_param(base_name='')
               
    # Get atlas in dwi space and atlas labels
    atlas = bids_strc_reg.get_path('atlas_in_dwi.nii.gz')
    atlas_labels = pd.read_csv(
        glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*atlas.label'))[0],
        sep=r'\s+', skiprows=14, header=None, names=['IDX', 'R', 'G', 'B', 'A', 'VIS', 'MSH', 'LABEL'], quotechar='"',)  
    bids_strc_reg_TPM  = create_bids_structure(subj=subj, sess=1, datatype='registration', description=cfg['atlas_TPM']+'_To_'+'allDelta-allb', root=data_path , 
                                 folderlevel='derivatives', workingdir=analysis_folder)
    bids_strc_reg_TPM.set_param(base_name='')
    TPMs = [bids_strc_reg_TPM.get_path('atlas_TPM_GM_in_dwi.nii.gz'), 
           bids_strc_reg_TPM.get_path('atlas_TPM_WM_in_dwi.nii.gz'),
           bids_strc_reg_TPM.get_path('atlas_TPM_CSF_in_dwi.nii.gz')]
    

    # Loop through ROIs
    ROI_ctr = 0
    for ROI in cfg['ROIs']:
        
        mask_indexes = create_ROI_mask(atlas, atlas_labels, TPMs, ROI, bids_strc_reg)
        
        # Loop through each parameter outputed in the model
        pattern_ctr = 0
        for pattern in patterns:
            
            # Load parameter data
            parameter_data = nib.load(glob.glob(os.path.join(output_path, pattern))[0]).get_fdata()
            
            # Mask the parameter with the ROI
            param_masked = parameter_data * 1 * mask_indexes
            data = param_masked.reshape(param_masked.shape[0]*param_masked.shape[1]*param_masked.shape[2], 1)
            data = data[~(np.isnan(data).any(axis=1) | (data == 0).any(axis=1))]
            
            # Store the mean of the masked data
            Data[k,ROI_ctr, pattern_ctr] = np.nanmean(data) if len(data) > 0 else np.nan
    
            pattern_ctr += 1
        
        ROI_ctr += 1
   
    k += 1

# Plot
num_patterns = Data.shape[2]  
num_ROIs = Data.shape[1]  
bar_width = 0.4  
fig, axes = plt.subplots(num_patterns, 1, figsize=(10, 4 * num_patterns), sharex=True)

for param_idx in range(num_patterns):  
    ax = axes[param_idx]  

    x = np.arange(num_ROIs)  

    # Get data for Standard and Random 
    standard_values = Data[0, :, param_idx]  
    random_values = Data[1, :, param_idx]

    ax.bar(x - bar_width/2, standard_values, width=bar_width, color='blue', label='MP-PCA' if param_idx == 0 else "")
    ax.bar(x + bar_width/2, random_values, width=bar_width, color='red', label='t-MP-PCA' if param_idx == 0 else "")

    # Finish
    paramname = patterns[param_idx][1:-1]
    ax.set_ylabel(f'{paramname}', fontsize=12, fontweight='bold')
    if param_idx == 1:
        ax.set_ylim(1, 4)
    elif param_idx == 2:
        ax.set_ylim(0.5, 1.5)
    elif param_idx == 3:
        ax.set_ylim(0, 1)

    ax.grid(True)

# Add ROI labels to the x-axis
axes[-1].set_xticks(x)
axes[-1].set_xticklabels(cfg['ROIs'], rotation=45, ha='right')

# Add legend
axes[0].legend()
plt.subplots_adjust(wspace=0.05,hspace=0.05, top=0.95, bottom=0.1, left=0.1, right=0.90)  
plt.suptitle('Parameter Values Across ROIs', fontsize=14, fontweight='bold')
plt.savefig(os.path.join(save_path,'Denoising_MPPCAvstMPPCA_NexiEstimates.png'))


# ########################## PLOT DECAYS - delete ##########################       
# cfg                         = {}
# cfg['atlas']                = 'Atlas_WHS_v4'
# cfg['atlas_TPM']            = 'TPM_C57Bl6'
# cfg['common_folder']        = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','common')
# ROI = 'CC'

# # Atlas
# bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description='Atlas_WHS_v4_To_allDelta-allb', root=data_path, 
#                              folderlevel='derivatives', workingdir='analysis_MPPCA')
# bids_strc_reg.set_param(base_name='')
# atlas_labels = pd.read_csv(glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*atlas.label'))[0], sep=r'\s+', skiprows=14, header=None,
#                                            names=['IDX', 'R', 'G', 'B', 'A', 'VIS', 'MSH', 'LABEL'], quotechar='"')
# bids_strc_reg_TPM  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas_TPM']+'_To_'+'allDelta-allb', root=data_path , 
#                                             folderlevel='derivatives', workingdir='analysis_MPPCA')
# bids_strc_reg_TPM.set_param(base_name='')
# TPMs = [bids_strc_reg_TPM.get_path(f'atlas_TPM_{tissue}_in_dwi.nii.gz') for tissue in ['GM', 'WM', 'CSF']]

# mask_indexes = create_ROI_mask(atlas, atlas_labels, TPMs, ROI, bids_strc_reg)


# data_dwi = vol_dwi[mask_indexes > 0,:]  # Select only voxels inside the ROI
# data_MPPCA = vol_MPPCA[mask_indexes > 0,:]  # Select only voxels inside the ROI
# data_tMPPCA = vol_tMPPCA[mask_indexes > 0,:]  # Select only voxels inside the ROI

    
# # Violin plot
# plt.close('all')
# selected_voxels = range(11, 12)  # Example voxel range
# data_records = []

# for idx in selected_voxels:
#     for vol, label in zip([data_dwi, data_MPPCA, data_tMPPCA], ['original', 'MPPCA', 'tMPPCA']):
#         signals = vol[idx, :]
#         for bval, signal in zip(bvals, signals):
#             data_records.append({
#                 'bval': int(bval),
#                 'signal': signal,
#                 'method': label
#             })

# # Convert to DataFrame for easier plotting
# df = pd.DataFrame(data_records)

# # Plot violin plot
# plt.figure(figsize=(12, 5))
# ax = plt.gca()
# for i, method in enumerate(['original', 'MPPCA', 'tMPPCA']):
#     subset = df[df['method'] == method]
#     parts = ax.violinplot(
#         [subset[subset['bval'] == b]['signal'] for b in sorted(set(bvals))],
#         positions=np.arange(len(set(bvals))) + i * 0.25 - 0.25,  # Offset violins
#         widths=0.2,
#         showmeans=False,
#         showextrema=False,
#         showmedians=False
#     )
#     for pc in parts['bodies']:
#         pc.set_facecolor(['black', 'red', 'blue'][i])
#         pc.set_edgecolor('black')
#         pc.set_alpha(0.6)

# ax.set_xticks(np.arange(len(set(bvals))))
# ax.set_xticklabels(sorted(set(bvals)))
# ax.set_xlabel('b-value')
# ax.set_ylabel('Signal')
# ax.set_title('Signal Decay Violin Plot per b-value')
# ax.legend(['original', 'MPPCA', 'tMPPCA'])
# plt.grid(True)
# plt.tight_layout()
# plt.show()

# from matplotlib.patches import Patch

# # Define legend handles manually
# legend_handles = [
#     Patch(facecolor='black', edgecolor='black', alpha=0.6, label='original'),
#     Patch(facecolor='red', edgecolor='black', alpha=0.6, label='MPPCA'),
#     Patch(facecolor='blue', edgecolor='black', alpha=0.6,  label='tMPPCA')
# ]

# ax.legend(handles=legend_handles, loc='upper right')


#### DELETE ####  delete

# Do pwd_avg
# dwi_MPPCA       = bids_strc_MPPCA.get_path('dwi.nii.gz')
# dwi_dn_MPPCA    = bids_strc_MPPCA.get_path('dwi_dn.nii.gz')
# dwi_dn_tMPPCA   = bids_strc_tMPPCA.get_path('dwi_dn.nii.gz')
# bvalsNom        = bids_strc_MPPCA.get_path('bvalsNom.txt')
# bvalsEff        = bids_strc_MPPCA.get_path('bvalsEff.txt')


# bids_strc_MPPCA.set_param(workingdir='analysis_MPPCA',description='pwd_dwi')
# create_directory(bids_strc_MPPCA.get_path())
# calculate_pwd_avg(dwi_MPPCA, bvalsNom, bvalsEff, bids_strc_MPPCA.get_path(), np.nan)
# file_orig  =  find_files_with_pattern(bids_strc_MPPCA,'pwd_avg_norm.nii.gz')[0]


# bids_strc_MPPCA.set_param(workingdir='analysis_MPPCA',description='pwd_dwi_dn')
# create_directory(bids_strc_MPPCA.get_path())
# calculate_pwd_avg(dwi_dn_MPPCA, bvalsNom, bvalsEff, bids_strc_MPPCA.get_path(), np.nan)
# fileMPPCA  =  find_files_with_pattern(bids_strc_MPPCA,'pwd_avg_norm.nii.gz')[0]

# bids_strc_tMPPCA.set_param(workingdir='analysis_tMPPCA',description='pwd_dwi_dn')
# create_directory(bids_strc_tMPPCA.get_path())
# calculate_pwd_avg(dwi_dn_tMPPCA, bvalsNom, bvalsEff, bids_strc_tMPPCA.get_path(), np.nan)
# filetMPPCA  =  find_files_with_pattern(bids_strc_tMPPCA,'pwd_avg_norm.nii.gz')[0]
#
# Get data
# vol_dwi  = nib.load(file_orig).get_fdata()
# vol_MPPCA  = nib.load(fileMPPCA).get_fdata()
# vol_tMPPCA  = nib.load(filetMPPCA).get_fdata()
# bvals_pwd = read_numeric_txt(find_files_with_pattern(bids_strc_tMPPCA,'bvalsNom_avg.txt')[0])[0]
