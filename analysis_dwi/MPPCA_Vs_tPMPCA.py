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


#### PREPARE DATA ####  
plt.close('all')
data_path   = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','dMRI_dMRSI_Pilot_20250428')
save_path = os.path.join(data_path,'results')
create_directory(save_path)
bvals_pershell =  [4, 12, 16, 24, 30, 40]

bids_strc_MPPCA = create_bids_structure(subj='sub-01', sess=1, datatype='dwi', root=data_path, description= 'allDelta-allb',
                                folderlevel='derivatives', workingdir='preprocessed_MPPCA')
bids_strc_tMPPCA = create_bids_structure(subj='sub-01', sess=1, datatype='dwi', root=data_path, description= 'allDelta-allb',
                                folderlevel='derivatives', workingdir='preprocessed_tMPPCA')

# Load images
file0  =  bids_strc_MPPCA.get_path('dwi.nii.gz')
file1  =  bids_strc_MPPCA.get_path('dwi_dn.nii.gz')
file2  =  bids_strc_tMPPCA.get_path('dwi_dn.nii.gz')
mask_f =  bids_strc_MPPCA.get_path('mask.nii.gz')

vol0  = nib.load(file0).get_fdata()
vol1  = nib.load(file1).get_fdata()
vol2  = nib.load(file2).get_fdata()

mask  = nib.load(mask_f).get_fdata()

SNR1 = nib.load(bids_strc_MPPCA.get_path('dwi_snr.nii.gz')).get_fdata()
SNR2 = nib.load(bids_strc_tMPPCA.get_path('dwi_snr.nii.gz')).get_fdata()
for v in range(SNR1.shape[-1]):
    SNR1[:, :, :, v] = np.multiply(SNR1[:, :, :, v], mask)
for v in range(SNR2.shape[-1]):
    SNR2[:, :, :, v] = np.multiply(SNR2[:, :, :, v], mask)
 

# Order data by bvals
bvals_path = bids_strc_MPPCA.get_path('bvalsNom.txt')
with open(bvals_path, 'r') as file:
    data = file.read()
    bvals = [float(num) for num in data.split()]

bvals_ordered, bvals_indx = order_bvals(bvals)
bvals_unique  = np.unique(bvals_ordered)

vol0_order    = vol0[:,:,:,bvals_indx]
vol1_order    = vol1[:,:,:,bvals_indx]
vol2_order    = vol2[:,:,:,bvals_indx]
res1  = vol1_order-vol0_order
res2  = vol2_order-vol0_order
SNR1_order    = SNR1[:,:,:,bvals_indx]
SNR2_order    = SNR2[:,:,:,bvals_indx]

 


########################## PLOT SNR ##########################       

data_SNR1 = SNR1_order.reshape(SNR1.shape[0]*SNR1.shape[1]*SNR1.shape[2], SNR1.shape[3]);
data_SNR1[data_SNR1 == 0] = np.nan
SNR1_mean=np.nanmean(np.nanmean(data_SNR1[:,0:bvals_pershell[0]*3], axis=0))
data_SNR2 = SNR2_order.reshape(SNR2.shape[0]*SNR2.shape[1]*SNR2.shape[2], SNR2.shape[3]);
data_SNR2[data_SNR2 == 0] = np.nan
SNR2_mean =np.nanmean(np.nanmean(data_SNR2[:,0:bvals_pershell[0]*3], axis=0))

fig, ax = plt.subplots()
ax.plot(bvals_ordered, np.nanmean(data_SNR1, axis=0), 'bo', markersize=3,label='MPPCA')
ax.plot(bvals_ordered, np.nanmean(data_SNR2, axis=0), 'ro', markersize=3,label='t-MPPCA')
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

slice_to_plot = vol1.shape[1] // 2  

fig, axes = plt.subplots(len(vols_to_plot), 5, figsize=(10, 10))
fig.subplots_adjust(wspace=0.05,hspace=0.02, top=0.95, bottom=0.05, left=0.05, right=0.95)  
k = 0
for vol_idx in vols_to_plot:
    title = str(int(bvals_ordered[vol_idx]))
    
    vol0_s = imutils.rotate(vol0_order[:, slice_to_plot,:, vol_idx], angle=90)
    vol1_s = imutils.rotate(vol1_order[:, slice_to_plot,:, vol_idx], angle=90)
    vol2_s = imutils.rotate(vol2_order[:, slice_to_plot,:, vol_idx], angle=90)
    res1_s = imutils.rotate(res1[:, slice_to_plot,:, vol_idx], angle=90)
    res2_s = imutils.rotate(res2[:, slice_to_plot,:, vol_idx], angle=90)

    axes[k,0].imshow(vol0_s, cmap='gray',vmin=0, vmax=1e4)
    axes[k,0].set_ylabel(f'b-val {title}, vol {vol_idx}')
    axes[k,0].set_xticks([])
    axes[k,0].set_yticks([]) 
    
    axes[k,1].imshow(vol1_s, cmap='gray',vmin=0, vmax=1e4)
    axes[k,1].axis('off') 
    
    axes[k,2].imshow(res1_s, cmap='gray',vmin=-3e3, vmax=3e3)
    axes[k,2].axis('off') 
    
    axes[k,3].imshow(vol2_s, cmap='gray',vmin=0, vmax=1e4)
    axes[k,3].axis('off') 
    
    axes[k,4].imshow(res2_s, cmap='gray',vmin=-3e3, vmax=3e3)
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

# fig, axes = plt.subplots(2, 2, figsize=(10, 10))

# data = res1.reshape(res1.shape[0]*res1.shape[1]*res1.shape[2]*res1.shape[3], 1)
# data = data[~(np.isnan(data).any(axis=1) | (data == 0).any(axis=1))]
# data = np.squeeze(data)

# stat, p1 = shapiro(data)
# axes[0,0].hist(data,bins=200)
# axes[0,0].set_title('MP-PCA Res Hist')
# sm.qqplot(data, ax=axes[1,0],line='45', fit=True)
# axes[1,0].set_title('MP-PCA Res QQplot')

# data = res2.reshape(res2.shape[0]*res2.shape[1]*res2.shape[2]*res2.shape[3], 1)
# data = data[~(np.isnan(data).any(axis=1) | (data == 0).any(axis=1))]
# data = np.squeeze(data)

# stat, p2 = shapiro(data)
# axes[0,1].hist(data,bins=200)
# axes[0,1].set_title('t-MP-PCA Res Hist')
# sm.qqplot(data, ax=axes[1,1],line='45', fit=True)
# axes[1,1].set_title('t-MP-PCA Res QQplot')

# fig.suptitle('Residuals distribution') # or plt.suptitle('Main title')
# plt.savefig(os.path.join(save_path,'Denoising_MPPCAvstMPPCA2.png'))



########################## PLOT DECAYS ##########################       


# Atlas
bids_strc_reg  = create_bids_structure(subj='sub-01', sess=1, datatype='registration', description='Atlas_WHS_v4_To_allDelta-allb', root=data_path, 
                             folderlevel='derivatives', workingdir='analysis_MPPCA')
bids_strc_reg.set_param(base_name='')
atlas = bids_strc_reg.get_path('atlas_in_dwi.nii.gz')
atlas_labels = pd.read_csv(
     glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*atlas.label'))[0],
     sep=r'\s+',
     skiprows=14,  
     header=None,  
     names=['IDX', 'R', 'G', 'B', 'A', 'VIS', 'MSH', 'LABEL'], 
     quotechar='"',)  
mask_indexes = create_ROI_mask(atlas, atlas_labels, 'CC', bids_strc_reg)


vol0 = vol0[mask_indexes > 0,:]  # Select only voxels inside the ROI
vol1 = vol1[mask_indexes > 0,:]  # Select only voxels inside the ROI
vol2 = vol2[mask_indexes > 0,:]  # Select only voxels inside the ROI

    
# Violin plot
plt.close('all')
selected_voxels = range(11, 12)  # Example voxel range
data_records = []

for idx in selected_voxels:
    for method, vol, label in zip([data0, data1, data2], [vol0, vol1, vol2], ['original', 'MPPCA', 'tMPPCA']):
        signals = vol[idx, :]
        for bval, signal in zip(bvals, signals):
            data_records.append({
                'bval': int(bval),
                'signal': signal,
                'method': label
            })

# Convert to DataFrame for easier plotting
df = pd.DataFrame(data_records)

# Plot violin plot
plt.figure(figsize=(12, 5))
ax = plt.gca()
for i, method in enumerate(['original', 'MPPCA', 'tMPPCA']):
    subset = df[df['method'] == method]
    parts = ax.violinplot(
        [subset[subset['bval'] == b]['signal'] for b in sorted(set(bvals))],
        positions=np.arange(len(set(bvals))) + i * 0.25 - 0.25,  # Offset violins
        widths=0.2,
        showmeans=False,
        showextrema=False,
        showmedians=False
    )
    for pc in parts['bodies']:
        pc.set_facecolor(['black', 'red', 'blue'][i])
        pc.set_edgecolor('black')
        pc.set_alpha(0.6)

ax.set_xticks(np.arange(len(set(bvals))))
ax.set_xticklabels(sorted(set(bvals)))
ax.set_xlabel('b-value')
ax.set_ylabel('Signal')
ax.set_title('Signal Decay Violin Plot per b-value')
ax.legend(['original', 'MPPCA', 'tMPPCA'])
plt.grid(True)
plt.tight_layout()
plt.show()

from matplotlib.patches import Patch

# Define legend handles manually
legend_handles = [
    Patch(facecolor='black', edgecolor='black', alpha=0.6, label='original'),
    Patch(facecolor='red', edgecolor='black', alpha=0.6, label='MPPCA'),
    Patch(facecolor='blue', edgecolor='black', alpha=0.6,  label='tMPPCA')
]

ax.legend(handles=legend_handles, loc='upper right')

########################## COMPUTE RATIOS ##########################   
 

mask_dil  = nib.load(bids_strc_MPPCA.get_path('mask.nii.gz')).get_fdata()
inverted_mask = 1 - mask_dil
for i in range(0,vol0.shape[3]):
    vol0[:,:,:,i] =  vol0[:,:,:,i]*mask_dil
    vol1[:,:,:,i] =  vol1[:,:,:,i]*mask_dil
    vol2[:,:,:,i] =  vol2[:,:,:,i]*mask_dil

## ratio
indices_min = np.where(bvals == min(bvals[bvals>0]))[0]
indices_max = np.where(bvals == max(bvals[bvals>0]))[0]

sigma0 = np.mean(np.std(vol0[:,:,:,indices_min].reshape(vol0.shape[0]*vol0.shape[1]*vol0.shape[2],len(indices_min)),axis=0)) / np.mean(np.std(vol0[:,:,:,indices_max].reshape(vol0.shape[0]*vol0.shape[1]*vol0.shape[2],len(indices_max)),axis=0))
sigma1 = np.mean(np.std(vol1[:,:,:,indices_min].reshape(vol1.shape[0]*vol1.shape[1]*vol1.shape[2],len(indices_min)),axis=0)) / np.mean(np.std(vol1[:,:,:,indices_max].reshape(vol1.shape[0]*vol1.shape[1]*vol1.shape[2],len(indices_max)),axis=0))
sigma2 = np.mean(np.std(vol2[:,:,:,indices_min].reshape(vol2.shape[0]*vol2.shape[1]*vol2.shape[2],len(indices_min)),axis=0)) / np.mean(np.std(vol2[:,:,:,indices_max].reshape(vol2.shape[0]*vol2.shape[1]*vol2.shape[2],len(indices_max)),axis=0))



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
# file0  =  find_files_with_pattern(bids_strc_MPPCA,'pwd_avg_norm.nii.gz')[0]


# bids_strc_MPPCA.set_param(workingdir='analysis_MPPCA',description='pwd_dwi_dn')
# create_directory(bids_strc_MPPCA.get_path())
# calculate_pwd_avg(dwi_dn_MPPCA, bvalsNom, bvalsEff, bids_strc_MPPCA.get_path(), np.nan)
# file1  =  find_files_with_pattern(bids_strc_MPPCA,'pwd_avg_norm.nii.gz')[0]

# bids_strc_tMPPCA.set_param(workingdir='analysis_tMPPCA',description='pwd_dwi_dn')
# create_directory(bids_strc_tMPPCA.get_path())
# calculate_pwd_avg(dwi_dn_tMPPCA, bvalsNom, bvalsEff, bids_strc_tMPPCA.get_path(), np.nan)
# file2  =  find_files_with_pattern(bids_strc_tMPPCA,'pwd_avg_norm.nii.gz')[0]
#
# Get data
# vol0  = nib.load(file0).get_fdata()
# vol1  = nib.load(file1).get_fdata()
# vol2  = nib.load(file2).get_fdata()
# bvals_pwd = read_numeric_txt(find_files_with_pattern(bids_strc_tMPPCA,'bvalsNom_avg.txt')[0])[0]
