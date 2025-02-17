"""
Compare two images 
Last changed Feb 2025
@author: Rita O
"""

import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
import imutils

# Load images
file0 = "/home/localadmin/Documents/Rita/Data/dMRI_Pilot_20250207/derivatives/preprocessed_designer/sub-01/ses-01/dwi/allDelta-allb_veraart/sub-01_ses-01_allDelta-allb_dwi.nii.gz"  
file1 = "/home/localadmin/Documents/Rita/Data/dMRI_Pilot_20250207/derivatives/preprocessed_designer/sub-01/ses-01/dwi/allDelta-allb_veraart/sub-01_ses-01_allDelta-allb_dwi_dn.nii.gz"  
file2 = "/home/localadmin/Documents/Rita/Data/dMRI_Pilot_20250207/derivatives/preprocessed_designer/sub-01/ses-01/dwi/allDelta-allb_jespersen/sub-01_ses-01_allDelta-allb_dwi_dn.nii.gz"  



vol0  = nib.load(file0).get_fdata()
vol1  = nib.load(file1).get_fdata()
vol2  = nib.load(file2).get_fdata()

diff = np.abs(vol1 - vol2)

# Chose slice and vols to plot
slice_idx = vol1.shape[1] // 2  

bvals_path = '/home/localadmin/Documents/Rita/Data/dMRI_Pilot_20250207/derivatives/preprocessed_designer/sub-01/ses-01/dwi/allDelta-allb_veraart/sub-01_ses-01_allDelta-allb_bvalsNom.txt'
with open(bvals_path, 'r') as file:
    data = file.read()
    bvals = [float(num) for num in data.split()]


vols = [vol1.shape[3] - 5, vol1.shape[3] - 10, vol1.shape[3] - 50]
vols = [52, 109, 300, 377]

# Plot
fig, axes = plt.subplots(len(vols), 4, figsize=(10, 10))
fig.subplots_adjust(wspace=0.05,hspace=0.02, top=0.95, bottom=0.05, left=0.05, right=0.95)  
k = 0
for vol_idx in vols:
    title = str(int(bvals[vol_idx]))
    
    vol0_s = imutils.rotate(vol0[:, slice_idx,:, vol_idx], angle=90)
    vol1_s = imutils.rotate(vol1[:, slice_idx,:, vol_idx], angle=90)
    vol2_s = imutils.rotate(vol2[:, slice_idx,:, vol_idx], angle=90)
    diff_s = imutils.rotate(diff[:, slice_idx,:, vol_idx], angle=90)
    
    axes[k,0].imshow(vol0_s, cmap='gray',vmin=0, vmax=1e4)
    axes[k,0].set_ylabel(f'b-val {title}, vol {vol_idx}')
    axes[k,0].set_xticks([])
    axes[k,0].set_yticks([]) 
    
    axes[k,1].imshow(vol1_s, cmap='gray',vmin=0, vmax=1e4)
    axes[k,1].axis('off') 

    axes[k,2].imshow(vol2_s, cmap='gray',vmin=0, vmax=1e4)
    axes[k,2].axis('off') 
    
    axes[k,3].imshow(diff_s, cmap='hot',vmin=0, vmax=1e3)
    axes[k,3].axis('off') 
    
    if k ==0:
        axes[k,0].set_title('Original')
        axes[k,1].set_title('Veraart')
        axes[k,2].set_title('Jespersen')
        axes[k,3].set_title('Abs(Jespersen - Veraart)')
        
    k += 1
    
plt.show()
plt.savefig('/home/localadmin/Documents/Rita/Data/dMRI_Pilot_20250207/results/Denoising.png')
