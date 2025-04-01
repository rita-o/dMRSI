"""
Compare two images 
Last changed Feb 2025
@author: Rita O
"""

import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
import imutils
from scipy.stats import shapiro
import statsmodels.api as sm 

# Load images
file0  = "/home/localadmin/Bureau/Rita/Data/Test_Denoise/sub-01_ses-01_Delta_15_fwd_dwi.nii"  
file1  = "/home/localadmin/Bureau/Rita/Data/Test_Denoise/sub-01_ses-01_Delta_15_fwd_dwi_dn.nii"  
file2  = "/home/localadmin/Bureau/Rita/Data/Test_Denoise/sub-01_ses-01_Delta_15_fwd_dwi_dn2.nii"  
mask_f = "/home/localadmin/Documents/Rita/Data/dMRI_Pilot_20250207/derivatives/preprocessed_designer/sub-01/ses-01/dwi/allDelta-allb_veraart/sub-01_ses-01_allDelta-allb_mask.nii.gz"  


vol0  = nib.load(file0).get_fdata()
vol1  = nib.load(file1).get_fdata()
vol2  = nib.load(file2).get_fdata()

mask  = nib.load(mask_f).get_fdata()
res1  = vol1-vol0
res2  = vol2-vol0

for i in range(0,res1.shape[3]):
    res1[:,:,:,i] =  res1[:,:,:,i]*mask
    res2[:,:,:,i] =  res2[:,:,:,i]*mask


# Chose slice and vols to plot
slice_idx = vol1.shape[1] // 2  

bvals_path = '/home/localadmin/Documents/Rita/Data/dMRI_Pilot_20250207/derivatives/preprocessed_designer/sub-01/ses-01/dwi/allDelta-allb_veraart/sub-01_ses-01_allDelta-allb_bvalsNom.txt'
with open(bvals_path, 'r') as file:
    data = file.read()
    bvals = [float(num) for num in data.split()]


vols = [vol1.shape[3] - 5, vol1.shape[3] - 10, vol1.shape[3] - 50]
vols = [70, 100, 120]

# Plot
fig, axes = plt.subplots(len(vols), 5, figsize=(10, 10))
fig.subplots_adjust(wspace=0.05,hspace=0.02, top=0.95, bottom=0.05, left=0.05, right=0.95)  
k = 0
for vol_idx in vols:
    title = str(int(bvals[vol_idx]))
    
    vol0_s = imutils.rotate(vol0[:, slice_idx,:, vol_idx], angle=90)
    vol1_s = imutils.rotate(vol1[:, slice_idx,:, vol_idx], angle=90)
    vol2_s = imutils.rotate(vol2[:, slice_idx,:, vol_idx], angle=90)
    res1_s = imutils.rotate(res1[:, slice_idx,:, vol_idx], angle=90)
    res2_s = imutils.rotate(res2[:, slice_idx,:, vol_idx], angle=90)

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
        axes[k,1].set_title('5D')
        axes[k,2].set_title('5D Res')
        axes[k,3].set_title('4D')
        axes[k,4].set_title('4D Res')

        
    k += 1
    
plt.savefig('/home/localadmin/Documents/Rita/Data/dMRI_Pilot_20250207/results/Denoising_5Dvs4D.png')

# Check gaussian
fig, axes = plt.subplots(2, 2, figsize=(10, 10))

data = res1.reshape(res1.shape[0]*res1.shape[1]*res1.shape[2]*res1.shape[3], 1)
data = data[~(np.isnan(data).any(axis=1) | (data == 0).any(axis=1))]
data = np.squeeze(data)

stat, p1 = shapiro(data)
axes[0,0].hist(data,bins=200)
axes[0,0].set_title('5D Res Hist')
sm.qqplot(data, ax=axes[1,0],line='45', fit=True)
axes[1,0].set_title('5D Res QQplot')

data = res2.reshape(res2.shape[0]*res2.shape[1]*res2.shape[2]*res2.shape[3], 1)
data = data[~(np.isnan(data).any(axis=1) | (data == 0).any(axis=1))]
data = np.squeeze(data)

stat, p2 = shapiro(data)
axes[0,1].hist(data,bins=200)
axes[0,1].set_title('4D Res Hist')
sm.qqplot(data, ax=axes[1,1],line='45', fit=True)
axes[1,1].set_title('4D Res QQplot')

fig.suptitle('Residuals distribution') # or plt.suptitle('Main title')
plt.savefig('/home/localadmin/Documents/Rita/Data/dMRI_Pilot_20250207/results/Denoising_5Dvs4D_2.png')
