#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 14:22:29 2025

@author: localadmin
"""
import os
import sys
import matplotlib.pyplot as plt
import glob
import pandas as pd

plt.close('all');
os.system('clear')
os.system('cls')

############################## ADD CODE PATH ##############################
sys.path.append(os.path.join(os.path.expanduser('~'),  'Documents', 'Rita','Codes_GitHub','dMRSI','processing_dwi'))
sys.path.append(os.path.join(os.path.expanduser('~'),  'Documents', 'Rita','Codes_GitHub','dMRSI'))


from bids_structure import *
from custom_functions import *
import importlib

importlib.reload(sys.modules['custom_functions'])
importlib.reload(sys.modules['bids_structure'])

########################## DATA PATH AND SUBJECTS ##########################
cfg                         = {}
cfg['data_path']            = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','dMRI_Pilot_20250207')
cfg['prep_foldername']      = 'preprocessed'
cfg['analysis_foldername']  = 'analysis'
cfg['common_folder']        = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','common')
cfg['scan_list_name']       = 'ScanList.xlsx'
cfg['atlas']                = 'Atlas_WHS_v4'
cfg['ROIs_GM']              = ['hippocampus','M1','M2','S1','S2', 'V1', 'PL','CG', 'Thal', 'WB']
cfg['ROIs_WM']              = ['CC']

save_path                   = '/home/localadmin/Documents/Rita/Data/dMRI_Pilot_20250207/results'

################# HELPER FUNCTION ####################
def process_bids_data(bids_strc):
    outpath = bids_strc.get_path()
    
    bvals   = read_numeric_txt(os.path.join(outpath,'powderaverage.bval'))
    data    = nib.load(os.path.join(outpath,'powderaverage_dwi.nii.gz')).get_fdata()
    if not os.path.exists(os.path.join(outpath,'b0.nii.gz')):
        extract_vols(os.path.join(outpath,'powderaverage_dwi.nii.gz'), os.path.join(outpath,'b0.nii.gz'), 0, 1)
    b0      = nib.load(os.path.join(outpath,'b0.nii.gz')).get_fdata()

    mask      = nib.load(os.path.join(outpath,'Masked','updated_mask_masked.nii.gz')).get_fdata()
    
    # Apply mask to data
    for v in range(data.shape[-1]):
        data[..., v] *= mask
    
    return bvals, data, b0

################# PLOT Powder averaged signal decay ####################

# Define BIDS structures
bids_strc_standard = create_bids_structure(subj='sub-01', sess=1, datatype='dwi', root=cfg['data_path'] , 
             folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='Nexi')
bids_strc_random = create_bids_structure(subj='sub-01', sess=2, datatype='dwi', root=cfg['data_path'] , 
             folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='Nexi')

bvals_standard, data_standard, b0_standard = process_bids_data(bids_strc_standard)
bvals_random, data_random, b0_random = process_bids_data(bids_strc_random)

# Reshape and remove zeros (mask)
data_standard = data_standard.reshape(-1, data_standard.shape[3])
data_standard[data_standard == 0] = np.nan

data_random = data_random.reshape(-1, data_random.shape[3])
data_random[data_random == 0] = np.nan

# Plot
fig, ax = plt.subplots()
ax.plot(np.transpose(bvals_standard) * 1e3, np.nanmean(data_standard, axis=0), 'ro', markersize=3)
ax.plot(np.transpose(bvals_random) * 1e3, np.nanmean(data_random, axis=0), 'bo', markersize=3)
ax.set_xlabel('Nominal b-val', fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})
ax.set_ylabel('S/S0', fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})
ax.set_xlim(500, 7500)
ax.set_ylim(0, 0.7)
ax.grid(True)
ax.legend(['Linear order','Random order'])
plt.title('Powder Averaged Signal Decay')
plt.savefig(os.path.join(save_path, 'StdVsRandom_Signaldecay_pwd.png'), bbox_inches='tight', dpi=300)

################# REGISTRATION BTW SESSIONS ####################
    
bids_strc_reg  = create_bids_structure(subj='sub-01', sess=1, datatype='registration', description='sess1_in_sess2', root=cfg['data_path'], 
                            folderlevel='derivatives', workingdir=cfg['analysis_foldername'])

# register sess1 --> sess2
create_directory(bids_strc_reg.get_path())
antsreg_simple(bids_strc_random.get_path('b0.nii.gz'), # fixed
        bids_strc_standard.get_path('b0.nii.gz'),  # moving
        bids_strc_reg.get_path('sess1_2_sess2'))

Nexi_sess2 = glob.glob(os.path.join(bids_strc_random.get_path(),'Masked','*.nii.gz'))
                       
for file in Nexi_sess2:
    new_name = os.path.basename(file)
    
    # apply inverse transform to put sess2 in sess1
    ants_apply_transforms([file],  # input 
                    bids_strc_standard.get_path('b0.nii.gz'), # moving
                    [os.path.join(bids_strc_reg.get_path(),new_name.replace('.nii.gz','_in_sess1.nii.gz'))], # output
                    [bids_strc_reg.get_path('sess1_2_sess20GenericAffine.mat'), 1], # transform 1
                    bids_strc_reg.get_path('sess1_2_sess21InverseWarp.nii.gz'))   # transform 2

# apply inverse transform to put sess2 in sess1
ants_apply_transforms([bids_strc_random.get_path('b0.nii.gz')],  # input 
                bids_strc_standard.get_path('b0.nii.gz'), # moving
                [os.path.join(bids_strc_reg.get_path(),'b0_in_sess1.nii.gz')], # output
                [bids_strc_reg.get_path('sess1_2_sess20GenericAffine.mat'), 1], # transform 1
                bids_strc_reg.get_path('sess1_2_sess21InverseWarp.nii.gz'))   # transform 2

################# EXTRACT VALUES ####################

bids_strc_reg  = create_bids_structure(subj='sub-01', sess=1, datatype='registration', description='allDelta-allb', root=cfg['data_path'], 
                             folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
standard_path  = '/home/localadmin/Documents/Rita/Data/dMRI_Pilot_20250207/derivatives/analysis/sub-01/ses-01/dwi/Nexi/Masked/'
random_path    = '/home/localadmin/Documents/Rita/Data/dMRI_Pilot_20250207/derivatives/analysis/sub-01/ses-01/registration/sess1_in_sess2/'
 
# Get atlas in dwi space and atlas labels
atlas = bids_strc_reg.get_path('atlas_in_dwi.nii.gz')
atlas_labels = pd.read_csv(
    glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*atlas.label'))[0],
    sep=r'\s+',
    skiprows=14,  
    header=None,  
    names=['IDX', 'R', 'G', 'B', 'A', 'VIS', 'MSH', 'LABEL'], 
    quotechar='"',)  

patterns, lims = get_param_names_model('Nexi')

# Get values
k = 0
Data = np.zeros((2,len(cfg['ROIs_GM']), len(patterns)))  

# Loop through standard vs random sessions
for path_data_estimates in [standard_path,random_path]:
    
    # Loop through ROIs
    ROI_ctr = 0
    for ROI in cfg['ROIs_GM']:
        
        mask_indexes = create_ROI_mask(atlas, atlas_labels, ROI, bids_strc_reg)
        
        # Loop through each parameter outputed in the model
        pattern_ctr = 0
        for pattern in patterns:
            
            # Load parameter data
            parameter_data = nib.load(glob.glob(os.path.join(path_data_estimates, pattern))[0]).get_fdata()
            
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

    ax.bar(x - bar_width/2, standard_values, width=bar_width, color='red', label='Linear' if param_idx == 0 else "")
    ax.bar(x + bar_width/2, random_values, width=bar_width, color='blue', label='Random' if param_idx == 0 else "")

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
axes[-1].set_xticklabels(cfg['ROIs_GM'], rotation=45, ha='right')

# Add legend
axes[0].legend()
plt.subplots_adjust(wspace=0.05,hspace=0.05, top=0.95, bottom=0.1, left=0.1, right=0.90)  
plt.suptitle('Parameter Values Across ROIs', fontsize=14, fontweight='bold')
plt.savefig(os.path.join(save_path, 'StdVsRandom_NexiEstimates.png'), bbox_inches='tight', dpi=300)

