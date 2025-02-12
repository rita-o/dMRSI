#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 14:22:29 2025

@author: localadmin
"""
import os
import sys
import matplotlib.pyplot as plt

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

################# Compare random order b vals acquisition ####################

# create bids structure
bids_strc_standard = create_bids_structure(subj='sub-01', sess=1, datatype='dwi', root=cfg['data_path'] , 
             folderlevel='derivatives', workingdir=cfg['prep_foldername'],description='Delta_26_fwd')
bids_strc_random = create_bids_structure(subj='sub-01', sess=2, datatype='dwi', root=cfg['data_path'] , 
             folderlevel='derivatives', workingdir=cfg['prep_foldername'],description='Delta_26_fwd')

# calculate powder average standard
out_path = os.path.join(bids_strc_standard.get_path(),'pwd_avg')
create_directory(out_path)
calculate_pwd_avg(bids_strc_standard.get_path('dwi.nii.gz'),
                  bids_strc_standard.get_path('bvalsNom.txt'),
                  bids_strc_standard.get_path('bvalsEff.txt'),
                  out_path,
                  15)
bvals_standard = read_numeric_txt(os.path.join(out_path,bids_strc_standard.get_param('base_name')+'bvalsNom_avg.txt'))
bvals_standard_all = read_numeric_txt(bids_strc_standard.get_path('bvalsNom.txt'))
standard_all = nib.load(bids_strc_standard.get_path('dwi.nii.gz')).get_fdata()
standard = nib.load(os.path.join(out_path,bids_strc_standard.get_param('base_name')+'dwi_pwd_avg_norm.nii.gz')).get_fdata()

mask_standard = nib.load(bids_strc_standard.get_path('b0_avg_mask.nii.gz')).get_fdata()

# calculate poweder average random
out_path = os.path.join(bids_strc_random.get_path(),'pwd_avg')
create_directory(out_path)
calculate_pwd_avg(bids_strc_random.get_path('dwi.nii.gz'),
                  bids_strc_random.get_path('bvalsNom.txt'),
                  bids_strc_random.get_path('bvalsEff.txt'),
                  out_path,
                  15)
bvals_random = read_numeric_txt(os.path.join(out_path,bids_strc_random.get_param('base_name')+'bvalsNom_avg.txt'))
bvals_random_all = read_numeric_txt(bids_strc_random.get_path('bvalsNom.txt'))
random_all = nib.load(bids_strc_random.get_path('dwi.nii.gz')).get_fdata()
random = nib.load(os.path.join(out_path,bids_strc_random.get_param('base_name')+'dwi_pwd_avg_norm.nii.gz')).get_fdata()
mask_random = nib.load(bids_strc_random.get_path('b0_avg_mask.nii.gz')).get_fdata()

# mask in brain region
for v in range(standard.shape[-1]):
    standard[:, :, :, v] = np.multiply(standard[:, :, :, v], mask_standard)
    standard_all[:, :, :, v] = np.multiply(standard_all[:, :, :, v], mask_standard)

data_standard = standard.reshape(standard.shape[0]*standard.shape[1]*standard.shape[2], standard.shape[3]);
data_standard[data_standard == 0] = np.nan
data_standard_all = standard_all.reshape(standard_all.shape[0]*standard_all.shape[1]*standard_all.shape[2], standard_all.shape[3]);
data_standard_all[data_standard_all == 0] = np.nan

    
for v in range(random.shape[-1]):
    random[:, :, :, v] = np.multiply(random[:, :, :, v], mask_random)    
    random_all[:, :, :, v] = np.multiply(random_all[:, :, :, v], mask_random)    

data_random = random.reshape(random.shape[0]*random.shape[1]*random.shape[2], random.shape[3]);
data_random[data_random == 0] = np.nan
data_random_all = random_all.reshape(random_all.shape[0]*random_all.shape[1]*random_all.shape[2], random_all.shape[3]);
data_random_all[data_random_all == 0] = np.nan  

output_folder = os.path.join(cfg['data_path'],'results')
create_directory(output_folder)

# figure 1
fig, ax = plt.subplots()
ax.plot(np.transpose(bvals_random)*1e3, np.nanmean(data_random, axis=0), 'bo', markersize=3)
ax.plot(np.transpose(bvals_standard)*1e3, np.nanmean(data_standard, axis=0), 'ro', markersize=3)
ax.set_xlabel('Nominal b-val',
              fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})
ax.set_ylabel('S/S0', fontdict={
              'size': 12, 'weight': 'bold', 'style': 'italic'})
plt.xlim(500, 7500) 
plt.ylim(0, 0.7) 
ax.grid(True)
ax.legend(['Random order','Linear order'])
plt.title('Powder averaged signal decay')
plt.savefig(os.path.join(output_folder, 'Signaldecay_pwd.png'),
             bbox_inches='tight', dpi=300)
# figure 2
b0_random=np.mean(np.nanmean(data_random_all, axis=0)[0:4])
b0_standard=np.mean(np.nanmean(data_standard_all, axis=0)[0:4])

fig, ax = plt.subplots()
ax.plot(np.squeeze(np.transpose(bvals_random_all)), np.nanmean(data_random_all, axis=0)/b0_random, 'bo', markersize=3)
ax.plot(np.squeeze(np.transpose(bvals_standard_all)), np.nanmean(data_standard_all, axis=0)/b0_standard, 'ro', markersize=3)
ax.set_xlabel('Nominal b-val',
              fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})
ax.set_ylabel('S/S0', fontdict={
              'size': 12, 'weight': 'bold', 'style': 'italic'})
ax.grid(True)
ax.legend(['Random order','Linear order'])
plt.title('Signal decay for all bval')
plt.xlim(500, 7500) 
plt.ylim(0, 0.7) 
plt.savefig(os.path.join(output_folder, 'Signaldecay.png'),
             bbox_inches='tight', dpi=300)
