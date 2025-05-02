"""
Compare LTE vs STE signals
Last changed Feb 2025
@author: Rita O
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
from scipy.optimize import curve_fit

importlib.reload(sys.modules['custom_functions'])
importlib.reload(sys.modules['bids_structure'])

def linear_model(b, m, b_int):
    return m * b + b_int

########################## DATA PATH AND SUBJECTS ##########################
cfg                         = {}
cfg['data_path']            = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','dMRI_Pilot_20250207')
cfg['prep_foldername']      = 'preprocessed'
cfg['analysis_foldername']  = 'analysis'
cfg['common_folder']        = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','common')
cfg['scan_list_name']       = 'ScanList.xlsx'
cfg['atlas']                = 'Atlas_WHS_v4'
cfg['ROIs']                 = ['hippocampus','M1','M2','S1','S2', 'V1', 'PL','CG', 'Thal', 'CC']
cfg['ROIs']                 = ['CC','Thal','hippocampus','M1','CSF','cerebellum']
  
subj_list = ['sub-01']
import pandas as pd
import glob
import copy

output_folder = os.path.join(cfg['data_path'],'results')
create_directory(output_folder)

################# Compare random order b vals acquisition ####################

scan_list   = pd.read_excel(os.path.join(cfg['data_path'] , 'ScanList.xlsx'))
data_type_LTE ='Delta_38'

for subj in subj_list:
    
    print('Getting model estimates ' + subj + '...')

    # Extract data for subject
    subj_data    = scan_list[scan_list['newstudyName'] == subj].reset_index(drop=True)
    
    # List of acquisition sessions
    sess_list    = [x for x in list(subj_data['blockNo'].unique()) if not math.isnan(x)] # clean NaNs
    
    sess_list = [1] # rita
    
    for sess in sess_list:
        print('Getting model estimates on session ' + str(sess) + '...')
    
        
        ### Handle LTE data ###
        bids_LTE      = create_bids_structure(subj=subj, sess=1, datatype='dwi', root=cfg['data_path'] , 
                     folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description=f'pwd_avg_{data_type_LTE}')
        bvals_LTE = read_numeric_txt(find_files_with_pattern(bids_LTE,'bvalsNom')[0])
        S_S0_LTE  = nib.load(find_files_with_pattern(bids_LTE,'pwd_avg_norm.nii.gz')[0]).get_fdata()
        
        ### Get MD ###
        bids_temp      = create_bids_structure(subj=subj, sess=1, datatype='dwi', root=cfg['data_path'] , 
                     folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description=f'DTI_DKI_{data_type_LTE}')
        bids_temp.set_param(base_name='')
        MD = nib.load(find_files_with_pattern(bids_temp,'md_dki')[0]).get_fdata()
        FA = nib.load(find_files_with_pattern(bids_temp,'fa_dki')[0]).get_fdata()

        
        ###  Handle STE data and register them to LTE using previously computed transforms ### 
        bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                     folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='pwd_avg')
        bids_LTE_prep      = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=cfg['data_path'] , 
                      folderlevel='derivatives', workingdir=cfg['prep_foldername'],description=data_type_LTE+'_fwd')
        bids_STE_prep      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                      folderlevel='derivatives', workingdir=cfg['prep_foldername'],description='STE_fwd')
        bids_strc_reg_ste  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description='STE_To_LTE_'+data_type_LTE+'_fwd', root=cfg['data_path'], 
                                       folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
        bids_strc_reg_ste.set_param(base_name='')
        bvals_STE = read_numeric_txt(find_files_with_pattern(bids_STE,'bvalsNom')[0])[0]
        files = find_files_with_pattern(bids_STE,'norm_')
        files = [f for f in files if '_in_LTE' not in f]
        selected_files = []
        for bval in bvals_STE:
            bval_int = int(bval * 1000)  
            for f in files:
                filename = os.path.basename(f)
                if f"_{bval_int}" in filename:
                    selected_files.append(f)       
        S_S0_STE = np.zeros([S_S0_LTE.shape[0],S_S0_LTE.shape[1],S_S0_LTE.shape[2],bvals_STE.shape[0]]) 
        k=0
        for f in selected_files:
            ants_apply_transforms_simple([f],  # input
                                bids_LTE_prep.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'),# reference
                                [f.replace('.nii.gz','_in_LTE.nii.gz')],  # output
                                [bids_strc_reg_ste.get_path('STE2dwi0GenericAffine.mat'), 0])  # transform 1
            S_S0_STE[:, :, :, k]  = nib.load(f.replace('.nii.gz','_in_LTE.nii.gz')).get_fdata()
            k += 1
    
              
        ###  Get atlas in dwi space and atlas labels ### 
        bids_strc_reg  = create_bids_structure(subj=subj, sess=1, datatype='registration', description=cfg['atlas']+'_To_'+data_type_LTE+'_fwd', root=cfg['data_path'] , 
                                     folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
        bids_strc_reg.set_param(base_name='')
        atlas = bids_strc_reg.get_path('atlas_in_dwi.nii.gz')
        atlas_labels = pd.read_csv(
             glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*atlas.label'))[0],
             sep=r'\s+',
             skiprows=14,  
             header=None,  
             names=['IDX', 'R', 'G', 'B', 'A', 'VIS', 'MSH', 'LABEL'], 
             quotechar='"',)  

        
        fig, axs = plt.subplots(1, len( cfg['ROIs']), figsize=(12, 4))  
        fig.subplots_adjust(wspace=0.05,hspace=0.02, top=0.90, bottom=0.14, left=0.09, right=0.95)  
   
        
        ### PLOT Signal FOR EACH ROI ###
        k=0
        for ROI in cfg['ROIs']:
     
            print('Working on ' + ROI + '...')

            mask_indexes = create_ROI_mask(atlas, atlas_labels, ROI, bids_strc_reg)

            # Mask LTE data 
            masked_LTE = S_S0_LTE[mask_indexes > 0]  # Select only voxels inside the ROI
            masked_LTE = masked_LTE.reshape(-1, S_S0_LTE.shape[-1])  # Reshape to (n_voxels, n_bvals)
            
            # Remove bad voxels (NaN or 0 in any bval)
            good_voxels_mask = ~(np.isnan(masked_LTE).any(axis=1) | (masked_LTE == 0).any(axis=1))
            masked_LTE = masked_LTE[good_voxels_mask]
            
            # Average across good voxels
            filtered_LTE = np.nanmean(masked_LTE, axis=0)
            
            # Insert b0
            filtered_bvals = bvals_LTE
            filtered_bvals = np.insert(filtered_bvals, 0, 0)
            filtered_LTE = np.insert(filtered_LTE, 0, 1)
            
            # Define data to plot
            data = np.log(filtered_LTE)
            #data = filtered_LTE

            # Plot data
            axs[k].plot(np.transpose(filtered_bvals), data, 'ko', markersize=3)

            # Fit LTE 
            coeffs = np.polyfit(np.transpose(filtered_bvals), data, 3)
            model = np.poly1d(coeffs)
            b_fit = np.linspace(0,bvals_LTE.max(), 100)  # Smooth b-values 
            y_fit = model(b_fit)
            axs[k].plot(b_fit, y_fit, linestyle="--", color="black",label='_nolegend_')
    
            # Plot MD
            masked_MD = MD[mask_indexes > 0]
            masked_MD = masked_MD.reshape(-1, 1)  # Reshape to (n_voxels, n_bvals)
            good_voxels_mask = ~(np.isnan(masked_MD).any(axis=1) | (masked_MD == 0).any(axis=1))
            masked_MD = masked_MD[good_voxels_mask]
            masked_MD = np.nanmean(masked_MD, axis=0)
 
            ln_signal_exp = -b_fit*masked_MD
            axs[k].plot(b_fit, ln_signal_exp, linestyle="-", color="red",label='_nolegend_')

            # FA
            masked_FA = FA[mask_indexes > 0]
            masked_FA = masked_FA.reshape(-1, 1)  # Reshape to (n_voxels, n_bvals)
            good_voxels_mask = ~(np.isnan(masked_FA).any(axis=1) | (masked_FA == 0).any(axis=1))
            masked_FA = masked_FA[good_voxels_mask]
            mean_FA = np.round(np.nanmean(masked_FA, axis=0),2)[0]
            
            # Mask STE data 
            masked_STE = S_S0_STE[mask_indexes > 0]  # Select only voxels inside the ROI
            masked_STE = masked_STE.reshape(-1, S_S0_STE.shape[-1])  # Reshape to (n_voxels, n_bvals)
            
            # Remove bad voxels (NaN or 0 in any bval)
            good_voxels_mask = ~(np.isnan(masked_STE).any(axis=1) | (masked_STE == 0).any(axis=1))
            masked_STE = masked_STE[good_voxels_mask]
            
           
            # Average across good voxels
            filtered_STE = np.nanmean(masked_STE, axis=0)
            
            # Insert b0
            filtered_bvals = bvals_STE
            filtered_bvals = np.insert(filtered_bvals, 0, 0)
            filtered_STE = np.insert(filtered_STE, 0, 1)
            
            # Define data to plot
            data = np.log(filtered_STE)
            #data = filtered_STE

            # Fit STE 
            # popt, pcov = curve_fit(linear_model, np.transpose(filtered_bvals), data) 
            # b_fit = np.linspace(0, 3, 100)  # Smooth b-values 
            # signal_fit = linear_model(b_fit, *popt)
            coeffs = np.polyfit(np.transpose(filtered_bvals), data, 2)
            model = np.poly1d(coeffs)
            b_fit = np.linspace(0,bvals_STE.max(), 100)  # Smooth b-values 
            y_fit = model(b_fit)
            
            # if sess == 1:
            axs[k].plot(np.transpose(filtered_bvals), data, 'bo', markersize=3)
            #axs[k].plot(np.transpose(bvals_STE), (np.nanmean(STE, axis=0)), 'bo', markersize=3)
            axs[k].plot(b_fit, y_fit, linestyle="--", color="blue",label='_nolegend_')

            # if sess == 2:
            #     axs[k].plot(np.transpose(filtered_bvals), data, 'ro', markersize=3)
            #     #axs[k].plot(np.transpose(bvals_STE), (np.nanmean(STE, axis=0)), 'ro', markersize=3)
            #     axs[k].plot(b_fit, signal_fit, linestyle="--", color="red",label='_nolegend_')
    
            
            # Set axes
            if k==0:
                axs[k].set_ylabel(r'$ln(S / S_0)$', fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})  
                axs[k].legend(['LTE','STE'], loc='upper right',prop={'size': 6})
            axs[k].set_xlabel('b-val', fontdict={'size': 12, 'style': 'italic'})
            axs[k].grid(True)
            axs[k].set_title(f'{ROI} - FA: {mean_FA}')
            axs[k].set_xlim([-0.5, 3])
            axs[k].set_ylim([-1.8, 0.1])
            if k != 0:  # Remove y-axis labels for all but the first subplot
                axs[k].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
            #axs[k].set_xticks(np.unique(bvals_LTE))
            
         
            k += 1
 

    plt.savefig(os.path.join(output_folder,'STE_vs_LTE.png'))
    #plt.close(fig)
    
   


                

       