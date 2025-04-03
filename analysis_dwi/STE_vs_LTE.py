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
cfg['ROIs_GM']              = ['hippocampus','M1','M2','S1','S2', 'V1', 'PL','CG', 'Thal', 'CC']

subj_list = ['sub-01']
import pandas as pd
import glob
import copy

output_folder = os.path.join(cfg['data_path'],'results')
create_directory(output_folder)

################# Compare random order b vals acquisition ####################

scan_list   = pd.read_excel(os.path.join(cfg['data_path'] , 'ScanList.xlsx'))

for subj in subj_list:
    
    print('Getting model estimates ' + subj + '...')

    # Extract data for subject
    subj_data    = scan_list[scan_list['newstudyName'] == subj].reset_index(drop=True)
    
    # Handle dwi data
    bids_LTE      = create_bids_structure(subj=subj, sess=1, datatype='dwi', root=cfg['data_path'] , 
                 folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='Nexi')
    extract_vols(os.path.join(bids_LTE.get_path(),'powderaverage_dwi.nii.gz'), bids_LTE.get_path('b0.nii.gz'), 0, 1)

    # Get atlas in dwi space and atlas labels
    bids_strc_reg  = create_bids_structure(subj=subj, sess=1, datatype='registration', description='allDelta-allb', root=cfg['data_path'] , 
                                 folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
    atlas = bids_strc_reg.get_path('atlas_in_dwi.nii.gz')
    atlas_labels = pd.read_csv(
         glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*atlas.label'))[0],
         sep=r'\s+',
         skiprows=14,  
         header=None,  
         names=['IDX', 'R', 'G', 'B', 'A', 'VIS', 'MSH', 'LABEL'], 
         quotechar='"',)  
    
    bvals_LTE = read_numeric_txt(os.path.join(bids_LTE.get_path(),'powderaverage.bval'))
    S_S0_LTE  = nib.load(os.path.join(bids_LTE.get_path(),'powderaverage_dwi.nii.gz')).get_fdata()    
 
    # List of acquisition sessions
    sess_list    = [x for x in list(subj_data['blockNo'].unique()) if not math.isnan(x)] # clean NaNs
    
    fig, axs = plt.subplots(1, len( cfg['ROIs_GM']), figsize=(12, 4))  
    fig.subplots_adjust(wspace=0.05,hspace=0.02, top=0.90, bottom=0.14, left=0.09, right=0.95)  
    for sess in sess_list:

            # Create BIDS structures
            bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                         folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='pwd_avg')
            extract_vols(find_files_with_pattern(bids_STE,'pwd_avg_norm.nii.gz')[0], bids_STE.get_path('b0.nii.gz'), 0, 1)

            # REGISTRATION 
            # Register STE --> dwi
            # antsreg_simple(bids_dwi.get_path('b0.nii.gz'),  # fixed
            #         bids_STE.get_path('b0.nii.gz'),# moving
            #         os.path.join(out_path,'STE2dwi'))
    
            # # Apply inverse transform to put T2w in dwi space
            # ants_apply_transforms_simple([bids_STE.get_path('b0.nii.gz')],  # input
            #                     bids_dwi.get_path('b0.nii.gz'),# reference
            #                     [bids_STE.get_path('b0.nii.gz').replace('.nii.gz','_in_dwi.nii.gz')],  # output
            #                     os.path.join(out_path,'STE2dwi0GenericAffine.mat'))  # transform 1

           
            # load data
            bvals_STE = read_numeric_txt(find_files_with_pattern(bids_STE,'bvalsNom')[0])
            S_S0_STE  = nib.load(find_files_with_pattern(bids_STE,'pwd_avg_norm.nii.gz')[0]).get_fdata()
            
           
            ### PLOT Signal FOR EACH ROI ###
            k=0
            for ROI in cfg['ROIs_GM']:
 
                mask_indexes = create_ROI_mask(atlas, atlas_labels, ROI, bids_strc_reg)
                
                # LTE signal
                if sess == 1:
                    
                    # Get values
                    S_S0_masked = copy.deepcopy(S_S0_LTE)
                    for v in range(S_S0_masked.shape[-1]):
                        S_S0_masked[:, :, :, v] = np.multiply(S_S0_masked[:, :, :, v], mask_indexes)
                    LTE = S_S0_masked.reshape(S_S0_masked.shape[0]*S_S0_masked.shape[1]*S_S0_masked.shape[2], S_S0_masked.shape[3])
                    LTE = LTE[~(np.isnan(LTE).any(axis=1) | (LTE == 0).any(axis=1))]
                    mask = bvals_LTE < 2.5
                    filtered_bvals = bvals_LTE[mask]
                    filtered_bvals = np.insert(filtered_bvals, 0, 0)
                    filtered_LTE = np.nanmean(LTE, axis=0)
                    filtered_LTE = np.expand_dims(filtered_LTE, axis=0)  
                    filtered_LTE = filtered_LTE[mask]
                    filtered_LTE = np.insert(filtered_LTE, 0, 1)
                    
                    # Plot data
                    axs[k].plot(np.transpose(filtered_bvals), np.log(filtered_LTE), 'ko', markersize=3)

                    # Fit LTE 
                    m, b = np.polyfit(np.transpose(filtered_bvals), np.log(filtered_LTE), 1)
                    model = np.poly1d(np.polyfit(np.transpose(filtered_bvals), np.log(filtered_LTE), 1))
                    b_fit = np.linspace(0, 3, 100)  # Smooth b-values 
                    axs[k].plot(b_fit, m*b_fit+b, linestyle="--", color="black",label='_nolegend_')

      
                # STE signal
                
                # Get values
                S_S0_masked = copy.deepcopy(S_S0_STE)
                for v in range(S_S0_masked.shape[-1]):
                    S_S0_masked[:, :, :, v] = np.multiply(S_S0_masked[:, :, :, v], mask_indexes)
                STE = S_S0_masked.reshape(S_S0_masked.shape[0]*S_S0_masked.shape[1]*S_S0_masked.shape[2], S_S0_masked.shape[3])
                STE = STE[~(np.isnan(STE).any(axis=1) | (STE == 0).any(axis=1))]
                filtered_bvals = np.insert(bvals_STE, 0, 0)
                filtered_STE= np.nanmean(STE, axis=0)
                filtered_STE = np.expand_dims(filtered_STE, axis=0)  
                filtered_STE = np.insert(filtered_STE, 0, 1)
                
                
                # Fit STE 
                popt, pcov = curve_fit(linear_model, np.transpose(filtered_bvals), np.log(filtered_STE)) 
                b_fit = np.linspace(0, 3, 100)  # Smooth b-values 
                signal_fit = linear_model(b_fit, *popt)
                
                if sess == 1:
                    axs[k].plot(np.transpose(filtered_bvals), np.log(filtered_STE), 'bo', markersize=3)
                    #axs[k].plot(np.transpose(bvals_STE), (np.nanmean(STE, axis=0)), 'bo', markersize=3)
                    axs[k].plot(b_fit, signal_fit, linestyle="--", color="blue",label='_nolegend_')

                if sess == 2:
                    axs[k].plot(np.transpose(filtered_bvals), np.log(filtered_STE), 'ro', markersize=3)
                    #axs[k].plot(np.transpose(bvals_STE), (np.nanmean(STE, axis=0)), 'ro', markersize=3)
                    axs[k].plot(b_fit, signal_fit, linestyle="--", color="red",label='_nolegend_')

                
                # Set axes
                if k==0:
                    axs[k].set_ylabel(r'ln($S / S_0$)', fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})  
                    axs[k].legend(['LTE','STE - 16 avg','STE - 64 avg'], loc='upper right',prop={'size': 6})
                axs[k].set_xlabel('b-val', fontdict={'size': 12, 'style': 'italic'})
                axs[k].grid(True)
                axs[k].set_title(ROI)
                axs[k].set_xlim([-0.5, 2.6])
                axs[k].set_ylim([-1.6, 0.05])
                if k != 0:  # Remove y-axis labels for all but the first subplot
                    axs[k].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
                axs[k].set_xticks([ 1, 2])
                
             
                k += 1
        
           

    plt.savefig(os.path.join(output_folder,'STE_vs_LTE.png'))
    #plt.close(fig)
    
   


                

       