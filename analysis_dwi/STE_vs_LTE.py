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
    bids_dwi      = create_bids_structure(subj=subj, sess=1, datatype='dwi', root=cfg['data_path'] , 
                 folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='Nexi')
    extract_vols(os.path.join(bids_dwi.get_path(),'powderaverage_dwi.nii.gz'), bids_dwi.get_path('b0.nii.gz'), 0, 1)

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
    
    bvals_dwi = read_numeric_txt(os.path.join(bids_dwi.get_path(),'powderaverage.bval'))
    S_S0_dwi  = nib.load(os.path.join(bids_dwi.get_path(),'powderaverage_dwi.nii.gz')).get_fdata()    
 
    # List of acquisition sessions
    sess_list    = [x for x in list(subj_data['blockNo'].unique()) if not math.isnan(x)] # clean NaNs
    
    fig, axs = plt.subplots(1, len( cfg['ROIs_GM']), figsize=(12, 4))  
    fig.subplots_adjust(wspace=0.05,hspace=0.02, top=0.90, bottom=0.14, left=0.09, right=0.95)  
    for sess in sess_list:

            # Create BIDS structures
            bids_DOR_temp = create_bids_structure(subj=subj, sess=sess, datatype='dwi_DOR', root=cfg['data_path'] , 
                         folderlevel='derivatives', workingdir=cfg['prep_foldername'],description='fwd')
            bids_DOR      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_DOR', root=cfg['data_path'] , 
                         folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
          
            
            # Create pwd average of DOR even though is not neeeded, but like this its on the analysis folder (maybe just do a copy later)
            out_path = os.path.join(bids_DOR.get_path(),'pwd_avg')
            create_directory(out_path)
            calculate_pwd_avg(bids_DOR_temp.get_path('dwi.nii.gz'),
                              bids_DOR_temp.get_path('bvalsNom.txt'),
                              bids_DOR_temp.get_path('bvalsEff.txt'),
                              out_path,
                              10)
            bids_DOR.set_param(description='pwd_avg')
            
            extract_vols(find_files_with_pattern(bids_DOR,'pwd_avg_norm.nii.gz')[0], bids_DOR.get_path('b0.nii.gz'), 0, 1)

            # REGISTRATION 
            # Register dor --> dwi
            # antsreg_simple(bids_dwi.get_path('b0.nii.gz'),  # fixed
            #         bids_DOR.get_path('b0.nii.gz'),# moving
            #         os.path.join(out_path,'dor2dwi'))
    
            # # Apply inverse transform to put T2w in dwi space
            # ants_apply_transforms_simple([bids_DOR.get_path('b0.nii.gz')],  # input
            #                     bids_dwi.get_path('b0.nii.gz'),# reference
            #                     [bids_DOR.get_path('b0.nii.gz').replace('.nii.gz','_in_dwi.nii.gz')],  # output
            #                     os.path.join(out_path,'dor2dwi0GenericAffine.mat'))  # transform 1

           
            # load data
            bvals_DOR = read_numeric_txt(find_files_with_pattern(bids_DOR,'bvalsNom')[0])
            S_S0_DOR  = nib.load(find_files_with_pattern(bids_DOR,'pwd_avg_norm.nii.gz')[0]).get_fdata()
            
           
            ### PLOT Signal FOR EACH ROI ###
            k=0
            for ROI in cfg['ROIs_GM']:
 
                mask_indexes = create_ROI_mask(atlas, atlas_labels, ROI, bids_strc_reg)
                
                # LTE signal
                if sess == 1:
                    
                    # Plot data
                    S_S0_masked = copy.deepcopy(S_S0_dwi)
                    for v in range(S_S0_masked.shape[-1]):
                        S_S0_masked[:, :, :, v] = np.multiply(S_S0_masked[:, :, :, v], mask_indexes)
                    DWI = S_S0_masked.reshape(S_S0_masked.shape[0]*S_S0_masked.shape[1]*S_S0_masked.shape[2], S_S0_masked.shape[3])
                    DWI = DWI[~(np.isnan(DWI).any(axis=1) | (DWI == 0).any(axis=1))]
                    axs[k].plot(np.transpose(bvals_dwi), np.log(np.nanmean(DWI, axis=0)), 'ko', markersize=3)
                    #axs[k].plot(np.transpose(bvals_dwi), np.nanmean(DWI, axis=0), 'ko', markersize=3)

                    # Fit LTE 
                    mask = bvals_dwi < 2.5
                    filtered_bvals = bvals_dwi[mask]
                    filtered_DWI = np.nanmean(DWI, axis=0)
                    filtered_DWI = np.expand_dims(filtered_DWI, axis=0)  
                    filtered_DWI = filtered_DWI[mask]

                    #model = np.poly1d(np.polyfit(np.transpose(bvals_dwi)[:,0], np.log(np.array(np.nanmean(DWI, axis=0))), 1))
                    m, b = np.polyfit(np.transpose(filtered_bvals), np.log(filtered_DWI), 1)
                    model = np.poly1d(np.polyfit(np.transpose(filtered_bvals), np.log(filtered_DWI), 1))
                    b_fit = np.linspace(0, 3, 100)  # Smooth b-values 
                    # signal_fit = linear_model(b_fit, *popt)
                    #axs[k].plot(b_fit, signal_fit, linestyle="--", color="black")
                    #axs[k].plot(b_fit, model(b_fit), linestyle="--", color="black",label='_nolegend_')
                    axs[k].plot(b_fit, m*b_fit+b, linestyle="--", color="black",label='_nolegend_')

      
                # STE signal
                S_S0_masked = copy.deepcopy(S_S0_DOR)
                for v in range(S_S0_masked.shape[-1]):
                    S_S0_masked[:, :, :, v] = np.multiply(S_S0_masked[:, :, :, v], mask_indexes)
                DOR = S_S0_masked.reshape(S_S0_masked.shape[0]*S_S0_masked.shape[1]*S_S0_masked.shape[2], S_S0_masked.shape[3])
                DOR = DOR[~(np.isnan(DOR).any(axis=1) | (DOR == 0).any(axis=1))]
                
                # Fit STE 
                popt, pcov = curve_fit(linear_model, np.transpose(bvals_DOR)[:,0], np.log(np.array(np.nanmean(DOR, axis=0)))) 
                #popt, pcov = curve_fit(linear_model, np.transpose(bvals_DOR)[:,0], (np.array(np.nanmean(DOR, axis=0)))) 
                b_fit = np.linspace(0, 3, 100)  # Smooth b-values 
                signal_fit = linear_model(b_fit, *popt)
                
                if sess == 1:
                    axs[k].plot(np.transpose(bvals_DOR), np.log(np.nanmean(DOR, axis=0)), 'bo', markersize=3)
                    #axs[k].plot(np.transpose(bvals_DOR), (np.nanmean(DOR, axis=0)), 'bo', markersize=3)
                    axs[k].plot(b_fit, signal_fit, linestyle="--", color="blue",label='_nolegend_')

                if sess == 2:
                    axs[k].plot(np.transpose(bvals_DOR), np.log(np.nanmean(DOR, axis=0)), 'ro', markersize=3)
                    #axs[k].plot(np.transpose(bvals_DOR), (np.nanmean(DOR, axis=0)), 'ro', markersize=3)
                    axs[k].plot(b_fit, signal_fit, linestyle="--", color="red",label='_nolegend_')

                
                # Set axes
                if k==0:
                    axs[k].set_ylabel(r'ln($S / S_0$)', fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})  
                    axs[k].legend(['LTE','STE - 16 avg','STE - 64 avg'], loc='upper right',prop={'size': 6})
                axs[k].set_xlabel('b-val', fontdict={'size': 12, 'style': 'italic'})
                axs[k].grid(True)
                axs[k].set_title(ROI)
                axs[k].set_xlim([0, 2.5])
                axs[k].set_ylim([-1.5, 0])
                if k != 0:  # Remove y-axis labels for all but the first subplot
                    axs[k].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
                axs[k].set_xticks([ 1, 2])
                
               
                #axs[k].set_yscale("log")
                
              
                # Fit STE 
                #popt, pcov = curve_fit(linear_model, np.transpose(bvals_DOR), np.nanmean(DOR, axis=0))
                
                k += 1
            
            
    plt.savefig(os.path.join(output_folder,'STE_vs_LTE.png'))
    #plt.close(fig)


                

       