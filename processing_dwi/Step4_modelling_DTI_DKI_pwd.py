"""
Script to analyse WM microstructure by fitting models such as Nexi, Sandi, ...
It uses the Designer module installed in the Designer python environment with a docker.

Last changed Jan 2025
@author: Rita O
"""

import os
import sys
import pandas as pd
import platform
import math
import importlib, sys
from custom_functions import *
from bids_structure import *

def Step4_modelling_DTI_DKI_pwd(subj_list, cfg):
    
    data_path   = cfg['data_path']     
    scan_list   = pd.read_excel(os.path.join(data_path, 'ScanList.xlsx'))
    
    # Define path to docker
    docker_path = '/data'
    
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in subj_list:
        
        print('Modelling ' + subj + '...')
    
        # Extract data for subject
        subj_data      = scan_list[(scan_list['newstudyName'] == subj)].reset_index(drop=True)
          
        # Define path to docker
        docker_path = '/data'
        
        ######## SESSION-WISE OPERATIONS ########
        for sess in list(subj_data['blockNo'].unique()):
          
            filtered_data = subj_data[(subj_data['phaseDir'] == 'fwd') & (subj_data['blockNo'] == sess) & (subj_data['noBval'] > 1) & (subj_data['acqType'] == 'PGSE') & (subj_data['scanQA'] == 'ok')]
            Delta_list = filtered_data['diffTime'].unique()
            
            Delta_list = [int(min(filtered_data["diffTime"])), 
                          int(max(filtered_data["diffTime"]))] 
            
            ######## DELTA-WISE OPERATIONS ########
            for Delta in Delta_list:
               
               
                # Define bids structure for the processed data
                bids_strc_prep = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path, 
                                                folderlevel='derivatives', workingdir=cfg['prep_foldername'])
                bids_strc_prep.set_param(description='Delta_' + str(Delta) + '_fwd')

                
                ######## Run DTI and DKI ########  
                 
                # Define BIDS structure for the analysis data
                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype="dwi", description=f'DTI_DKI_Delta_{Delta}', root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
              
                # Make output folder 
                output_path = bids_strc_analysis.get_path()
                input_path = os.path.join(output_path,'inputs')
                create_directory(input_path)
                    
                # Copy necessary files for analysis and rename the path to the docker path
                dwi   = copy_files_BIDS(bids_strc_prep,input_path, 'dwi_dn_gc_ec.mif').replace(data_path,docker_path)
                mask  = copy_files_BIDS(bids_strc_prep,input_path,  'mask.nii.gz').replace(data_path,docker_path)
                out_folder   = output_path.replace(data_path, docker_path)
 
                # Run model
                call = [f'docker run -v {data_path}:/{docker_path} nyudiffusionmri/designer2:v2.0.10 tmi -DTI -DKI',
                        f'{dwi} {out_folder}'] #
         
                print(' '.join(call))
                os.system(' '.join(call))
                
                # Rename paths to local folder
                bids_strc_analysis.set_param(root=data_path)
                bids_strc_prep.set_param(root=data_path)
    
                output_path = bids_strc_analysis.get_path()
    
                # Put with the same header as original image because Designer always changes everything (rolling eyes intensively)
                for filename in os.listdir(output_path):
                    if filename.endswith(".nii"):
                        in_img = os.path.join(output_path, filename)
                        ref_img = mask.replace(docker_path,data_path)

                        call = [f'flirt',
                             f'-in  {in_img}',
                             f'-ref {ref_img}',
                             f'-out {in_img}',
                             f'-applyxfm -usesqform']
                        os.system(' '.join(call))
                        
                        os.system('rm ' + f'{in_img}')
                        os.system('gunzip ' + f'{in_img}' + '.gz')               
                        
                # Mask output with brain mask for better visualization
                for filename in os.listdir(output_path):
                    if filename.endswith(".nii"):
                        multiply_by_mask(os.path.join(output_path, filename), # filename input
                                         os.path.join(output_path,'Masked'), # output folder
                                         bids_strc_prep.get_path('mask.nii.gz')) # mask
                
                ######## Compute PWD ######## 
                
                # Create BIDS structures for LTE
                bids_LTE_temp = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=cfg['data_path'] , 
                             folderlevel='derivatives', workingdir=cfg['prep_foldername'],description=f'Delta_{Delta}_fwd')
                bids_LTE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=cfg['data_path'] , 
                             folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description=f'pwd_avg_Delta_{Delta}')
              
                
                # Create pwd average of LTE 
                create_directory(bids_LTE.get_path())
                calculate_pwd_avg(bids_LTE_temp.get_path('dwi_dn_gc_ec.nii.gz'),
                                  bids_LTE_temp.get_path('bvalsNom.txt'),
                                  bids_LTE_temp.get_path('bvalsEff.txt'),
                                  bids_LTE.get_path(),
                                  np.nan)
             
            ######## Compute PWD ######## 

            # Create BIDS structures for STE
            bids_STE_temp = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                          folderlevel='derivatives', workingdir=cfg['prep_foldername'],description='STE_fwd')
            bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                          folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='pwd_avg')
            if os.path.exists(bids_STE_temp.get_path('dwi_dn_gc_topup.nii.gz')):                
                  # Create pwd average of STE 
                  create_directory(bids_STE.get_path())
                  calculate_pwd_avg(bids_STE_temp.get_path('dwi_dn_gc_topup.nii.gz'),
                                    bids_STE_temp.get_path('bvalsNom.txt'),
                                    bids_STE_temp.get_path('bvalsEff.txt'),
                                    bids_STE.get_path(),
                                    np.nan)
                  
                ######## Compute Micro FA - wrong ######## 
 
                # # Load LTE data
                # bids_LTE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=cfg['data_path'] , 
                #              folderlevel='derivatives', workingdir=cfg['analysis_foldername'], description=f'pwd_avg_Delta_{Delta}') 
                # bvals_LTE = read_numeric_txt(find_files_with_pattern(bids_LTE,'bvalsNom')[0])
                # S_S0_LTE  = nib.load(find_files_with_pattern(bids_LTE,'pwd_avg_norm.nii.gz')[0]).get_fdata()
                
                # # Load STE data
                # bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                #              folderlevel='derivatives', workingdir=cfg['analysis_foldername'], description='pwd_avg')
                # bvals_STE = read_numeric_txt(find_files_with_pattern(bids_STE,'bvalsNom')[0])
                # S_S0_STE  = nib.load(find_files_with_pattern(bids_STE,'pwd_avg_norm.nii.gz')[0]).get_fdata()
                
                # # Find b-value that is common to both
                # common_bvalue = np.intersect1d(bvals_STE, bvals_LTE)[0] # [ms/um²]
                # vol_STE = int(np.where(bvals_STE == common_bvalue)[1][0])
                # vol_LTE = int(np.where(bvals_LTE == common_bvalue)[1][0])

                # # Get MD value
                # MD = nib.load(os.path.join(bids_strc_analysis.get_path(),'md_dki.nii')).get_fdata() # [um²/ms]
                
                # # Calculate micro FA
                # mask_data = nib.load(ref_img).get_fdata()
                # E_STE = np.squeeze(S_S0_STE[:,:,:,vol_STE]) #*mask_data
                # E_LTE = np.squeeze(S_S0_LTE[:,:,:,vol_LTE])#*mask_data
                # var_u = np.log(E_LTE/E_STE) * 2 * (common_bvalue)**(-2) * MD**(-2)
                # var_u_safe = np.copy(var_u)
                # var_u_safe[np.abs(var_u_safe) < 1e-6] = 1e-6
                # #var_u_safe_clean = np.nan_to_num(var_u_safe, nan=300)  # Replaces NaN with 1
                
                # microFA = np.sqrt(3 / 2) * (1 + (2/5)*(1/var_u_safe) )**-0.5
                # bad_FA = np.where(E_LTE-E_STE < 0) 
                
                # for i in range(0,len(bad_FA[0])):
                #     microFA[bad_FA[0][i],bad_FA[1][i],bad_FA[2][i]]=np.nan;
                
                # # Save image
                # affine = nib.load(ref_img).affine
                # img = nib.Nifti1Image(microFA, affine)
                # nib.save(img, os.path.join(bids_strc_analysis.get_path(),'microFA.nii'))

            
        
            