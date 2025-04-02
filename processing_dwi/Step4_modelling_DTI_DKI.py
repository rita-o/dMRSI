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

def Step4_modelling_DTI_DKI(subj_list, cfg):
    
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
            
            Delta_list = [15, 38] # change
            
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
                mask  = copy_files_BIDS(bids_strc_prep,input_path,  'b0_avg_mask.nii.gz').replace(data_path,docker_path)
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
                                         bids_strc_prep.get_path('b0_avg_mask.nii.gz')) # mask
                
                ######## Compute Micro FA ######## 
                
                # Load LTE data
                bids_dwi      = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=cfg['data_path'] , 
                             folderlevel='derivatives', workingdir=cfg['analysis_foldername'], description='Nexi')
                bvals_dwi = read_numeric_txt(os.path.join(bids_dwi.get_path(),'powderaverage.bval'))
                S_S0_dwi  = nib.load(os.path.join(bids_dwi.get_path(),'powderaverage_dwi.nii.gz')).get_fdata()    
                
                # Load STE data
                bids_DOR      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_DOR', root=cfg['data_path'] , 
                             folderlevel='derivatives', workingdir=cfg['analysis_foldername'], description='pwd_avg')
                bvals_DOR = read_numeric_txt(find_files_with_pattern(bids_DOR,'bvalsNom')[0])
                S_S0_DOR  = nib.load(find_files_with_pattern(bids_DOR,'pwd_avg_norm.nii.gz')[0]).get_fdata()
                
                # Find b-value that is common to both
                common_bvalue = np.intersect1d(bvals_DOR, bvals_dwi)[0] # [ms/um²]
                vol_DOR = int(np.where(bvals_DOR == common_bvalue)[1][0])
                vol_DWI = int(np.where(bvals_dwi == common_bvalue)[1][0])

                # Get MD value
                MD = nib.load(os.path.join(bids_strc_analysis.get_path(),'md_dki.nii')).get_fdata() # [um²/ms]
                
                # Calculate micro FA
                E_STE = np.squeeze(S_S0_DOR[:,:,:,vol_DOR])
                E_LTE = np.squeeze(S_S0_dwi[:,:,:,vol_DWI])
                var_u = np.log(E_LTE/E_STE) * 2 * (common_bvalue)**(-2) * MD**(-2)
                microFA = np.sqrt(3 / 2) * (1 + (2/5)*(1/var_u) )**-0.5
                bad_FA = np.where(E_LTE-E_STE < 0) 
                
                for i in range(0,len(bad_FA[0])):
                    microFA[bad_FA[0][i],bad_FA[1][i],bad_FA[2][i]]=np.nan;
                
                # Save image
                affine = nib.load(os.path.join(bids_strc_analysis.get_path(),'md_dki.nii')).affine
                img = nib.Nifti1Image(microFA, affine)
                nib.save(img, os.path.join(bids_strc_analysis.get_path(),'microFA.nii'))

            
        
            