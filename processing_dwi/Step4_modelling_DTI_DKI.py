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
    
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in subj_list:
        
        print('Modelling ' + subj + '...')
    
        # Extract data for subject
        subj_data      = scan_list[(scan_list['newstudyName'] == subj)].reset_index(drop=True)
          
        # Define path to docker
        docker_path = '/data'
        
        ######## SESSION-WISE OPERATIONS ########
        for sess in list(subj_data['blockNo'].unique()):
          
            for type in ['dwi','dwi_DOR']:
               
                # Define bids structure for the processed data
                bids_strc_prep = create_bids_structure(subj=subj, sess=sess, datatype=type, root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['prep_foldername'])
    
                if type=='dwi':
                    # By default get dataset with largest diffusion time
                    filtered_data = subj_data[(subj_data['phaseDir'] == 'fwd') & (subj_data['noBval'] > 1)]
                    ind_folder =filtered_data["diffTime"].idxmax()
                    bids_strc_prep.set_param(description='Delta_'+str(int(filtered_data['diffTime'][ind_folder]))+'_'+filtered_data['phaseDir'][ind_folder])
                    type_name = 'LTE'
                    dwi_filename = 'dwi_dn_gc_ec.mif'
                    mask_filename = 'b0_avg_mask.nii.gz'
                    
                else:
                    bids_strc_prep.set_param(description='fwd')
                    type_name = 'STE'
                    dwi_filename = 'dwi_dn_gc_topup.mif'
                    mask_filename = 'b0_mask.nii.gz'
                 
                ######## MODEL-WISE OPERATIONS ########     
                 
                # Define BIDS structure for the analysis data
                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype="dwi", description=f'DTI_DKI_{type_name}', root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
              
                # Make output folder 
                output_path = bids_strc_analysis.get_path()
                input_path = os.path.join(output_path,'inputs')
                create_directory(input_path)
                    
                # Copy necessary files for analysis and rename the path to the docker path
                dwi   = copy_files_BIDS(bids_strc_prep,input_path, dwi_filename).replace(data_path,docker_path)
                mask  = copy_files_BIDS(bids_strc_prep,input_path, mask_filename).replace(data_path,docker_path)
                out_folder   = output_path.replace(data_path,'/data')
                
                # Run model
                call = [f'docker run -v {data_path}:/data nyudiffusionmri/designer2:v2.0.10 tmi -DTI -DKI',
                        f'{dwi} {out_folder}'] #
         
                print(' '.join(call))
                os.system(' '.join(call))
                
                # Calculate microFA
                affine = nib.load(os.path.join(bids_strc_analysis.get_path(),'md_dki.nii')).affine
                MD = nib.load(os.path.join(bids_strc_analysis.get_path(),'md_dki.nii')).get_fdata()
                MK = nib.load(os.path.join(bids_strc_analysis.get_path(),'mk_dki.nii')).get_fdata()
                vADC = (MK*MD**2)/3
                microFA = np.sqrt(3 / 2) * (1 + (MD**2)/(5/2*vADC) )**-0.5
                img = nib.Nifti1Image(microFA, affine)
                nib.save(img, os.path.join(bids_strc_analysis.get_path(),'microFA.nii'))
    
    
                # Rename paths to local folder
                bids_strc_analysis.set_param(root=data_path)
                bids_strc_prep.set_param(root=data_path)
    
                output_path = bids_strc_analysis.get_path()
    
                # Mask output with brain mask for better visualization
                for filename in os.listdir(output_path):
                    if filename.endswith(".nii"):
                        multiply_by_mask(os.path.join(output_path, filename), # filename input
                                         os.path.join(output_path,'Masked'), # output folder
                                         bids_strc_prep.get_path('b0_avg_mask.nii.gz')) # mask
            
        
            