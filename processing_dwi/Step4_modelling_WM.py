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

def Step4_modelling_WM(subj_list, cfg):
    
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
          
            # Define bids structure for the processed data
            bids_strc_prep = create_bids_structure(subj=subj, sess=sess, datatype="dwi", root=data_path, 
                                        folderlevel='derivatives', workingdir=cfg['prep_foldername'])


            # By default get dataset with largest diffusion time
            filtered_data = subj_data[(subj_data['phaseDir'] == 'fwd') & (subj_data['noBval'] > 1)]
            ind_folder =filtered_data["diffTime"].idxmax()
            bids_strc_prep.set_param(description='Delta_'+str(int(filtered_data['diffTime'][ind_folder]))+'_'+filtered_data['phaseDir'][ind_folder])
               
            
            ######## MODEL-WISE OPERATIONS ########
            for model in cfg['model_list_WM']:
                    
                # Define BIDS structure for the analysis data
                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype="dwi", description=model, root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
              
                # Make output folder 
                output_path = bids_strc_analysis.get_path()
                input_path = os.path.join(output_path,'inputs')
                create_directory(input_path)
                
                # Copy necessary files for analysis and rename the path to the docker path
                dwi   = copy_files_BIDS(bids_strc_prep,input_path,'dwi_dn_gc_ec.mif').replace(data_path,docker_path)
                mask  = copy_files_BIDS(bids_strc_prep,input_path,'b0_avg_mask.nii.gz').replace(data_path,docker_path)
                sigma = copy_files_BIDS(bids_strc_prep,input_path,'dwi_dn_sigma.nii.gz').replace(data_path,docker_path)
                #b0 = copy_files_BIDS(bids_strc_prep,input_path,'b0_dn_gc_ec_avg.nii.gz')
                
                # RUN MODEL ESTIMATE
                estim_SMI_designer(dwi,
                                   mask, 
                                   sigma,
                                   output_path.replace(data_path,docker_path), data_path)
                
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
                
            
            