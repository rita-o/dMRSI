"""
Script to analyse GM microstructure by fitting models such as Nexi, Sandi, ...
It uses the SwissKnife module installed in the SwissKnife python environment.

Last changed Jan 2025
@author: Rita O
"""

from graymatter_swissknife import estimate_model
import os
import sys
import pandas as pd
import platform
import math
import importlib, sys
from custom_functions import *
from bids_structure import *

def Step4_modelling_GM(subj_list, cfg):
    
    data_path   = cfg['data_path']     
    scan_list   = pd.read_excel(os.path.join(data_path, 'ScanList.xlsx'))
    
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in subj_list:
        
        print('Modelling ' + subj + '...')
    
        # Extract data for subject
        subj_data      = scan_list[(scan_list['newstudyName'] == subj)].reset_index(drop=True)
        
        ######## SESSION-WISE OPERATIONS ########
        for sess in list(subj_data['blockNo'].unique()) :
          

            ######## MODEL-WISE OPERATIONS ########
            for model in cfg['model_list_GM']:
                
                if model=='Nexi':
                    data_used = 'allDelta-allb'
                elif model=='Sandi':
                    filtered_data = subj_data[(subj_data['acqType'] == 'PGSE') & (subj_data['phaseDir'] == 'fwd') & (subj_data['blockNo'] == sess) & (subj_data['noBval'] > 1)]
                    ind_folder = getattr(filtered_data["diffTime"], 'idxmin')()
                    data_used = 'Delta_'+str(int(filtered_data['diffTime'][ind_folder]))+'_fwd'  
                
                # Define bids structure 
                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description=model)
                bids_strc_prep = create_bids_structure(subj=subj, sess=sess, datatype="dwi", description=data_used, root=data_path, 
                                          folderlevel='derivatives', workingdir=cfg['prep_foldername'])
                
                # Make output folder 
                output_path = bids_strc_analysis.get_path()
                input_path = os.path.join(output_path,'inputs')
                create_directory(input_path)

                # Copy necessary files for analysis
                dwi         = copy_files_BIDS(bids_strc_prep,input_path,'dwi_dn_gc_ec.nii.gz')
                big_delta   = copy_files_BIDS(bids_strc_prep,input_path,'DiffTime.txt')
                small_delta = copy_files_BIDS(bids_strc_prep,input_path,'DiffDuration.txt')
                bvals       = copy_files_BIDS(bids_strc_prep,input_path,'bvalsNom.txt')
                #b0        = copy_files_BIDS(bids_strc_prep,input_path,'b0_dn_gc_ec_avg.nii.gz') # for atlas

                # Get diffusion duration (assumes the same value for all acquisitions)
                small_delta = np.loadtxt(small_delta)[0]
         
                # Modify units of bvals for NEXI          
                new_bvals = bvals.replace('.txt','_units.txt')
                modify_units_bvals(bvals, new_bvals )
        
                # Copy necessary files for analysis 
                if model=='Nexi':
                    bids_strc_lowb = create_bids_structure(subj=subj, sess=sess, datatype="dwi", description="allDelta-lowb", root=data_path, 
                                                folderlevel='derivatives', workingdir=cfg['prep_foldername'])
                    sigma     = copy_files_BIDS(bids_strc_lowb,input_path,'dwi_dn_sigma.nii.gz')
                elif model=='Sandi':
                    sigma     = copy_files_BIDS(bids_strc_prep,input_path,'dwi_dn_sigma.nii.gz')
   

                # Estimate model
                estimate_model(
                    model,
                    dwi,
                    new_bvals,
                    big_delta,
                    small_delta,
                    sigma,
                    output_path
                )
                
                # Mask output for better visualization
                for filename in os.listdir(output_path):
                    if filename.endswith(".nii.gz"):
                        multiply_by_mask(os.path.join(output_path, filename), # filename input
                                         os.path.join(output_path,'Masked'), # output folder
                                         bids_strc_prep.get_path('mask.nii.gz')) # mask


