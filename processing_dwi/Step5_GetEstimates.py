"""
Script to retreive model estimates within regions of interest.
It does not use a particular python environment.

* not finished yet * 
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
import glob
import numpy as np
import SimpleITK as sitk
import numpy.ma as ma

def Step5_GetEstimates(data_path, subj_list, cfg):
    
    
    data_path   = cfg['data_path']     
    scan_list   = pd.read_excel(os.path.join(data_path, 'ScanList.xlsx'))
    atlas_path  = cfg['common_folder']
                 
    ########################## SUBJECT-WISE OPERATIONS ##########################
    for subj in subj_list:
        
        print('Getting model estimates ' + subj + '...')
    
        # Extract data for subject
        subj_data    = scan_list[scan_list['newstudyName'] == subj].reset_index(drop=True)
        
        # List of acquisition sessions
        sess_list    = [x for x in list(subj_data['blockNo'].unique()) if not math.isnan(x)] # clean NaNs
        
        ###### SESSION-WISE OPERATIONS 
        for sess in list(subj_data['blockNo'].unique()):
            
         
            for model in cfg['model_list']:
                
                # Define BIDS structure for the analysis data
                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype=model, description='inputs', root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                # Atlas - not finished
                if not os.path.exists(bids_strc_analysis.get_path('template_in_dwi.nii.gz')):
                    atlas    = glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*atlas.nii.gz'))[0]
                    template = glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*T2s.nii.gz'))[0]
                    b0       = find_files_with_pattern(bids_strc_analysis,'b0_dn_gc_ec_avg.nii.gz')[0]
                    b0_bc    = b0.replace('.nii.gz','_bc.nii.gz')
                    N4_unbias(b0,b0_bc)
                    output_path = os.path.dirname(b0_bc)
                    antsreg(template,b0_bc, os.path.join(output_path,'atlas2dwi'))
                    ants_apply_transforms([template,   atlas],
                                          b0_bc,
                                          [os.path.join(output_path,'template_in_dwi.nii.gz'),os.path.join(output_path,'atlas_in_dwi.nii.gz')],
                                           os.path.join(output_path,'atlas2dwi1Warp.nii.gz'),
                                           os.path.join(output_path,'atlas2dwi0GenericAffine.mat'))  
                
                
              
                atlas_path = os.path.join(output_path,'template_in_dwi.nii.gz')
                atlas = nib.load(os.path.join(output_path,'template_in_dwi.nii.gz')).get_fdata()

                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype=model, root=data_path, 
                             folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                output_path = bids_strc_analysis.get_path()
                
                atlas_labels = pd.read_csv(
                    glob.glob(os.path.join(cfg['common_folder'], cfg['atlas'], '*atlas.label'))[0],
                    delim_whitespace=True,  # Use whitespace as the delimiter
                    skiprows=14,  # Skip the header lines (modify as per actual file structure)
                    header=None,  # No column headers in the file
                    names=['IDX', 'R', 'G', 'B', 'A', 'VIS', 'MSH', 'LABEL'],  # Assign column names
                    quotechar='"',  # Handle quoted strings for the LABEL column
                )
   

                for ROI in cfg['ROIs']:
                    
                    ind_list=atlas_labels["LABEL"].str.find(ROI)
                    match_idx = ind_list[ind_list != -1].index
                    match_idx= atlas_labels["IDX"][match_idx].to_numpy()
                    mask_indexes = (atlas==match_idx)
                                      # List of specific patterns to match
                    patterns = ["*de.nii.gz", "*t.nii.gz"]
                    
                    # Loop through each specified pattern and process the first matching file
                    for pattern in patterns:
                        parameter_data = nib.load(glob.glob(os.path.join(output_path, pattern))[0]).get_fdata()  
                        param_masked = parameter_data * mask_indexes
                        
                        
                        data = np.multiply(mask_data.get_fdata(),)
                   


              
                        df = pd.read_table(atlas_labels,skiprows=13,sep='\s{3}')
