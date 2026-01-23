#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 14:28:52 2026

@author: localadmin
"""
import os
import sys
import pandas as pd
import math
from custom_functions import *
from bids_structure import *
from pathlib import Path

def Step1_preproc(subj_list, cfg):
    
    
    code_path = cfg['code_path2']
    
    data_path       = cfg['data_path']     
    scan_list       = pd.read_excel(os.path.join(data_path, cfg['scan_list_name'] ))
    
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in  subj_list:
    
        print('Processing dmrs of subject ' + subj + '...')
    
        # Extract data for this subject
        subj_data      = scan_list[(scan_list['study_name'] == subj)].reset_index(drop=True)
        
        # Generate paths and convert data 
        raw_path        = os.path.join(data_path, 'raw_data', list(subj_data['raw_data_folder'].unique())[0]) 
        preproc_path    = os.path.join(data_path, 'derivatives', cfg['prep_foldername'], subj,'dmrs')

        # Extract data for subject
        subj_data = scan_list[
                            (scan_list['study_name'] == subj) &
                            (scan_list['acqType'] == 'SPECIAL') &
                            (scan_list['analyse'] == 'y')
                             ].reset_index(drop=True)

        # List of acquisition sessions
        sess_list    = [x for x in list(subj_data['sessNo'].unique()) if not math.isnan(x)] # clean NaNs
    
        ######## SESSION-WISE OPERATIONS ########
        for sess in list(subj_data['sessNo'].unique()) :
            
            print(f"Processing session {sess} ...")
            
            bids_strc = create_bids_structure(subj=subj, sess=sess, datatype="dmrs", root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['prep_foldername'])
            # Loop for different Mixing Times
            TM_list             = np.unique(subj_data['TM'].astype(int).tolist())

            for TM in TM_list:
                            
                print(f"Processing mixing time: {TM} ms ...")
                
                # Get the scan numbers for the water reference. Assumes there is only one
                subj_data_water = subj_data[
                                    (subj_data['sessNo'] == sess) &
                                    (subj_data['phaseDir'] == 'water') &
                                    (subj_data['TM'] == TM) 
                                     ].reset_index(drop=True)
                water_reference_sequence_number =  subj_data_water['scanNo'].iloc[0]
                
                # Get the scan numbers for the metabolite data 
                subj_data_metab = subj_data[
                                    (subj_data['sessNo'] == sess) &
                                    (subj_data['phaseDir'] == 'metab') &
                                    (subj_data['TM'] == TM)
                                     ].reset_index(drop=True)
                metab_sequence_number = subj_data_metab['scanNo'].dropna().str.split(',').explode().astype(int).tolist()
    
                # Define parameters to input in matlab   
                input_path        = raw_path
                output_path       = os.path.join(bids_strc.get_path(),f'TM_{str(TM)}')
                create_directory(output_path)
                scan_list = metab_sequence_number
                scan_list_format_matlab = "[" + " ".join(map(str, scan_list)) + "]"
                coil_type         = cfg['coil_type'] 
                toolbox_path      = os.path.join(cfg['code_path2'],'dSPECIAL_matlab_codes_Toi')
                LCMpath           = cfg['LC_model']
    
                basis_sets_folder = cfg['basis_set']
                #basis_set         = str(next(Path(basis_sets_folder).glob(f"*TM{TM}*")))
                basis_set         = '/home/localadmin/Documents/Rita/Data/common/mrs_basis_sets/Basis_Set_dSPECIAL_differentTM/9p4T_Toi_dSPECIAL_TM150_20260106.BASIS'
                    
                # Matlab command
                matlab_cmd = (
                    "try, "
                    f"addpath('{code_path}'); "
                    f"processing_dmrs_matlab("
                    f"'{input_path}', "
                    f"'{output_path}', "
                    f"'{scan_list_format_matlab}', "
                    f"'{coil_type}', "
                    f"'{toolbox_path}', "
                    f"'{basis_set}', "
                    f"'{LCMpath}'); "
                    "catch ME, disp(getReport(ME)); exit(1); "
                    "end; exit(0);"
                )
                
                cmd = [
                    "matlab", "-nodisplay", "-nosplash", "-nodesktop",
                    "-r", matlab_cmd
                ]
                
                print("\nMATLAB command:")
                print(" ".join(cmd))
                print()
                subprocess.run(cmd)
