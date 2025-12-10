#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 20 14:30:54 2025

@author: localadmin
"""

import os
import sys
import pandas as pd
import platform
import math
import shutil
import importlib, sys

#import my modules
cfg_path = sys.argv[1] 
config_file = os.path.join(cfg_path, '.config.json')
import json
with open(config_file, 'r') as f:
    cfg = json.load(f)

sys.path.append(cfg['code_path'])
sys.path.append(cfg['code_path2'])
from custom_functions import *
from bids_structure import *

def Step2_raw2nii2bids_human(subj_list,cfg):
    
    data_path       = cfg['data_path']   

    scan_list       = pd.read_excel(os.path.join(data_path, cfg['scan_list_name']))
    
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in  subj_list:
    
        print('Converting to nifti ' + subj + '...')
    
        # Extract data for this subject
        subj_data      = scan_list[(scan_list['study_name'] == subj)].reset_index(drop=True)
        
        # Generate paths and convert data from DICOM to NIFTI
        raw_path        = os.path.join(data_path, 'raw_data', list(subj_data['raw_data_folder'].unique())[0]) 
        nifti_path      = os.path.join(data_path, 'nifti_data', 'sorted', subj)
        create_directory(nifti_path)
        
        
        ######## SESSION-WISE OPERATIONS ########
        for sess in list(subj_data['sessNo'].unique()) :
            
            bids_strc = create_bids_structure(subj=subj, sess=sess, datatype="dwi", root=data_path, 
                                        folderlevel='nifti_data', workingdir='sorted')
            
            # Index of scans for this session 
            study_indx  = subj_data.index[subj_data['sessNo'] == sess].tolist()


            ###### SCAN-WISE OPERATIONS ######
            for scn_ctr in study_indx:
                 
                 # Scan folder number
                 scan_no = subj_data['scanNo'][scn_ctr]
                 
                 # Convert indidual folders
                 for file in os.listdir(raw_path):
                    if '0'+str(scan_no) + '-' in file:
                        dir_to_convert = os.path.join(raw_path,file)
                        if subj_data['analyse'][scn_ctr] =='y' and subj_data['acqType'][scn_ctr]=='PGSE':
                            delta_val = int(subj_data['diffTime'][scn_ctr])
                            dir_to_save = os.path.join(nifti_path,f"ses-{sess:02}",'dwi','Delta_'+str(delta_val)+'_'+subj_data['phaseDir'][scn_ctr])
                        else:
                            dir_to_save = os.path.join(nifti_path,f"ses-{sess:02}",'anat')
                        create_directory(dir_to_save)
                        call = [f'dcm2niix -ba n -o {dir_to_save} {dir_to_convert}']
                        os.system(' '.join(call))
                        
                       
                        # Organize
                        for filename in os.listdir(dir_to_save):
                            if subj_data['analyse'][scn_ctr] =='y' and subj_data['acqType'][scn_ctr]=='PGSE':
                                 bids_strc = create_bids_structure(subj=subj, sess=sess, datatype="dwi", root=data_path,description='Delta_'+str(delta_val)+'_'+subj_data['phaseDir'][scn_ctr], 
                                                            folderlevel='nifti_data', workingdir='sorted')
                                
                                 if 'nii' in filename:
                                     os.rename(os.path.join(dir_to_save,filename), bids_strc.get_path('dwi.nii'))
                                     gzip_file( bids_strc.get_path('dwi.nii'))
                                     os.remove( bids_strc.get_path('dwi.nii'))
                                 elif 'bval' in filename:
                                     os.rename(os.path.join(dir_to_save,filename),  bids_strc.get_path('bvalsNom.txt'))
                                     shutil.copy2(bids_strc.get_path('bvalsNom.txt'), bids_strc.get_path('bvalsEff.txt'))
                                 elif 'bvec' in filename:
                                     os.rename(os.path.join(dir_to_save,filename), bids_strc.get_path('bvecs.txt'))
                                 elif 'json' in filename:
                                     os.rename(os.path.join(dir_to_save,filename), bids_strc.get_path('dwi.json'))
                            elif subj_data['analyse'][scn_ctr] =='y' and subj_data['acqType'][scn_ctr]=='T1W':
                                 bids_strc = create_bids_structure(subj=subj, sess=sess, datatype="anat", root=data_path, 
                                                           folderlevel='nifti_data', workingdir='sorted')
                                 if 'nii' in filename:
                                     os.rename(os.path.join(dir_to_save,filename), bids_strc.get_path('T1w.nii'))
                                     gzip_file(bids_strc.get_path('T1w.nii'))
                                     os.remove(bids_strc.get_path('T1w.nii'))                           
                                 elif 'json' in filename:
                                     os.rename(os.path.join(dir_to_save,filename), bids_strc.get_path('T1w.json'))
                    
      

           

if __name__ == "__main__":
    import json
    import sys
    import os

    cfg_data_path = str(sys.argv[1])
    with open(os.path.join(cfg_data_path, '.config.json')) as f:
        cfg = json.load(f)

    # Add code paths
    sys.path.append(cfg['code_path2'])

    from Step2_raw2nii2bids_human import *

    subj_list = cfg['subj_list']
    Step2_raw2nii2bids_human(subj_list, cfg)