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

    scan_list       = pd.read_excel(os.path.join(data_path, 'ScanList.xlsx'))
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in  subj_list:
    
        print('Converting to nifti ' + subj + '...')
    
        # Extract data for this subject
        subj_data      = scan_list[(scan_list['newstudyName'] == subj)].reset_index(drop=True)
        
        # Generate paths and convert data from DICOM to NIFTI
        raw_path        = os.path.join(data_path, 'raw_data', list(subj_data['studyName'].unique())[0]) 
        nifti_path      = os.path.join(data_path, 'nifti_data', 'unsorted', subj)
        if not os.path.exists(nifti_path):
           os.makedirs(nifti_path)
        
        call = [f'dcm2niix -ba n -o {nifti_path} {raw_path}']
        os.system(' '.join(call))
        
        # Convert from NIFTI to BIDS
        nifti_path2      = os.path.join(data_path, 'nifti_data', 'sorted')
        if not os.path.exists(nifti_path2):
           os.makedirs(nifti_path2)
        call = [f' niix2bids -i {nifti_path} -o {nifti_path2}']
        os.system(' '.join(call))
        
        # Remove folders
        for name in os.listdir(nifti_path2):
            if not name.startswith("sub"):
                path = os.path.join(nifti_path2, name)
                if os.path.isfile(path):
                    os.remove(path)
                elif os.path.isdir(path):
                    shutil.rmtree(path)
    
        # Rename folder
        for name in os.listdir(nifti_path2):
            path = os.path.join(nifti_path2, name)
            if name.startswith("sub") and os.path.isdir(path):
                new_path = os.path.join(nifti_path2, subj)
                os.rename(path, new_path)
                
        # Organize in folders
        pattern = re.compile(r'D(\d{2})')
        for sess in list(subj_data['blockNo'].unique()) :
            
            bids_strc = create_bids_structure(subj=subj, sess=sess, datatype="dwi", root=data_path, 
                                        folderlevel='nifti_data', workingdir='sorted')
            
            nifti_path3  =  os.path.join(nifti_path2, subj, f'ses-0{sess}')
            os.rename(os.path.join(nifti_path2, subj, f'ses-{sess}'), nifti_path3)
            
            nifti_path3  = os.path.join(nifti_path3 ,'dwi')
            for filename in os.listdir(nifti_path3):
                if 'run-1' in filename:
                    match = pattern.search(filename)
                    if match:
                        delta_val = match.group(1)
                        bids_strc.set_param(datatype='dwi',description='Delta_'+str(delta_val)+'_fwd')

                        folder_name = f"Delta_{delta_val}_fwd"
                        dest_folder = os.path.join(nifti_path3, folder_name)
                        os.makedirs(dest_folder, exist_ok=True)
            
                        src = os.path.join(nifti_path3, filename)
                        dst = os.path.join(dest_folder, filename)
                        if 'nii' in filename:
                            shutil.copy2(src, bids_strc.get_path('dwi.nii'))
                            gzip_file(bids_strc.get_path('dwi.nii'))
                            os.remove(bids_strc.get_path('dwi.nii'))
                        elif 'bval' in filename:
                                shutil.copy2(src, bids_strc.get_path('bvalsNom.txt'))
                                shutil.copy2(src, bids_strc.get_path('bvalsEff.txt'))
                        elif 'bvec' in filename:
                                shutil.copy2(src, bids_strc.get_path('bvecs.txt'))
                        elif 'json' in filename:
                                shutil.copy2(src, bids_strc.get_path('dwi.json'))
                        print(f"Copied: {src} -> {dst}")
                        
            # reverse phasing
            pattern = re.compile(r'PA')
            for filename in os.listdir(nifti_path3):
                if pattern.search(filename):
                    src = os.path.join(nifti_path3, filename)
                    bids_strc.set_param(datatype='dwi',description='Delta_26_'+'rev')
                    os.makedirs( bids_strc.get_path(), exist_ok=True)

                    if 'nii' in filename:
                        shutil.copy2(src, bids_strc.get_path('dwi.nii'))
                        gzip_file(bids_strc.get_path('dwi.nii'))
                        os.remove(bids_strc.get_path('dwi.nii'))
                    elif 'bval' in filename:
                            shutil.copy2(src, bids_strc.get_path('bvalsNom.txt'))
                            shutil.copy2(src, bids_strc.get_path('bvalsEff.txt'))
                    elif 'bvec' in filename:
                            shutil.copy2(src, bids_strc.get_path('bvecs.txt'))
                    elif 'json' in filename:
                            shutil.copy2(src, bids_strc.get_path('dwi.json'))
            
                    print(f"Copied: {src} -> {dst}")
                                
            # Remove folders
            for name in os.listdir(nifti_path3):
                 if not name.startswith("Delta") and not name.startswith("rev"):
                     path = os.path.join(nifti_path3, name)
                     if os.path.isfile(path):
                         os.remove(path)
                     elif os.path.isdir(path):
                         shutil.rmtree(path)
            # Anat 
            nifti_path3  = os.path.join(data_path, 'nifti_data', 'sorted' ,subj, f'ses-0{sess}','anat')
            bids_strc.set_param(datatype='anat',description='')
            for filename in os.listdir(bids_strc.get_path()):
                if 'run-1' in filename:
                        if 'nii' in filename:
                            os.rename(os.path.join(nifti_path3,filename),  bids_strc.get_path('T1w.nii'))
                            gzip_file(bids_strc.get_path('T1w.nii'))
                            os.remove(bids_strc.get_path('T1w.nii'))
                            extract_vols(bids_strc.get_path('T1w.nii.gz'),bids_strc.get_path('T1w.nii.gz'),0,1)
                        elif 'json' in filename:
                            os.rename(os.path.join(nifti_path3,filename),  bids_strc.get_path('T1w.json'))
                if not 'run-1' in filename:
                         os.remove(os.path.join(nifti_path3,filename))

           

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