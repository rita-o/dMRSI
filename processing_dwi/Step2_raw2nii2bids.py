import os
import sys
import pandas as pd
import platform
import math
from custom_functions import *
from bids_structure import *

def Step2_raw2nii2bids(subj_list,cfg):
    
    data_path       = cfg['data_path']     
    scan_list       = pd.read_excel(os.path.join(data_path, 'ScanList.xlsx'))
    
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in  subj_list:
    
        print('Converting to nifti ' + subj + '...')
    
        # Extract data by study
        subj_data      = scan_list[(scan_list['newstudyName'] == subj)].reset_index(drop=True)
        
        # List of scans for each study
        study_scanNo  = list(subj_data['scanNo'])
       
        # Generate paths and convert data 
        raw_path        = os.path.join(data_path, 'raw_data', list(subj_data['studyName'].unique())[0]) 
        nifti_path      = os.path.join(data_path, 'nifti_data', subj)
        create_directory(nifti_path)
        raw_to_nifti(raw_path, nifti_path)
    
        ######## SESSION-WISE OPERATIONS ########
        for sess in list(subj_data['blockNo'].unique()) :
            
            bids_strc = create_bids_structure(subj=subj, sess=sess, datatype="dwi", root=data_path, 
                                        folderlevel='derivatives', workingdir=cfg['prep_foldername'])
            
            ###### SCAN-WISE OPERATIONS ######
            for scn_ctr in range (len(study_scanNo)):
                
                
                if subj_data['scanQA'][scn_ctr] == 'ok' and subj_data['acqType'][scn_ctr] == 'PGSE':
        
                    # Get paths and directories
                    method_path   = os.path.join(raw_path,str(study_scanNo[scn_ctr]), 'method')
                    nii_path    = os.path.join(nifti_path,str(study_scanNo[scn_ctr]) + '_1_' + subj_data['acqSeq'][scn_ctr])
                    bids_strc.set_param(datatype='dwi',description='E'+str(study_scanNo[scn_ctr]))
        
                    # Transfer files
                    create_directory(bids_strc.get_path())
                    copy_file([os.path.join(nii_path, '1.nii.gz')], [bids_strc.get_path('dwi.nii.gz')])
                    extract_methods(method_path, bids_strc, 'diff')
                    plot_bvals(bids_strc)
        
    
                elif subj_data['scanQA'][scn_ctr] == 'ok' and subj_data['acqType'][scn_ctr] == 'T2W':
        
                    # Get paths and directories
                    nii_path    = os.path.join(nifti_path,str(study_scanNo[scn_ctr]) + '_1_' + subj_data['acqSeq'][scn_ctr])
                    bids_strc.set_param(datatype='anat',description=None)

                    # Transfer files
                    create_directory(bids_strc.get_path())
                    copy_file([os.path.join(nii_path, '1.nii.gz')], [bids_strc.get_path('T2w.nii.gz')])
        
                elif subj_data['scanQA'][scn_ctr] == 'ok' and subj_data['acqType'][scn_ctr] == 'B0': # probably not working
        
                    # Get paths and directories
                    nii_path    = os.path.join(nifti_path,str(study_scanNo[scn_ctr])  + '_1_' + subj_data['acqSeq'][scn_ctr])
                    bids_strc.set_param(datatype='B0',description=None)

                    # Transfer files
                    create_directory(bids_strc.get_path())
                    copy_file([os.path.join(nii_path, '1.nii.gz')], [bids_strc.get_path('B0.nii.gz')])
        
            #shutil.rmtree(os.path.join(data_path, 'niftiData'))