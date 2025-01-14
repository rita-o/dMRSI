import os
import sys
import pandas as pd
import platform
import math
from custom_functions import *
from bids_structure import *

def Step2_correct_orientation(subj_list,cfg):
    
    data_path       = cfg['data_path']     
    scan_list       = pd.read_excel(os.path.join(data_path, 'ScanList.xlsx'))
    
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in  subj_list:
    
        print('Coreecting of orientation of ' + subj + '...')
    
        # Extract data by study
        subj_data      = scan_list[(scan_list['newstudyName'] == subj)].reset_index(drop=True)

        # List of scans for each study
        study_scanNo     = list(subj_data['scanNo'])
        new_orientation  = list(subj_data['Notes'])
    
        ######## SESSION-WISE OPERATIONS ########
        for sess in list(subj_data['blockNo'].unique()) :
            
            bids_strc = create_bids_structure(subj=subj, sess=sess, datatype="dwi", root=data_path, 
                                        folderlevel='derivatives', workingdir=cfg['prep_foldername'])
          
            ###### SCAN-WISE OPERATIONS ######
            for scn_ctr in range (len(study_scanNo)):
                
            
                if subj_data['scanQA'][scn_ctr] == 'ok' and subj_data['acqType'][scn_ctr] == 'PGSE':
                    bids_strc.set_param(datatype='dwi',description='E'+str(study_scanNo[scn_ctr]))
                   
                    reorient_nifit(bids_strc.get_path('dwi.nii.gz'), new_orientation[scn_ctr])
        
                elif subj_data['scanQA'][scn_ctr] == 'ok' and subj_data['acqType'][scn_ctr] == 'T2W':
        
                    # Get paths and directories
                    bids_strc.set_param(datatype='anat',description=None)

                    
                    reorient_nifit(bids_strc.get_path('T2w.nii.gz'), new_orientation[scn_ctr])
                
                