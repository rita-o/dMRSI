"""
Script to correct orientation labels of the nifties that come out of Bruker system
because they are usually not correct (for now)
It does not use a particular python environment

Last changed Jan 2025
@author: Rita O
"""

import os
import sys
import pandas as pd
import platform
import math
from custom_functions import *
from bids_structure import *

def Step2_correct_orientation(subj_list,cfg):
    
    data_path       = cfg['data_path']     
    scan_list       = pd.read_excel(os.path.join(data_path, cfg['scan_list_name'] ))
    
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in  subj_list:
    
        print('Correcting of orientation of ' + subj + '...')
    
        # Extract data for this subject
        subj_data      = scan_list[(scan_list['newstudyName'] == subj)].reset_index(drop=True)
    
        ######## SESSION-WISE OPERATIONS ########
        for sess in list(subj_data['blockNo'].unique()) :
            
            bids_strc = create_bids_structure(subj=subj, sess=sess, datatype="dwi", root=data_path, 
                                        folderlevel='nifti_data', workingdir='sorted')
          
            # Index of scans for this session 
            study_indx  = subj_data.index[subj_data['blockNo'] == sess].tolist()
            
            ###### SCAN-WISE OPERATIONS ######
            for scn_ctr in study_indx:
                
                # reorient images and save in respective path
                if subj_data['scanQA'][scn_ctr] == 'ok' and subj_data['acqType'][scn_ctr] == 'PGSE':
                    bids_strc.set_param(datatype='dwi',description='Delta_'+str(int(subj_data['diffTime'][scn_ctr]))+'_'+subj_data['phaseDir'][scn_ctr])
                    reorient_nifit(bids_strc.get_path('dwi.nii.gz'), subj_data['Notes'][scn_ctr])
        
                elif subj_data['scanQA'][scn_ctr] == 'ok' and subj_data['acqType'][scn_ctr] == 'T2W':
                    bids_strc.set_param(datatype='anat',description=None)
                    reorient_nifit(bids_strc.get_path('T2w.nii.gz'), subj_data['Notes'][scn_ctr])
                
                elif subj_data['scanQA'][scn_ctr] == 'ok' and subj_data['acqType'][scn_ctr] == 'STE':
                    bids_strc.set_param(datatype='dwi_STE',description='STE_'+subj_data['phaseDir'][scn_ctr])
                    reorient_nifit(bids_strc.get_path('dwi.nii.gz'), subj_data['Notes'][scn_ctr])