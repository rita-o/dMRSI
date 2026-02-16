"""
Script to correct orientation labels of the nifties that come out of Bruker system
because they are usually not correct (for now)
It does not use a particular python environment

Last changed Jan 2025
@author: Rita O
"""

import os
import pandas as pd
from dmri_dmrs_toolbox.misc.bids_structure import create_bids_structure
from dmri_dmrs_toolbox.misc.custom_functions import reorient_nifit

def Step2_correct_orientation(cfg):
    
    data_path       = cfg['data_path']     
    scan_list       = pd.read_excel(os.path.join(data_path, cfg['scan_list_name'] ))
    
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in cfg['subj_list']:
    
        print('Correcting of orientation of ' + subj + '...')
    
        # Extract data for this subject
        subj_data      = scan_list[(scan_list['study_name'] == subj)].reset_index(drop=True)
    
        ######## SESSION-WISE OPERATIONS ########
        for sess in list(subj_data['sessNo'].unique()) :
            
            bids_strc = create_bids_structure(subj=subj, sess=sess, datatype="dwi", root=data_path, 
                                        folderlevel='nifti_data', workingdir='sorted')
          
            # Index of scans for this session 
            study_indx  = subj_data.index[subj_data['sessNo'] == sess].tolist()
            
            ###### SCAN-WISE OPERATIONS ######
            for scn_ctr in study_indx:
                
                # reorient images and save in respective path
                if subj_data['analyse'][scn_ctr] =='y' and subj_data['acqType'][scn_ctr] == 'PGSE':
                    bids_strc.set_param(datatype='dwi',description='Delta_'+str(int(subj_data['diffTime'][scn_ctr]))+'_'+subj_data['phaseDir'][scn_ctr])
                    reorient_nifit(bids_strc.get_path('dwi.nii.gz'), subj_data['Reorient'][scn_ctr],cfg)
        
                elif subj_data['analyse'][scn_ctr] =='y' and subj_data['acqType'][scn_ctr] == 'T2W':
                    bids_strc.set_param(datatype='anat',description=None)
                    reorient_nifit(bids_strc.get_path('T2w.nii.gz'), subj_data['Reorient'][scn_ctr],cfg)
                
                elif subj_data['analyse'][scn_ctr] =='y' and subj_data['acqType'][scn_ctr] == 'STE':
                    bids_strc.set_param(datatype='dwi_STE',description='STE_'+subj_data['phaseDir'][scn_ctr])
                    reorient_nifit(bids_strc.get_path('dwi.nii.gz'), subj_data['Reorient'][scn_ctr],cfg)