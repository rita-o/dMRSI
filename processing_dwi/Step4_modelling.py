"""
Script to analyse GM microstructure by fitting models such as Nexi, Sandi, ...
It uses the SwissKnife module installed in the SwissKnife python environment.

Last changed Jan 2025
@author: Rita O
"""
import sys
import os
import json
#from graymatter_swissknife import estimate_model
import pandas as pd
import platform
import math
import importlib, sys

# import my modules
# cfg_path = sys.argv[1] 
# config_file = os.path.join(cfg_path, '.config.json')
# import json
# with open(config_file, 'r') as f:
#     cfg = json.load(f)

# sys.path.append(cfg['code_path'])
# sys.path.append(os.path.join(cfg['code_path'], 'processing_dwi'))
from custom_functions import *
from bids_structure import *

def Step4_modelling(subj_list, cfg):
    
    data_path   = cfg['data_path']     
    scan_list   = pd.read_excel(os.path.join(data_path, 'ScanList.xlsx'))
    cfg['model_list'] = cfg['model_list_GM'] + cfg['model_list_WM']
        
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in subj_list:
        
        print('Modelling ' + subj + '...')
    
        # Extract data for subject
        subj_data      = scan_list[(scan_list['newstudyName'] == subj)].reset_index(drop=True)
        
        ######## SESSION-WISE OPERATIONS ########
        for sess in list(subj_data['blockNo'].unique()) :
          
            print('Working on session ' + str(sess) + '...')

            ######## MODEL-WISE OPERATIONS ########
            for model in cfg['model_list']:
                
                print('Working with ' + model + '...')

                if model=='Nexi':
                    data_used = 'allDelta-allb'
                elif model=='Sandi': # lowest diff time
                    filtered_data = subj_data[(subj_data['acqType'] == 'PGSE') & (subj_data['phaseDir'] == 'fwd') & (subj_data['blockNo'] == sess) & (subj_data['noBval'] > 1)]
                    ind_folder = getattr(filtered_data["diffTime"], 'idxmin')()
                    data_used = 'Delta_'+str(int(filtered_data['diffTime'][ind_folder]))+'_fwd'  
                elif model=='SMI' or model=='SMI_wSTE': # largest diff time
                    filtered_data = subj_data[(subj_data['acqType'] == 'PGSE') & (subj_data['phaseDir'] == 'fwd') & (subj_data['blockNo'] == sess) & (subj_data['noBval'] > 1)]
                    ind_folder = getattr(filtered_data["diffTime"], 'idxmax')()
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
                sigma       = copy_files_BIDS(bids_strc_prep,input_path,'dwi_dn_sigma.nii.gz')
                mask        = copy_files_BIDS(bids_strc_prep,input_path,'mask.nii.gz')

                # Get diffusion duration (assumes the same value for all acquisitions)
                #small_delta = np.loadtxt(small_delta)[0]
         
                # Modify units of bvals for NEXI          
                new_bvals = bvals.replace('.txt','_units.txt')
                modify_units_bvals(bvals, new_bvals )
        
                # Copy necessary files for analysis 
                if model=='Nexi':
                    bids_strc_lowb = create_bids_structure(subj=subj, sess=sess, datatype="dwi", description="allDelta-lowb", root=data_path, 
                                                folderlevel='derivatives', workingdir=cfg['prep_foldername'])
                    sigma     = copy_files_BIDS(bids_strc_lowb,input_path,'dwi_dn_sigma.nii.gz')
                elif model=='SMI':
                    input_file =  copy_files_BIDS(bids_strc_prep,input_path,'dwi_dn_gc_ec.mif').replace(data_path,docker_path)
                    others = '' 
                elif model=='SMI_wSTE':
                    bids_strc_STE = create_bids_structure(subj=subj, sess=sess, datatype="dwi_STE", root=data_path, 
                                                folderlevel='derivatives', workingdir=cfg['prep_foldername'],description='STE_fwd')
                    STE        = copy_files_BIDS(bids_strc_STE,input_path,'dwi_dn_gc_topup.mif').replace(data_path,docker_path)
                    LTE        = copy_files_BIDS(bids_strc_prep,input_path,'dwi_dn_gc_ec.mif').replace(data_path,docker_path)
                    input_file = LTE + ',' + STE
                    others     = '-echo_time 51,51 -bshape 1,0 -compartments EAS,IAS -debug'
           
                # Run SwissKnife models
                if model=='Nexi' or model=='Sandi':
                    
                    # Choose arguments for Nexi
                    args = [model, 
                            output_path, 
                            dwi,  
                            new_bvals,  
                            big_delta,  
                            small_delta, 
                            sigma,
                            '--debug']
        
                    # Run script
                    command = ["conda", "run", "-n", "SwissKnife", "python", os.path.join(cfg['code_path'], 'auxiliar_modelling.py')] + args  
                    subprocess.run(command, check=True)
        
                    
                # Run Designer models
                elif model=='SMI' or model=='SMI_wSTE':  
                    estimate_SMI_designer(input_file,
                                       mask.replace(data_path,docker_path), 
                                       sigma.replace(data_path,docker_path),
                                       output_path.replace(data_path,docker_path), 
                                       data_path.replace(data_path,docker_path),
                                       others)
                    
                    # Choose arguments for Nexi
                    args = [model, 
                            output_path.replace(data_path,docker_path), 
                            input_file,  
                            mask.replace(data_path,docker_path),  
                            sigma.replace(data_path,docker_path), 
                            data_path.replace(data_path,docker_path),
                            others,
                            '--debug']
                    
                    # Run script
                    command = ["conda", "run", "-n", "base", "python", os.path.join(cfg['code_path'], 'auxiliar_modelling.py')] + args  
                    subprocess.run(command, check=True)
                    
                    # # put ma~
                    # bids_strc_analysis.set_param(root=data_path)
                    # bids_strc_prep.set_param(root=data_path)

                    # output_path = bids_strc_analysis.get_path()
                    
                    
                    
                
                # Mask output for better visualization
                for filename in os.listdir(output_path):
                    if filename.endswith(".nii.gz"):
                        multiply_by_mask(os.path.join(output_path, filename), # filename input
                                         os.path.join(output_path,'Masked'), # output folder
                                         bids_strc_prep.get_path('mask.nii.gz')) # mask




# if __name__ == "__main__":
#     import json
#     import sys
#     import os

#     cfg_data_path = str(sys.argv[1])
#     with open(os.path.join(cfg_data_path, '.config.json')) as f:
#         cfg = json.load(f)

#     # Add code paths
#     sys.path.append(os.path.join(cfg['code_path'], 'processing_dwi'))

#     from Step4_modelling_GM import *

#     subj_list = cfg['subj_list']
#     Step4_modelling_GM(subj_list, cfg)
    
 