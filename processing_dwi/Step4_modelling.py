"""
Script to analyse GM microstructure by fitting models such as Nexi, Sandi, ...
It uses the SwissKnife module installed in the SwissKnife python environment.

Last changed Jan 2025
@author: Rita O
"""
import sys
import os
import json
import pandas as pd
import platform
import math
import importlib, sys
from custom_functions import *
from bids_structure import *

def Step4_modelling(subj_list, cfg):
    
    data_path   = cfg['data_path']     
    scan_list   = pd.read_excel(os.path.join(data_path, cfg['scan_list_name'] ))
    cfg['model_list'] = cfg['model_list_GM'] + cfg['model_list_WM']
    # Define path to docker
    docker_path = '/data' 
    
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in subj_list:
        
        print('Modelling ' + subj + '...')
    
        # Extract data for subject
        subj_data      = scan_list[(scan_list['study_name'] == subj)].reset_index(drop=True)
        subj_data      = subj_data[subj_data['analyse'] == 'y']

        ######## SESSION-WISE OPERATIONS ########
        for sess in list(subj_data['sessNo'].unique()) :
          
            print('Working on session ' + str(sess) + '...')
            
            ########################## DELTA-WISE OPERATIONS ##########################      
            filtered_data = subj_data[(subj_data['phaseDir'] == 'fwd') & (subj_data['sessNo'] == sess) & (subj_data['noBval'] > 1) & (subj_data['acqType'] == 'PGSE')]
            Delta_list = filtered_data['diffTime'].unique().astype(int).tolist()
            
            for Delta in Delta_list:

                # Define bids structure for the processed data
                bids_strc_prep = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path, 
                                                folderlevel='derivatives', workingdir=cfg['prep_foldername'])
                #bids_strc_prep.set_param(description='Delta_' + str(Delta) + '_fwd') # repricated: done on each Delta processed separately
                bids_strc_prep.set_param(description='allDelta-allb/Delta_' + str(Delta)) # new: done on each Delta processed together ('allDelta')

                ######## Run DTI and DKI ########  
                 
                # Define BIDS structure for the analysis data
                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype="dwi", description=f'DTI_DKI_Delta_{Delta}', root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
              
                # Make output folder 
                output_path = bids_strc_analysis.get_path()
                
                # Decide which model to run
                run_DTIDKI = 1
                bvals = np.unique(read_numeric_txt(get_file_in_folder(bids_strc_prep, '*bvalsNom.txt'))[0])
                has_low  = np.any(bvals <= 1000)
                has_high = np.any(bvals > 1000)
                
                do_model = []
                if has_high:  # at least one shell above 1000
                    do_model.append('-DKI')
                if has_low:  # all bvals <= 1000
                    do_model.append('-DTI')
                if has_high==0 and has_low==0:  # no valid shells found
                    run_DTIDKI = 0
                    print('No valid shells found for neither DTI nor DKI!')
                do_model = ' '.join(do_model)

                # Just run model if it doesn't exist on the folder yet
                if not os.path.exists(output_path) or cfg['redo_modelling']:
                    
                    input_path = os.path.join(output_path,'inputs')
                    create_directory(input_path)
                    
                    # Designer option:
                    # Copy necessary files for analysis and rename the path to the docker path. New: done on each Delta processed together ('allDelta')
                    dwi   = copy_files_BIDS(bids_strc_prep,input_path, 'dwi_dn_gc_ec.mif').replace(data_path,docker_path)
                    mask  = copy_files_BIDS(bids_strc_prep,input_path,  'mask.nii.gz').replace(data_path,docker_path)
                    out_folder   = output_path.replace(data_path, docker_path)

                    # Run model
                    if run_DTIDKI==1:
                        if cfg['is_alive']=='ex_vivo':
                            extra = '-maxb 7'
                        else:
                            extra = ''
                            
                        call = [f'docker run -v {data_path}:/{docker_path} nyudiffusionmri/designer2:v2.0.13 tmi {do_model} {extra}',
                                        f'{dwi} {out_folder} -fit_constraints 0,1,0'] #-fit_constraints 0,1,0 -akc_outliers
                    
                        print(' '.join(call))
                        os.system(' '.join(call))
                        
                        # # Rename paths to local folder
                        bids_strc_analysis.set_param(root=data_path)
                        bids_strc_prep.set_param(root=data_path)
                    
                        output_path = bids_strc_analysis.get_path()
                         
                   
                if run_DTIDKI==1:
                    # Mask output with brain mask for better visualization
                    for filename in os.listdir(output_path):
                        if filename.endswith(".nii"):
                            multiply_by_mask(os.path.join(output_path, filename), # filename input
                                             os.path.join(output_path,'Output_masked'), # output folder
                                                     os.path.join(output_path,'inputs', 'mask.nii.gz')) # mask
                            
                    # Plot summary in dwi space
                    bids_strc_prep.set_param(description='allDelta-allb') # new: done on each Delta processed together ('allDelta')
                    plot_summary_params_model(os.path.join(output_path,'Output_masked'), 'DTI_DKI', cfg,bids_strc_prep.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'))
                    
                    # Register to anat space
                    bids_strc_anat = create_bids_structure(subj=subj, sess=sess, datatype="anat", root=data_path, 
                                              folderlevel='derivatives', workingdir=cfg['prep_foldername'])
                    register_outputfits_to_anat(os.path.join(output_path,'Output_masked'),
                                                os.path.join(output_path,'Output_in_anat'),
                                            'DTI_DKI',cfg, bids_strc_anat, bids_strc_prep)
                
                    # Plot summary plot in anat space
                    if cfg['subject_type']=='organoid':
                             bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas']  +'-To-'+cfg['anat_format'], root=data_path, 
                                                         folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                             bids_strc_reg.set_param(base_name='')
                             atlas=bids_strc_reg.get_path(f"atlas_in_{cfg['anat_format']}.nii.gz")
                             atlas_labels = prepare_atlas_labels(cfg['atlas'], glob.glob(os.path.join(bids_strc_anat.get_path(), '*label*'))[0])
                             mask = create_ROI_mask(atlas, atlas_labels, [], 'organoids', cfg['tpm_thr'], bids_strc_reg)
                             plot_summary_params_model(os.path.join(output_path,'Output_in_anat'), 'DTI_DKI', cfg, bids_strc_anat.get_path(f'{cfg['anat_format']}_bc_brain.nii.gz'), bids_strc_reg.get_path('mask_organoids.nii.gz'))
    
                    else:
                        plot_summary_params_model(os.path.join(output_path,'Output_in_anat'), 'DTI_DKI', cfg, bids_strc_anat.get_path(f'{cfg['anat_format']}_bc_brain.nii.gz'))

                
                ######## Compute PWD for LTE data ######## 
                
                # Create BIDS structures for LTE
                bids_LTE_temp = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=cfg['data_path'] , 
                             folderlevel='derivatives', workingdir=cfg['prep_foldername'],description=f'allDelta-allb/Delta_{Delta}')
                bids_LTE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=cfg['data_path'] , 
                             folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description=f'pwd_avg_Delta_{Delta}')
              
                # Create pwd average of LTE, just if it doesn't exist yet
                if not os.path.exists(bids_LTE.get_path()) or cfg['redo_modelling']:
                    create_directory(bids_LTE.get_path())
                    calculate_pwd_avg(get_file_in_folder(bids_LTE_temp,'*dwi_dn_gc_ec.nii.gz'),
                                      get_file_in_folder(bids_LTE_temp,'*bvalsNom.txt'),
                                      get_file_in_folder(bids_LTE_temp,'*bvalsEff.txt'),
                                      bids_LTE.get_path(),
                                      np.nan)
             

            ######## Compute PWD for STE data (if exists) ######## 
            
            # Create BIDS structures for STE
            bids_STE_temp = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                          folderlevel='derivatives', workingdir=cfg['prep_foldername'],description='STE_fwd')
            bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                          folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='pwd_avg')
            
            # Create pwd average of STE, just if it doesn't exist yet
            if os.path.exists(bids_STE_temp.get_path('dwi_dn_gc_topup.nii.gz')):                
                 
                  if not os.path.exists(bids_STE.get_path()) or cfg['redo_modelling']:
                      create_directory(bids_STE.get_path())
                      calculate_pwd_avg(bids_STE_temp.get_path('dwi_dn_gc_topup.nii.gz'),
                                        bids_STE_temp.get_path('bvalsNom.txt'),
                                        bids_STE_temp.get_path('bvalsEff.txt'),
                                        bids_STE.get_path(),
                                        np.nan)
                      
            ######## Compute PWD for STE data registered to LTE (if exists) ######## 
             
            # Create BIDS structures for STE
            bids_STE_temp = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                           folderlevel='derivatives', workingdir=cfg['prep_foldername'],description='STE_fwd')
            bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                           folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='pwd_avg_in_LTE')
            bids_STE_reg      = create_bids_structure(subj=subj, sess=sess, datatype='registration', root=cfg['data_path'] , 
                          folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='STE-To-LTE_'+ "allDelta-allb")
            bids_STE_reg.set_param(base_name='')
            
            # Create pwd average of STE, just if it doesn't exist yet
            if get_file_in_folder(bids_STE_reg,'*dn_gc_topup.nii.gz'):                
                  
                   if not os.path.exists(bids_STE.get_path()) or cfg['redo_modelling']:
                       create_directory(bids_STE.get_path())
                       calculate_pwd_avg(get_file_in_folder(bids_STE_reg,'*dn_gc_topup.nii.gz'),
                                         bids_STE_temp.get_path('bvalsNom.txt'),
                                         bids_STE_temp.get_path('bvalsEff.txt'),
                                         bids_STE.get_path(),
                                         np.nan)
                          
            
            ######## Compute MicroFA if data exists ########  
            if os.path.exists(bids_STE_temp.get_path('dwi_dn_gc_topup.nii.gz')):   
               
                # 1. Computing microFA in matlab 
                      # Define BIDS structure
                # bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                #               folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='microFA_matlab')
                # output_path = bids_STE.get_path()
                # bids_STE_reg      = create_bids_structure(subj=subj, sess=sess, datatype='registration', root=cfg['data_path'] , 
                #               folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='STE-To-LTE_'+ "allDelta-allb")
                # bids_STE_reg.set_param(base_name='')
                # bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                #               folderlevel='derivatives', workingdir=cfg['prep_foldername'],description='STE_fwd')
                # bids_LTE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=cfg['data_path'] , 
                #              folderlevel='derivatives', workingdir=cfg['prep_foldername'],description=f"allDelta-allb/Delta_{cfg['LTEDelta_for_microFA']}")
                # header        =  get_file_in_folder(bids_LTE,'*mask.nii.gz')
                
                     # Just run model if it doesn't exist on the folder yet
                # if not os.path.exists(output_path) or cfg['redo_modelling']:
                #     mdm_matlab(bids_LTE, bids_STE, bids_STE_reg, header, output_path, cfg['code_path2'], cfg['toolboxes'],low_b=True)

                # 2. Computing microFA in python
                    # Define BIDS structure
                bids_LTE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=cfg['data_path'] , 
                             folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description=f"pwd_avg_Delta_{cfg['LTEDelta_for_microFA']}")
                bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                               folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='pwd_avg_in_LTE')
                bids_LTE_aux  = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=cfg['data_path'] , 
                             folderlevel='derivatives', workingdir=cfg['prep_foldername'],description=f"allDelta-allb/Delta_{cfg['LTEDelta_for_microFA']}")
                mask          =  get_file_in_folder(bids_LTE_aux,'*mask.nii.gz')
                bids_STE_out      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                              folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='microFA')
                output_path = bids_STE_out.get_path()
                
                    # Just run model if it doesn't exist on the folder yet
                if not os.path.exists(output_path) or cfg['redo_modelling']:
                    compute_micro_FA(bids_LTE, bids_STE, mask, output_path)
                
                  
            ########################## MODEL-WISE OPERATIONS ##########################       
            for model in cfg['model_list']:
                
                print('Working with ' + model + '...')

                if 'Nexi' in model or 'Smex' in model or 'Sandix' in model:  
                    data_used = 'allDelta-allb'
                elif 'Sandi' in model: # lowest diff time
                    filtered_data = subj_data[(subj_data['acqType'] == 'PGSE') & (subj_data['phaseDir'] == 'fwd') & (subj_data['sessNo'] == sess) & (subj_data['noBval'] > 1)]
                    ind_folder = getattr(filtered_data["diffTime"], 'idxmin')()
                    #data_used = 'Delta_'+str(int(filtered_data['diffTime'][ind_folder]))+'_fwd'   # depricated
                    data_used = 'allDelta-allb/Delta_'+str(int(filtered_data['diffTime'][ind_folder]))   # new
                elif model=='SMI' or model=='SMI_wSTE': # largest diff time
                    filtered_data = subj_data[(subj_data['acqType'] == 'PGSE') & (subj_data['phaseDir'] == 'fwd') & (subj_data['sessNo'] == sess) & (subj_data['noBval'] > 1)]
                    ind_folder = getattr(filtered_data["diffTime"], 'idxmax')()
                    #data_used = 'Delta_'+str(int(filtered_data['diffTime'][ind_folder]))+'_fwd'   # depricated
                    data_used = 'allDelta-allb/Delta_'+str(int(filtered_data['diffTime'][ind_folder]))   # new

                print(f'Using data: {data_used}')
                # Define bids structure 
                bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype='dwi', root=data_path, 
                                            folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description=model)
                bids_strc_prep = create_bids_structure(subj=subj, sess=sess, datatype="dwi", description=data_used, root=data_path, 
                                          folderlevel='derivatives', workingdir=cfg['prep_foldername'])
                
                # Define paths   
                output_path = bids_strc_analysis.get_path()
                input_path = os.path.join(output_path,'inputs')

                # Just run model if it doesn't exist on the folder yet
                if not os.path.exists(output_path) or cfg['redo_modelling']:
                   
                    # Make output folder 
                    create_directory(input_path)

                    # Copy necessary files for analysis
                    dwi         = copy_files_BIDS(bids_strc_prep,input_path,'dwi_dn_gc_ec.nii.gz')
                    big_delta   = copy_files_BIDS(bids_strc_prep,input_path,'DiffTime.txt')
                    small_delta = copy_files_BIDS(bids_strc_prep,input_path,'DiffDuration.txt')
                    bvals       = copy_files_BIDS(bids_strc_prep,input_path,'bvalsNom.txt')
                    sigma       = copy_files_BIDS(bids_strc_prep,input_path,'dwi_dn_sigma.nii.gz')
                    mask        = copy_files_BIDS(bids_strc_prep,input_path,'mask_dil.nii.gz')
        
                    # Get diffusion duration (assumes the same value for all acquisitions)
                    #small_delta = np.loadtxt(small_delta)[0]
             
                    # Modify units of bvals for NEXI          
                    new_bvals = bvals.replace('.txt','_units.txt')
                    modify_units_bvals(bvals, new_bvals )
            
                    # Copy necessary files for analysis 
                    if 'Nexi'in model or 'Sandi' in model or 'Smex' in model or 'Sandix' in model:  
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
                    if 'Nexi'in model or 'Sandi' in model or 'Smex' in model or 'Sandix' in model:  
                        
                        # Define arguments 
                        args = [model, 
                                output_path, 
                                dwi,  
                                new_bvals, 
                                big_delta,  
                                small_delta, 
                                sigma,
                                mask,
                                cfg['is_alive'],
                                '--debug']
                                                  
                        # Run script
                        if "wmicroFA" in model:
                            bids_STE      = create_bids_structure(subj=subj, sess=sess, datatype='dwi_STE', root=cfg['data_path'] , 
                                          folderlevel='derivatives', workingdir=cfg['analysis_foldername'],description='microFA')
                            uFA = os.path.join(bids_STE.get_path(),'Uaniso.nii')
                            args.insert(-1, uFA)  
                            command = ["conda", "run", "-n", "SwissKnife_exp", "python", os.path.join(cfg['code_path'], 'auxiliar_modelling.py')] + args  
                            subprocess.run(command, check=True)
                            
                        elif "mrsinformed" in model:
                            bids_mrs      = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=f'dmrs-to-allDelta-allb', root=cfg['data_path'], 
                                           folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                            args[7] = get_file_in_folder(bids_mrs,'*voxel_mrs.nii.gz')

                            # mrs_radius_s = os.path.join(bids_mrs.get_path(),'mrs_radius_s.nii')
                            args.insert(-1, '10.5')  # sub-01
                            args.insert(-1, '10')  # sub-01
                            command = ["conda", "run", "-n", "SwissKnife_exp", "python", os.path.join(cfg['code_path'], 'auxiliar_modelling.py')] + args  
                            subprocess.run(command, check=True)
                        else:
                            command = ["conda", "run", "-n", "SwissKnife", "python", os.path.join(cfg['code_path'], 'auxiliar_modelling.py')] + args  
                            subprocess.run(command, check=True)
            
                        
                    # Run Designer models
                    elif model=='SMI' or model=='SMI_wSTE':  
                        # Define arguments 
                        args = [model, 
                                output_path.replace(data_path,docker_path), 
                                input_file,  
                                mask.replace(data_path,docker_path),  
                                sigma.replace(data_path,docker_path), 
                                data_path,
                                others]
                    
                        # Run script
                        command = ["conda", "run", "-n", "base", "python", os.path.join(cfg['code_path'], 'auxiliar_modelling.py')] + args  
                        subprocess.run(command, check=True)
                     

                # Mask output for better visualization
                patterns, lims, maximums = get_param_names_model(model,cfg['is_alive'])
                for filename in os.listdir(output_path):
                    if any(fnmatch.fnmatch(filename, pattern) for pattern in patterns):
                        multiply_by_mask(os.path.join(output_path, filename), # filename input
                                         os.path.join(output_path,'Output_masked'), # output folder
                                         os.path.join(output_path,'inputs', 'mask_dil.nii.gz')) # mask
                # Plot summary plot in dwi space
                bids_strc_prep.set_param(description='allDelta-allb') # new: done on each Delta processed together ('allDelta')
                plot_summary_params_model(os.path.join(output_path,'Output_masked'), model, cfg, bids_strc_prep.get_path('b0_dn_gc_ec_avg_bc_brain.nii.gz'))
                
                # Register to anat space
                bids_strc_anat = create_bids_structure(subj=subj, sess=sess, datatype="anat", root=data_path, 
                                          folderlevel='derivatives', workingdir=cfg['prep_foldername'])                
                register_outputfits_to_anat(os.path.join(output_path,'Output_masked'),
                                            os.path.join(output_path,'Output_in_anat'),
                                            model,cfg, bids_strc_anat, bids_strc_prep)
                
                # Plot summary plot in anat space
                if cfg['subject_type']=='organoid':
                         bids_strc_reg  = create_bids_structure(subj=subj, sess=sess, datatype='registration', description=cfg['atlas']  +'-To-'+cfg['anat_format'], root=data_path, 
                                                     folderlevel='derivatives', workingdir=cfg['analysis_foldername'])
                         bids_strc_reg.set_param(base_name='')
                         atlas=bids_strc_reg.get_path(f"atlas_in_{cfg['anat_format']}.nii.gz")
                         atlas_labels = prepare_atlas_labels(cfg['atlas'], glob.glob(os.path.join(bids_strc_anat.get_path(), '*label*'))[0])
                         mask = create_ROI_mask(atlas, atlas_labels, [], 'organoids', cfg['tpm_thr'], bids_strc_reg)
                         plot_summary_params_model(os.path.join(output_path,'Output_in_anat'), model, cfg, bids_strc_anat.get_path(f'{cfg['anat_format']}_bc_brain.nii.gz'), bids_strc_reg.get_path('mask_organoids.nii.gz'))

                else:
                    plot_summary_params_model(os.path.join(output_path,'Output_in_anat'), model, cfg, bids_strc_anat.get_path(f'{cfg['anat_format']}_bc_brain.nii.gz'))

           
            
            


