import os
import sys
import pandas as pd
import platform
import math
import importlib, sys
import numpy as np
from custom_functions import *
from bids_structure import *
import copy
import glob

def Step3_preproc(subj_list, cfg):
    
    data_path   = cfg['data_path']     
    scan_list   = pd.read_excel(os.path.join(data_path, 'ScanList.xlsx'))
    topupon     = cfg['do_topup'] 
    
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in subj_list:
    
        print('Preprocessing ' + subj + '...')
    
        # Extract data for subject
        subj_data      = scan_list[(scan_list['newstudyName'] == subj)].reset_index(drop=True)
        
        # List of scans for subject
        study_scanNo    = list(subj_data['scanNo'])
        
        # List of acquisition sessions
        sess_list    = [x for x in list(subj_data['blockNo'].unique()) if not math.isnan(x)] # clean NaNs
    
        ######## SESSION-WISE OPERATIONS ########
        for sess in sess_list:
            
            ###### ANAT FILE OPERATIONS ######
            bids_strc_anat = create_bids_structure(subj=subj, sess=sess, datatype="anat", root=data_path, 
                                        folderlevel='derivatives', workingdir=cfg['prep_foldername'])
            
            # BET
            if not os.path.exists(bids_strc_anat.get_path('T2w_brain.nii.gz')):
                #brain_extract_BREX(bids_strc_anat.get_path('T2w.nii.gz'))
                brain_extract_RATS(bids_strc_anat.get_path('T2w.nii.gz'))
                #make_mask(bids_strc_anat.get_path('T2w_brain.nii.gz'), bids_strc_anat.get_path('T2w_brain_mask.nii.gz'), 100)                
                QA_brain_extract(bids_strc_anat.get_path('T2w.nii.gz'),os.path.join(bids_strc_anat.get_path(),'QA_brain_extract'))
            
            ###### DIFF SCAN-WISE OPERATIONS ######
            bids_strc = create_bids_structure(subj=subj, sess=sess, datatype="dwi", root=data_path, 
                                        folderlevel='derivatives', workingdir=cfg['prep_foldername'])
            
            # Generate paths for fwd and rev acquisition types
            paths_dwi_fwd=[]; paths_b0_fwd=[]; paths_dwi_rev=[]; paths_b0_rev=[]
            paths_bvals_fwd =[];paths_bvecs_fwd =[]; paths_bvals_rev =[]; paths_bvecs_rev =[]
            dwi_indices = np.where(np.logical_and(np.array(subj_data['acqType']) == 'PGSE',np.array(subj_data['scanQA']) == 'ok'))[0] 
            for scn_ctr in dwi_indices:    
                bids_strc.set_param(description='E'+str(study_scanNo[scn_ctr]))
               
                if subj_data['phaseDir'][scn_ctr] == 'fwd' :
                    target_lists = (paths_dwi_fwd, paths_b0_fwd, paths_bvals_fwd, paths_bvecs_fwd)             
                elif  subj_data['phaseDir'][scn_ctr] == 'rev':
                    target_lists = (paths_dwi_rev, paths_b0_rev, paths_bvals_rev, paths_bvecs_rev)
                                  
                for lst, suffix in zip(target_lists, ['dwi.nii.gz', 'b0.nii.gz', 'bvalsNom.txt', 'bvecs.txt']):
                    lst.append(bids_strc.get_path(suffix))
           
            
            # Extract b0's for each scan and register brain mask with avg b0 
            masks_paths = []
            for kk in range(len(paths_dwi_fwd)):
                
                # AVG B0
                if not os.path.exists(paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg.nii.gz')):
                    extract_b0(    [paths_dwi_fwd[kk]], \
                                    [paths_dwi_fwd[kk].replace('dwi.nii.gz', 'bvecs.txt')], \
                                    [paths_dwi_fwd[kk].replace('dwi.nii.gz', 'bvalsNom.txt')], \
                                    [paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0.nii.gz')]    )
                    make_avg(3, [paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0.nii.gz')], [paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg.nii.gz')])
                
                # REGISTER  
                if not os.path.exists(paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg_mask.nii.gz')): 
                    antsreg(bids_strc_anat.get_path('T2w.nii.gz'), paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg.nii.gz'), paths_dwi_fwd[kk].replace('dwi.nii.gz', 'anat2dwi'))
                    ants_apply_transforms([bids_strc_anat.get_path('T2w.nii.gz'),bids_strc_anat.get_path('T2w_brain.nii.gz')],
                                          paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg.nii.gz'),
                                          [paths_dwi_fwd[kk].replace('dwi.nii.gz', 'T2w_in_dwi.nii.gz'),paths_dwi_fwd[kk].replace('dwi.nii.gz', 'T2w_brain_in_dwi.nii.gz')],
                                          paths_dwi_fwd[kk].replace('dwi.nii.gz', 'anat2dwi0GenericAffine.mat'),
                                          paths_dwi_fwd[kk].replace('dwi.nii.gz', 'anat2dwi1Warp.nii.gz'))                  
                    QA_reg(paths_dwi_fwd[kk].replace('dwi.nii.gz', 'T2w_brain_in_dwi.nii.gz'),paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg.nii.gz'),os.path.join(os.path.dirname(paths_dwi_fwd[kk]), 'QA_reg'))
                    make_mask(paths_dwi_fwd[kk].replace('dwi.nii.gz', 'T2w_brain_in_dwi.nii.gz'), paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg_mask.nii.gz'), 5e3)                
                
                masks_paths.append(paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg_mask.nii.gz'))

                # option 2: make brain mask from avg_b0 based on intensity values
                #threshold_image(paths_dwi_fwd[kk].replace('dwi.nii.gz', 'b0_avg.nii.gz'), paths_dwi_fwd[kk].replace('dwi.nii.gz', 'mask.nii.gz'), 4e3, 1e6)
                #masks_paths.append(paths_dwi_fwd[kk].replace('dwi.nii.gz', 'mask.nii.gz'))

            for kk in range(len(paths_dwi_rev)):
     
                # AVG B0
                if not os.path.exists(paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg.nii.gz')):
                    extract_b0(     [paths_dwi_rev[kk]], \
                                    [paths_dwi_rev[kk].replace('dwi.nii.gz', 'bvecs.txt')], \
                                    [paths_dwi_rev[kk].replace('dwi.nii.gz', 'bvalsNom.txt')], \
                                    [paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0.nii.gz')]    )
                    make_avg(3, [paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0.nii.gz')], [paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg.nii.gz')])
                 
                # REGISTER  
                if not os.path.exists(paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg_mask.nii.gz')): 
                    antsreg(bids_strc_anat.get_path('T2w.nii.gz'), paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg.nii.gz'), paths_dwi_rev[kk].replace('dwi.nii.gz', 'anat2dwi'))
                    ants_apply_transforms([bids_strc_anat.get_path('T2w.nii.gz'),bids_strc_anat.get_path('T2w_brain.nii.gz')],
                                          paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg.nii.gz'),
                                          [paths_dwi_rev[kk].replace('dwi.nii.gz', 'T2w_in_dwi.nii.gz'),paths_dwi_rev[kk].replace('dwi.nii.gz', 'T2w_brain_in_dwi.nii.gz')],
                                          paths_dwi_rev[kk].replace('dwi.nii.gz', 'anat2dwi0GenericAffine.mat'),
                                          paths_dwi_rev[kk].replace('dwi.nii.gz', 'anat2dwi1Warp.nii.gz'))                  
                    make_mask(paths_dwi_rev[kk].replace('dwi.nii.gz', 'T2w_brain_in_dwi.nii.gz'), paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg_mask.nii.gz'), 5e3)                
                
                masks_paths.append(paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg_mask.nii.gz'))

                # option 2: make brain mask from avg_b0 based on intensity values
                #threshold_image(paths_dwi_rev[kk].replace('dwi.nii.gz', 'b0_avg.nii.gz'), paths_dwi_rev[kk].replace('dwi.nii.gz', 'mask.nii.gz'), 4e3, 1e6)
                #masks_paths.append(paths_dwi_rev[kk].replace('dwi.nii.gz', 'mask.nii.gz'))
    
            ###### COMBINED SCANS OPERATIONS #####           
            if cfg['preproc_type'] == 'combined':
              
                # Generate combined output path
                bids_strc.set_param(description='allE-allb')
                create_directory(bids_strc.get_path())
                
                # COMBINE DATA
                if not os.path.exists(bids_strc.get_path('dwi.nii.gz')):
                    concat_niftis(paths_dwi_fwd,  bids_strc.get_path('dwi.nii.gz'), 'all')
                    concat_files(paths_bvals_fwd, bids_strc.get_path('bvalsNom.txt'))
                    concat_files(paths_bvecs_fwd, bids_strc.get_path('bvecs.txt')) 
                    concat_param(np.array(subj_data['diffTime'][subj_data['phaseDir'] == 'fwd']),paths_bvals_fwd,bids_strc.get_path('DiffTime.txt'))
                    concat_param(np.array(subj_data['diffDuration'][subj_data['phaseDir'] == 'fwd']),paths_bvals_fwd,bids_strc.get_path('DiffDuration.txt'))
    
                    # Extract only low b values
                    o_bids_strc = copy.deepcopy(bids_strc)
                    bids_strc.set_param(description='allE-lowb')
                    create_directory(bids_strc.get_path())
 
    
                    desiredbvals = [0, 1000]
                    old_dataset = {"dwi":   o_bids_strc.get_path('dwi.nii.gz'), 
                                    "bvals": o_bids_strc.get_path('bvalsNom.txt'),
                                    "bvecs": o_bids_strc.get_path('bvecs.txt')}
                    new_dataset = {"dwi":   bids_strc.get_path('dwi.nii.gz'), 
                                    "bvals": bids_strc.get_path('bvalsNom.txt'),
                                    "bvecs": bids_strc.get_path('bvecs.txt')}              
                   
                    dwi_extract(old_dataset, new_dataset,','.join(map(str, desiredbvals)))
                    param_extract(o_bids_strc.get_path('DiffTime.txt'), o_bids_strc.get_path('bvalsNom.txt'), desiredbvals, bids_strc.get_path('DiffTime.txt'))
                    param_extract(o_bids_strc.get_path('DiffDuration.txt'), o_bids_strc.get_path('bvalsNom.txt'), desiredbvals, bids_strc.get_path('DiffDuration.txt'))

                    # DENOISE low b values
                    if not os.path.exists(bids_strc.get_path('dwi_dn.nii.gz')):
                        denoise_vols_default_kernel(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('dwi_dn_sigma.nii.gz'))
                        #denoise_vols(bids_strc.get_path('dwi.nii.gz'), '11,11,11',bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('dwi_dn_sigma.nii.gz'))
                        calc_noise_floor(bids_strc.get_path('dwi_dn_sigma.nii.gz'), bids_strc.get_path('dwi_nf.nii.gz') ) # added by rita
                        calc_snr(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn_sigma.nii.gz'),bids_strc.get_path('dwi_snr.nii.gz'))
                        output_path = bids_strc.get_path();
                        QA_plotbvecs(bids_strc.get_path('bvecs.txt'), os.path.join(output_path, 'QA_acquisition'))
                        QA_denoise(bids_strc.get_path('dwi_dn_res.nii.gz'), os.path.join(output_path, 'QA_denoise'),-1000,1000)
                        QA_denoise(bids_strc.get_path('dwi_dn_sigma.nii.gz'), os.path.join(output_path, 'QA_denoise'),0,600)

                    # Generate combined output path
                    bids_strc.set_param(description='allE-allb')

            elif cfg['preproc_type'] == 'individual':
                
                # if individual, chose the folder with longest diffusion time
                filtered_data = subj_data[(subj_data['phaseDir'] == 'fwd') & (subj_data['noBval'] > 1)]
                ind_folder = filtered_data.loc[filtered_data["diffTime"].idxmax()]["scanNo"]
                bids_strc.set_param(description='E'+str(ind_folder))
                           
            # Create deformed mask
            union_niftis(masks_paths, bids_strc.get_path('mask_before_ec.nii.gz'))
            filter_clusters_by_size(bids_strc.get_path('mask_before_ec.nii.gz'), bids_strc.get_path('mask_before_ec.nii.gz'), 200)
            dilate_im(bids_strc.get_path('mask_before_ec.nii.gz'), bids_strc.get_path('mask_before_ec.nii.gz'), '1.5')

            # DENOISE
            if not os.path.exists(bids_strc.get_path('dwi_dn.nii.gz')):
                denoise_vols_default_kernel(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('dwi_dn_sigma.nii.gz'))
                #denoise_vols(bids_strc.get_path('dwi.nii.gz'), '11,11,11',bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('dwi_dn_sigma.nii.gz'))
                calc_noise_floor(bids_strc.get_path('dwi_dn_sigma.nii.gz'), bids_strc.get_path('dwi_nf.nii.gz') ) # added by rita
                calc_snr(bids_strc.get_path('dwi.nii.gz'), bids_strc.get_path('dwi_dn_sigma.nii.gz'),bids_strc.get_path('dwi_snr.nii.gz'))

            # GIBBS UNRINGING
            if not os.path.exists(bids_strc.get_path('dwi_dn_gc.nii.gz')):
                gibbs_corr(bids_strc.get_path('dwi_dn.nii.gz'), bids_strc.get_path('dwi_dn_gc.nii.gz'))
    
            # TOPUP
            if not os.path.exists(bids_strc.get_path('dwi_dn_gc_topup.nii.gz')) and topupon:
                topup_routine(paths_b0_fwd, paths_b0_rev, bids_strc.get_path('dwi_dn_gc.nii.gz'), bids_strc, os.path.join(data_path,'common', 'mycnf_fmri.cnf'))
            
            # EDDY
            if not os.path.exists(bids_strc.get_path('dwi_dn_gc_ec.nii.gz')):
                eddy_routine(bids_strc.get_path('dwi_dn_gc.nii.gz'), bids_strc.get_path('dwi_dn_gc_ec.nii.gz'), bids_strc.get_path('mask_before_ec.nii.gz'), bids_strc.get_path('bvalsNom.txt'), bids_strc.get_path('bvecs.txt'), topupon) # added the corr
                # convert_to_scans(os.path.join(output_path, base_name + dwistr + '_dn_gc_ec.nii.gz'), paths_dwi_fwd, '_gc_topup_eddy')
                extract_b0([bids_strc.get_path('dwi_dn_gc_ec.nii.gz')], [bids_strc.get_path('dwi_dn_gc_ec.eddy_rotated_bvecs')], \
                            [bids_strc.get_path('bvalsNom.txt')], [bids_strc.get_path('b0_dn_gc_ec.nii.gz')]) # check this image for motion
            
            # Generate non-deformed masks
            if not os.path.exists(bids_strc.get_path('mask.nii.gz')):

                make_avg(3, [bids_strc.get_path('b0_dn_gc_ec.nii.gz')], [bids_strc.get_path('b0_dn_gc_ec_avg.nii.gz')])
                #threshold_image(bids_strc.get_path('b0_dn_gc_ec_avg.nii.gz'), bids_strc.get_path('mask.nii.gz'), 1e4, 9e4)
                antsreg(bids_strc_anat.get_path('T2w.nii.gz'), bids_strc.get_path('b0_dn_gc_ec_avg.nii.gz'), bids_strc.get_path('anat2dwi'))
                ants_apply_transforms([bids_strc_anat.get_path('T2w.nii.gz'),bids_strc_anat.get_path('T2w_brain.nii.gz')],
                                      bids_strc.get_path('b0_dn_gc_ec_avg.nii.gz'),
                                      [bids_strc.get_path('T2w_in_dwi.nii.gz'),bids_strc.get_path('T2w_brain_in_dwi.nii.gz')],
                                      bids_strc.get_path('anat2dwi0GenericAffine.mat'),
                                      bids_strc.get_path('anat2dwi1Warp.nii.gz'))                  
                make_mask(bids_strc.get_path('T2w_brain_in_dwi.nii.gz'), bids_strc.get_path('mask.nii.gz'), 5e3)                
                #filter_clusters_by_size(bids_strc.get_path('mask.nii.gz'), bids_strc.get_path('mask.nii.gz'), 200)
                dilate_im(bids_strc.get_path('mask.nii.gz'), bids_strc.get_path('mask_dil.nii.gz'), '1')
    
    
            # Convert to mif in case
            nifti_to_mif(bids_strc.get_path('dwi_dn_gc_ec.nii.gz'), bids_strc.get_path('dwi_dn_gc_ec.eddy_rotated_bvecs'), bids_strc.get_path('bvalsNom.txt'), bids_strc.get_path('dwi_dn_gc_ec.mif'))

                  
            # Quality analysis
            output_path = bids_strc.get_path();
            QA_DTI_fit(bids_strc.get_path('dwi_dn_gc.nii.gz'), bids_strc.get_path('bvalsNom.txt'), bids_strc.get_path('bvecs.txt'), bids_strc.get_path('mask.nii.gz'), os.path.join(output_path, 'QA_dti_before_eddy'))
            QA_DTI_fit(bids_strc.get_path('dwi_dn_gc_ec.nii.gz'), bids_strc.get_path('bvalsNom.txt'), bids_strc.get_path('dwi_dn_gc_ec.eddy_rotated_bvecs'), bids_strc.get_path('mask.nii.gz'),os.path.join(output_path, 'QA_dti_after_eddy'))
            QA_plotSNR(bids_strc.get_path('dwi_snr.nii.gz'), bids_strc.get_path('dwi_nf.nii.gz'), bids_strc.get_path('mask.nii.gz'), bids_strc.get_path('bvalsNom.txt'),os.path.join(output_path, 'QA_acquisition'))
            QA_plotbvecs(bids_strc.get_path('bvecs.txt'), os.path.join(output_path, 'QA_acquisition'))
            QA_denoise(bids_strc.get_path('dwi_dn_res.nii.gz'), os.path.join(output_path, 'QA_denoise'),-1000,1000)
            QA_denoise(bids_strc.get_path('dwi_dn_sigma.nii.gz'), os.path.join(output_path, 'QA_denoise'),0,600)
            QA_eddy(bids_strc.get_path('mask.nii.gz'),bids_strc.get_path('mask_dil.nii.gz'), bids_strc.get_path('dwi_dn_gc.nii.gz'), bids_strc.get_path('dwi_dn_gc_ec.nii.gz'), os.path.join(output_path, 'QA_eddy'),bids_strc.get_path('bvalsNom.txt'),bids_strc)
                
                
            # show_waiting_window('Check masks!')
            scans_fwd = []; scans_rev =[];  paths_dwi_fwd = []; paths_dwi_rev =[]; paths_b0_fwd =[]; paths_b0_rev =[]; output_path =[] 
