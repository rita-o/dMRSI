#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 26 16:10:30 2025

@author: localadmin
"""

import os
import sys
import matplotlib.pyplot as plt

plt.close('all');
os.system('clear')

############################## ADD CODE PATH ##############################
dmrsi_path = os.path.join(os.path.expanduser('~'),  'Documents','Rita','Codes_GitHub','dMRSI')
sys.path.append(dmrsi_path)
sys.path.append(os.path.join(dmrsi_path,'processing_dmrs'))
sys.path.append(os.path.join(dmrsi_path,'common','nifti_mrs_from_raw'))


import importlib, sys
from custom_functions import *
from bids_structure import *
importlib.reload(sys.modules['custom_functions'])
importlib.reload(sys.modules['bids_structure'])
from Step0_convert_brukerraw2niimrs import *
from Step1_Fitting import *
from scipy.io import savemat
from collections import defaultdict
from scipy.io import loadmat


########################## DATA PATH AND SUBJECTS ##########################
subj_list = ['sub-01']#['sub-01','sub-02','sub-03']#

cfg                         = {}
cfg['data_path']            = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','dMRI_dMRS_Pilot_20250424')
#cfg['data_path']            = os.path.join('/media','localadmin','DATA','data','20250424')
cfg['fixed_phase_shift']    = 245
cfg['prep_foldername']      = 'preprocessed'
cfg['coil_combination_method'] = 'Bruker method' # 'FSL MRS'
cfg['analysis_foldername']  = 'analysis'
cfg['common_folder']        = os.path.join(dmrsi_path,'common')
cfg['basis_filename']       = os.path.join(cfg['common_folder'], 'mrs_basis')
cfg['diffusion_times']      = 'all'
cfg['atlas']                = 'Atlas_WHS_v4'
cfg['diffusion_models']     =  ['exp']#['biexp','callaghan','dki','exp'] # #[]
cfg['ppm_lim']              = [.2, 4.3]
cfg['baseline']             = 'poly, 0' # 'spline, moderate' #
cfg['save_fit']             = True
cfg['model']                = 'free_shift'
cfg['scan_list_name']       = 'ScanList_cut.xlsx'


data_path   = cfg['data_path']
scan_list   = pd.read_excel(os.path.join(data_path, cfg['scan_list_name']))

######## SUBJECT-WISE OPERATIONS ########
for subj in subj_list:
    
    print('Convert MRS of ' + subj + '...')

    # Extract data for subject
    subj_data      = scan_list[(scan_list['newstudyName'] == subj)].reset_index(drop=True)
        
    # List of acquisition sessions
    sess_list    = [x for x in list(subj_data['blockNo'].unique()) if not math.isnan(x)] # clean NaNs
    
    ######## SESSION-WISE OPERATIONS ########
    for sess in sess_list:


        # Index of scans for this session 
        study_indx = subj_data.index[
            (subj_data['blockNo'] == sess) & 
            (subj_data['acqType'] == 'SPECIAL') & 
            (subj_data['phaseDir'] == 'metab')
        ].tolist()       
        # subj_data_short = subj_data[
        #     (subj_data['blockNo'] == sess) &
        #     (subj_data['acqType'] == 'SPECIAL') &
        #     (subj_data['phaseDir'] == 'metab')
        # ].reset_index(drop=True)
        

        # Load FID data
        bids_strc = create_bids_structure(subj=subj, sess=sess, datatype='dmrs', root=data_path,
                                          folderlevel='derivatives', workingdir='preprocessed',description="combined")   
        dn_mat = loadmat(bids_strc.get_path('dmrs_allDelta-allb_FID_denoised.mat'))
        data_dn = dn_mat['data']  # shape: (num_bvals, num_deltas, num_avgs, num_fid_points)
        
        special_rows = subj_data[(subj_data['acqType'] == 'SPECIAL') & (subj_data['phaseDir'] == 'metab')]        
        bvals = np.unique((special_rows['bvals']))
        bvals = bvals[~np.isnan(bvals)]  
        Deltas = np.unique((special_rows['diffTime']))
        Deltas = Deltas[~np.isnan(Deltas)]  

        
        for i,b in enumerate(bvals):
            
            for j,d in enumerate(Deltas):

                # Read data
                match = special_rows[
                   (special_rows['bvals'] == b) & 
                   (special_rows['diffTime'] == d)]
               
                scan_no = match['scanNo'].values[0]
       
                bids_strc = create_bids_structure(subj=subj, sess=sess, datatype='dmrs', root=data_path,
                                                  folderlevel='derivatives', workingdir='preprocessed',description=f"seq-{scan_no}")   
                
                data_filename  = bids_strc.get_path('dmrs_comb_cc_shift_align_denoised.nii.gz')
                
                
                data_ex       = mrs_io.read_FID(bids_strc.get_path('dmrs_odd_cc_shift_align.nii.gz'))
                data_ex = data_ex.mrs()
                
                for shot in range(0,len(data_ex)):
                    data_ex[shot].FID = data_dn[i, j, shot, :]
                

                # Format data: (points, averages)
                fid_list = [shot.FID for shot in data_ex]  # list of 1D arrays

                # Convert to 2D array: (points, averages)
                dmrs_list = np.stack(fid_list, axis=1)  # shape: (fid_points, num_averages)


                #dmrs_list = np.stack(data_ex).T

                dmrs_list = dmrs_list.reshape(((1, 1, 1,) + dmrs_list.shape))
                
                # Create NIfTI-MRS object (make sure these values are correct for your acquisition)
                nifti_obj = ntools.create_nmrs.gen_nifti_mrs(
                    dmrs_list,
                    2.0e-04,             # in seconds (check your acquisition!)
                    4.0026e+08,       # in Hz (check your scanner frequency!)
                    dim_tags=['DIM_USER_0', None, None]
                )
                

                # nifti_obj.set_dim_tag(
                #    'DIM_USER_0',
                #    'DIM_USER_0',
                #    'b-value increment',
                #    {'b_value': {'Value': [b], 'Description': 'b-value in ms.μm^-2'}})

                 
                nifti_obj.save(data_filename)

    
        # # Now organize into 4D numpy array
        # all_bvals = sorted(set(subj_data['bvals'][study_indx]))
        # all_deltas = sorted(set(subj_data['diffTime'][study_indx]))

        # num_bvals = len(all_bvals)
        # num_deltas = len(all_deltas)
        
        # # Assume each FID is same length and each scan has same number of averages
        # example_key = next(iter(organized_data))
        # num_avgs = len(FID_list)
        # num_fid_points = len(FID_list[0])  # j.FID.shape[0]

        # # Initialize array
        # result_array = np.zeros((num_bvals, num_deltas, num_avgs, num_fid_points), dtype=np.complex64)
        # for (bval, Delta), avg_list in organized_data.items():
        #     b_idx = all_bvals.index(bval)
        #     d_idx = all_deltas.index(Delta)
        #     result_array[b_idx, d_idx, :, :] = np.squeeze(avg_list)
         
        # fig, ax = plt.subplots(1, 1)
        # delta_ctr=2
        # ax.plot(data_mrs[0].ppmAxis,np.fft.fft(np.mean(result_array[0,delta_ctr,:,:],axis=0).real),label='bvalindx=0')
        # ax.plot(data_mrs[0].ppmAxis,np.fft.fft(np.mean(result_array[1,delta_ctr,:,:],axis=0).real),label='bvalindx=1')
        # ax.plot(data_mrs[0].ppmAxis,np.fft.fft(np.mean(result_array[2,delta_ctr,:,:],axis=0).real),label='bvalindx=2')
        # ax.plot(data_mrs[0].ppmAxis,np.fft.fft(np.mean(result_array[3,delta_ctr,:,:],axis=0).real),label='bvalindx=3')
        # ax.plot(data_mrs[0].ppmAxis,np.fft.fft(np.mean(result_array[4,delta_ctr,:,:],axis=0).real),label='bvalindx=4')
        # ax.legend(loc='upper right',fontsize=8)
        # plt.show()
         
        # mrs_info = {
        #         'ppmAxis': data_mrs[0].ppmAxis,
        #         'centralFreq': data_mrs[0].centralFrequency,
        #         'dwelltime': data_mrs[0].dwellTime,
        #         'bvals': all_bvals,
        #         'deltas': all_deltas
            
        #     }

        # bids_strc.set_param(description='combined')
        # savemat(bids_strc.get_path(f'dmrs_allDelta-allb_FID.mat'), {'data': result_array})
        # savemat(bids_strc.get_path(f'dmrs_allDelta-allb_info.mat'), {'data_info': mrs_info})
