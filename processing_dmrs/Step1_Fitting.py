#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to fit the dMRS data. There are two options:
 - simple 1D fitting of one b-value
 - dynamic 2D fitting 
 
 
Last changed Jan 2025
@author: Rita O
"""
import os.path

import fsl_mrs.dynamic as dyn
import numpy as np
import pickle
from fsl_mrs.utils import mrs_io
from fsl_mrs.utils import plotting as splot
from fsl_mrs.utils.misc import parse_metab_groups
from fsl_mrs.utils import fitting
import array as arr
import matplotlib.pyplot as plt
from custom_functions import *
from bids_structure import *
import pandas as pd
from fsl_mrs.utils import report
import datetime
from fsl_mrs.core import nifti_mrs as ntools
import pandas as pd
import math


def Step1_Fitting(subj_list, cfg):
    
    # Initial definitions 
    path_to_data    = cfg['data_path']
    scan_list       = pd.read_excel(os.path.join(path_to_data, 'ScanList.xlsx'))

    ######## SUBJECT-WISE OPERATIONS ########
    for subj in subj_list:
        
        print('Fitting MRS of ' + subj + '...')
    
        # Extract data for subject
        subj_data      = scan_list[(scan_list['newstudyName'] == subj)].reset_index(drop=True)

        ######## SESSION-WISE OPERATIONS ########
        for sess in list(subj_data['blockNo'].unique()):

            basis_filename = cfg['basis_filename']
            
            # Get the scan numbers for the water reference. Assumes there is only one
            water_reference_sequence_number = subj_data.loc[
                    (subj_data['acqType'] == 'SPECIAL') &
                    (subj_data['blockNo'] == sess) &
                    (subj_data['phaseDir'] == 'water'),
                    'scanNo'
                ].iloc[0]
            # Get the scan numbers for the metabolite data 
            metab_sequence_numbers =  subj_data.loc[
                    (subj_data['acqType'] == 'SPECIAL') &
                    (subj_data['blockNo'] == sess) &
                    (subj_data['phaseDir'] == 'metab'),
                    'scanNo'
                ].tolist()

            for seq_no in metab_sequence_numbers:
                print(f'Sequence {seq_no}')
                # Read data
                bids_strc = create_bids_structure(subj=subj, sess=sess, datatype='dmrs', root=path_to_data,
                                                  folderlevel='derivatives', workingdir='preprocessed', description=f"seq-{seq_no}")
                data_filename  = bids_strc.get_path('dmrs_processed.nii.gz')
                data           = mrs_io.read_FID(data_filename)
                data_to_fit      = data.mrs(basis_file=basis_filename)



                ## 1. Simple fit - first b value - option 1 (using python fsl)

                # # Create output path
                bids_strc.set_param(workingdir=cfg['analysis_foldername'], description=f'{seq_no}_single_fit')
                out_path    = bids_strc.get_path()
                if not os.path.exists(bids_strc.get_path()):
                    os.makedirs(bids_strc.get_path())

                # Run the fitting
                Fitargs = {'ppmlim': cfg['ppm_lim'],
                           'method': 'Newton',
                           'baseline': cfg['baseline'],
                           'metab_groups': parse_metab_groups(data_to_fit,  ['Mac']),
                           'model': cfg['model'] ,
                           }
                           #'x0': }

                data_to_fit.processForFitting()

                res = fitting.fit_FSLModel(data_to_fit,**Fitargs)

                # Save and build report
                create_directory(out_path)
                splot.plot_fit(data_to_fit, res, out=os.path.join(out_path,'single_fit.png'))
                report.create_svs_report(
                    data_to_fit,
                    res,
                    fidfile=' ',
                    filename=os.path.join(out_path,'report.html'),
                    h2ofile=' ',
                    basisfile=basis_filename,
                    date=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))

                for result in ['summary', 'concentrations', 'qc', 'parameters', 'concentrations-mh','parameters-mh']:
                    res.to_file(os.path.join(out_path, result + '.txt'), result)

                call= [f'fsl_mrs',
                           f'--data {data_filename}',
                           f'--basis {basis_filename}',
                           #f'--lorentzian',
                           f'--metab_groups "Mac"',
                           f'--output {out_path}',
                           f'--report',
                           f'--overwrite',
                           f'--free_shift',
                           f'--baseline "{cfg['baseline']}"',
                       ]

                # print(' '.join(call))
                # os.system(' '.join(call))

            ## 2. Dynamic fit
            print('Dynamic fitting...')
            bids_strc = create_bids_structure(subj=subj, sess=sess, datatype='dmrs', root=path_to_data,
                                              folderlevel='derivatives', workingdir='preprocessed', description='combined')
            diffusion_times = []
            for dmrs in os.listdir(bids_strc.get_path()):
                if 'TD' in dmrs:
                    diffusion_times.append(dmrs.split('TD_')[-1].split('_dmrs.nii.gz')[0])

            for diffusion_time in diffusion_times:
               
                data_filename = bids_strc.get_path(f'TD_{diffusion_time}_dmrs.nii.gz')
                data = mrs_io.read_FID(data_filename)
                dmrs_list = data.mrs(basis_file=basis_filename)
                bvals = data.hdr_ext._dim_info[0]['hdr']['b_value']['Value']

                bvals_filename = bids_strc.get_path('bvals')

                for diffusion_model in cfg['diffusion_models']:
                    # Create output path
                    bids_strc.set_param(workingdir=cfg['analysis_foldername'], description='dyn_fit_'+diffusion_model+'_TD_'+str(diffusion_time))
                    out_path    = bids_strc.get_path()
                    if not os.path.exists(out_path):
                        os.makedirs(out_path)

                    # Create FSL MRS config file path
                    mrs_dyn_config_filename = os.path.join(cfg['common_folder'],'mrs_dyn_config_multi.py')

                    # Check that the basis has the right phase/frequency convention
                    for mrs in dmrs_list:
                        mrs.check_Basis(repair=True)

                    # write b vals file
                    with open(bvals_filename, 'w') as file:
                        for element in bvals:
                            file.write(f"{element} ")


                    Fitargs = {'ppmlim': cfg['ppm_lim'],
                               'baseline': cfg['baseline'],
                               'metab_groups': parse_metab_groups(dmrs_list[0], 'Mac'),
                               'model': cfg['model'] }

                    for mrs in dmrs_list:
                        mrs.processForFitting()

                    dobj = dyn.dynMRS(
                            dmrs_list,
                            bvals,
                            config_file=mrs_dyn_config_filename,
                            rescale=True, # apparently has no impact on results oO
                            **Fitargs)

                    init = dobj.initialise(verbose=True)
                    dres = dobj.fit(init=init, verbose = True)

                    splot.plotly_dynMRS(dmrs_list, dres.reslist, dobj.time_var)

                    # Save and build report
                    create_directory(out_path)
                    dobj.save(out_path)#, save_dyn_obj=args.full_save)

                    splot.plot_fit(data_to_fit, res, out=os.path.join(out_path,f'TD_{diffusion_time}_dyn_fit.png'))
                    report.create_dynmrs_report(
                        dres,
                        fidfile=data_filename,
                        filename=os.path.join(out_path,f'TD_{diffusion_time}_report.html'),
                        basisfile=basis_filename,
                        configfile=mrs_dyn_config_filename,
                        tvarfiles=bvals_filename,
                        date=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))


                    # Save predicted FID
                    if cfg['save_fit']:
                        # Get the predicted fit from the results list.
                        pred_data = np.stack([reslist.pred for reslist in dres.reslist]).T

                        # Reapply the scaling factor to ensure prediction has the same overall scale
                        pred_data /= dmrs_list[0].scaling['FID']
                        # Shape as SVS data
                        pred_data = pred_data.reshape((1, 1, 1) + pred_data.shape)

                        # Create NIfTI-MRS
                        from fsl_mrs.core.nifti_mrs import create_nmrs
                        # If this is going to be merged don't worry about getting the affine right.
                        affine = data.voxToWorldMat
                        pred = create_nmrs.gen_nifti_mrs(
                            pred_data,
                            data.dwelltime,
                            data.spectrometer_frequency[0],
                            nucleus=data.nucleus[0],
                            dim_tags=data.dim_tags,
                            affine=affine)
                        pred.save(os.path.join(out_path , f'TD_{diffusion_time}_fit.nii.gz'))
