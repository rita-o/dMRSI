#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to fit the dMRS data. There are two options:
 - simple 1D fitting of one b-value
 - dynamic 2D fitting 
 
 
Last changed Jan 2025
@author: Rita O
"""

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
    
    data_path   = cfg['data_path']
    scan_list   = pd.read_excel(os.path.join(data_path, 'ScanList.xlsx'))
    
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in subj_list:
        
        print('Fitting MRS of ' + subj + '...')
    
        # Extract data for subject
        subj_data      = scan_list[(scan_list['newstudyName'] == subj)].reset_index(drop=True)

        # List of acquisition sessions
        sess_list    = [x for x in list(subj_data['blockNo'].unique()) if not math.isnan(x)] # clean NaNs
    
        ######## SESSION-WISE OPERATIONS ########
        for sess in sess_list:
            
            # Read data
            bids_strc = create_bids_structure(subj=subj, sess=sess, datatype='dmrs', root=data_path,
                                                        folderlevel='derivatives', workingdir='preprocessed')
            basis_filename = os.path.join(cfg['common_folder'],'mrs_basis')
            #dyn_filename   = os.path.join(cfg['common_folder'],'mrs_dyn_config_multi.py')

            data_filename  = bids_strc.get_path('dmrs.nii.gz')
            data           = mrs_io.read_FID(data_filename)
            dmrs_list      = data.mrs(basis_file=basis_filename)
            bvals          = data.dynamic_hdr_vals()[-1].flatten().astype(float)

            bvals_filename = bids_strc.get_path('bvals')
            new_filename   = bids_strc.get_path('dmrs_lowb.nii.gz')


            ## 1. Simple fit - first b value - option 1 (using python fsl)
            
            # # Create output path
            bids_strc.set_param(workingdir=cfg['analysis_foldername'], description='single_fit')
            out_path    = bids_strc.get_path()

            # Run the fitting
            data_to_fit = dmrs_list[0]

            # Fitargs = {'ppmlim': (0.2, 5.0),
            #            'baseline_order': 1,
            #            'metab_groups': parse_metab_groups(data_to_fit, ['Ala', 'Asc','Asp', 'bHB','Cr', 'GABA', 'Glc', 'Gln', 'Glu', 'GPC','GSH','Ins', 'Lac', 'NAA', 'NAAG', 'PCho','PCr', 'PE', 'Scyllo', 'Tau']),
            #             'model': 'lorentzian'}
            Fitargs = {#'ppmlim': (0.2, 4.2),
                       #'method': 'Newton',
                       #'baseline_order': 4,
                       'metab_groups': parse_metab_groups(data_to_fit,  ['Mac']),
                       'model': 'lorentzian'}

            data_to_fit.processForFitting()#ppmlim=Fitargs['ppmlim']) # very important point!! If it's not done things go wrong

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


            ## 1. Simple fit - first b value - option 2 (using fsl mrs command line)
            
            # Create output path
            bids_strc.set_param(workingdir=cfg['analysis_foldername'], description='single_fit')
            out_path    = bids_strc.get_path()

            data_to_fit, _ = ntools.split(data, dimension='DIM_USER_0',index_or_indices=0)
            data_to_fit.set_dim_tag('DIM_USER_0', 'DIM_DYN')
            data_to_fit.save(new_filename)

            call= [f'fsl_mrs',
                       f'--data {new_filename}',
                       f'--basis {basis_filename}',
                       f'--lorentzian',
                       f'--metab_groups "Mac"',
                       f'--output {out_path}',
                       f'--report',
                       f'--overwrite']

            print(' '.join(call))
            os.system(' '.join(call))

            ## 2. Dynamic fit - option 2 (using fsl mrs command line)

            for diffusion_model in cfg['diffusion_models']:
                # Create output path
                bids_strc.set_param(workingdir=cfg['analysis_foldername'], description='dyn_fit_'+diffusion_model)
                out_path    = bids_strc.get_path()

                # Create FSL MRS config file path
                mrs_dyn_config_filename = os.path.join(cfg['common_folder'],'mrs_dyn_config_'+diffusion_model+'.py')

                # Check that the basis has the right phase/frequency convention
                for mrs in dmrs_list:
                    mrs.check_Basis(repair=True)

                # write b vals file
                with open(bvals_filename, 'w') as file:
                    for element in bvals:
                        file.write(f"{element} ")


                call= [f'fsl_dynmrs',
                        f'--data {data_filename}',
                        f'--basis {basis_filename}',
                        f'--dyn_config {mrs_dyn_config_filename}',
                        f'--time_variables {bvals_filename}',
                        f'--lorentzian',
                        f'--baseline_order 1',
                        f'--metab_groups "Mac"',
                        f'--output {out_path}',
                        f'--report',
                        f'--overwrite']

                #print(' '.join(call))
                #os.system(' '.join(call))


                # 2. Dynamic fit - option 1 (using python fsl)
                # to be improved

                # Create output path
                bids_strc.set_param(workingdir=cfg['analysis_foldername'], description='dyn_fit_'+diffusion_model)
                out_path    = bids_strc.get_path()

                # Check that the basis has the right phase/frequency convention
                for mrs in dmrs_list:
                    mrs.check_Basis(repair=True)

                Fitargs = {'ppmlim': (0.2, 5.0),
                           'baseline_order': 1,
                           'metab_groups': parse_metab_groups(dmrs_list, 'Mac'),
                            'model': 'lorentzian'}
                dobj = dyn.dynMRS(
                        dmrs_list,
                        bvals,
                        config_file={mrs_dyn_config_filename},
                        rescale=True,
                        **Fitargs)


                dres = dobj.fit()
               # _ = dres.plot_mapped()
                splot.plotly_dynMRS(dmrs_list, dres.reslist, dobj.time_var)
