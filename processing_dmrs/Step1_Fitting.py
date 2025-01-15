#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 15:05:10 2025

@author: localadmin
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

def Step1_Fitting(subj_list, cfg):
    
    
    data_path   = cfg['data_path']     
    scan_list   = pd.read_excel(os.path.join(data_path, 'ScanList.xlsx'))
    
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in subj_list:
        
        print('Fitting MRS of ' + subj + '...')
    
        # Extract data for subject
        subj_data      = scan_list[(scan_list['newstudyName'] == subj)].reset_index(drop=True)
          
        ######## SESSION-WISE OPERATIONS ########
        for sess in 1:
            
            bids_strc = create_bids_structure(subj='sub-01', sess=1, datatype='dmrs', root=data_path, 
                                                        folderlevel='derivatives', workingdir='preprocessed')
            
            data = mrs_io.read_FID(bids_strc.get_path('dmrs.nii.gz'))
            dmrs_list = data.mrs(basis_file=os.path.join(cfg['common_folder'],'mrs_basis'))
            bvals = data.dynamic_hdr_vals()[-1].flatten().astype(float)


            for b in range(len(bvals)):
                a = dmrs_list[0].FID
                a = np.fft.fft(a); 
                plt.figure()
                plt.plot(a)
                plt.show()
                 
 
            ## 1. Simple fit - one b value
            data_to_fit = dmrs_list[0]
            Fitargs = {'ppmlim': (0.2, 5.0),
                       'baseline_order': 1,
                       'metab_groups': parse_metab_groups(data_to_fit, ['Ala', 'Asc','Asp', 'bHB','Cr', 'GABA', 'Glc', 'Gln', 'Glu', 'GPC','GSH','Ins', 'Lac', 'Mac', 'NAA', 'NAAG', 'PCho','PCr', 'PE', 'Scyllo', 'Tau']),
                       'model': 'lorentzian'}
            Fitargs = {'ppmlim': (0.2, 5.0),
                       'baseline_order': 1,
                       'metab_groups': parse_metab_groups(data_to_fit, ['NAA']),
                       'model': 'lorentzian'}
            
            # Run the fitting and save
            res = fitting.fit_FSLModel(data_to_fit,**Fitargs)
            bids_strc.set_param(workingdir='analysis', description='single_fit')
            out_path = bids_strc.get_path()
            create_directory(out_path)
            splot.plot_fit(data_to_fit, res, out=os.path.join(out_path,'single_fit.png'))


            ## 2. Dynamic fit
            # Check that the basis has the right phase/frequency convention
            for mrs in dmrs_list:
                mrs.check_Basis(repair=True)
                
            
            
            dobj = dyn.dynMRS(
                    dmrs_list,
                    bvals,
                    config_file='/home/localadmin/Documents/Rita/Data/common/mrs_dyn_config.py',
                    rescale=True,
                    **Fitargs)
            
            
            dres = dobj.fit()
            _ = dres.plot_mapped()
splot.plotly_dynMRS(dmrs_list, dres.reslist, dobj.time_var)
