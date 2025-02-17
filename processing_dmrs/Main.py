#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main script to analyse dMRS data

Last changed Jan 2025
@author: Rita, Malte
"""
import os
import sys
import matplotlib.pyplot as plt

plt.close('all');
os.system('clear')
os.system('cls')

############################## ADD CODE PATH ##############################
dmrsi_path = os.path.join(os.path.expanduser('~'), 'Documents','Projects','dMRS_starting_data_cristina','dMRSI')
sys.path.append(dmrsi_path)
sys.path.append(os.path.join(dmrsi_path,'processing_dmrs'))


import importlib, sys
from custom_functions import *
from bids_structure import *
importlib.reload(sys.modules['custom_functions'])
importlib.reload(sys.modules['bids_structure'])
from Step1_Fitting import *

########################## DATA PATH AND SUBJECTS ##########################
subj_list = ['sub-01']#['sub-01','sub-02','sub-03']

cfg                         = {}
cfg['data_path']            = os.path.join(os.path.expanduser('~'), 'Documents','Projects','dMRS_starting_data_cristina','CristinasTestData')
cfg['prep_foldername']      = 'preprocessed'
cfg['analysis_foldername']  = 'analysis_all_models'
cfg['common_folder']        = os.path.join(dmrsi_path,'common')
cfg['scan_list_name']       = 'ScanList.xlsx'
cfg['atlas']                = 'Atlas_WHS_v4'
cfg['diffusion_models']     = []# ['callaghan']
cfg['ppm_lim']              = [0.2, 4.3]
cfg['baseline_order']       = 8

#### STEP 1. Fitting of data >>> Use fsl_mrs env
Step1_Fitting(subj_list, cfg)
