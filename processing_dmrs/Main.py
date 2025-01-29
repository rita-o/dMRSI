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
sys.path.append(os.path.join(os.path.expanduser('~'),  'Documents', 'Projects','dMRS_starting_data_cristina','dMRSI','processing_dmrs'))
sys.path.append(os.path.join(os.path.expanduser('~'),  'Documents', 'Projects','dMRS_starting_data_cristina','dMRSI'))


import importlib, sys
from custom_functions import *
from bids_structure import *
importlib.reload(sys.modules['custom_functions'])
importlib.reload(sys.modules['bids_structure'])
from Step1_Fitting import *

########################## DATA PATH AND SUBJECTS ##########################
subj_list = ['sub-01','sub-02','sub-03']

cfg                         = {}
cfg['data_path']            = os.path.join(os.path.expanduser('~'), 'Documents','Projects','dMRS_starting_data_cristina','CristinasTestData')
cfg['prep_foldername']      = 'preprocessed'
cfg['analysis_foldername']  = 'analysis'
cfg['common_folder']        = os.path.join(os.path.expanduser('~'), 'Documents','Projects','dMRS_starting_data_cristina','common')



#### STEP 1. Fitting of data >>> Use fsl_mrs env
Step1_Fitting(subj_list, cfg)
