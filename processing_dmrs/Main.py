#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 14:41:37 2024

@author: localadmin
"""

import os
import sys
import matplotlib.pyplot as plt


plt.close('all');
os.system('clear')
os.system('cls')

# from custom_functions import QA_denoise
# importlib.reload(sys.modules['custom_functions'])
# from custom_functions import QA_denoise


############################## ADD CODE PATH ##############################
sys.path.append(os.path.join(os.path.expanduser('~'),  'Documents', 'Rita','Codes'))
sys.path.append(os.path.join(os.path.expanduser('~'),  'Documents', 'Rita','Codes','processing_dmrs'))

import importlib, sys
from custom_functions import *
from bids_structure import *
importlib.reload(sys.modules['custom_functions'])
importlib.reload(sys.modules['bids_structure'])


########################## DATA PATH AND SUBJECTS ##########################
subj_list = ['sub-01','sub-02','sub-03','sub-04']
subj_list = ['sub-01']

cfg                         = {}
cfg['data_path']            = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','CristinasTestData')
cfg['prep_foldername']      = 'preprocessed2'
cfg['analysis_foldername']  = 'analysis'
cfg['common_folder']        = os.path.join(os.path.expanduser('~'), 'Documents','Rita','Data','common')



#### STEP 0. CONVERT DATA >>> Use Base env
Step1_fill_study_excel(cfg)   ## Do once or if new data is added to the excel study file