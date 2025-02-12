#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Subscript to do step 4 of dMRI processing to be run in the SwissKnife env
Last changed Jan 2025
@author: Malte O
"""
import os
import sys
import json

sys.path.append(os.path.join(os.path.expanduser('~'),  'Documents', 'Projects','dMRS_starting_data_cristina','dMRSI','processing_dwi'))
sys.path.append(os.path.join(os.path.expanduser('~'),  'Documents', 'Projects','dMRS_starting_data_cristina','dMRSI'))
from Step4_modelling_GM import *
from Step4_modelling_WM import *

os.system('clear')

if __name__ == "__main__":
    cfg_data_path = str(sys.argv[1])
    f = open(os.path.join(cfg_data_path, '.config.json'))
    cfg = json.load(f)
    f.close()

    subj_list = cfg['subj_list']

    Step4_modelling_GM(subj_list,cfg) ### Do more than once if needed
    Step4_modelling_WM(subj_list,cfg) ### Do more than once if needed
