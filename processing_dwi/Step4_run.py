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


os.system('clear')

if __name__ == "__main__":
    cfg_data_path = str(sys.argv[1])
    f = open(os.path.join(cfg_data_path, '.config.json'))
    cfg = json.load(f)
    f.close()

    sys.path.append(cfg['code_path'])
    
    from Step4_modelling_GM import *
    from Step4_modelling_WM import *

    from bids_structure import *
    from custom_functions import *

    subj_list = cfg['subj_list']

    Step4_modelling_GM(subj_list,cfg) ### Do more than once if needed
    Step4_modelling_WM(subj_list,cfg) ### Do more than once if needed
