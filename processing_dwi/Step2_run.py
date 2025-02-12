#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Subscript to do step 2 of dMRI processing to be run in the dicomifier env
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
    
    from Step2_raw2nii2bids import *
    from Step2_correct_orientation import *
    from bids_structure import *
    from custom_functions import *
    
    subj_list = cfg['subj_list']

    Step2_raw2nii2bids(subj_list, cfg)  # Do once for subject
    Step2_correct_orientation(subj_list, cfg)  # Do once for subject
