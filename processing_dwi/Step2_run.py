#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Subscript to do step 2 of dMRI processing to be run in the dicomifier env
Last changed Jan 2025
@author: Malte O
"""
import os
import sys
sys.path.append(os.path.join(os.path.expanduser('~'),  'Documents', 'Projects','dMRS_starting_data_cristina','dMRSI','processing_dwi'))
sys.path.append(os.path.join(os.path.expanduser('~'),  'Documents', 'Projects','dMRS_starting_data_cristina','dMRSI'))

from Step2_raw2nii2bids import *
from Step2_correct_orientation import *

import json

os.system('clear')

if __name__ == "__main__":
    cfg_data_path = str(sys.argv[1])
    f = open(os.path.join(cfg_data_path, 'config.json'))
    cfg = json.load(f)
    f.close()

    f = open(os.path.join(cfg_data_path, 'subj_list.json'))
    subj_list = json.load(f)
    f.close()

    Step2_raw2nii2bids(subj_list, cfg)  # Do once for subject
    Step2_correct_orientation(subj_list, cfg)  # Do once for subject
