#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  1 15:14:45 2025

@author: localadmin
"""
import sys
import numpy as np
import os

def Run_model():
   
    # Get arguments passed to the script
    model       = sys.argv[1] 

    if model == 'Nexi' or model =='Sandi' or model =='Smex':
       
        from graymatter_swissknife import estimate_model
        out_path    = sys.argv[2]
        dwi_path    = sys.argv[3]  
        bvals_path  = sys.argv[4]
        Delta_path  = sys.argv[5]
        delta_path  = sys.argv[6]
        sigma_path  = sys.argv[7]
        mask_path   = sys.argv[8]
        extra       = sys.argv[9]

        debug       = '--debug' in sys.argv  # Set debug flag if passed

        if extra=='ex_vivo':
            param_lims=np.array(([1, 150], (0, 2) , (0, 2), [0.1, 0.9]))
        else:
            param_lims=None
            
        estimate_model(
            model,
            dwi_path,
            bvals_path,
            Delta_path,
            np.loadtxt(delta_path)[0],
            sigma_path,
            out_path,
            mask_path=mask_path,
            adjust_parameter_limits=param_lims,
            debug=debug
        )
        
    elif model =='SMI' or model=='SMI_wSTE':
         out_path    = sys.argv[2]
         dwi_path    = sys.argv[3]  
         mask_path   = sys.argv[4]
         sigma_path  = sys.argv[5]
         data_path   = sys.argv[6]
         others      = sys.argv[7]
         
         call = [f'docker run -v {data_path}:/data nyudiffusionmri/designer2:v2.0.12 tmi -SMI',
                f'{others}',
                f'-sigma {sigma_path}',
                f'-mask {mask_path}',
                f'{dwi_path}',
                f'{out_path}']
    
         print(' '.join(call))
         os.system(' '.join(call))
        
         call = [f'docker run -v {data_path}:/data nyudiffusionmri/designer2:v2.0.10 chmod -R 777 {out_path}']
         print(' '.join(call))
         os.system(' '.join(call))
    
              
    
if __name__ == "__main__":
    Run_model()
    