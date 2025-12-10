#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  1 15:14:45 2025

@author: localadmin
"""
import sys
import subprocess

def RodentSkullStripping():
    
    # Get arguments passed to the script
    input_file  = sys.argv[1]

    # Build the rbm command
    cmd = ["rbm", input_file, input_file.replace('.nii.gz', '_brain_mask.nii.gz')]

    # Run command 
    subprocess.run(cmd, check=True)
    

if __name__ == "__main__":
    RodentSkullStripping()
