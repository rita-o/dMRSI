"""
Script to fill in an excel with information from th study and imaging data
It does not use a particular python environment.

Last changed Jan 2025
@author: Rita O
"""

import pandas as pd
import platform
import sys
import os
import re


def Step1_fill_study_excel(cfg):
    
    data_path       = cfg['data_path']      
    file_path       = os.path.join(data_path, 'ScanList.xlsx')
    list_methods    = pd.read_excel(file_path)
    
    for ii in range (list_methods.shape[0]):
    
        scan_path   = os.path.join(data_path, 'raw_data', list_methods['studyName'][ii], str(list_methods['scanNo'][ii]))
    
        # Extract the name of the sequence  
        with open(os.path.join(scan_path, 'acqp'), 'r') as f:
            for line in f:
                if '##$ACQ_protocol_name=' in line: # rita, before was '##$Method='
                    #seq_name = line.split(':')[1][:-2]
                    match=re.search(r'\<(.*?)\>',next(f)) # rita
                    seq_name=match.group(1) # rita
                    if seq_name == 'RARE':
                        list_methods.at[ii, 'acqSeq'] = 'T2_Turbo' + seq_name
                    elif seq_name == 'FieldMap':
                        list_methods.at[ii, 'acqSeq'] = 'B0Map-ADJ_B0MAP'
                    else: # for ihMT and the diffusion sequence the sequence name is ok                   
                        list_methods.at[ii, 'acqSeq'] = seq_name
        
        # Extract methods for the dMRI sequence
        if seq_name.find('DTI_EPI') != -1: # for the diffusion data # rita before was DtiEpi
    
            with open(os.path.join(scan_path, 'method'), 'r') as f: 
                for line in f:
    
                    if '##TITLE=' in line:
                        pv_version = line.split(',')[1]
                        list_methods.at[ii, 'PV'] = pv_version.strip()
    
                    if '##$PVM_DwGradSep=( 1 )' in line:
                        diff_time = next(f)
                        list_methods.at[ii, 'diffTime'] = float(diff_time)
    
                    if '##$EPI_YN_ForwardReverseMode=' in line:
                        rev_opt = line.split('=')[1]
                        if rev_opt.strip() == 'No':
                            list_methods.at[ii, 'phaseDir'] = 'fwd'
                        elif rev_opt.strip() == 'Yes':
                            list_methods.at[ii, 'phaseDir'] = 'rev'
    
                    if '##$PVM_DwNDiffDir=' in line:
                        no_dirs = line.split('=')[1]
                        list_methods.at[ii, 'noDirs'] = int(no_dirs)
    
                    if '##$PVM_DwNDiffExpEach=' in line:
                        no_bval = line.split('=')[1]
                        list_methods.at[ii, 'noBval'] = int(no_bval)
    
                    if '##$PVM_DwAoImages=' in line:
                        no_b0 = line.split('=')[1]
                        list_methods.at[ii, 'noB0'] = int(no_b0)
    
                    if '##$PVM_NRepetitions=' in line:
                        no_rep = line.split('=')[1]
                        list_methods.at[ii, 'NR'] = int(no_rep)
    
                    if '##$PVM_NAverages=' in line:
                        no_avg = line.split('=')[1]
                        list_methods.at[ii, 'NA'] = int(no_avg)
    
                    if '##$PVM_DummyScans=' in line:
                        no_dummy = line.split('=')[1]
                        list_methods.at[ii, 'noDummy'] = int(no_dummy)
    
                    if '##$PVM_EpiEchoSpacing=' in line:
                        SE = line.split('=')[1]
                        list_methods.at[ii, 'SE'] = float(SE)
    
                    if '##$PVM_EpiNEchoes=' in line:
                        n_echoes = line.split('=')[1]
                        list_methods.at[ii, 'EPIfact'] = int(n_echoes)
                    
                    if '##$PVM_DwGradDur=( 1 )' in line:
                        diff_dur = next(f)
                        list_methods.at[ii, 'diffDuration'] = float(diff_dur)
       
    
        else: # for all the other sequences, extract only the PV version
    
            with open(os.path.join(scan_path, 'method'), 'r') as f: 
                for line in f:
    
                    if '##TITLE=' in line:
                        pv_version = line.split(',')[1]
                        list_methods.at[ii, 'PV'] = pv_version.strip()
    
    df              = pd.DataFrame(list_methods)
    df.to_excel(file_path, index=False) 