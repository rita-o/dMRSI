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
    has_microfa = ('wmicroFA' in model)
    has_mrs     = ('mrsinformed' in model) 
    model_clean  = model.split('_')[0] if (has_microfa or has_mrs) else model 
 
    if model_clean == 'Nexi' or model_clean =='Sandi' or model_clean =='Smex' or model_clean =='Sandix':
       
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

        # parameter limits  
        param_lims = None
        if 'Nexi' in model_clean:
            if extra == 'ex_vivo':
                param_lims = np.array(([1, 150], (0, 2), (0, 2), [0.1, 0.9]))
            elif extra == 'in_vivo':
                param_lims = np.array([[1, 150], [0.1, 3.5], [0.1, 3.5], [0.1, 0.9]])
        
        elif 'Sandix' in model_clean:
            if extra == 'ex_vivo':
                param_lims = np.array([[1, 150], [0.1, 2], [0.1, 2], [0.05, 0.95], [1, 30], [0.05, 0.5]])
            elif extra == 'in_vivo':
                param_lims = np.array([[1, 150], [0.1, 3.5], [0.1, 3.5], [0.05, 0.95], [1, 30], [0.05, 0.5]])
        
        elif 'Sandi' in model_clean:
            if extra == 'ex_vivo':
                param_lims = np.array([[0.1, 2], [0.1, 2], [0.05, 0.95], [1, 30], [0.05, 0.5]])
            elif extra == 'in_vivo':
                param_lims = np.array([[0.1, 3.5], [0.1, 3.5], [0.05, 0.95], [1, 30], [0.05, 0.5]])


        small_delta_val = float(np.loadtxt(delta_path)[0])
        est_kwargs = dict(
            model_name=model_clean,                
            dwi_path=dwi_path,
            bvals_path=bvals_path,
            delta_path=Delta_path,
            small_delta=small_delta_val,
            lowb_noisemap_path=sigma_path,
            out_path=out_path,
            mask_path=mask_path,
            adjust_parameter_limits=param_lims,
            debug=debug
            )
            
        # Only add uA_path if needed; 
        if has_microfa:
            uFA = sys.argv[10]
            est_kwargs['uA_path'] = uFA
            
        if has_mrs:
            mrs_radius_s =  float(sys.argv[10])
            est_kwargs['mrs_radius_s'] = mrs_radius_s
           
        # if has_mrs:
        #     est_kwargs['prior_rs_target'] = float(sys.argv[10])
        #     est_kwargs['prior_rs_sigma']  = float(sys.argv[11])
 
        # Estimate model
        estimate_model(**est_kwargs)
        
        # Calculate fneurite and fsoma from f and fs
        if 'Sandi' in model_clean or 'Sandix' in model_clean:
            import glob as glob
            import nibabel as nib
                
            map_f = glob.glob(os.path.join(out_path, "*f.nii.gz"))
            base_name = os.path.basename(map_f[0]).split('f.nii.gz')[0]
            img_f  = nib.load(map_f[0])
            hdr = img_f.header.copy()
            hdr.set_data_dtype(np.float32)
        
            map_f = nib.load(map_f[0]).get_fdata()
            
            map_fs = glob.glob(os.path.join(out_path, "*fs.nii.gz"))
            map_fs = nib.load(map_fs[0]).get_fdata()
        
            f_neurite = map_f*(1-map_fs)
            f_soma = map_f*map_fs
        
            ni_neurite = nib.Nifti1Image(f_neurite.astype(np.float32), affine=img_f.affine, header=hdr)
            ni_soma    = nib.Nifti1Image(f_soma.astype(np.float32),    affine=img_f.affine, header=hdr)
        
            nib.save(ni_neurite, os.path.join(out_path, f"{base_name}fneurite.nii.gz"))
            nib.save(ni_soma, os.path.join(out_path, f"{base_name}fsoma.nii.gz"))

      
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
         
    #elif model =='DTI_DKI_dipy' :
        
        # from dipy.core.gradients import gradient_table
        # import nibabel as nib
        # import dipy.reconst.dti as dti
        # import dipy.reconst.dki as dki
    
        # # Load data
        # out_path    = sys.argv[2]
        # dwi_path    = sys.argv[3]  
        # bvals_path  = sys.argv[4]
        # bvecs_path  = sys.argv[5]
        # mask_path   = sys.argv[6]
        # bvecs=np.loadtxt(bvecs_path)
        # bvals=np.loadtxt(bvals_path)
        # img = nib.load(dwi_path)
        # data = img.get_fdata()
        # mask_data = nib.load(mask_path).get_fdata() > 0
                
        # # Prepare gradients
        # mask_b = bvals <= 3000
        # gtab = gradient_table(bvals=bvals[mask_b], bvecs=bvecs[:,mask_b])
        
        # # Fit
        # print('Building DKI model...')
        # dkimodel = dki.DiffusionKurtosisModel(gtab,fit_method='OLS')
        # print('Fitting DKI...')
        # dkifit = dkimodel.fit(data[:,:,:,mask_b],mask=mask_data)

        # # Save maps
        # MD = dkifit.md*1e3 # units um^2/ms
        # MD_img = nib.Nifti1Image(MD.astype("float32"), affine=img.affine, header=img.header)
        # nib.save(MD_img, os.path.join(out_path,'MD.nii'))
        
        # MK = dkifit.mk()   # no units
        # MK_img = nib.Nifti1Image(MK.astype("float32"), affine=img.affine, header=img.header)
        # nib.save(MK_img, os.path.join(out_path,'MK.nii'))
        
        # FA = dkifit.fa   # no units
        # FA_img = nib.Nifti1Image(FA.astype("float32"), affine=img.affine, header=img.header)
        # nib.save(FA_img, os.path.join(out_path,'FA.nii'))
        
    
              
    
if __name__ == "__main__":
    Run_model()
    