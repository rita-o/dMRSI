import os
import re
import shutil
import numpy as np
import nibabel as nib
# import nilearn.image as nlrn
os.environ.pop("MPLBACKEND", None)
import matplotlib.pyplot as plt
import nibabel as nib
from more_itertools import locate
import sys
import tkinter as tk
from tkinter import ttk
from scipy.ndimage import label, find_objects
import shutil
import distinctipy
import subprocess
import json
import math
from scipy.ndimage import binary_opening, label
from itertools import groupby
import gzip
import tempfile
from scipy import stats
import pandas as pd
import glob as glob
import fnmatch

##### FILES AND SYSTEM OPERATIONS #####

def gunzip_file(input_path):
    output_path = input_path.replace('.nii.gz','.nii')  # removes the .gz extension

    with gzip.open(input_path, 'rb') as f_in:
        with open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return output_path

def gzip_file(input_path):
    output_path = input_path.replace('.nii','.nii.gz')  # removes the .gz extension

    with open(input_path, 'rb') as f_in:
        with gzip.open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return output_path

def generate_paths(main_path, list_scanNo):

    list_paths = []
    for items in list_scanNo:
        list_paths.append(os.path.join(main_path, str(items)))
    return (dict(zip(list_scanNo, list_paths)))


# def gen_bids_path(suffix):

#     if bids_strc.get("description"):
#         description = bids_strc["description"]
#         base_name = f'{bids_strc["subject"]}_{
#             bids_strc["session"]}_{description}_'
#     else:
#         description = ""  # Default to an empty string if description is None or empty
#         base_name = f'{bids_strc["subject"]}_{bids_strc["session"]}_'

#     out_path = os.path.join(bids_strc["root"],
#                             bids_strc["subject"],
#                             bids_strc["session"],
#                             bids_strc["datatype"],
#                             *(description,) if description else (),
#                             base_name + suffix)
#     return out_path


def extract_paths(dict_paths, list_scans):

    extracted_paths = []
    for key in list_scans:
        if key != 0:
            extracted_paths.append(dict_paths.get(key))
    return extracted_paths


def show_waiting_window(message):

    # Create the root window
    root = tk.Tk()
    root.title(message)

    # Set the size of the window
    root.geometry("300x100")

    # Add a label to the window
    label = ttk.Label(root, text="Please wait...")
    label.pack(pady=10)

    # Add an OK button to the window
    ok_button = ttk.Button(root, text="OK", command=root.destroy)
    ok_button.pack(pady=10)

    # This makes the window appear in the center of the screen
    root.eval('tk::PlaceWindow . center')

    # Run the application
    root.mainloop()


def replace_string(input_list_string, string_to_replace, replacing_string):

    output_list_string = []
    for ii in range(len(input_list_string)):
        output_list_string.append(input_list_string[ii].replace(
            string_to_replace, replacing_string))
    return output_list_string


# def create_directory(path):

#     if not os.path.exists(path):
#         os.makedirs(path)

def delete_directory_sudo(path):
    
    call = [f' sudo rm -r {path}']
    os.system(' '.join(call))
    
def copy_file(source_paths, destination_paths):

    for source, destination in zip(source_paths, destination_paths):
        shutil.copyfile(source, destination)

def copy_files_BIDS(bids_strc, output_path, filename):

     out_f  = os.path.join(output_path, filename)
     in_f   =  get_file_in_folder(bids_strc, '*'+filename)
     copy_files([in_f],[out_f])
    
     # new_path = os.path.join(output_path,bids_strc.get_param("base_name")+filename)
     # copy_file([bids_strc.get_path(filename)], 
     #           [new_path])
     
     return out_f

def get_file_in_folder(bids_strc, pattern):
    file = glob.glob(os.path.join(bids_strc.get_path(), pattern))[0]
    return file
    
def get_subdir(directory_path):
    subdirectories = [name for name in os.listdir(
        directory_path) if os.path.isdir(os.path.join(directory_path, name))]
    if len(subdirectories) > 0:
        return os.path.join(directory_path, subdirectories[0])
    return None


def get_filename(file_path):

    file_name = os.path.basename(file_path)
    file_name_no_extension = file_name.split('.')[0]

    return file_name_no_extension


def get_files_with_extension(folder_path, extension):
    files = []
    for file_name in os.listdir(folder_path):
        # Create the full path to the file
        file_path = os.path.join(folder_path, file_name)
        # Check if the path points to a file and has the desired extension
        if os.path.isfile(file_path) and file_name.endswith(extension):
            files.append(file_path)
    return files


def remove_file(input_path):

    call = [f'rm',
            f'{input_path}']

    os.system(' '.join(call))


def concat_files(list_files, output_path):

    with open(output_path, 'w') as fout:
        for lines in zip(*[open(fname).readlines() for fname in list_files]):
            lines = [k.rstrip('\n') for k in lines]
            fout.write(' '.join(lines) + '\n')

def unconcat_files(file_path, delta_path):
    from collections import defaultdict
    
    # Load and clean delta values
    deltas = [int(float(x)) for x in read_numeric_txt(delta_path)[0]]
    
    # Load and clean file values
    with open(file_path, 'r') as f:
        values = [[float(x) for x in line.strip().split()] for line in f if line.strip()]
    
    # Group indices by delta
    delta_groups = defaultdict(list)
    for idx, delta in enumerate(deltas):
        delta_groups[delta].append(idx)

    # Save each group to a new NIfTI file
    for delta, indices in delta_groups.items():

        delta_folder = os.path.join(os.path.join(os.path.dirname(file_path)), f'Delta_{delta}')
        os.makedirs(delta_folder, exist_ok=True)
        
        name = os.path.basename(file_path).replace('allDelta-allb', f'Delta_{delta}')
        output_path = os.path.join(delta_folder, f'{name}')
        if '.eddy_rotated_bvecs' in output_path:
            #output_path = output_path.replace('.eddy_rotated_bvecs', '_bvecs_rotated.txt')
            output_path = re.sub(r'(Delta_\d+_).*$', r'\1bvecsRotated.txt', output_path)
            
        extracted_lines = []
        for line in values:
            new_values = [line[i] for i in indices]
            extracted_lines.append(' '.join(map(str, new_values)))

        # Write to output
        with open(output_path, 'w') as f:
            f.write('\n'.join(extracted_lines))
        print(f"Saved: {output_path}")
        
def write_txt(values, output_path, mode):
    ''' mode is 'w' for write and 'a' for append '''
    with open(output_path, mode) as f:
        f.write(' '.join(map(str, values)))
    f.close()


def copy_files(source_paths, destination_paths):

    for source, destination in zip(source_paths, destination_paths):
        shutil.copyfile(source, destination)


def read_numeric_txt(input_path):

    with open(input_path, 'r') as file:
        data = file.read()
        values = [float(num) for num in data.split()]

    values = np.array(values)[np.newaxis, :]
    return values


def remove_folder(folder_path):

    call = [f'rm -r {folder_path}']
    os.system(' '.join(call))

def update_cfg(cfg):
    
    ordered_steps = [
        'redo_all',
        'redo_bet_anat',
        'redo_b0_extract',
        'redo_merge_dwi',
        'redo_denoise',
        'redo_gibbs',
        'redo_topup',
        'redo_eddy',
        'redo_final_mask']

    update_flag = False
    for step in ordered_steps:
        if cfg[step] == 1:
            update_flag = True
        if update_flag:
            cfg[step] = 1

    with open(cfg['data_path'] + '/.config.json', 'w') as f:
        json.dump(cfg, f)

    return cfg

def create_composite_page(images, titles=None, heading=None, subtext=None, page_size=(1240, 1754), margin=60,bottom_margin=180):
    from PIL import Image, ImageDraw, ImageFont
    """Compose a page with up to 2 images and an optional heading and individual image titles."""
    img_w, img_h = (page_size[0] - 2 * margin, (page_size[1] - 4 * margin) // 2)
    page = Image.new('RGB', page_size, color='white')
    draw = ImageDraw.Draw(page)

    try:
        font_heading = ImageFont.truetype("DejaVuSans-Bold.ttf", 40)
        font_title = ImageFont.truetype("DejaVuSans.ttf", 18)
        font_subtext = ImageFont.truetype("DejaVuSans.ttf", 22)
    except:
        font_heading = font_title = ImageFont.load_default()

    y = margin

    # Folder-level heading (top of first page per group)
    if heading:
        draw.text((margin, y), heading, fill='black', font=font_heading)
        y += 50

    if subtext:
        draw.text((margin, y), subtext, fill='gray', font=font_subtext)
        y += 35
        
    # Draw each image with its title above it
    for i, img in enumerate(images):
        title = titles[i] if titles else f"Image {i+1}"

        # Draw title
        draw.text((margin, y), title, fill='black', font=font_title)
        y += 35  # space after title

        # Resize and paste image
        resized = img.resize((img_w, img_h))
        page.paste(resized, (margin, y))
        y += img_h + margin

    return page

def make_summary_pdf(base_path, output_pdf):
    from PIL import Image
    folder_metadata = {
        'QA_acquisition': {
            'heading': 'Acquisition',
            'subtext': 'Check noise floor is lower than signal, specially for high bvals. \n Check bvecs is sampled on the whole sphere.'
        },
        'QA_denoise': {
            'heading': 'Denoising',
            'subtext': 'Check residuals do not have structure, check sigma map is uniform'
        },
        'QA_gc': {
            'heading': 'Gibbs unringing',
            'subtext': 'Check that rings were removed'
        },
        'QA_topup': {
            'heading': 'Topup',
            'subtext': 'Check that the image is corrected for distortions'
        },
        'QA_eddy': {
            'heading': 'Eddy',
            'subtext': 'Check registration quality and motion correction'
        },
        'QA_dti_before_eddy': {
            'heading': 'DTI before eddy',
            'subtext': ''
        },
        'QA_dti_after_eddy': {
            'heading': 'DTI after eddy',
            'subtext': 'Tensor reconstruction after eddy correction; should match pre-eddy'
        },
        'QA_mask': {
            'heading': 'Masks',
            'subtext': 'Check that mask after preprocessing is good anf fits well the brain'
        }
    }



    pages = []

    for folder_name, meta in folder_metadata.items():
        folder_path = os.path.join(base_path, folder_name)
        png_paths = sorted(glob.glob(os.path.join(folder_path, '*.png')))
        if not png_paths:
            continue
    
        imgs = [Image.open(p).convert('RGB') for p in png_paths]
        titles = [os.path.basename(p) for p in png_paths]
    
        for i in range(0, len(imgs), 2):
            chunk_imgs = imgs[i:i+2]
            chunk_titles = titles[i:i+2]
            page = create_composite_page(
                chunk_imgs,
                titles=chunk_titles,
                heading=meta['heading'] if i == 0 else None,
                subtext=meta['subtext'] if i == 0 else None
            )
            pages.append(page)
    
        if not pages:
            raise ValueError("No PNGs found in QA folders.")

    first_page, rest = pages[0], pages[1:]
    first_page.save(output_pdf, save_all=True, append_images=rest)
    print(f"Saved PDF to: {output_pdf}")

    
##### TOPUP #####


def topup_routine(dwi_path, bids_strc, topupcfg_path):

    topup_input_files = create_topup_input_files(bids_strc, topupcfg_path)
    do_topup(topup_input_files)
    apply_topup(topup_input_files, dwi_path, bids_strc)


def create_topup_input_files(bids_strc, topupcfg_path):

    topup_input_files = {}

    # Forward and reverse b0 images
    concat_niftis([bids_strc.get_path('b0_fwd.nii.gz'), bids_strc.get_path('b0_rev.nii.gz')],
                  bids_strc.get_path('b0_fwd_rev.nii.gz'), 'all')

    topup_input_files['b0_fwd_rev'] = bids_strc.get_path('b0_fwd_rev.nii.gz')

    if any(dim <= 10 for dim in nib.load(bids_strc.get_path('b0_fwd_rev.nii.gz')).shape[:3]):
        print('Your data is a slab and that is not good for topup and eddy, we are padding it with zeros')
        data = bids_strc.get_path('b0_fwd_rev.nii.gz')
        data_pad = data.replace('.nii.gz','_padded.nii.gz')
        pad_image(data, data_pad)
        topup_input_files['b0_fwd_rev'] = data_pad

    # Acqp file
    nii_fwd = nib.load(bids_strc.get_path('b0_fwd.nii.gz')).get_fdata()
    nii_rev = nib.load(bids_strc.get_path('b0_rev.nii.gz')).get_fdata()
    
    # ensure there is a 4th dimention
    nii_fwd = np.expand_dims(nii_fwd, axis=-1) if nii_fwd.ndim == 3 else nii_fwd
    nii_rev = np.expand_dims(nii_rev, axis=-1) if nii_rev.ndim == 3 else nii_rev

    no_fwd = nii_fwd.shape[-1]
    no_rev = nii_rev.shape[-1]

    for ii in range(no_fwd):
        if ii == 0:
            write_txt([0, 1, 0, 0.01, '\n'], bids_strc.get_path('acqp_topup.txt'), 'w')  # 0.087
        else:
            write_txt([0, 1, 0, 0.01, '\n'], bids_strc.get_path('acqp_topup.txt'), 'a')
    for ii in range(no_rev):
        write_txt([0, -1, 0, 0.01, '\n'],
                  bids_strc.get_path('acqp_topup.txt'), 'a')

    topup_input_files['acqp'] = bids_strc.get_path('acqp_topup.txt')

    # Configuration file - rita
    copy_files([topupcfg_path], [bids_strc.get_path('mycnf_topup.cnf')])
    topup_input_files['config'] = bids_strc.get_path('mycnf_topup.cnf')
    topup_input_files['output_pref'] = bids_strc.get_path('b0_topup')

    return topup_input_files


def do_topup(topup_input_files):  # {config} rita

    fwd_rev_nii = topup_input_files['b0_fwd_rev']
    acqp = topup_input_files['acqp']
    config = topup_input_files['config']
    output_pref = topup_input_files['output_pref']

    call = [f'topup',
            f'--imain={fwd_rev_nii}',
            f'--datain={acqp}',
            f'--config={config}',
            f'--out={output_pref}']

    print(' '.join(call))
    os.system(' '.join(call))
    

def apply_topup(topup_input_files, dwi_path, bids_strc):

    imain = dwi_path
    datain = topup_input_files['acqp']
    inindex = 1
    topup = bids_strc.get_path('b0_topup_fieldcoef')
    out = dwi_path.replace('.nii.gz', '_topup.nii.gz')

    if any(dim <= 10 for dim in nib.load(imain).shape[:3]):
         print('Your data is a slab and that is not good for topup and eddy, we are padding it with zeros')
         pad_image(imain, imain.replace('.nii.gz','_padded.nii.gz'))
         imain = imain.replace('.nii.gz','_padded.nii.gz')
     
    call = [f'applytopup ',
            f'--imain={imain}',
            f'--datain={datain}',
            f'--inindex={inindex}',
            f'--topup={topup}',
            f'--method=jac',
            f'--out={out}']

    print(' '.join(call))
    os.system(' '.join(call))
    
    if 'padded' in topup_input_files['b0_fwd_rev']:
        copy_file([out],[out.replace('.nii.gz','_padded.nii.gz')])
        unpad_image(out,out)


##### EDDY #####


def eddy_routine(dwi_path, out_path, mask_path, bvals_path, bvecs_path, topupon):

    eddy_input_files = create_eddy_input_files(
        dwi_path, out_path, mask_path, bvals_path, bvecs_path, topupon)
    do_eddy(eddy_input_files)


def create_eddy_input_files(dwi_path, out_path, mask_path, bvals_path, bvecs_path, topupon):

    eddy_input_files = {}

    basename = get_filename(bvals_path)
    basename = basename.replace('bvalsNom', '')

    nii_dwi = nib.load(dwi_path).get_fdata()
    no_dwi = nii_dwi.shape[-1]

    # Acqp file
    write_txt([0, 1, 0, 0.01], dwi_path.replace(get_filename(dwi_path) + '.nii.gz', basename +
              # last number is "dwell time" multiplied by "number of PE steps - 1"
                                                'acqp_eddy.txt'), 'w')
    eddy_input_files['acqp'] = dwi_path.replace(
        get_filename(dwi_path) + '.nii.gz', basename + 'acqp_eddy.txt')
    # Index file
    write_txt(['1'] * no_dwi, dwi_path.replace(get_filename(dwi_path) +
              '.nii.gz', basename + 'index_eddy.txt'), 'w')
    eddy_input_files['index'] = dwi_path.replace(
        get_filename(dwi_path) + '.nii.gz', basename + 'index_eddy.txt')
    # Bvals
    eddy_input_files['bvals'] = bvals_path
    # Bvecs
    eddy_input_files['bvecs'] = bvecs_path
    # dwi
    eddy_input_files['dwi'] = dwi_path

    # Mask
    if not os.path.exists(mask_path):
        mask_array = np.ones(nii_dwi.shape[0:3])
        array_to_nii(dwi_path, mask_array, dwi_path.replace(
            get_filename(dwi_path), basename + 'mask_before_ec'))
        eddy_input_files['mask'] = dwi_path.replace(
            get_filename(dwi_path), basename + 'mask_before_ec')
    else:
        eddy_input_files['mask'] = mask_path

    # adapted by rita to use eddy Quad as there is a .nii.gz uncessary
    eddy_input_files['out'] = out_path.replace('.nii.gz', '')

    # Check topup files exist
    if os.path.exists(dwi_path.replace(get_filename(dwi_path), basename + 'b0_topup_fieldcoef')) and topupon:
        eddy_input_files['topup'] = dwi_path.replace(
            get_filename(dwi_path) + '.nii.gz', basename + 'b0_topup')
    elif not os.path.exists(dwi_path.replace(get_filename(dwi_path), basename + 'b0_topup_fieldcoef')) and topupon:
        print(f"Error: The file at '{dwi_path.replace(basename, 'b0_topup_fieldcoef')}' does not exist.")
        sys.exit(1)

    return eddy_input_files


def do_eddy(eddy_input_files):  # rita addes repol and slm linear

    acqp = eddy_input_files['acqp']
    index = eddy_input_files['index']
    bvals = eddy_input_files['bvals']
    bvecs = eddy_input_files['bvecs']
    mask = eddy_input_files['mask']
    output = eddy_input_files['out']
    dwi = eddy_input_files['dwi']
        
    # If the image is like a slab (one dimension is very short), pad the image with zeros along 
    # that dimension otherwise eddy has problems and crashes.
    # I dont understand really deeply the cause of this problem but this seems
    # to be a good workaround
    if any(dim <= 10 for dim in nib.load(mask).shape[:3]):
        dwi_pad = dwi.replace('.nii.gz','_padded.nii.gz')
        mask_pad = mask.replace('.nii.gz','_padded.nii.gz')
        pad_image(mask, mask_pad)
        dwi = dwi_pad
        mask = mask_pad
        output_orig = output
        output = output+ '_padded'
        have_to_unpad = True
    else:
        have_to_unpad = False
      
    # Write call to eddy
    call = [f'eddy_cuda10.2',
            f'--imain={dwi}',
            f'--mask={mask}',
            f'--index={index}',
            f'--acqp={acqp}',
            f'--bvecs={bvecs}',
            f'--bvals={bvals}', \
            # f'--slm=linear',  #if data not acquired all sphere
            f'--out={output}', \
            f'--data_is_shelled --verbose']

    # If there is topup data, use it
    if eddy_input_files.get('topup'):
        topup = eddy_input_files['topup']
        call.insert(6, f'--topup={topup}')

    # Run eddy
    print(' '.join(call))
    os.system(' '.join(call))
    
    # Unpad the images for the next steps
    if have_to_unpad:
        unpad_image(output + '.nii.gz', output_orig +'.nii.gz')
        copy_files([output + '.eddy_rotated_bvecs'], [output_orig + '.eddy_rotated_bvecs'])
        

def pad_image(img_path, out_path, pad_before=5, pad_after=5):
    img = nib.load(img_path)
    data = img.get_fdata()
    affine = img.affine
    header = img.header
 
    shape = data.shape
    spatial_shape = shape[:3]
 
    # Determine the smallest dimension index (0=X, 1=Y, 2=Z)
    smallest_axis = np.argmin(spatial_shape)
 
    # Construct padding for each axis
    if data.ndim == 3:
        padding = [(0, 0), (0, 0), (0, 0)]
    elif data.ndim == 4:
        padding = [(0, 0), (0, 0), (0, 0), (0, 0)]
    else:
        raise ValueError(f"Unsupported image dimension: {data.ndim}")
 
    # Apply padding only to the smallest axis
    padding[smallest_axis] = (pad_before, pad_after)
 
    padded_data = np.pad(data, tuple(padding), mode='constant', constant_values=0)
    padded_img = nib.Nifti1Image(padded_data, affine, header)
    nib.save(padded_img, out_path)
      
    
def unpad_image(img_path, out_path, pad_before=5, pad_after=5):
    img = nib.load(img_path)
    data = img.get_fdata()
    affine = img.affine
    header = img.header

    shape = data.shape
    spatial_shape = shape[:3]

    # Determine the smallest spatial axis (0=X, 1=Y, 2=Z)
    smallest_axis = np.argmin(spatial_shape)

    # Define slicing to remove padding
    slices = [slice(None)] * data.ndim  # e.g., [slice(None), slice(None), slice(None)] or 4D
    slices[smallest_axis] = slice(pad_before, spatial_shape[smallest_axis] - pad_after)

    unpadded_data = data[tuple(slices)]
    unpadded_img = nib.Nifti1Image(unpadded_data, affine, header)
    nib.save(unpadded_img, out_path)
    
##### QA #####


def QA_DTI_fit(nifti_path, bvals_path, bvecs_path, mask_path, output_path):

    create_directory(output_path)
    out_base = os.path.join(output_path, 'dti')

    call = [f'dtifit',
            f'--data={nifti_path}',
            f'--mask={mask_path}',
            f'--bvecs={bvecs_path}',
            f'--bvals={bvals_path}',
            f'--out={out_base}',]

    os.system(' '.join(call))
    
    FA = os.path.join(output_path,'dti_FA.nii.gz')
    V1 = os.path.join(output_path,'dti_V1.nii.gz')
    png_path = os.path.join(output_path,'V1.png')
    
    # Define threhsold intensity and voxels to plot
    dim1    = int(np.ceil(nib.load(FA).shape[0]/2))
    dim2    = int(np.ceil(nib.load(FA).shape[1]/2))
    dim3    = int(np.ceil(nib.load(FA).shape[2]/2))
    
    call = [f'fsleyes render --hideCursor --hidex --hidez  --voxelLoc {dim1} {dim2} {dim3} ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30',
            f'--outfile {png_path}',
            f'{FA}',
            f'{V1} --overlayType rgbvector '
            f'--modulateImage {FA}']
    
    print(' '.join(call))
    os.system(' '.join(call))

def QA_brain_extract(anat_path,output_path,anat_format):
    
    create_directory(output_path)
    
    img     = nib.load(anat_path)
    slicee  = int(np.ceil(img.shape[1]/2))
    dim1    = int(np.ceil(img.shape[0]/2))
    dim3    = int(np.ceil(img.shape[2]/2))
    maxint  = int(round(np.ceil(np.max(img.get_fdata())),1))

    
    png_path = os.path.join(output_path, f'{anat_format}.png')
    call = [f'fsleyes render --hideCursor --hidex --hidez  --voxelLoc {dim1} {slicee} {dim3} ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{anat_path}',
            f'-dr 0 {maxint} ',]

    print(' '.join(call))
    os.system(' '.join(call))
    
    
    anat_brain_path = anat_path.replace('.nii.gz','_brain.nii.gz')
    png_path = os.path.join(output_path, f'{anat_format}_brain.png')
    call = [f'fsleyes render --hideCursor --hidex --hidez  --voxelLoc {dim1} {slicee} {dim3} ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{anat_brain_path}',
            f'-dr 0 {maxint} ',]

    print(' '.join(call))
    os.system(' '.join(call))
    
    # make countour mask
    anat_brain_mask_path     = anat_path.replace('.nii.gz','_brain_mask.nii.gz')
    anat_brain_mask_dil_path = anat_brain_mask_path.replace('.nii.gz','_dil.nii.gz');
    dilate_im(anat_brain_mask_path, anat_brain_mask_dil_path, '1')
    countour_path = anat_brain_path.replace('.nii.gz', '_contour.nii.gz')
    
    call = [f'fslmaths',
            f'{anat_brain_mask_dil_path}',
            f'-add',
            f'{anat_brain_mask_path}',
            f'-uthr 1',
            f'{countour_path}']
    print(' '.join(call))
    os.system(' '.join(call))


    png_path = os.path.join(output_path, f'{anat_format}_with_{anat_format}brain.png')
    call = [f'fsleyes render --hideCursor --hidex --hidez  --voxelLoc {dim1} {slicee} {dim3} ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{anat_path}',
            f'-dr 0 {maxint} ',
            f'{countour_path}']

    print(' '.join(call))
    os.system(' '.join(call))
  

def QA_mask(dwi, mask, mask_dil, mask_orig, output_path):
    
    create_directory(output_path)
    
    # Define threhsold intensity and voxels to plot
    dim1    = int(np.ceil(nib.load(mask).shape[0]/2))
    dim2    = int(np.ceil(nib.load(mask).shape[1]/2))
    dim3    = int(np.ceil(nib.load(mask).shape[2]/2))
    
    create_directory(output_path)
    
    png_path = os.path.join(output_path, 'dwi_mask_afterproc.png')
    call = [f'fsleyes render --hideCursor --hidex --voxelLoc {dim1} {dim2} {dim3} ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{mask} ',
            f'-dr 0 1']

    print(' '.join(call))
    os.system(' '.join(call))
    
    png_path = os.path.join(output_path, 'dwi_mask_beforeproc.png')
    call = [f'fsleyes render --hideCursor --hidex   --voxelLoc {dim1} {dim2} {dim3} ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{mask_orig} ',
            f'-dr 0 1']

    print(' '.join(call))
    os.system(' '.join(call))
    
  
    # make countour mask
    countour = os.path.basename(mask).replace('.nii.gz','_countour.nii.gz')
    call = [f'fslmaths',
            f'{mask_dil}',
            f'-add',
            f'{mask}',
            f'-uthr 1',
            f'{countour}']
    print(' '.join(call))
    os.system(' '.join(call))


    png_path = os.path.join(output_path, 'dwi_with_mask_afterproc.png')
    maxint  = int(round(np.ceil(np.max(nib.load(dwi).get_fdata())),1))
    call = [f'fsleyes render --hideCursor --hidex  --voxelLoc {dim1} {dim2} {dim3} ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{dwi}',
            f'-dr 0 {maxint} ',
            f'{countour}']

    print(' '.join(call))
    os.system(' '.join(call))
    
def QA_denoise(bids_strc, res, sigma, output_path):
    
    res_path    = bids_strc.get_path(res)
    sigma_path  = bids_strc.get_path(sigma)
    
    # Define threhsold intensity and voxels to plot
    maxintsigma = int(round(0.9 * np.max(np.max(nib.load(sigma_path).get_fdata()))))
    maxintres   = int(round(0.9 * np.max(np.max(nib.load(res_path).get_fdata()))))
    dim1    = int(np.ceil(nib.load(res_path).shape[0]/2))
    dim2    = int(np.ceil(nib.load(res_path).shape[1]/2))
    dim3    = int(np.ceil(nib.load(res_path).shape[2]/2))
    
    create_directory(output_path)
    
    png_path = os.path.join(output_path, res.replace('.nii.gz','.png'))
    call = [f'fsleyes render --hideCursor --hidex --hidez  --voxelLoc {dim1} {dim2} {dim3} ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{res_path} ',
            f'-dr -{maxintres} {maxintres}']

    print(' '.join(call))
    os.system(' '.join(call))
    
    png_path = os.path.join(output_path, sigma.replace('.nii.gz','.png'))
    call = [f'fsleyes render --hideCursor --hidex --hidez  --voxelLoc {dim1} {dim2} {dim3} ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{sigma_path} --cmap red-yellow',
            f'-dr 0 {maxintsigma}']

    print(' '.join(call))
    os.system(' '.join(call))
    

def QA_topup(bids_strc, before, after, output_path):

    
    before_path = bids_strc.get_path(before)
    after_path = bids_strc.get_path(after)
    
    # Define threhsold intensity and voxels to plot
    maxint  = int(round(0.8 * np.max(np.max(nib.load(before_path).get_fdata()))))
    dim1    = int(np.ceil(nib.load(before_path).shape[0]/2))
    dim2    = int(np.ceil(nib.load(before_path).shape[1]/2))
    dim3    = int(np.ceil(nib.load(before_path).shape[2]/2))
    
    create_directory(output_path)
    
    png_path = os.path.join(output_path, before.replace('.nii.gz','.png'))
    call = [f'fsleyes render --hideCursor --hidex  --voxelLoc {dim1} {dim2} {dim3}',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{before_path} ',
            f'-dr 0 {maxint} ']

    print(' '.join(call))
    os.system(' '.join(call))
    
    png_path = os.path.join(output_path, after.replace('.nii.gz','.png'))
    call = [f'fsleyes render --hideCursor --hidex  --voxelLoc {dim1} {dim2} {dim3}',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{after_path} ',
            f'-dr  0 {maxint} ']

    print(' '.join(call))
    os.system(' '.join(call))
    
    
def QA_gc(bids_strc, before, after, output_path):
    
    before_path = bids_strc.get_path(before)
    after_path = bids_strc.get_path(after)
    
    maxint = int(round(0.8* np.max(np.max(nib.load(before_path).get_fdata()))))
    dim1    = int(np.ceil(nib.load(before_path).shape[0]/2))
    dim2    = int(np.ceil(nib.load(before_path).shape[1]/2))
    dim3    = int(np.ceil(nib.load(before_path).shape[2]/2))
    
    create_directory(output_path)
    
    png_path = os.path.join(output_path, before.replace('.nii.gz','.png'))
    call = [f'fsleyes render --hideCursor --hidex   --voxelLoc {dim1} {dim2} {dim3} ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{before_path}',
            f'-dr 0 {maxint} ']

    print(' '.join(call))
    os.system(' '.join(call))
    
    png_path = os.path.join(output_path, after.replace('.nii.gz','.png'))
    call = [f'fsleyes render --hideCursor --hidex --voxelLoc  {dim1} {dim2} {dim3}',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{after_path} ',
            f'-dr 0 {maxint} ']

    print(' '.join(call))
    os.system(' '.join(call))
     
   

def QA_eddy(mask_path, mask_dil_path, dwi_path, dwi_ec_path, output_path, bvals_path,bids_strc):

    create_directory(output_path)

    # make countour mask
    countour_path = mask_path.replace('.nii.gz', '_bo_contour.nii.gz')
    call = [f'fslmaths',
            f'{mask_dil_path}',
            f'-add',
            f'{mask_path}',
            f'-uthr 1',
            f'{countour_path}']
    print(' '.join(call))
    os.system(' '.join(call))

    # retreive bvals = b0 position and chose those volumes
    bvals = np.loadtxt(bvals_path)
    b0s = np.where(bvals == 0)
    maxint = int(round(0.8* np.max(np.max(nib.load(dwi_path).get_fdata()))))
    dim1    = int(np.ceil(nib.load(mask_path).shape[0]/2))
    dim2    = int(np.ceil(nib.load(mask_path).shape[1]/2))
    dim3    = int(np.ceil(nib.load(mask_path).shape[2]/2))

    for ii in np.linspace(1, len(b0s[0])-1, num=3, dtype='int'):
        volume = b0s[0][ii]

        # plot dwi before eddy
        png_path = os.path.join(
            output_path, 'nodifcontour_v' + str(volume) + '_dwi.png')
        call = [f'fsleyes render --hideCursor --hidex  --voxelLoc {dim1} {dim2} {dim3} ',
                f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
                f'--outfile {png_path}',
                f'{dwi_path}',
                f'-v {volume}',
                f'-dr 0 {maxint} ',
                f'{countour_path}']

        print(' '.join(call))
        os.system(' '.join(call))

        # plot dwi after eddy
        png_path = os.path.join(
            output_path, 'nodifcontour_v' + str(volume) + '_eddy.png')
        call = [f'fsleyes render --hideCursor --hidex --voxelLoc {dim1} {dim2} {dim3}',
                f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
                f'--outfile {png_path}',
                f'{dwi_ec_path}',
                f'-v {volume}',
                f'-dr 0 {maxint}',
                f'{countour_path}']
        print(' '.join(call))
        os.system(' '.join(call))
          
    output_qc = os.path.join(output_path,'qc')
    if os.path.exists(output_qc):
        shutil.rmtree(output_qc)
    
    indx_eddy = bids_strc.get_path('index_eddy.txt')
    acqc_eddy = bids_strc.get_path('acqp_eddy.txt')
    suffix    = dwi_ec_path.replace('.nii.gz', '')

    call = [f'eddy_quad {suffix}', \
            f'-idx {indx_eddy} -par {acqc_eddy}', \
            f'-m {mask_path} -b {bvals_path}', \
            f'-o {output_qc}']
    print(' '.join(call))
    os.system(' '.join(call))
        

def QA_reg(moving_path, fixed_path, output_path):

    create_directory(output_path)
    os.chdir(output_path)
    call = [f'slicer ',
            f'{fixed_path}',
            f'{moving_path}',
            f'-s 2 -y 0.35 slb.png -y 0.45 slc.png -y 0.55 sld.png -y 0.65 sle.png']

    os.system(' '.join(call))


    out_image = os.path.join(output_path,'Mov2Fixed.png')
    call = [f'pngappend slb.png + slc.png + sld.png + sle.png',
            f'{out_image}; rm sl*.png']
    os.system(' '.join(call))


def plot_bvals(bids_strc):

    bvals_nom = read_numeric_txt(bids_strc.get_path('bvalsNom.txt'))
    bvals_eff = read_numeric_txt(bids_strc.get_path('bvalsEff.txt'))
    bvals_avg = calc_avg_bval(bvals_nom, bvals_eff)

    bvals_nom_unique = np.unique(bvals_nom)

    # QA plot: effective bvals as a function of nominal bvals
    fig, ax = plt.subplots()
    # Plot the data and save figure
    ax.plot(bvals_nom, bvals_eff, 'bo', markersize=8)
    # Plot only the mean values
    ax.plot(bvals_nom_unique, bvals_avg, 'ro', markersize=8)
    # Annotate the mean values
    for ii in range(len(bvals_avg)):
        ax.annotate(str(round(bvals_avg[ii], 1)), xy=(bvals_nom_unique[ii], bvals_avg[ii]),
                    xytext=(bvals_nom_unique[ii] - 100.0, bvals_avg[ii] + 50.0))
    # Set axes
    ax.set_xlabel('Nominal b-val',
                  fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})
    ax.set_ylabel('Effective b-val',
                  fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})
    ax.grid(True)
    plt.savefig(bids_strc.get_path('bvals_Eff_Nom.png'))
                #bbox_inches='tight', dpi=300)
    plt.close()


def QA_plotbvecs(bvec_path, bval_path, output_path):

    create_directory(output_path)

    # Load the bvecs file (assuming it is a whitespace-delimited text file)
    bvecs = np.loadtxt(bvec_path)
    bvals = np.loadtxt(bval_path)
    unique_bvals = np.unique(bvals)
    
    ndir_per_shell = []
    for i in range(len(unique_bvals)):
        a=list(bvals)
        ndir_per_shell.append(a.count(unique_bvals[i]))
    
    # normalize to the outter shell
    bvals_norm = np.multiply(bvals,1/max(bvals))
    bvecs_norm = np.multiply(bvecs,np.sqrt(bvals_norm))
    
    # create colors
    color_list = distinctipy.get_colors(len(unique_bvals),pastel_factor=0.5)
    
    # Create a 3D plot
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the data
    k=0
    for b in unique_bvals:
        bvecs_to_plot = bvecs_norm[:,bvals==b]
        ax.scatter(bvecs_to_plot[0, :], bvecs_to_plot[1, :], bvecs_to_plot[2, :], color=color_list[k], marker='*')
        k=k+1

    # Set axis limits
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])

    # Set axis properties for better visualization
    ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio
    plt.savefig(os.path.join(output_path, 'QA_bvecs.png'),
                bbox_inches='tight', dpi=300)


def QA_plotSNR(bids_strc, dwi, snr, dwi_sigma, mask, bvals, output_path):

    # Load data
    bvals_nom  = read_numeric_txt(bids_strc.get_path(bvals))
    mask       = nib.load(bids_strc.get_path(mask)).get_fdata()

    create_directory(output_path)

    # Compute noise floor
    nf        = nib.load(bids_strc.get_path(dwi_sigma)).get_fdata()*np.sqrt(np.pi/2)
    nf_masked = nf * mask
    nf_masked = nf.reshape(nf.shape[0]*nf.shape[1]*nf.shape[2], 1)
    nf_masked = nf_masked[~(np.isnan(nf_masked).any(axis=1) | (nf_masked == 0).any(axis=1))]

    # Plot SNR, masked by whole brain mask
    SNR = nib.load(bids_strc.get_path(snr)).get_fdata()
    for v in range(SNR.shape[-1]):
        SNR[:, :, :, v] = np.multiply(SNR[:, :, :, v], mask)

    data = SNR.reshape(SNR.shape[0]*SNR.shape[1]*SNR.shape[2], SNR.shape[3]);
    data[data == 0] = np.nan

    fig, ax = plt.subplots()
    ax.plot(np.transpose(bvals_nom), np.nanmean(data, axis=0), 'bo', markersize=3)
    ax.set_xlabel('Nominal b-val',
                  fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})
    ax.set_ylabel('Mean SNR', fontdict={
                  'size': 12, 'weight': 'bold', 'style': 'italic'})
    ax.grid(True)
    plt.savefig(os.path.join(output_path, 'QA_SNR.png'),
                bbox_inches='tight', dpi=300)
    
    
    # Plot S/S0, masked by whole brain mask
    dwi = nib.load(bids_strc.get_path(dwi)).get_fdata()
    
    for v in range(dwi.shape[-1]):
        dwi[:, :, :, v] = np.multiply(dwi[:, :, :, v], mask)

    data = dwi.reshape(dwi.shape[0]*dwi.shape[1]*dwi.shape[2], dwi.shape[3]);
    data[data == 0] = np.nan

    fig, ax = plt.subplots()
    ax.plot(np.transpose(bvals_nom), np.nanmean(data, axis=0), 'bo', markersize=3)
    ax.plot(np.transpose(bvals_nom), np.repeat(np.nanmean(nf_masked), np.transpose(bvals_nom).shape[0]))
    ax.set_xlabel('Nominal b-val',
                  fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})
    ax.set_ylabel('Mean S and noise floor', fontdict={
                  'size': 12, 'weight': 'bold', 'style': 'italic'})
    ax.grid(True)
    plt.savefig(os.path.join(output_path, 'QA_S_nf.png'),
                bbox_inches='tight', dpi=300)
    
def QA_ROIs(roi_paths, template, output_path):
   
    out_image = os.path.join(output_path, 'ROIs.png')
    output_txt_path = os.path.join(output_path, 'mask_ROIs.txt')
    output_mask_path = os.path.join(output_path, 'mask_ROIs.nii.gz')
    
    # Make mask all ROIs
    ref_img = nib.load(roi_paths[0])
    data_shape = ref_img.shape
    affine = ref_img.affine
    header = ref_img.header
    merged_mask = np.zeros(data_shape, dtype=np.uint8)
    label = 1
    label_mapping = []
    for roi_path in roi_paths:
        if 'WB' in roi_path:
            continue  # skip unwanted masks
    
        roi_img = nib.load(roi_path)
        roi_data = roi_img.get_fdata()
    
        # Binarize the ROI and label it (overwrite label only where it's 0)
        merged_mask[(roi_data > 0) & (merged_mask == 0)] = label
        roi_name = os.path.splitext(os.path.basename(roi_path))[0]
        label_mapping.append(f"{label}: {roi_name}")
        label += 1
        
    with open(output_txt_path, 'w') as f:
        f.write("\n".join(label_mapping))
    
    # Save merged mask
    merged_mask_img = nib.Nifti1Image(merged_mask, affine, header)
    nib.save(merged_mask_img, output_mask_path)


    # Make fsl eyes call
    call = [
        'fsleyes', 'render',
        '--scene', 'lightbox',
        '--outfile', out_image ,
        '--hideCursor',
        '--hidex',
        '--hidez',
        '--xcentre', '0', '0',
        '--ycentre', '0', '0',
        '--zcentre', '0', '0',
        '--labelSize', '30',
        '--zrange','0.38','0.65',
        '--sliceSpacing','0.02',
        '--sampleSlices','centre',
        template,                
        '--name', 'Base',
    ]
    call.extend([
                os.path.join(output_path, 'mask_ROIs.nii.gz'),
                '--name', 'ROI',
                '--overlayType', 'volume',
                '--cmap', 'HSV',
                '--displayRange', '0', '12',
                '--volume', '0',
                '--alpha', '70.0'
            ])

    
    print(' '.join(call))
    subprocess.run(call)
    
def plot_summary_params_model(output_path, model, cfg, template_path=None, countour_path=None):
    
    
   # Make colorbar
   import matplotlib.colors as mcolors
   import matplotlib.cm as cm
   import imutils

   jet = cm.get_cmap('jet', 256)
   if model == 'Nexi' or model =='Smex':
       jet_colors = jet(np.linspace(0, 1, 256))
       fade_len = 20
       fade = np.linspace(0, 1, fade_len).reshape(-1, 1)
       jet_colors[:fade_len, :3] *= fade  # Keep alpha (4th channel) unchanged
       jet_colors[-fade_len:, :3] *= fade[::-1]  # Reverse fade for the end
       custom_jet_black = mcolors.ListedColormap(jet_colors)
   else:
       jet_colors = jet(np.linspace(0, 1, 256))
       fade_len = 1
       fade = np.linspace(0, 1, fade_len).reshape(-1, 1)
       jet_colors[:fade_len, :3] *= fade  # Keep alpha (4th channel) unchanged
       custom_jet_black = mcolors.ListedColormap(jet_colors)

   
   # Get model parameter names and display ranges
   patterns, lims, maximums = get_param_names_model(model,cfg['is_alive'])
   
   # Load contour data
   if countour_path is not None:
       countour_data = nib.load(countour_path).get_fdata()
       if cfg['subject_type'] =='rat' :
            slicee = int(np.ceil(nib.load(countour_path).shape[1]/2))
            contour = imutils.rotate(countour_data[:,slicee, :], angle=90)
            fact = int((max(contour.shape) - min(contour.shape)) / 2)
            if  fact != 0:
                contour = contour[fact:-fact, :]
       elif cfg['subject_type'] =='human' :
            slicee = int(np.ceil(nib.load(countour_path).shape[2]/2))
            contour = imutils.rotate(countour_data[:,:, slicee], angle=90)
       elif cfg['subject_type'] =='organoid':
            slicee = int(np.ceil(nib.load(countour_path).shape[1]/2))
            contour = imutils.rotate(countour_data[:,slicee,:], angle=90)
       contour[np.isnan(contour)] = 0

   # Create subplot grid
   n_params = len(patterns)
   n_rows = 1 if n_params <= 4 else 2
   n_cols = math.ceil(n_params / n_rows)

   if template_path is not None:
       n_cols = math.ceil((n_params +1 )/ n_rows)
           
   fig, axs = plt.subplots(n_rows, n_cols, figsize=(12, 3.6))
   axs = axs.flatten()
   fig.subplots_adjust(wspace=0.05, hspace=0.11, top=0.92, bottom=0.1, left=0.05, right=0.95)
   
   # Display each parameter slice
   for ax, pattern, lim in zip(axs, patterns, lims):
       
       # Load file
       matched_file = glob.glob(os.path.join(output_path, pattern))
       if pattern=='*sandi*f*':
           matched_file = [
                f for f in matched_file
                if 'fs' not in os.path.basename(f).lower()] 
       param_path = matched_file[0]
       param_data = nib.load(param_path).get_fdata()
   
       # Extract and process middle slice
       if cfg['subject_type'] =='rat' :
           slicee = int(np.ceil(nib.load(param_path).shape[1]/2))
           img = imutils.rotate(param_data[:,slicee, :], angle=90)
           fact = int((max(img.shape) - min(img.shape)) / 2)
           if  fact != 0:
               img = img[fact:-fact, :]
           img[np.isnan(img)] = 0
       elif cfg['subject_type'] =='human' :
           slicee = int(np.ceil(nib.load(param_path).shape[2]/2))
           img = imutils.rotate(param_data[:,:, slicee], angle=90)
           img[np.isnan(img)] = 0
       elif cfg['subject_type'] =='organoid':
           slicee = int(np.ceil(nib.load(param_path).shape[1]/2))
           img = imutils.rotate(param_data[:,slicee,:], angle=90)
           img[np.isnan(img)] = 0
   
       # Show slice
       #if model in ['Nexi', 'Smex']:
         #  pattern = original_pattern
       im = ax.imshow(img, cmap=custom_jet_black, vmin=lim[0], vmax=lim[1])
       cleaned_pattern = re.sub(model, '', pattern, flags=re.IGNORECASE)
       cleaned_pattern = cleaned_pattern.replace('[^s]', '')
       cleaned_pattern = re.sub(r'\*{2,}', '*', cleaned_pattern).strip('*')
       ax.set_title(cleaned_pattern)
       ax.axis('off')
       
       # overlay countour
       if countour_path is not None:
           ax.contour(contour, levels=[0.5], colors='white', linewidths=1)
           #overlay = np.ma.masked_where(contour == 0, contour)
           #ax.imshow(overlay, cmap='gray', alpha=0.4)
   
       # Add colorbar
       cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.02)
       cbar.set_ticks([lim[0], lim[1]])
       cbar.ax.set_xticklabels([lim[0], lim[1]], rotation=-45)
       cbar.ax.tick_params(labelsize=10)
      
     
   # Add template
   if template_path is not None:
       template_data = nib.load(template_path).get_fdata()
       if cfg['subject_type'] =='rat' :
           slicee = int(np.ceil(nib.load(template_path).shape[1]/2))
           img = imutils.rotate(template_data[:,slicee, :], angle=90)
           fact = int((max(img.shape) - min(img.shape)) / 2)
           if  fact != 0:
               img = img[fact:-fact, :]
           img[np.isnan(img)] = 0
           maxint = int(np.round(0.9*np.ceil(np.max(img))))
       elif cfg['subject_type'] =='human' :
           slicee = int(np.ceil(nib.load(template_path).shape[2]/2))
           img = imutils.rotate(template_data[:,:, slicee], angle=90)
           img[np.isnan(img)] = 0
           maxint = int(np.round(0.9*np.ceil(np.max(img))))
       elif cfg['subject_type'] =='organoid':
           slicee = int(np.ceil(nib.load(template_path).shape[1]/2))
           img = imutils.rotate(template_data[:,slicee,:], angle=90)
           img[np.isnan(img)] = 0
           maxint = int(np.round(0.9*np.ceil(np.max(img))))

    
       # Show slice
       ax = axs[-1]  # next axis after params
       im = ax.imshow(img, cmap='gray')
       ax.set_title('template')
       ax.axis('off')
       
       # overlay countour
       if countour_path is not None:
           ax.contour(contour, levels=[0.5], colors='white', linewidths=1)
    
       # Add colorbar
       cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.02)
       cbar.set_ticks([0, maxint])
       cbar.ax.set_xticklabels([], rotation=-45)
       cbar.ax.tick_params(labelsize=10)
       
   # Hide unused axes
   n_used = n_params + (1 if template_path is not None else 0)
   if len(axs) > n_used:
        for ax in axs[n_used:]:
            ax.set_visible(False)

   
   # Save and close figure
  # plt.tight_layout(rect=[0, 0, 1, 1])
   plt.savefig(os.path.join(output_path, f'{model}_summary.png'))
   plt.close(fig)
    
##### BVALS #####


def modify_units_bvals(path_bvals, new_path_bvals):

    with open(new_path_bvals, 'w') as f:
        f.write(' '.join(map(str, np.loadtxt(path_bvals)*1e-3)))
    f.close()


def param_extract(param_path, bvals_path, desiredbvals, output_path):

    param = np.loadtxt(param_path)
    bvals = np.loadtxt(bvals_path)
    new_param = param[np.in1d(bvals, desiredbvals)]

    with open(output_path, 'w') as f:
        f.write(' '.join(map(str, new_param)))
    f.close()


def prep_bval_file(input_path, output_path):

    #  Step 1: Load the file with numbers (assuming each number is space-separated)
    with open(input_path, 'r') as file:
        # Read the entire line and split it into a list of strings
        numbers = file.readline().split()

    # Step 2: Process the numbers (divide by 1000 and round to 1 decimal)
    processed_numbers = [round(float(num) / 1000, 1) for num in numbers]

    # Step 3: Save the processed numbers back into a new text file
    with open(output_path, 'w') as file:
        # Join the numbers back into a string with space separation and write to the file
        file.write(' '.join(map(str, processed_numbers)))


def calc_avg_bval(bvals_nom, bvals_eff):

    # Calculate average effective values for each nominal shell value
    unique_numbers = np.unique(bvals_nom)
    means = []
    for unique_number in unique_numbers:
        corresponding_indices = np.where(bvals_nom == unique_number)
        corresponding_values = bvals_eff[corresponding_indices]
        mean_value = np.mean(corresponding_values)
        means.append(mean_value)

    return means

##### STUDY OPERATIONS #####


def concat_param(excel_param, path_bvals, output_path):

    array_param = []
    for k in range(0, len(excel_param)):
        array_param.extend(
            np.repeat(excel_param[k], len(np.loadtxt(path_bvals[k]))))

        with open(output_path, 'w') as f:
            f.write(' '.join(map(str, array_param)))
        f.close()


def extract_scan_no(scan_list, scan_idx, study_scanMode, paths_fwd, paths_rev):

    if scan_list['scanQA'][scan_idx] == 'ok' and scan_list['acqType'][scan_idx] == study_scanMode:
        # Create list with the scan numbers
        if scan_list['phaseDir'][scan_idx] == 'fwd':
            if scan_list['preprocOrder'][scan_idx] == 'first':
                paths_fwd.insert(0, scan_list['scanNo'][scan_idx])
            else:
                paths_fwd.append(scan_list['scanNo'][scan_idx])
        if scan_list['phaseDir'][scan_idx] == 'rev':
            paths_rev.append(scan_list['scanNo'][scan_idx])
    return paths_fwd, paths_rev


def extract_methods(methods_in, bids_strc, acqp, cfg=None):

    if acqp == 'PGSE':

        bmatrix = []
        bvals_eff = []
        bvals_nom = []
        dirs = []
        no_dirs_pershell = []
        with open(methods_in, 'r') as f:
            for line in f:
                if '##TITLE=' in line:
                    PV_version = (line.split(' ')).pop().strip()
                if '##$PVM_NRepetitions=' in line:
                    no_rep = int(line.split('=')[1])
                if '##$PVM_DwNDiffExpEach=' in line:
                    noShells = int(line.split('=')[1])
                if 'PVM_DwNDiffExp=' in line:
                    no_dwi = int(line.split('=')[1])
                if '##$PVM_DwDir=' in line:
                    values = re.findall(r'\d+', line.split('=')[1])
                    dims_dirs = [int(x) for x in values]
                    line = next(f)
                    while not '##$PVM_DwDgSwitch=' in line:
                        dirs.extend([float(x) for x in line.split()])
                        line = next(f)
                if '##$PVM_DwBMat=' in line:
                    values = re.findall(r'\d+', line.split('=')[1])
                    dims_bmat = [int(x) for x in values]
                    line = next(f)
                    while not '##$PVM_DwBMatImag' in line:
                        bmatrix.extend([float(x) for x in line.split()])
                        line = next(f)
                if '##$PVM_DwAoImages=' in line:
                    no_b0 =  int(line.split('=')[1])
                if '##$PVM_DwShellNDir=' in line: 
                    next_line = next(f).strip()
                    no_dirs_pershell = [int(x) for x in next_line.split()]
                if 'PVM_DwBvalEach=' in line:
                    line = next(f)
                    while not '##$PVM_' in line:
                        bvals_nom.extend([float(x) for x in line.split()])
                        line = next(f)
                if '##$PVM_DwEffBval=' in line:
                    line = next(f)
                    while not '##$PVM_' in line:
                        bvals_eff.extend([float(x) for x in line.split()])
                        line = next(f)

        if PV_version == 'V1.1':

            # Reshape and duplicate initial dirs matrix for all the shells
            dirs = np.array(dirs).reshape(
                dims_dirs[1], dims_dirs[0], order='F')
            dirs = np.repeat(dirs, noShells, axis=1)
            no_dirs = dims_dirs[0]
            no_b0 = dims_bmat[0] - dirs.shape[1]

            # Pad with zeros for b=0 images
            dirs = np.concatenate((np.zeros((3, no_b0)), dirs), axis=1)
            # Repeat directions for each Nrep
            dirs = np.tile(dirs, no_rep)

            # Create nominal bvals
            bvals_nom = [x for i in range(no_dirs) for x in bvals_nom]
            bvals_nom = np.concatenate(
                (np.zeros((no_b0, 1)), np.array(bvals_nom).reshape(-1, 1))).T
            bvals_nom = np.tile(bvals_nom, no_rep)

            # Format effective bvals
            bvals_eff = np.tile(np.array(bvals_eff)[:, np.newaxis].T, no_rep)

            # Format bmatrix
            bmatrix = np.array(bmatrix).reshape(
                dims_bmat[1]*dims_bmat[2], dims_bmat[0], order='F')
            bmatrix = np.tile(bmatrix, no_rep)

            # Save all files
            np.savetxt(bids_strc.get_path('bvalsEff.txt'), bvals_eff, delimiter=' ', fmt='%.1f')
            np.savetxt(bids_strc.get_path('bvalsNom.txt'), bvals_nom, delimiter=' ', fmt='%.1f')
            np.savetxt(bids_strc.get_path('bvecs.txt'), dirs, delimiter=' ', fmt='%.16f')
            np.savetxt(bids_strc.get_path('bmatrix.txt'), bmatrix, delimiter=' ', fmt='%.16f')

        elif PV_version == 'V3.5': 
            # bvals,and direction are saved differently from previous version 
            # direction number is already the total number of directions in a multi-shell protocol, so no need to replicate for all shells

            # Reshape for 3 x number of dirs
            dirs                = np.array(dirs).reshape(dims_dirs[1], dims_dirs[0], order='F')
            no_dirs             = dims_dirs[0]
            
            # Pad with zeros for b=0 images
            dirs                = np.concatenate((np.zeros((3, no_b0)), dirs), axis=1)
            
            # Create nominal bvals
            bvals_nom_t = []
            bvals_nom_t = [value for value, count in zip(bvals_nom, no_dirs_pershell) for _ in range(count)]
            bvals_nom_t = np.concatenate(
                (np.zeros((no_b0, 1)), np.array(bvals_nom_t).reshape(-1, 1))).T
          
            # Format bmatrix
            bmatrix             = np.array(bmatrix).reshape(dims_bmat[1]*dims_bmat[2], dims_bmat[0], order='F')
            #bmatrix             = np.tile(bmatrix, no_rep)
            
            # Format effective bvals
            bvals_eff           = np.tile(np.array(bvals_eff)[:,np.newaxis].T, no_rep)
            
            # Save all files
            np.savetxt(bids_strc.get_path('bvalsEff.txt'), bvals_eff, delimiter=' ', fmt='%.1f')
            np.savetxt(bids_strc.get_path('bvalsNom.txt'), bvals_nom_t, delimiter=' ', fmt='%.1f')
            np.savetxt(bids_strc.get_path('bvecs.txt'), dirs, delimiter=' ', fmt='%.16f')
            np.savetxt(bids_strc.get_path('bmatrix.txt'), bmatrix, delimiter=' ', fmt='%.16f')
            
        else:
            raise ValueError('ParaVision version not recognized!')
    
    elif acqp == 'STE':

        bvals_nom  =  np.loadtxt(os.path.join(cfg['common_folder'],'STE_bvalsNom.txt'))
        bvals_eff  =  np.loadtxt(os.path.join(cfg['common_folder'],'STE_bvalsEff.txt'))
        bvecs_fake =  np.loadtxt(os.path.join(cfg['common_folder'],'STE_bvecs_fake.txt'))


        np.savetxt(bids_strc.get_path('bvalsNom.txt'), bvals_nom[np.newaxis, :], delimiter=' ', fmt='%.1f')
        np.savetxt(bids_strc.get_path('bvalsEff.txt'), bvals_eff[np.newaxis, :], delimiter=' ', fmt='%.1f')
        np.savetxt(bids_strc.get_path('bvecs_fake.txt'), bvecs_fake, delimiter=' ', fmt='%.1f')

def order_bvals(bvals):
      
    # 1. Create groups
    blocks = []
    start_idx = 0
    for value, group in groupby(enumerate(bvals), key=lambda x: x[1]):
        group = list(group)
        indices = [i for i, _ in group]
        blocks.append({
            'bval': value,
            'indices': indices,
            'values': [bvals[i] for i in indices]
        })
    
    # 2: Sort blocks by b-value
    sorted_blocks = sorted(blocks, key=lambda b: b['bval'])
    
    # 3: Reassemble
    sorted_bvals = []
    sorted_indices = []
    
    for block in sorted_blocks:
        sorted_bvals.extend(block['values'])
        sorted_indices.extend(block['indices'])

    return sorted_bvals, sorted_indices

##### IMAGE OPERATIONS - HANDLING #####

def binary_op(input_path1, input_path2, operation, output_path):

    call = [f'fslmaths',
            f'{input_path1}',
            f'{operation}',
            f'{input_path2}'
            f' {output_path}']

    print(' '.join(call))
    os.system(' '.join(call))


def image_operations(input_image1, input_image2, operation, output_image):

    # to choose between 'm' - multiply
    #                   '+' - add
    #                   '-' - substract
    #                   '/' - divide
    #                   '^2' - power

    nii_1_shape = nib.load(input_image1).shape
    nii_2_shape = nib.load(input_image2).shape

    if len(nii_1_shape) != len(nii_2_shape):
        raise ValueError('Images need to have the same shape!')
    else:
        nii_ndims = len(nii_1_shape)

    call = [f'ImageMath',
            f'{nii_ndims}',
            f'{output_image}',
            f'{operation}',
            f'{input_image1}',
            f'{input_image2}']

    os.system(' '.join(call))


def erode_im(input_path, output_path, sigma):

    nii_shape = nib.load(input_path).shape
    nii_ndims = len(nii_shape)

    call = [f'ImageMath ',
            f'{nii_ndims}',
            f'{output_path}',
            f'ME',
            f'{input_path}',
            f'{sigma}']

    os.system(' '.join(call))

def erode_im_fsl(input_path, output_path):

    call = [f'fslmaths ',
            f'{input_path}',
            f'-ero -ero -ero -ero',
            f'{output_path}']

    os.system(' '.join(call))
  
    
def fsl_mult(input_path1, input_path2, output_path):
    
    call = [f'fslmaths',
            f'{input_path1}',
            f'-mul {input_path2}',
            f'{output_path}']
    os.system(' '.join(call))
    
def dilate_im(input_path, output_path, sigma):

    nii_shape = nib.load(input_path).shape
    nii_ndims = len(nii_shape)

    call = [f'ImageMath ',
            f'{nii_ndims}',
            f'{output_path}',
            f'MD',
            f'{input_path}',
            f'{sigma}']

    os.system(' '.join(call))


def disassemble_4D(input_path, output_prefix):

    call = [f'ImageMath 4',
            f'{output_prefix}',
            f'TimeSeriesDisassemble',
            f'{input_path}']

    os.system(' '.join(call))


def filter_clusters_by_size(input_path, output_path, min_size):
    # Load the binary image
    img = nib.load(input_path)
    data = img.get_fdata()

    # Label connected clusters of 1s
    labeled_data, num_features = label(data)

    # Find objects in the labeled array
    objects = find_objects(labeled_data)

    # Initialize an empty array to store the filtered data
    filtered_data = np.zeros_like(data)

    # Loop through each object and check its size
    for i, obj in enumerate(objects):
        cluster = (labeled_data[obj] == (i + 1))
        cluster_size = cluster.sum()

        # If the cluster size is greater than or equal to the min_size, keep it
        if cluster_size >= min_size:
            filtered_data[obj][cluster] = 1

    # Create a new Nifti1Image
    filtered_img = nib.Nifti1Image(filtered_data, img.affine, img.header)

    # Save the new image
    nib.save(filtered_img, output_path)


def threshold_image(input_image, output_image, thr_low, thr_high):  # rita

    nii_shape = nib.load(input_image).shape
    nii_ndims = len(nii_shape)

    call = [f'ThresholdImage',
            f'{nii_ndims}',
            f'{input_image}',
            f'{output_image}',
            f'{thr_low}',
            f'{thr_high}']

    os.system(' '.join(call))


def dilate_image(input_image, output_image):

    call = [f'fslmaths',
            f'{input_image}',
            f'-dilD',
            # add -kernel flag in order to modify the box kernel / by default is 3x3x3
            f'{output_image}']

    os.system(' '.join(call))


def make_mask(input_path, output_path, val):
    call = [f'fslmaths',
            f'{input_path}',
            f'-thr {val} -bin',
            f'{output_path}']

    os.system(' '.join(call))
    
    # Fill in holes by creating mask
    call = [f'fslmaths',
            f'{output_path} -fillh'  ,
            f'{output_path}']
    
    print(' '.join(call))
    os.system(' '.join(call))

def fsl_reorient(input_path):
    call = [f'fslreorient2std',
            f'{input_path}',
            f'{input_path}']

    os.system(' '.join(call))

def multiply_by_mask(input_path, output_path, mask_path):
    
    create_directory(output_path)
    if input_path.endswith('.nii.gz'):
        output_path = os.path.join(output_path,os.path.basename(input_path).replace('.nii.gz', '_masked.nii.gz'))
    elif input_path.endswith('.nii'):
        output_path = os.path.join(output_path,os.path.basename(input_path).replace('.nii', '_masked.nii.gz'))
            
    call = [f'fslmaths',
            f'{input_path}',
            f'-mul {mask_path}',
            f'{output_path}']

    os.system(' '.join(call))

def create_countour(mask):
    
    dilate_im(mask, mask.replace('.nii.gz','_dil.nii.gz'), '1')

    # make countour mask
    countour_path = mask.replace('.nii.gz', '_bo_contour.nii.gz')
    call = [f'fslmaths',
            f"{mask.replace('.nii.gz','_dil.nii.gz')}",
            f'-add',
            f'{mask}',
            f'-uthr 1',
            f"{mask.replace('.nii.gz', '_contour.nii.gz')}"]
    print(' '.join(call))
    os.system(' '.join(call))
               
    
def create_inverse_mask(mask_path, brain_mask, output_path):
    
    create_directory(output_path)
    if mask_path.endswith('.nii.gz'):
        output_path = os.path.join(output_path,os.path.basename(mask_path).replace('_mask.nii.gz', '_inv_mask.nii.gz'))
    elif mask_path.endswith('.nii'):
        output_path = os.path.join(output_path,os.path.basename(mask_path).replace('_mask.nii', '_inv_mask.nii.gz'))
            
    call = [f'fslmaths',
            f'{brain_mask}',
            f'-sub {mask_path}',
            f'{output_path}']

    os.system(' '.join(call))
    
##### IMAGE OPERATIONS - PROCESSING #####

def extract_b0(dwi_path, bvec_path, bval_path, output_path):

    for ii in range(len(dwi_path)):

        dwi = dwi_path[ii]
        bvec = bvec_path[ii]
        bval = bval_path[ii]
        output = output_path[ii]

        call = [f'dwiextract',
                f'-bzero',
                f'-fslgrad {bvec} {bval}',
                f'{dwi}',
                f'{output}',
                f'-force']

        os.system(' '.join(call))

def extract_vols(dwi_path, outputpath, volstart, volend):
    
    call = [f'fslroi {dwi_path}', 
            f'{outputpath}',
            f'{volstart} {volend}']
    
    os.system(' '.join(call))



def dwi_extract(old_dataset, new_dataset, bvals_list,):

    old_dwi = old_dataset["dwi"]
    old_bvals = old_dataset["bvals"]
    old_bvecs = old_dataset["bvecs"]
    new_dwi = new_dataset["dwi"]
    new_bvals = new_dataset["bvals"]
    new_bvecs = new_dataset["bvecs"]

    call = [f'dwiextract',
            f'{old_dwi}',
            f'{new_dwi}',
            f'-shells {bvals_list}',
            f'-fslgrad {old_bvecs} {old_bvals}',
            f'-export_grad_fsl {new_bvecs} {new_bvals} -force']

    print(' '.join(call))
    os.system(' '.join(call))


def denoise_vols_default_kernel(input_path, output_path, noise_path):
    """
    Function that denoises data with mrtrix
   
    Args:
        input_path (str)       : path where the original data is
        output_path (str)      : path where to save the denoised data
        noise_path  (str)      : path where to save the sigma map generated
      
    Returns:
        none
    """
    
    # denoise data
    call = [f'dwidenoise',
            f'{input_path}',
            f'{output_path}',
            f'-noise {noise_path}',
            f'-debug',
            f'-force -info']
    os.system(' '.join(call))

    # calculate residuals
    res_path = output_path.replace('.nii.gz', '_res.nii.gz')
    call = [f'mrcalc',
            f'{input_path}',
            f'{output_path}',
            f'-subtract',
            f'{res_path} -force']
    os.system(' '.join(call))
        
    # calculate sigma map by hand, different definition as implemented in mrtrix
    sigma_path  = output_path.replace('.nii.gz','_sigma2.nii.gz')
    res         = nib.load(res_path).get_fdata()
    template    = nib.load(res_path)
    sigma       = np.std(res,3)
    sigma_img   = nib.Nifti1Image(sigma, affine=template.affine, header=template.header)
    nib.save(sigma_img,  sigma_path) 
    

def denoise_designer(input_path, bvecs, bvals, output_path, data_path, algorithm_name):
    """
    Function that denoises data in Designer Toolbox implemented in Docker
   
    Args:
        input_path (str)       : path where the original data is
        bvecs (str)            : path to the bvecs file
        bvals (str)            : path to the bvals file
        output_path (str)      : path where to save the denoised data
        data_path  (str)       : path to the main folder of analysis to mount the data folder in Docker
        algorithm_name (str)   : name of the algortihm to use for the denoising. 
                                    Options are: veraart, jespersen
            
    Returns:
        none
    """
     
    # calculate kernel size
    num_vols  = len(read_numeric_txt(bvals).T)
    N = math.ceil(num_vols ** (1/3))  
    if N % 2 == 0:
        N += 1  
        
    # convert to mif
    nifti_to_mif(input_path, bvecs, bvals, input_path.replace('.nii.gz','.mif'))

    # run denoising
    docker_path  = '/data'
    input_path   = input_path.replace(data_path,docker_path)
    input_path   = input_path.replace('.nii.gz','.mif')
    output_path  = output_path.replace(data_path,docker_path)
    output_path  = output_path.replace('.nii.gz','.mif')

    call = [f'docker run -v {data_path}:/data nyudiffusionmri/designer2:v2.0.10 designer -denoise',
            f'{input_path}',
            f'{output_path} -pf 0.75 -pe_dir i -algorithm {algorithm_name} -extent {N},{N},{N} -debug']
    
    os.system(' '.join(call))
    print(' '.join(call))

    # convert back to nii.gz and put with the same header as original image
    output_path  = output_path.replace(docker_path,data_path)
    nifti_to_mif(output_path, output_path.replace('.mif','.bvec'), output_path.replace('.mif','.bval'), output_path.replace('.mif','.nii.gz'))
    input_path   = input_path.replace(docker_path,data_path)
    input_path   = input_path.replace('.mif','.nii.gz')
    output_path  = output_path.replace('.mif','.nii.gz')
    call = [f'flirt',
         f'-in  {output_path}',
         f'-ref {input_path}',
         f'-out {output_path}',
         f'-applyxfm -usesqform']
    os.system(' '.join(call))

    # calculate residuals
    res_path = output_path.replace('.nii.gz','_res.nii.gz')
    call     = [f'mrcalc',
            f'{input_path}',
            f'{output_path}',
            f'-subtract',
            f'{res_path} -force']
    os.system(' '.join(call))
    
    # calculate sigma map
    sigma_path  = output_path.replace('.nii.gz','_sigma.nii.gz')
    res         = nib.load(res_path).get_fdata()
    template    = nib.load(res_path)
    sigma       = np.std(res,3)
    sigma_img   = nib.Nifti1Image(sigma, affine=template.affine, header=template.header)
    nib.save(sigma_img,  sigma_path) 
    

def denoise_matlab(input_path, output_path, delta_path, code_path, toolbox_path, dn_type):
    """
    Function that denoises data in matlab

    Args:
        input_path (str)       : path where the original data is 
        out_path (str)         : path where to save the denoised data. 
        delta_path (str)       : path to where the diffusion times associated with each volume in input_path is
        code_path  (str)       : path to where the matlab code that does the denoising is.
        toolbox_path  (str)    : path to where the toolboxes for matlab are. 
                    They are added to the path in Matlab directly
        dn_type(str)           : string referring to the type of denoising. Options are: MPPCA, tMPPCA-4D, tMPPCA-5D

    Returns:
        none
    """

    # unzip file because spm doesn't like .gz
    input_path = gunzip_file(input_path)

    # replace .nii.gz by just .nii
    output_path = output_path.replace('.nii.gz','.nii')
    
    # read delta times and calculate how many volumes per diffusion time exists
    Deltas = read_numeric_txt(delta_path)
    unique_deltas, counts = np.unique(Deltas, return_counts=True)

    # calculate kernel size
    num_vols  = nib.load(input_path).shape[-1]
    #num_vols = counts[0]
    N = math.ceil(num_vols ** (1/3))  
    if N % 2 == 0:
        N += 1  
     
    # matlab command
    matlab_cmd = (
        f"try, "
        f"addpath('{code_path}'); "
        f"denoise_in_matlab('{input_path}', '{output_path}' ,'{counts}','{N}', '{toolbox_path}', '{dn_type}'); "
        f"catch, exit(1), end, exit(0)"
    )
    cmd = [
        "matlab", "-nodisplay", "-nosplash", "-nodesktop",
        "-r", matlab_cmd
    ]

    subprocess.run(cmd)

    # gzip file
    output_path = gzip_file(output_path)
    
    # calculate residuals
    input_path = input_path.replace('.nii','.nii.gz')
    res_path = output_path.replace('.nii.gz','_res.nii.gz')
    call     = [f'mrcalc',
            f'{input_path}',
            f'{output_path}',
            f'-subtract',
            f'{res_path} -force']
    os.system(' '.join(call))
    
    gzip_file(output_path.replace('.nii.gz','_sigma.nii'))
    
    # calculate sigma map
    # sigma_path  = output_path.replace('.nii.gz','_sigma.nii.gz')
    # res         = nib.load(res_path).get_fdata()
    # template    = nib.load(res_path)
    # sigma       = np.std(res,3)
    # sigma_img   = nib.Nifti1Image(sigma, affine=template.affine, header=template.header)
    # nib.save(sigma_img,  sigma_path) 

# def denoise_vols(input_path, kernel_size, ouput_path, noise_path):

#     call = [f'dwidenoise',
#             f'{input_path}',
#             f'{ouput_path}',
#             f'-extent {kernel_size}',
#             f'-noise {noise_path}',
#             f'-debug',
#             f'-force']
#     os.system(' '.join(call))


# def denoise_vols_mask(input_path, kernel_size, mask_path, ouput_path, noise_path):

#     call = [f'dwidenoise',
#             f'{input_path}',
#             f'{ouput_path}',
#             f'-extent {kernel_size}',
#             f'-noise {noise_path}',
#             f'-mask {mask_path}',
#             f'-debug',
#             f'-force']

#     os.system(' '.join(call))


def denoise_img(input_path, dim, ouput_path):
    """
    Function that denoises data with Ants tools
   
    Args:
        input_path (str)       : path where the original data is
        dim (int)              : dimension of the dataset
        output_path (str)      : path where to save the denoised data
      
    Returns:
        none
    """
    
    call = [f'DenoiseImage',
            f'-d {dim}',
            f'-i {input_path}',
            f'-o {ouput_path}',
            f'-v 1']

    os.system(' '.join(call))


def extract_bvals(dwi_path, bval_nom_path, bval_eff_path, bvec_path, bvals_list, output_pref):

    bval_nom_name = get_filename(bval_nom_path)
    bval_eff_name = get_filename(bval_eff_path)
    bvec_name = get_filename(bvec_path)

    # Extract dwi's corresponding to bvals_list
    dwi_extract(dwi_path, bval_nom_path, bvec_path, bvals_list, output_pref)

    bvals_nom = read_numeric_txt(bval_nom_path)
    bvals_eff = read_numeric_txt(bval_eff_path)

    # Extract b-vals nom and eff corresponding to bvals_list
    bvals_thr = list(map(int, bvals_list.split(',')))
    idx_bvals_thr = np.where(np.isin(bvals_nom.flatten(), bvals_thr))

    bvals_nom_thr = [bvals_nom.flatten()[index]
                     for index in list(idx_bvals_thr)]
    bvals_eff_thr = [bvals_eff.flatten()[index]
                     for index in list(idx_bvals_thr)]

    flat_bvals_nom_thr = np.concatenate(bvals_nom_thr).tolist() if isinstance(
        bvals_nom_thr[0], np.ndarray) else bvals_nom_thr
    write_txt(list(flat_bvals_nom_thr), bval_nom_path.replace(
        bval_nom_name, bval_nom_name + output_pref), 'w')

    flat_bvals_eff_thr = np.concatenate(bvals_eff_thr).tolist() if isinstance(
        bvals_eff_thr[0], np.ndarray) else bvals_eff_thr
    write_txt(list(flat_bvals_eff_thr), bval_eff_path.replace(
        bval_eff_name, bval_eff_name + output_pref), 'w')

    with open(bvec_path, 'r') as file:
        lines = file.readlines()

    # Extract columns based on indices
    extracted_lines = []
    for line in lines:
        # Split the line into columns (assuming space-separated values)
        columns = line.strip().split()
        # Extract the desired columns
        extracted_columns = [columns[i] for i in list(idx_bvals_thr[0])]
        # Join the extracted columns back into a string
        extracted_line = ' '.join(extracted_columns)
        extracted_lines.append(extracted_line)

    # Write the extracted lines to the output file
    with open(bvec_path.replace(bvec_name, bvec_name + output_pref), 'w') as file:
        for line in extracted_lines:
            file.write(line + '\n')
            
# def estim_DTI_DKI_designer(input_mif, 
#                            mask_path, output_path, data_path):

#     call = [f'docker run -v {data_path}:/data nyudiffusionmri/designer2:v2.0.10 tmi -SMI',
#             f'{input_mif}',
#             f'{output_path}',
#             f'-mask {mask_path}']

#     print(' '.join(call))
#     os.system(' '.join(call))

def mdm_matlab(bids_LTE, bids_STE, bids_STE_reg, header, output_path, code_path, toolbox_path, low_b=False):
    """
    Function that denoises data in matlab

    Args:
        input_path (str)       : path where the original data is 
        out_path (str)         : path where to save the denoised data. 
        delta_path (str)       : path to where the diffusion times associated with each volume in input_path is
        code_path  (str)       : path to where the matlab code that does the denoising is.
        toolbox_path  (str)    : path to where the toolboxes for matlab are. 
                    They are added to the path in Matlab directly
        low_b  (Bollean)       : True if the user wnats to use only the low bvalues to compute the microFa.
                    Useful for regions like the CSF where at high b-values the noise floor is reached and then microFA doesn't make sense

    Returns:
        none
    """
    
    # unzip file because spm doesn't like .gz
    STE_path  = bids_STE_reg.get_path('STE_in_LTE_dn_gc_topup.nii.gz')
    LTE_path  = get_file_in_folder(bids_LTE,'*dwi_dn_gc_ec.nii.gz')
    LTE_path = gunzip_file(LTE_path)
    STE_path = gunzip_file(STE_path)
    header   = gunzip_file(header)

    # replace .nii.gz by just .nii
    output_path = output_path.replace('.nii.gz','.nii')
    create_directory(output_path)
    
    # prepare files
    copy_files([LTE_path, STE_path, header], [os.path.join(output_path,'LTE.nii'),os.path.join(output_path,'STE.nii'),os.path.join(output_path,'mask.nii')])
    copy_files([get_file_in_folder(bids_LTE,'*bvalsNom.txt'), get_file_in_folder(bids_LTE,'*bvecsRotated.txt')], [os.path.join(output_path,'LTE.bval'),os.path.join(output_path,'LTE.bvec')])
    copy_files([bids_STE.get_path('bvalsNom.txt'), bids_STE.get_path('bvecs_fake.txt')], [os.path.join(output_path,'STE.bval'),os.path.join(output_path,'STE.bvec')])

    if low_b==True:
        # extract low b vals
        bvals_LTE = read_numeric_txt(os.path.join(output_path,'LTE.bval'))
        desiredbvals = np.unique(bvals_LTE[bvals_LTE<=2000])
        old_dataset = {"dwi":  os.path.join(output_path,'LTE.nii'), 
                        "bvals": os.path.join(output_path,'LTE.bval'),
                        "bvecs": os.path.join(output_path,'LTE.bvec')}
        new_dataset = {"dwi":   os.path.join(output_path,'LTE.nii'), 
                        "bvals": os.path.join(output_path,'LTE.bval'),
                        "bvecs": os.path.join(output_path,'LTE.bvec')}              
        dwi_extract(old_dataset, new_dataset,','.join(map(str, desiredbvals)))
    
    # new paths
    LTE_path = os.path.join(output_path,'LTE.nii')
    STE_path = os.path.join(output_path,'STE.nii')

    # matlab command
    matlab_cmd = (
        f"try, "
        f"addpath('{code_path}'); "
        f"calculate_microFA('{LTE_path}', '{STE_path}' ,'{header}','{output_path}', '{toolbox_path}'); "
        f"catch, exit(1), end, exit(0)"
    )
    cmd = [
        "matlab", "-nodisplay", "-nosplash", "-nodesktop",
        "-r", matlab_cmd
    ]
    
    subprocess.run(cmd)

# def estimate_SMI_designer(input_mif, mask_path, sigma_path, output_path, data_path, others):

#     call = [f'docker run -v {data_path}:/data nyudiffusionmri/designer2:v2.0.12 tmi -SMI',
#             f'{others}',
#             f'-sigma {sigma_path}',
#             f'-mask {mask_path}',
#             f'{input_mif}',
#             f'{output_path}']

#     print(' '.join(call))
#     os.system(' '.join(call))
    
#     call = [f'docker run -v {data_path}:/data nyudiffusionmri/designer2:v2.0.10 chmod -R 777 {output_path}']
#     print(' '.join(call))
#     os.system(' '.join(call))
    
def calculate_pwd_avg(dwi_path, bval_nom_path, bval_eff_path, output_path, diff_time):

    bvals_nom = read_numeric_txt(bval_nom_path)
    bvals_nom_unique, idx = np.unique(bvals_nom, return_index=True)
    sorted_bvals_nom_unique = bvals_nom_unique[np.argsort(idx)]

    bvals_eff = read_numeric_txt(bval_eff_path)
    bvals_eff = bvals_eff.flatten()
    bval_eff_avg = []

    dwi_img = nib.load(dwi_path)
    dwi_array = dwi_img.get_fdata()

    dwi_name = get_filename(dwi_path)
    bval_nom_name = get_filename(bval_nom_path)
    bval_eff_name = get_filename(bval_eff_path)

    # Calculate average b0
    b0_idx = np.where(bvals_nom == 0.0)[1]
    b0_vols = dwi_array[:, :, :, b0_idx]
    b0_avg = np.mean(b0_vols, axis=3)

    dwi_avg = np.empty(tuple(
        dwi_img.shape[0:3]) + (len(sorted_bvals_nom_unique[sorted_bvals_nom_unique != 0]),))
    dwi_avg_norm = np.empty(tuple(
        dwi_img.shape[0:3]) + (len(sorted_bvals_nom_unique[sorted_bvals_nom_unique != 0]),))

    for ii in range(len(sorted_bvals_nom_unique)):

        bval_idx = np.where(bvals_nom == sorted_bvals_nom_unique[ii])[1]
        dwi_vols = dwi_array[:, :, :, bval_idx]
        dwi_temp = np.mean(dwi_vols, axis=3)
        dwi_temp_norm = dwi_temp / b0_avg

        array_to_nii(dwi_img, dwi_temp, os.path.join(output_path, dwi_name +
                     '_' + str(sorted_bvals_nom_unique[ii].astype(int)) + '.nii.gz'))
        array_to_nii(dwi_img, dwi_temp_norm, os.path.join(output_path, dwi_name +
                     '_norm_' + str(sorted_bvals_nom_unique[ii].astype(int)) + '.nii.gz'))

        write_txt(sorted_bvals_nom_unique[sorted_bvals_nom_unique != 0] /
                  1e3, os.path.join(output_path, bval_nom_name + '_avg.txt'), 'w')
        write_txt([str(diff_time)] * len(sorted_bvals_nom_unique[sorted_bvals_nom_unique != 0]),
                  os.path.join(output_path, 'td.txt'), 'w')

        if sorted_bvals_nom_unique[ii] != 0:

            dwi_avg[:, :, :, ii-1] = dwi_temp
            dwi_avg_norm[:, :, :, ii-1] = dwi_temp_norm

            bval_eff_avg.append(
                round(np.mean(bvals_eff[np.array(bval_idx)])/1e3, 2))

            array_to_nii(dwi_img, dwi_avg, os.path.join(
                output_path, dwi_name + '_pwd_avg.nii.gz'))
            array_to_nii(dwi_img, dwi_avg_norm, os.path.join(
                output_path, dwi_name + '_pwd_avg_norm.nii.gz'))

        write_txt(bval_eff_avg, os.path.join(
            output_path, bval_eff_name + '_avg.txt'), 'w')


def N4_unbias(input_path, output_path):

    size = nib.load(input_path).shape
    dim = len(size)

    if dim == 2 or dim == 3 or dim == 4:
        call = [f'N4BiasFieldCorrection',
                f'-d {dim} ',
                f'-i {input_path}',
                f'-o {output_path}',
                f'-c [50x50x50x50,0.00001]',
                f'-v 1']
        os.system(' '.join(call))

    else:
        raise ValueError('Image dimension not supported!')


def gibbs_corr(input_path, output_path):

    call = [f'mrdegibbs ',
            f'{input_path}',
            f'{output_path}',
            f'-force']

    os.system(' '.join(call))


def make_avg(dim, input_path, output_path):  # rita

    for ii in range(len(input_path)):

        input = input_path[ii]
        output = output_path[ii]

        call = [f'antsMotionCorr',
                f'-d {dim}',
                f'-a {input}',
                f'-o {output}']

        os.system(' '.join(call))


def calc_snr(input_path1, input_path2, output_path):
    
    binary_op(input_path1, input_path2, '-div', output_path)

def brain_extract_BET(input_path):
    call = [f'bet',
            f'{input_path}',
            f'{input_path.replace(".nii.gz", "_brain.nii.gz")} -R -m']
    
    print(' '.join(call))
    os.system(' '.join(call))
    
def brain_extract_RATS(input_path, anat_thr):

    # Segment with ANTS to create brain mask
    RATS_path = 'RATS_MM'
    
    call = [f'{RATS_path}',
            f'{input_path}',
            f'{input_path.replace(".nii.gz", "_brain_mask.nii.gz")}',
            f'-t 1500']
    print(' '.join(call))
    os.system(' '.join(call))  
    
    # Use brain mask to get just the T2w brain image
    binary_op(input_path,input_path.replace(".nii.gz", "_brain_mask.nii.gz"), '-mul', input_path.replace(".nii.gz", "_brain.nii.gz"))
    
    # Apply extra threshold on intensity
    call = [f'fslmaths',
            f'{input_path.replace(".nii.gz", "_brain.nii.gz")}',
            f'-thr {anat_thr}', # 4000, 2100
            f'{input_path.replace(".nii.gz", "_brain.nii.gz")}']
    
    print(' '.join(call))
    os.system(' '.join(call))
    
    # Fill in holes by creating mask
    call = [f'fslmaths',
            f'{input_path.replace(".nii.gz", "_brain.nii.gz")} -fillh'  ,
            f'{input_path.replace(".nii.gz", "_brain_mask.nii.gz")}']
    
    print(' '.join(call))
    os.system(' '.join(call))
    
    # Extra clean of mask
    img     = nib.load(input_path.replace(".nii.gz", "_brain_mask.nii.gz"))
    data    = img.get_fdata()
    cleaned = binary_opening(data, structure=np.ones((3,3,3)))
    labeled, n = label(cleaned) # Keep only the largest connected component
    largest_component = (labeled == np.argmax(np.bincount(labeled.flat)[1:]) + 1)
    nib.save(nib.Nifti1Image(largest_component.astype(np.uint8), img.affine), input_path.replace(".nii.gz", "_brain_mask.nii.gz"))
    
    # Use new clean brain mask (there is never too many masks xD ) to get just the T2w brain image
    binary_op(input_path,input_path.replace(".nii.gz", "_brain_mask.nii.gz"), '-mul', input_path.replace(".nii.gz", "_brain.nii.gz"))

def brain_extract_organoids(input_path,val):
    
    #Make brain mask
    make_mask(input_path, input_path.replace(".nii.gz", "_brain_mask.nii.gz"), val)
    
    # Use brain mask to get just the T2w brain image
    binary_op(input_path,input_path.replace(".nii.gz", "_brain_mask.nii.gz"), '-mul', input_path.replace(".nii.gz", "_brain.nii.gz"))
    
    
def brain_extract_BREX(input_path,BREX_path):
    out_path = os.path.dirname(input_path)
    os.chdir(out_path)
    copy_files([os.path.join(BREX_path, 'atlasBREX.sh'),
                os.path.join(BREX_path, 'templates',
                             'mouse_T2w', 'b_mouse_T2w.nii.gz'),
                os.path.join(BREX_path, 'templates', 'mouse_T2w', 'nb_mouse_T2w.nii.gz')],
               [os.path.join(out_path, 'atlasBREX.sh'),
                os.path.join(out_path, 'b_mouse_T2w.nii.gz'),
                os.path.join(out_path, 'nb_mouse_T2w.nii.gz')])
   
    # antsreg(os.path.join(out_path, 'b_mouse_T2w.nii.gz'),input_path, input_path.replace('.nii.gz', 'template2anat'))
    # ants_apply_transforms([os.path.join(out_path, 'b_mouse_T2w.nii.gz'),os.path.join(out_path, 'nb_mouse_T2w.nii.gz')],
    #                               input_path,
    #                               [os.path.join(out_path, 'b_mouse_T2w_in_anat.nii.gz'),os.path.join(out_path, 'nb_mouse_T2w_in_anat.nii.gz')],
    #                               input_path.replace('.nii.gz', 'template2anat0GenericAffine.mat'),
    #                               input_path.replace('.nii.gz', 'template2anat1Warp.nii.gz'))                  
          

    call = [f'bash {os.path.join(out_path, "atlasBREX.sh")}',
            f'-b {os.path.join(out_path, "b_mouse_T2w.nii.gz")}',
            f'-nb {os.path.join(out_path, "nb_mouse_T2w.nii.gz")}',
            f'-h {input_path} -f 0.01 ']
    

    os.system(' '.join(call))
    
    remove_file(os.path.join(out_path, 'atlasBREX.sh'))
    remove_file(os.path.join(out_path, 'b_mouse_T2w.nii.gz'))
    remove_file(os.path.join(out_path, 'nb_mouse_T2w.nii.gz'))


def register_outputfits_to_anat(output_path, new_output_path,model,cfg, bids_strc_anat,bids_strc_prep):
     # auxiliar function to register output model fits to anatomical space
 
     anat_format = cfg['anat_format']
     create_directory(new_output_path)
     patterns, lims, maximums = get_param_names_model(model,cfg['is_alive'])
    
     for filename in os.listdir(output_path):
         if  any(fnmatch.fnmatch(filename, pattern) for pattern in patterns):
             
             if cfg['subject_type']=='organoid':
                  # pad image temporarily for registration
                  pad_image(os.path.join(output_path, filename),os.path.join(output_path, filename))
                  pad_image(bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'),bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'))
    
                  # Apply inverse transform to put template anat in dwi
                  ants_apply_transforms([os.path.join(output_path, filename)],  # input 
                                       bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'), # reference
                                       [os.path.join(new_output_path,filename.replace('.nii','_padded.nii'))], # output
                                       [bids_strc_prep.get_path(f'dwiafterpreproc2{anat_format}0GenericAffine.mat'), 0], # transform 1
                                       bids_strc_prep.get_path(f'dwiafterpreproc2{anat_format}1Warp.nii.gz')) # transform 2
    
                  # unpad the images previousy padded
                  unpad_image(os.path.join(output_path, filename),os.path.join(output_path, filename))
                  unpad_image(bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'),bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'))
                  unpad_image(os.path.join(new_output_path,filename.replace('.nii','_padded.nii')),os.path.join(new_output_path, filename))
                  remove_file(os.path.join(new_output_path,filename.replace('.nii','_padded.nii')))
                              
             else:
                 # Apply inverse transform to put template anat in dwi
                 ants_apply_transforms_simple([os.path.join(output_path, filename)],  # input 
                                      bids_strc_anat.get_path(f'{anat_format}_bc_brain.nii.gz'), # reference
                                      [os.path.join(new_output_path,filename)], # output
                                      [bids_strc_prep.get_path(f'dwiafterpreproc2{anat_format}0GenericAffine.mat'), 0]) # transform 1
        
        
##### ATLAS HANDLE #####


def prepare_atlas(atlas_name, atlas_folder, atlas_type):
    import glob

    if atlas_name== 'Atlas_WHS_v4' and atlas_type=='atlas':
        
        # Define atlas 
        atlas      = glob.glob(os.path.join(atlas_folder, atlas_name, '*atlas.nii.gz'))[0]
        template   = glob.glob(os.path.join(atlas_folder, atlas_name, '*template_brain.nii.gz'))[0]
       
        # remove extra regions to make it more like brain only
        img  = nib.load(atlas)
        data_a = img.get_fdata()
        data_a[data_a == 42] = 0
        data_a[data_a == 41] = 0
        data_a[data_a == 45] = 0
        data_a[data_a == 76] = 0
        img  = nib.load(template)
        data_t = img.get_fdata()           
        mask = (data_a != 0).astype(np.uint8)
        data_t = data_t*mask
        nib.save(nib.Nifti1Image(data_a, img.affine), atlas.replace('.nii.gz', '_crop.nii.gz'))
        nib.save(nib.Nifti1Image(data_t, img.affine), template.replace('.nii.gz', '_crop.nii.gz'))
        
        for image in (atlas,template):
         
         # Crop template/atlas - otherwise too much data to register
         img  = nib.load(image.replace('.nii.gz', '_crop.nii.gz'))
         data = img.get_fdata()
         masked_data = np.zeros_like(data)
         masked_data[:, 230:840, :] = data[:, 230:840, :]
         nib.save(nib.Nifti1Image(masked_data, img.affine), image.replace('.nii.gz', '_crop.nii.gz'))
         
         # Downsample template/atlas to avoid segmentation faults
         input_img = nib.load(image.replace('.nii.gz', '_crop.nii.gz'))
         if image==atlas:
             resampled_img = nib.processing.resample_to_output(input_img, [0.05, 0.5, 0.05],order=0)
         elif image==template:
             resampled_img = nib.processing.resample_to_output(input_img, [0.05, 0.5, 0.05])
         nib.save(resampled_img,  image.replace('.nii.gz', '_crop_lowres.nii.gz')) 
       
        # Define atlas 
        atlas      = atlas.replace('.nii.gz', '_crop_lowres.nii.gz')
        template   = template.replace('.nii.gz', '_crop_lowres.nii.gz')
        
    elif atlas_name == 'TPM_C57Bl6'  and atlas_type=='TPM':
    
        # Define TPM 
        atlas    = glob.glob(os.path.join(atlas_folder, atlas_name, '*TPM_C57Bl6_n30.nii'))[0]
        template = glob.glob(os.path.join(atlas_folder, atlas_name, '*C57Bl6_T2_n10_template_brain.nii'))[0]
        
        for image in (atlas,template):
            
            # Correct scale
            img = nib.load(image)
            data = img.get_fdata()
            affine = img.affine.copy()
            affine[:3, :3] /= 10  # Correct for the scale factor
            corrected_img = nib.Nifti1Image(data, affine, img.header)
            
            # Save the rescaled image
            nib.save(corrected_img, image.replace('.nii', '_rescaled.nii'))
            
            # Separate TPMs into different files
            if image==atlas:
                
                # get path of rescaled image
                input_path = image.replace('.nii', '_rescaled.nii')
                out_path = image.replace('.nii', '_rescaled_vol_')
                
                call = [f'fslsplit',
                        f'{input_path}',
                        f'{out_path} -t']
                print(' '.join(call))
    
                os.system(' '.join(call))
    
    
            #     # Resample each 3D volume of the 4D TPM atlas
            #     data = input_img.get_fdata()
            #     affine = input_img.affine
            #     header = input_img.header
        
            #     resampled_volumes = []
            #     for dim in range(data.shape[-1]):
            #        vol_3d = nib.Nifti1Image(data[..., dim], input_img.affine)
            #        resampled_vol = nib.processing.resample_to_output(vol_3d, voxel_sizes=[0.08, 0.5, 0.08], order=3)
            #        nib.save(resampled_img, image.replace('.nii', '_rescaled_lowres.nii'))
    
            #        resampled_volumes.append(resampled_vol.get_fdata())
        
            #     resampled_data = np.stack(resampled_volumes, axis=-1)
            #     resampled_img = nib.Nifti1Image(resampled_data, resampled_vol.affine, resampled_vol.header)
        
            # elif image==template:
            #     # Resample the 3D template directly
            #     resampled_img = nibabel.processing.resample_to_output(input_img, [0.08, 0.5, 0.08])
         
            # # Save final lowres version
            # nib.save(vol_3d, image.replace('.nii', '_rescaled_lowres.nii'))
           
        # Define TPM 
        atlas      = atlas.replace('.nii', '_rescaled.nii')
        template   = template.replace('.nii', '_rescaled.nii')
  
    # If no adjustments need to be made in atlas it just reads the files
    elif atlas_type=='atlas':
   
       # Define atlas 
       atlas      = glob.glob(os.path.join(atlas_folder, atlas_name, '*atlas*'))[0]
       template   = glob.glob(os.path.join(atlas_folder, atlas_name, '*template_brain*'))[0]
   
    # If no adjustments need to be made in atlas it just reads the files
    elif atlas_type=='TPM':
     
        # Define TPM 
        atlas_list    = glob.glob(os.path.join(atlas_folder, atlas_name, '*TPM*'))
        atlas = [a for a in atlas_list if 'vol' not in os.path.basename(a)][0]
        template = glob.glob(os.path.join(atlas_folder,atlas_name, '*template_brain*'))[0]
        
        # Separate TPMs into different files
        out_path = atlas.replace('.nii', '_vol_')     
        call = [f'fslsplit',
                f'{atlas}',
                f'{out_path} -t']
        print(' '.join(call))
        os.system(' '.join(call))
      
    return atlas, template

def prepare_atlas_labels(atlas_name, atlas_label_path):
    import glob
    import random
    
    #atlas_label_path = glob.glob(os.path.join(atlas_folder, atlas_name, '*label*'))[0]

    # Handle Atlas_WHS_v4 which have the labels in a specific formar 
    if atlas_name== 'Atlas_WHS_v4' or atlas_label_path.endswith('.label'):
        atlas_labels = pd.read_csv(atlas_label_path, sep=r'\s+', skiprows=14, header=None,
                                   names=['IDX', 'R', 'G', 'B', 'A', 'VIS', 'MSH', 'LABEL'], quotechar='"')
     
    # Handle DKT-style .txt file
    elif atlas_label_path.endswith('.txt') and atlas_name== 'Atlas_DKT':
        with open(atlas_label_path, 'r') as f:
            content = f.read()
    
        # Extract label entries like [1002, "left caudal anterior cingulate"]
        matches = re.findall(r'\[\s*(\d+)\s*,\s*"([^"]+)"\s*\]', content)
    
        labels = []
        for match in matches:
            index = int(match[0])
            name = match[1].strip()
    
            # Assign random RGB colors
            R, G, B = random.randint(0, 255), random.randint(0, 255), random.randint(0, 255)
            A, VIS, MSH = 1, 1, 0
    
            labels.append({
                'IDX': index,
                'R': R,
                'G': G,
                'B': B,
                'A': A,
                'VIS': VIS,
                'MSH': MSH,
                'LABEL': name
            })
    
        atlas_labels = pd.DataFrame(labels)
        atlas_labels.sort_values(by='IDX', inplace=True)
    
    
    # If none of the previous, assumes there is an xml file    
    else:
        import xml.etree.ElementTree as ET

        # Load the XML file
        tree = ET.parse(atlas_label_path)
        root = tree.getroot()
        
        # Extract label data
        labels = []
        for label in root.findall(".//label"):
            
            if 'index' in label.attrib:
                # Format 1: attribute-based
                index = int(label.attrib['index'])
                name = label.text.strip()
            else:
                # Format 2: element-based
                index_elem = label.find('index')
                name_elem = label.find('name')
                if index_elem is not None:
                    index = int(index_elem.text.strip())
                else:
                    continue  # Skip if no index
                name = name_elem.text.strip() if name_elem is not None else ''
        
            # Assign random RGB colors
            R = random.randint(0, 255)
            G = random.randint(0, 255)
            B = random.randint(0, 255)
        
            # Constants for alpha, visibility, mesh
            A = 1
            VIS = 1
            MSH = 0
        
            # for some reason this atlas have the labels shifted by 1 value, so we need to correct
            if 'Atlas_Juelich' in atlas_name: 
               index = index +1
               
            labels.append({
                'IDX': index,
                'R': R,
                'G': G,
                'B': B,
                'A': A,
                'VIS': VIS,
                'MSH': MSH,
                'LABEL': name
            })
        
        # Convert to DataFrame
        df = pd.DataFrame(labels)
        
        # Optional: Sort by index
        df.sort_values(by='IDX', inplace=True)

        atlas_labels = df


    return atlas_labels

def make_atlas_manual_organoid(organoid_mask,folder,label_path):
    
    atlas = organoid_mask.replace('_mask.nii.gz', '_atlas.nii.gz')
    
    candidate_masks = glob.glob(os.path.join(folder, '*organoid?_mask.nii.gz'))
    pattern = re.compile(r'organoid([A-Z])_mask\.nii\.gz$')
    organoid_masks = [f for f in candidate_masks if pattern.search(os.path.basename(f))]

    temp_files = []
    label_lines = []
    color_base = (255, 105, 180)  # Hot pink base

    for idx, mask_path in enumerate(organoid_masks, start=1):
        temp_file = mask_path.replace('.nii.gz', f'_temp.nii.gz')
        temp_files.append(temp_file)

        # Multiply mask by its unique label index
        call = ['fslmaths', mask_path, '-mul', str(idx), temp_file]
        os.system(' '.join(call))

        # Add a label line
        r, g, b = color_base  # You could vary color if needed
        match = re.search(r'(organoid[A-Za-z])_mask\.nii\.gz$', os.path.basename(mask_path))
        label = match.group(1) if match else 'unknown'
        label_lines.append(f'   {idx}   {r}  {g}  {b}        1  1  0    "{label}"\n')

    # Merge all the labeled masks
    merged_command = ['fslmaths', temp_files[0]]
    for temp_file in temp_files[1:]:
        merged_command += ['-add', temp_file]
    merged_command.append(atlas)
    os.system(' '.join(merged_command))

    # Remove temporary files
    for temp_file in temp_files:
        os.remove(temp_file)

    # Write the label file
    header = """\
    ################################################
    # ITK-SnAP Label Description File
    # File format: 
    # IDX   -R-  -G-  -B-  -A--  VIS MSH  LABEL
    # Fields: 
    #    IDX:   Zero-based index 
    #    -R-:   Red color component (0..255)
    #    -G-:   Green color component (0..255)
    #    -B-:   Blue color component (0..255)
    #    -A-:   Label transparency (0.00 .. 1.00)
    #    VIS:   Label visibility (0 or 1)
    #    MSH:   Label mesh visibility (0 or 1)
    #  LABEL:   Label description 
    ################################################
    """

    with open(label_path, 'w', encoding='utf-8') as f:
        f.write(header + '\n')
        f.writelines(label_lines)
        
   
        
#def make_atlas_label_organoid(organoid_mask, non_organoid_mask,label_path):
    
    # temp_non_organoid = non_organoid_mask.replace('.nii.gz', '_temp.nii.gz')
    # call = ['fslmaths',
    #          f'{non_organoid_mask}',
    #          '-mul', '2',
    #          f'{temp_non_organoid}'
    #          ]
  
    # os.system(' '.join(call))
    
    
    # atlas = organoid_mask.replace('_mask.nii.gz', '_atlas.nii.gz')
    # call = [f'fslmaths',
    #         f'{organoid_mask}',
    #         f'-add {temp_non_organoid}',
    #         f'{atlas}']

    # os.system(' '.join(call))
    # remove_file(temp_non_organoid)
    
    # header = """\
    #             ################################################
    #             # ITK-SnAP Label Description File
    #             # File format: 
    #             # IDX   -R-  -G-  -B-  -A--  VIS MSH  LABEL
    #             # Fields: 
    #             #    IDX:   Zero-based index 
    #             #    -R-:   Red color component (0..255)
    #             #    -G-:   Green color component (0..255)
    #             #    -B-:   Blue color component (0..255)
    #             #    -A-:   Label transparency (0.00 .. 1.00)
    #             #    VIS:   Label visibility (0 or 1)
    #             #    MSH:   Label mesh visibility (0 or 1)
    #             #  LABEL:   Label description 
    #             ################################################
    #             """

    # label_line = ['   1   255  105  180        1  1  0    "organoids"\n'
    #               '   2   255  105  180        1  1  0    "medium"\n']

    # with open(label_path, "w", encoding="utf-8") as f:
    #     f.write(header)
    #     f.write("\n")
    #     f.writelines(label_line)        
    
def create_ROI_mask(atlas, atlas_labels, TPMs, ROI, tpm_thr, bids_strc_reg):
 
     # Define tpms and threshold at 0.9
     atlas_shape = nib.load(atlas).shape
     
     # GM
     tpm_GM = [f for f in TPMs if 'GM' in f and f and os.path.exists(f)]
     tmp_GM = nib.load(tpm_GM[0]).get_fdata() > tpm_thr if tpm_GM else np.ones(atlas_shape, dtype=bool)
    
     # WM
     tpm_WM = [f for f in TPMs if 'WM' in f and f and os.path.exists(f)]
     tmp_WM = nib.load(tpm_WM[0]).get_fdata() > tpm_thr if tpm_WM else np.ones(atlas_shape, dtype=bool)
    
     # CSF
     tpm_CSF = [f for f in TPMs if 'CSF' in f and f and os.path.exists(f)]
     tmp_CSF = nib.load(tpm_CSF[0]).get_fdata() > tpm_thr if tpm_CSF else np.ones(atlas_shape, dtype=bool)
 

     if 'Atlas_WHS_v4' in atlas:
         # Define ROI labels for each ROI asked
         roi_definitions = {
            'hippocampus': ['hippocampus'],
            'M1': ['Primary motor area'],
            'M2': ['Secondary motor area'],
            'S1': ['Primary somatosensory'],
            'S2': ['Secondary somatosensory'],
            'V1': ['Primary visual'],
            'PL': ['Prelimbic'],
            'CG': ['Cingulate'],
            'CC': ['corpus callosum'],
            'Thal': ['thalamic nucleus'],
            'CSF': ['Ventricular'],
            'Cereb GM': ['erebellum'],
            'Cereb WM': ['erebellum'],
            'putamen': ['putamen'],
            'insula': ['insula'],
            'cornu': ['ornu'],
            'WB': ['whole brain']
        } 
     elif 'Atlas_Juelich' in atlas:
        
         roi_definitions = {
            'hippocampus': ['ippocampus'],
            'V1': ['V1'],
            'V2': ['V2'],
            'premotor': ['Premotor'],
            'parietal': ['parietal'],
            'S1': ['Primary somatosensory'],
            'S2': ['Secondary somatosensory'],
            'M1': ['Primary motor'],
            'Broca': ['Broca'],
            'CC': ['callosal'],
            'WB': ['whole brain']
        } 
     elif 'Atlas_Neuromorphometrics' in atlas:
        
         roi_definitions = {
            'frontal': ['frontal cortex', 'frontal gyrus'],
            'precentral': ['precentral gyrus'],
            'postcentral': ['postcentral gyrus'],
            'occipital': ['occipital gyrus'],
            'parietal': ['parietal lobule'],
            'temporal': ['temporal lobe','temporal gyrus'],
            'WB': ['whole brain']
        } 
         
     elif 'Atlas_DKT' in atlas:
        
         roi_definitions = {
            'frontal': ['frontal'],
            'precentral': ['precentral'],
            'postcentral': ['postcentral'],
            'cuneus': [' cuneus'],
            'occipital': ['occipital'],
            'parietal': ['parietal'],
            'temporal': ['temporal'],
            'WB': ['whole brain']
        } 
         
     elif 'anat_space_organoids' in atlas:
         
          roi_definitions = {
             'organoids': ['organoid'],
             'medium': ['medium'],

         } 
          
     else:
         print('>> Error: you did not defined which regions corresponds to each ROI. Please add a condition for that atlas in create_ROI_mask')
        
     # Load the atlas data
     template = nib.load(atlas)
     atlas_data = template.get_fdata()
      
     # Find matching indices for ROI
     if ROI=='WB':
         if  'Atlas_WHS_v4' in atlas:
            label_mask = ~atlas_labels["LABEL"].str.contains("olfactory", regex=True)
            match_idx = atlas_labels.loc[label_mask & (atlas_labels["IDX"] != 0), "IDX"].to_numpy()
         else:
            match_idx = atlas_labels["IDX"].to_numpy()[atlas_labels["IDX"].to_numpy() != 0]
     elif ROI == 'Thal':
         #  exclude "hypothalamic"
         ind_list = atlas_labels["LABEL"].str.contains('thalamic nucleus', regex=True) & ~atlas_labels["LABEL"].str.contains('Hypothalamic', regex=True) & ~atlas_labels["LABEL"].str.contains('unspecified', regex=True)
         match_idx = atlas_labels["IDX"][ind_list].to_numpy()
     else:
         ind_list = atlas_labels["LABEL"].str.contains("|".join(roi_definitions[ROI]), regex=True)
         match_idx = atlas_labels["IDX"][ind_list].to_numpy()
     
     # Create the mask for the ROI
     mask_indexes = np.isin(atlas_data, match_idx)
     masked_data = (mask_indexes > 0).astype(np.uint8)
    
     # Save the mask
     masked_img = nib.Nifti1Image(masked_data, affine=template.affine, header=template.header)
     nib.save(masked_img, bids_strc_reg.get_path(f'mask_{ROI}.nii.gz'))
    
     # Multiply by tissue probability map
     if ROI=='CC': # is WM
         masked_data = masked_data*tmp_WM
     elif ROI=='CSF':
         masked_data = masked_data*tmp_CSF
     elif ROI=='Cereb GM':
         masked_data = masked_data*tmp_GM 
     elif ROI=='Cereb WM':
         masked_data = masked_data*tmp_WM 
     elif ROI=='WB':
         masked_data = masked_data* (tmp_GM | tmp_WM | tmp_CSF)
     else:
         masked_data = masked_data*tmp_GM

     # Save the mask
     masked_img = nib.Nifti1Image(masked_data, affine=template.affine, header=template.header)
     nib.save(masked_img, bids_strc_reg.get_path(f'mask_{ROI}.nii.gz'))


     return masked_data
 
def create_ROI_mask_fromindx(atlas, atlas_labels, indx, bids_strc_reg):
    
    # Load the atlas data
    template = nib.load(atlas)
    atlas_data = template.get_fdata()
    
    mask_indexes = np.isin(atlas_data, indx)
    masked_data = (mask_indexes > 0).astype(np.uint8)
    

    return masked_data
    
##### NIFTI HANDLE #####


def nifti_to_mif(nifti_path, bvecs_path, bvals_path, mif_path):

    call = [f'mrconvert',
            f'-fslgrad {bvecs_path} {bvals_path} ',
            f'{nifti_path} {mif_path} -force']

    os.system(' '.join(call))


def reorient_nifit(file_path, new_orient):

    call = [f'fslorient -deleteorient {file_path}']
    os.system(' '.join(call))

    call = [f'fslswapdim {file_path} {new_orient} {file_path}']
    os.system(' '.join(call))

    call = [f'fslorient -setqformcode 1 {file_path}']
    os.system(' '.join(call))


def union_niftis(file_paths, output_path):

    # Load the NIfTI files
    images = [nib.load(path) for path in file_paths]

    # Get the data from the NIfTI files
    image_data = [img.get_fdata() for img in images]

    # Check that all images have the same shape
    shape = image_data[0].shape
    if not all(img.shape == shape for img in image_data):
        raise ValueError("All input images must have the same dimensions")

    # Perform the union operation (logical OR) on the images
    union_data = np.zeros(shape)
    for data in image_data:
        union_data = np.logical_or(union_data, data)

    # Convert the boolean mask back to the original data type
    union_data = union_data.astype(image_data[0].dtype)

    # Create a new NIfTI image with the union data
    union_img = nib.Nifti1Image(union_data, images[0].affine, images[0].header)

    # Save the resulting image
    nib.save(union_img, output_path)


def apply_mrconvert(input_path, output_path, axis_no, axis_low, axis_high):

    call = [f'mrconvert',
            f'{input_path}',
            f'{output_path}',
            f'-coord {str(axis_no) + " " + str(axis_low) + ":" + str(axis_high)}',
            f'-force']

    os.system(' '.join(call))


def convert_to_scans(input_path, init_paths, ext):

    input_img = nib.load(input_path).get_fdata()
    output_paths = replace_string(init_paths, '.nii.gz', ext + '.nii.gz')

    for jj in range(len(output_paths)):

        ref_nii = nib.load(init_paths[jj])
        ref_shape = ref_nii.shape
        no_vols = ref_shape[-1]

        temp_nifti = input_img[:, :, :, :no_vols]
        input_img = np.delete(input_img, np.s_[:no_vols], axis=-1)

        array_to_nii(ref_nii, temp_nifti, output_paths[jj])

    if input_img.shape[-1] != 0:
        raise ValueError("There is remaining data in the initial nifti!")

def concat_niftis(list_niftis, output_path, opt):
    ''' Concatenates volumes from a list of NIfTI files.
        If 'all' is selected, all volumes are concatenated.
        If an integer is provided, it selects up to that number of volumes from each NIfTI. '''

    nifti_objs = [nib.load(n) for n in list_niftis]
    volumes = []

    for nifti in nifti_objs:
        img_data = nifti.get_fdata()
        
        # Ensure all images are at least 4D (expand if necessary)
        if img_data.ndim == 3:
            img_data = np.expand_dims(img_data, axis=-1)

        if opt == 'all':
            volumes.append(img_data)
        elif isinstance(opt, int):
            volumes.append(img_data[:, :, :, :opt])
        else:
            raise ValueError("Please enter either 'all' or an integer number of volumes!")

    combined_nifti = np.concatenate(volumes, axis=3)
    array_to_nii(nifti_objs[0], combined_nifti, output_path)

def unconcat_niftis(nifti_path, delta_path):
    from collections import defaultdict

    # Load the 4D NIfTI file
    img = nib.load(nifti_path)
    data = img.get_fdata()
    affine = img.affine
    header = img.header

    # Load delta values
    deltas = [int(x) for x in read_numeric_txt(delta_path)[0]]
     
    if len(deltas) != data.shape[3]:
        raise ValueError("Number of delta values does not match number of 3D volumes in NIfTI file.")

    # Group indices by delta value
    delta_groups = defaultdict(list)
    for idx, delta in enumerate(deltas):
        delta_groups[int(float(delta))].append(idx)

    # Write separate NIfTI files for each delta group
    for delta, indices in delta_groups.items():
        sub_data = data[..., indices]
        new_img = nib.Nifti1Image(sub_data, affine, header)
        
        # Make folder per delta
        delta_folder = os.path.join(os.path.join(os.path.dirname(nifti_path)), f'Delta_{delta}')
        os.makedirs(delta_folder, exist_ok=True)
        
        name = os.path.basename(nifti_path).replace('allDelta-allb', f'Delta_{delta}')
        output_path = os.path.join(delta_folder, f'{name}')
        nib.save(new_img, output_path)
        print(f"Saved: {output_path}")
        
        bvecs_path = re.sub(r'(Delta_\d+_).*$', r'\1bvecsRotated.txt', output_path)
        bvals_path = re.sub(r'(Delta_\d+_).*$', r'\1bvalsNom.txt', output_path)

        nifti_to_mif(output_path, bvecs_path, bvals_path, output_path.replace('.nii.gz','.mif'))
        print(f"Converted to mif: {output_path}")
        
        mask_path         = nifti_path.replace('dwi_dn_gc_ec.nii.gz','mask.nii.gz')
        new_mask_path     = os.path.join(delta_folder,'mask.nii.gz')
        mask_dil_path     = nifti_path.replace('dwi_dn_gc_ec.nii.gz','mask_dil.nii.gz')
        new_mask_dil_path = os.path.join(delta_folder,'mask_dil.nii.gz')
        copy_files([mask_path,mask_dil_path],[new_mask_path,new_mask_dil_path])
        print(f"Copyied mask files: {output_path}")

        sigma_path         = nifti_path.replace('dwi_dn_gc_ec.nii.gz','dwi_dn_sigma.nii.gz')
        new_sigma_path     = os.path.join(delta_folder,'dwi_dn_sigma.nii.gz')
        copy_files([sigma_path],[new_sigma_path])
        print(f"Copyied sigma files: {output_path}")

        

def array_to_nii(in_img, in_array, out_img):

    if isinstance(in_img, nib.nifti1.Nifti1Image):
        img = in_img
    else:
        img = nib.load(in_img)

    aff, hdr = img.affine, img.header
    newimg = nib.Nifti1Image(in_array, affine=aff, header=hdr)
    nib.save(newimg, out_img)


def raw_to_nifti(input_path, output_path):

    if not os.listdir(output_path):

        # Convert data
        call = [os.environ['HOME']+f'/anaconda3/envs/Dicomifier/bin/dicomifier ',
                f'to-nifti',
                f'-z',
                f'{input_path}',
                f'{output_path}']
        os.system(' '.join(call))

        # Retrieve directory in which the converted data is
        data_dir_lvl1 = get_subdir(output_path)
        data_dir_lvl2 = get_subdir(data_dir_lvl1)

        # Move all files and subdirectories to the main scan path
        destination_path = os.path.abspath(
            os.path.join(data_dir_lvl2, '..', '..'))
        for item in os.listdir(data_dir_lvl2):
            item_path = os.path.join(data_dir_lvl2, item)
            destination_item_path = os.path.join(destination_path, item)
            shutil.move(item_path, destination_item_path)
        # Remove the empty source_path
        os.rmdir(data_dir_lvl2)
        os.rmdir(data_dir_lvl1)

    else:
        print('Data was already converted!')

##### REGISTRATION #####

def antsreg_full(fixed_path, moving_path, out_transform, lesion_mask_moving=None):
 
    out_im = out_transform + '.nii.gz'
    
    call = [f'antsRegistration -d 3 --interpolation Linear',
            f'--winsorize-image-intensities [0.01,0.99] --use-histogram-matching 1 ',
            f'--initial-moving-transform [{fixed_path}, {moving_path},1]',
            f'--transform Rigid[0.1] --convergence [1000x500x250x0,1e-6,10] --shrink-factors 12x8x4x1 --smoothing-sigmas 5x4x3x1vox ',
            f'--metric MI[{fixed_path}, {moving_path},1,32,Regular,0.25]',
            #f'--metric CC[{fixed_path}, {moving_path},0.5,4]',
            f'--transform Affine[0.15] --convergence [1000x500x250x0,1e-6,10] --shrink-factors 12x8x4x1 --smoothing-sigmas 5x4x3x1vox ',
            f'--metric MI[{fixed_path}, {moving_path},1,32,Regular,0.25]', \
            #f'--metric CC[{fixed_path}, {moving_path},0.5,4]' ,\
            f'--transform SyN[0.05,3,0] --convergence [200x100x50x20,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox ', \
            f'--metric MI[{fixed_path}, {moving_path},0.5,32,Random,0.25]' ,\
            f'--metric CC[{fixed_path}, {moving_path},0.5,4]']
           # f'-o [{out_transform},{out_im}] ']
        
    if lesion_mask_moving:
        call += ['-x', f'[,{lesion_mask_moving}]']

    call.append(f'-o [{out_transform},{out_im}] ')

    print(' '.join(call))
    os.system(' '.join(call))

# def antsreg(fixed_path, moving_path, out_transform):

#     out_im = out_transform + '.nii.gz'
    
#     call = [f'antsRegistration -d 3 --interpolation Linear',
#             f'--winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0 ',
#             f'--initial-moving-transform [{fixed_path}, {moving_path},1]',
#             f'--transform Rigid[0.1] --convergence [1000x500x250x0,1e-7,10] --shrink-factors 12x8x4x1 --smoothing-sigmas 5x4x3x1vox ',
#             f'--metric MI[{fixed_path}, {moving_path},0.5,32,Regular,0.25]',
#             #f'--metric CC[{fixed_path}, {moving_path},0.5,4]',
#             f'--transform Affine[0.15] --convergence [1000x500x250x0,1e-7,10] --shrink-factors 12x8x4x1 --smoothing-sigmas 5x4x3x1vox ',
#             f'--metric MI[{fixed_path}, {moving_path},1.25,32,Random,0.25]', \
#             f'--metric CC[{fixed_path}, {moving_path},0.5,4]' ,\
#             f'--transform SyN[0.1,4,0] --convergence [100x70x50x20,1e-7,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox ', \
#             f'--metric MI[{fixed_path}, {moving_path},1.25,32,Random,0.25]' ,\
#             f'--metric CC[{fixed_path}, {moving_path},1,4]', \
#             f'-o [{out_transform},{out_im}] ']

#     print(' '.join(call))
#     os.system(' '.join(call))
    
def antsreg_simple(fixed_path, moving_path, out_transform, lesion_mask_moving=None):

    out_im = out_transform + '.nii.gz'
    
    call = [f'antsRegistration -d 3 --interpolation Linear',
            f'--winsorize-image-intensities [0.005,0.995] --use-histogram-matching 1 ',
            f'--initial-moving-transform [{fixed_path}, {moving_path},1]',
            f'--transform Rigid[0.1] --convergence [1000x500x250x0,1e-6,10] --shrink-factors 12x8x4x1 --smoothing-sigmas 5x4x3x1vox ',
            f'--metric MI[{fixed_path}, {moving_path},1,32,Regular,0.25]',
            #f'--metric CC[{fixed_path}, {moving_path},0.5,4]',
            #f'--transform Affine[0.15] --convergence [1000x500x250x0,1e-6,10] --shrink-factors 12x8x4x1 --smoothing-sigmas 5x4x3x1vox ',
            #f'--metric MI[{fixed_path}, {moving_path},1,32,Regular,0.25]', \
            # f'--metric CC[{fixed_path}, {moving_path},0.5,4]' ,\
            #f'--transform SyN[0.1,4,0] --convergence [100x70x50x20,1e-7,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox ', \
            #f'--metric MI[{fixed_path}, {moving_path},1.25,32,Random,0.25]' ,\
            #f'--metric CC[{fixed_path}, {moving_path},1,4]', \
            ]
        
    if lesion_mask_moving:
        call += ['-x', f'[,{lesion_mask_moving}]']
  
    call.append(f'-o [{out_transform},{out_im}] ')
  
    print(' '.join(call))
    os.system(' '.join(call))
      
def antsreg_Affine(fixed_path, moving_path, out_transform):

    out_im = out_transform + '.nii.gz'
    
    call = [f'antsRegistration -d 3 --interpolation Linear',
            f'--winsorize-image-intensities [0.005,0.995] --use-histogram-matching 1 ',
            f'--initial-moving-transform [{fixed_path}, {moving_path},1]',
            f'--transform Rigid[0.1] --convergence [1000x500x250x0,1e-6,10] --shrink-factors 12x8x4x1 --smoothing-sigmas 5x4x3x1vox ',
            f'--metric MI[{fixed_path}, {moving_path},1,32,Regular,0.25]',
            #f'--metric CC[{fixed_path}, {moving_path},0.5,4]',
            f'--transform Affine[0.15] --convergence [1000x500x250x0,1e-6,10] --shrink-factors 12x8x4x1 --smoothing-sigmas 5x4x3x1vox ',
            f'--metric MI[{fixed_path}, {moving_path},1,32,Regular,0.25]', \
            # f'--metric CC[{fixed_path}, {moving_path},0.5,4]' ,\
            #f'--transform SyN[0.1,4,0] --convergence [100x70x50x20,1e-7,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox ', \
            #f'--metric MI[{fixed_path}, {moving_path},1.25,32,Random,0.25]' ,\
            #f'--metric CC[{fixed_path}, {moving_path},1,4]', \
            f'-o [{out_transform},{out_im}] ']

    print(' '.join(call))
    os.system(' '.join(call))
    
def antsreg_atlas(fixed_path, moving_path, out_transform):

    out_im = out_transform + '.nii.gz'

    #moving_path = new_path
    call = [f'antsRegistration -d 3 --interpolation Linear',
            f'--winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0 ',
            f'--initial-moving-transform [{fixed_path}, {moving_path},1]',
            f'--transform Rigid[0.1] --convergence [1000x500x250x0,1e-7,10] --shrink-factors 12x8x4x1 --smoothing-sigmas 5x4x3x1vox ',
            f'--metric MI[{fixed_path}, {moving_path},0.5,32,Regular,0.25]',
            #f'--metric CC[{fixed_path}, {moving_path},0.5,4]',
            f'--transform Affine[0.15] --convergence [1000x500x250x0,1e-7,10] --shrink-factors 12x8x4x1 --smoothing-sigmas 5x4x3x1vox ',
            f'--metric MI[{fixed_path}, {moving_path},1.25,32,Random,0.25]', \
            # f'--metric CC[{fixed_path}, {moving_path},0.5,4]' ,\
            f'--transform SyN[0.1,4,0] --convergence [100x70x50x20,1e-7,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox ', \
            f'--metric MI[{fixed_path}, {moving_path},1.25,32,Random,0.25]' ,\
            f'--metric CC[{fixed_path}, {moving_path},0.5,4]', \
            f'-o [{out_transform},{out_im}] --verbose ']

    print(' '.join(call))
    os.system(' '.join(call))
    
def antsreg_syn(fixed_path, moving_path, output_prefix, transformation):

    call = [f'antsRegistrationSyN.sh -d 3',
            f'-f {fixed_path}',
            f'-m {moving_path}',
            f'-o {output_prefix}',
            f'-t {transformation}']

    print(' '.join(call))
    os.system(' '.join(call))


        

def ants_apply_transforms_simple(input_path, ref_path, output_path, transf_1, extra=None):
    
    for ii in range(len(input_path)):
        input_temp = input_path[ii]
        output_temp = output_path[ii]

        call = [
            'antsApplyTransforms',
            '-d 3',
            f'-i {input_temp}', \
            f'-r {ref_path}', \
            f'-t {transf_1}', \
            f'-o {output_temp}']

        # Add more arguments if exists
        if extra is not None:
            if isinstance(extra, str):
               call.extend(extra.split())  # split the string into a list
            else:
               call.extend(extra)


        # Always verbose
        call.append('--verbose')

        # Run the command
        print(' '.join(call))
        os.system(' '.join(call))
        
def ants_apply_transforms_simple_4D(input_path, ref_path, output_path, transf_1, extra=None):
    
    for ii in range(len(input_path)):
        input_temp = input_path[ii]
        output_temp = output_path[ii]
    
        img_4d = nib.load(input_temp)
        data_4d = img_4d.get_fdata()
        affine = img_4d.affine
        header = img_4d.header

        # Ensure it's 4D
        if data_4d.ndim != 4:
            raise ValueError(f"Expected 4D image, got shape {data_4d.shape}")
    
        transformed_volumes = []
        for vol_idx in range(data_4d.shape[3]):
            vol_data = data_4d[..., vol_idx]
    
            with tempfile.NamedTemporaryFile(suffix=".nii.gz", delete=False) as temp_infile:
                temp_in = temp_infile.name
                nib.Nifti1Image(vol_data, affine, header).to_filename(temp_in)
    
            with tempfile.NamedTemporaryFile(suffix=".nii.gz", delete=False) as temp_outfile:
                temp_out = temp_outfile.name
    
            dim = 3  # Since each volume is 3D
            call = [
                'antsApplyTransforms',
                f'-d {dim}',
                f'-i {temp_in}',
                f'-r {ref_path}',
                f'-t {transf_1}',
                f'-o {temp_out}'
            ]
    
            if extra is not None:
                if isinstance(extra, str):
                    call.extend(extra.split())
                else:
                    call.extend(extra)
    
            call.append('--verbose')
            print(' '.join(call))
            os.system(' '.join(call))
    
            # Load transformed volume
            transformed_vol = nib.load(temp_out).get_fdata()
            transformed_volumes.append(transformed_vol)
    
            # Cleanup temp files
            os.remove(temp_in)
            os.remove(temp_out)
    
        # Stack back into 4D
        transformed_4d = np.stack(transformed_volumes, axis=3)
        nib.Nifti1Image(transformed_4d, affine, header).to_filename(output_temp)
    
    

def ants_apply_transforms(input_path, ref_path, output_path, transf_1, transf_2, extra=None):  # input_type

    for ii in range(len(input_path)):

        input_temp = input_path[ii]
        output_temp = output_path[ii]

        call = [f'antsApplyTransforms',
                f'-d 3', \
                # f'-e {input_type}', \
                f'-i {input_temp}', \
                f'-r {ref_path}', \
                f'-t {transf_1}', \
                f'-t {transf_2}', \
                f'-o {output_temp}']
        
        # Add more arguments if exists
        if extra is not None:
            if isinstance(extra, str):
               call.extend(extra.split())  # split the string into a list
            else:
               call.extend(extra)


        # Always verbose
        call.append('--verbose')

        # Run the command
        print(' '.join(call))
        os.system(' '.join(call))

##### MODELS #####

def linear_model(b, m, b_int):
    return m * b + b_int

def get_param_names_model(model, is_alive):
    
    if model=='Nexi':
        if is_alive=='ex_vivo':
            patterns = ["*nexi*t_ex*", "*nexi*di*","*nexi*de*","*nexi*f*"]
            lims     = [(0, 50), (0, 2), (0, 2),  (0, 0.4)]
            maximums = np.array([[1, 80], [0.1, 2], [0, 2], [0.1, 0.9]])
        else:
            patterns = ["*nexi*t_ex*", "*nexi*di*","*nexi*de*","*nexi*f*"]
            lims     = [(0, 100), (0, 3.5), (0, 3.5),  (0, 1)]
            maximums = np.array([[1, 80], [0.1, 3.5], [0.1, 3.5], [0.1, 0.9]])
    
    elif model=='Smex':
        if is_alive=='ex_vivo':
            patterns = ["*smex*t_ex*", "*smex*di*","*smex*de*","*smex*f*"]
            lims     = [(0, 50), (0, 2), (0, 2),  (0, 0.4)]
            maximums = np.array([[1, 80], [0.1, 2], [0, 2], [0.1, 0.9]])
        else:
            patterns = ["*smex*t_ex*", "*smex*di*","*smex*de*","*smex*f*"]
            lims     = [(0, 80), (0, 3.5), (0, 2),  (0, 0.85)]
            maximums = np.array([[1, 80], [0.1, 3.5], [0.1, 3.5], [0.1, 0.9]])
    
    elif model=='Sandi':
        if is_alive=='ex_vivo':
            patterns = ["*sandi*rs*","*sandi*fs*", "*sandi*di*","*sandi*de*","*sandi*f*"]
            lims = [(0, 25), (0,0.5), (0, 2), (0, 2),  (0, 0.9)]
            maximums = np.full((len(patterns), 2), np.inf)
            maximums[:, 0] = -np.inf  
        else:
            patterns = ["*sandi*rs*","*sandi*fs*", "*sandi*di*","*sandi*de*","*sandi*f*"]
            lims = [(0, 25), (0,0.5), (0, 3.5), (0, 3.5),  (0, 0.9)]
            maximums = np.full((len(patterns), 2), np.inf)
            maximums[:, 0] = -np.inf  
        
    elif model=='Sandix':
        if is_alive=='ex_vivo':
             patterns = ["*sandix*t_ex*", "*sandix*di*","*sandix*de*","*sandix*f*","*sandix*rs*","*sandix*fs*"]
             lims     = [(0, 50), (0, 2), (0, 2),  (0, 0.4), (0, 25), (0,0.5)]
             maximums = np.full((len(patterns), 2), np.inf)
             maximums[:, 0] = -np.inf 
        else:
            patterns = ["*sandix*t_ex*", "*sandix*di*","*sandix*de*","*sandix*f*","*sandix*rs*","*sandix*fs*"]
            lims     = [(0, 100), (0, 3.5), (0, 3.5),  (0, 1), (0, 25), (0,0.5)]
            maximums = np.full((len(patterns), 2), np.inf)
            maximums[:, 0] = -np.inf 
        
    elif model=='SMI':
        patterns = ["*Da*", "*DePar*", "*DePerp*", "*f*", "*fw*", "*p2*", "*p4*"]
        lims = [(0, 4), (0, 4), (0, 4),  (0, 0.85), (0, 3), (0, 0.5), (0,0.5)]
        maximums = np.full((len(patterns), 2), np.inf)
        maximums[:, 0] = -np.inf  

    elif model=='DTI_DKI':
        if is_alive=='ex_vivo':
            patterns = ['*md_dki*','*mk_dki*','*fa_dki*']
            lims = [(0.5, 1.5), (0.2, 0.8), (0, 0.2)]
            #maximums = np.full((len(patterns), 2), np.inf)
            #maximums[:, 0] = -np.inf 
            maximums = np.array([[0, 5], [0, 50], [0, 50]])
        else:
            patterns = ['*md_dki*','*mk_dki*','*fa_dki*']
            lims = [(0, 2), (0, 2), (0, 1)]
            #maximums = np.full((len(patterns), 2), np.inf)
            #maximums[:, 0] = -np.inf 
            maximums = np.array([[0, 3], [0, 3], [0, 1]])
            
    elif model=='Micro_FA':
            patterns = ['*microFA*','*MD*']
            lims = [(0, 1), (0, 3)]
            maximums = np.array([[0, 1], [0, 3]])

                
    
    return patterns, lims, maximums

def run_script_in_conda_environment(script_path,env_name):
    subprocess.run(f"""conda init
        source ~/.bashrc
        source activate base
        conda activate """+env_name+f"""
        python """+script_path,
            shell=True, executable='/bin/bash', check=True)


def create_mrs_dyn_config(diffusion_model, path, cfg):
    """
    Creates the config file needed for dynamic fitting in FSL MRS.
    diffusion_model: chosen diffusion model S(b), currently implemted:
        'callaghan': Callaghan model of randomly oriented sticks.
        'dki': diffusion kurtosis imaging signal representation.
        'biexp': biexponential model provided by FSL MRS.
    path: path where to store the config file.
    cfg: dictionary with the configuration parameters.
    """
    f = open(cfg['common_folder']+'/mrs_dyn_param_' + diffusion_model + '.py', 'r')
    model_parametrization = f.read()
    f = open(cfg['common_folder']+'/mrs_dyn_models.py', 'r')
    models = f.read()
    configtxt = model_parametrization + '\n\n'+models
    f = open(path, 'w')
    f.write(configtxt)
    f.close()

def create_directory(directory_name):
     try:
         if not os.path.exists(directory_name):
            os.makedirs(directory_name)
     except FileExistsError:
         print(f"Warn: Directory '{directory_name}' already exists.")
     except PermissionError:
         print(f"Err: Permission denied: Unable to create '{directory_name}'.")
     except Exception as e:
        print(f"An error occurred: {e}")
