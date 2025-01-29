import os
import re
import shutil
import numpy as np
import nibabel as nib
# import nilearn.image as nlrn
import matplotlib.pyplot as plt
import nibabel as nib
from more_itertools import locate
import sys
import tkinter as tk
from tkinter import ttk
from scipy.ndimage import label, find_objects
import shutil
import distinctipy

##### FILES AND SYSTEM OPERATIONS #####

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


def create_directory(path):

    if not os.path.exists(path):
        os.makedirs(path)

def delete_directory_sudo(path):
    
    call = [f' sudo rm -r {path}']
    os.system(' '.join(call))
    
def copy_file(source_paths, destination_paths):

    for source, destination in zip(source_paths, destination_paths):
        shutil.copyfile(source, destination)

def copy_files_BIDS(bids_strc, output_path, filename):

     new_path = os.path.join(output_path,bids_strc.get_param("base_name")+filename)
     copy_file([bids_strc.get_path(filename)], 
               [new_path])
     
     return new_path


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



##### TOPUP #####


def topup_routine(path_fwd, path_rev, dwi_path, bids_strc, topupcfg_path):

    topup_input_files = create_topup_input_files(
        path_fwd, path_rev, bids_strc, topupcfg_path)
    do_topup(topup_input_files)
    apply_topup(topup_input_files, dwi_path, bids_strc)


def create_topup_input_files(paths_fwd, paths_rev, bids_strc, topupcfg_path):

    topup_input_files = {}

    # Forward and reverse b0 images
    concat_niftis(paths_fwd, bids_strc.get_path('b0_fwd.nii.gz'), 1)
    concat_niftis(paths_rev, bids_strc.get_path('b0_rev.nii.gz'), 1) # assumes only one B0 value was collected in rev direction
    concat_niftis([bids_strc.get_path('b0_fwd.nii.gz'), bids_strc.get_path('b0_rev.nii.gz')],
                  bids_strc.get_path('b0_fwd_rev.nii.gz'), 'all')

    topup_input_files['b0_fwd_rev'] = bids_strc.get_path('b0_fwd_rev.nii.gz')

    # Acqp file
    nii_fwd = nib.load(bids_strc.get_path('b0_fwd.nii.gz')).get_fdata()
    nii_rev = nib.load(bids_strc.get_path('b0_rev.nii.gz')).get_fdata()
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
    copy_files([topupcfg_path], [bids_strc.get_path('mycnf_fmri.cnf')])
    topup_input_files['config'] = bids_strc.get_path('mycnf_fmri.cnf')
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

    call = [f'applytopup ',
            f'--imain={imain}',
            f'--datain={datain}',
            f'--inindex={inindex}',
            f'--topup={topup}',
            f'--method=jac',
            f'--out={out}']

    print(' '.join(call))
    os.system(' '.join(call))

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

    call = [f'eddy_cuda10.2',
            f'--imain={dwi}',
            f'--mask={mask}',
            f'--index={index}',
            f'--acqp={acqp}',
            f'--bvecs={bvecs}',
            f'--bvals={bvals}', \
            f'--slm=linear', \
            f'--out={output}', \
            f'--data_is_shelled --verbose']

    if eddy_input_files.get('topup'):
        topup = eddy_input_files['topup']
        call.insert(6, f'--topup={topup}')

    print(' '.join(call))
    os.system(' '.join(call))

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
    
    call = [f'fsleyes render --hideCursor --hidex --hidez  --voxelLoc 60 12 50 ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30',
            f'--outfile {png_path}',
            f'{FA}',
            f'{V1} --overlayType rgbvector '
            f'--modulateImage {FA}']
    
    print(' '.join(call))
    os.system(' '.join(call))

def QA_brain_extract(anat_path,output_path):
    
    create_directory(output_path)
    
    png_path = os.path.join(output_path, 'T2W.png')
    call = [f'fsleyes render --hideCursor --hidex --hidez  --voxelLoc 60 13 50 ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{anat_path}',
            f'-dr 0 40000 ',]

    print(' '.join(call))
    os.system(' '.join(call))
    
    
    anat_brain_path = anat_path.replace('.nii.gz','_brain.nii.gz')
    png_path = os.path.join(output_path, 'T2W_brain.png')
    call = [f'fsleyes render --hideCursor --hidex --hidez  --voxelLoc 60 13 50 ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{anat_brain_path}',
            f'-dr 0 40000 ',]

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


    png_path = os.path.join(output_path, 'T2W_with_T2wbrain.png')
    call = [f'fsleyes render --hideCursor --hidex --hidez  --voxelLoc 60 13 50 ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{anat_path}',
            f'-dr 0 40000 ',
            f'{countour_path}']

    print(' '.join(call))
    os.system(' '.join(call))
  

def QA_denoise(bids_strc, res, sigma, output_path):
    
    res_path = bids_strc.get_path(res)
    sigma_path = bids_strc.get_path(sigma)

    create_directory(output_path)
    
    png_path = os.path.join(output_path, res.replace('.nii.gz','.png'))
    call = [f'fsleyes render --hideCursor --hidex --hidez  --voxelLoc 60 12 50 ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{res_path} ',
            f'-dr -1000 1000']

    print(' '.join(call))
    os.system(' '.join(call))
    
    png_path = os.path.join(output_path, sigma.replace('.nii.gz','.png'))
    call = [f'fsleyes render --hideCursor --hidex --hidez  --voxelLoc 60 12 50 ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{sigma_path} --cmap red-yellow',
            f'-dr 0 600']

    print(' '.join(call))
    os.system(' '.join(call))
    

def QA_topup(bids_strc, before, after, output_path):
    
    before_path = bids_strc.get_path(before)
    after_path = bids_strc.get_path(after)
    
    create_directory(output_path)
    
    png_path = os.path.join(output_path, before.replace('.nii.gz','.png'))
    call = [f'fsleyes render --hideCursor --hidex --hidez  --voxelLoc 60 10 50 ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{before_path} ',
            f'-dr -1000 100000 ']

    print(' '.join(call))
    os.system(' '.join(call))
    
    png_path = os.path.join(output_path, after.replace('.nii.gz','.png'))
    call = [f'fsleyes render --hideCursor --hidex --hidez  --voxelLoc 60 10 50 ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{after_path} ',
            f'-dr -1000 100000 ']

    print(' '.join(call))
    os.system(' '.join(call))
    
    
def QA_gc(bids_strc, before, after, output_path):
    
    before_path = bids_strc.get_path(before)
    after_path = bids_strc.get_path(after)
    
    create_directory(output_path)
    
    png_path = os.path.join(output_path, before.replace('.nii.gz','.png'))
    call = [f'fsleyes render --hideCursor --hidex --hidez  --voxelLoc 60 10 50 ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{before_path}',
            f'-dr -1000 100000 ']

    print(' '.join(call))
    os.system(' '.join(call))
    
    png_path = os.path.join(output_path, after.replace('.nii.gz','.png'))
    call = [f'fsleyes render --hideCursor --hidex --hidez  --voxelLoc 60 10 50 ',
            f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
            f'--outfile {png_path}',
            f'{after_path} ',
            f'-dr -1000 100000 ']

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
    b0s = np.where(bvals == 1000)

    for ii in np.linspace(1, len(b0s[0])-1, num=3, dtype='int'):
        volume = b0s[0][ii]

        # plot dwi before eddy
        png_path = os.path.join(
            output_path, 'nodifcontour_v' + str(volume) + '_dwi.png')
        call = [f'fsleyes render --hideCursor --hidex --hidez  --voxelLoc 60 12 50 ',
                f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
                f'--outfile {png_path}',
                f'{dwi_path}',
                f'-v {volume}',
                f'-dr 0 40000 ',
                f'{countour_path}']

        print(' '.join(call))
        os.system(' '.join(call))

        # plot dwi after eddy
        png_path = os.path.join(
            output_path, 'nodifcontour_v' + str(volume) + '_eddy.png')
        call = [f'fsleyes render --hideCursor --hidex --hidez  --voxelLoc 60 12 50 ',
                f'--xcentre -0 0 --ycentre -0 0 --zcentre -0 0 --labelSize 30 ',
                f'--outfile {png_path}',
                f'{dwi_ec_path}',
                f'-v {volume}',
                f'-dr 0 40000',
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


def QA_plotSNR(bids_strc, snr_path, nf_path, mask_path, bvals_path, output_path):

    snr_path = bids_strc.get_path(snr_path)
    nf_path = bids_strc.get_path(nf_path)
    mask_path = bids_strc.get_path(mask_path)
    bvals_path = bids_strc.get_path(bvals_path)

    create_directory(output_path)

    mask = nib.load(mask_path).get_fdata()
    nf = np.multiply(nib.load(nf_path).get_fdata(), mask)
    SNR = nib.load(snr_path).get_fdata()

    for v in range(SNR.shape[-1]):
        SNR[:, :, :, v] = np.multiply(SNR[:, :, :, v], mask)

    data = SNR.reshape(SNR.shape[0]*SNR.shape[1]*SNR.shape[2], SNR.shape[3]);
    data[data == 0] = np.nan

    # filename    = get_filename(snr_path)
    bvals_nom = read_numeric_txt(bvals_path)

    fig, ax = plt.subplots()
    ax.plot(np.transpose(bvals_nom), np.nanmean(data, axis=0), 'bo', markersize=3)

    ax.plot(np.transpose(bvals_nom), np.repeat(
        np.mean(nf), np.transpose(bvals_nom).shape[0]))

    # Set axes
    ax.set_xlabel('Nominal b-val',
                  fontdict={'size': 12, 'weight': 'bold', 'style': 'italic'})
    ax.set_ylabel('Mean SNR', fontdict={
                  'size': 12, 'weight': 'bold', 'style': 'italic'})
    ax.grid(True)
    plt.savefig(os.path.join(output_path, 'QA_SNR.png'),
                bbox_inches='tight', dpi=300)

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


def extract_methods(methods_in, bids_strc, acqp):

    if acqp == 'diff':

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


def dilate_im(input_path, output_path, sigma):

    nii_shape = nib.load(input_path).shape
    nii_ndims = len(nii_shape)

    call = [f'/home/localadmin/SOFTWARES/ants-2.5.3/bin/ImageMath ',
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

    call = [f'/home/localadmin/SOFTWARES/ants-2.5.3/bin/ThresholdImage',
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

    call = [f'dwidenoise',
            f'{input_path}',
            f'{output_path}',
            f'-noise {noise_path}',
            f'-debug',
            f'-force -info']
    os.system(' '.join(call))

    res_path = output_path.replace('.nii.gz', '_res.nii.gz')
    call = [f'mrcalc',
            f'{input_path}',
            f'{output_path}',
            f'-subtract',
            f'{res_path}']
    os.system(' '.join(call))


def denoise_vols(input_path, kernel_size, ouput_path, noise_path):

    call = [f'dwidenoise',
            f'{input_path}',
            f'{ouput_path}',
            f'-extent {kernel_size}',
            f'-noise {noise_path}',
            f'-debug',
            f'-force']
    os.system(' '.join(call))


def denoise_vols_mask(input_path, kernel_size, mask_path, ouput_path, noise_path):

    call = [f'dwidenoise',
            f'{input_path}',
            f'{ouput_path}',
            f'-extent {kernel_size}',
            f'-noise {noise_path}',
            f'-mask {mask_path}',
            f'-debug',
            f'-force']

    os.system(' '.join(call))


def denoise_img(input_path, dim, ouput_path):

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
            
def estim_DTI_DKI_designer(input_mif, 
                           mask_path, output_path, data_path):

    call = [f'docker run -v {data_path}:/data nyudiffusionmri/designer2:v2.0.10 tmi -SMI',
            f'{input_mif}',
            f'{output_path}',
            f'-mask {mask_path}']

    print(' '.join(call))
    os.system(' '.join(call))

def denoise_designer(input_path, output_path, data_path):

    docker_path = '/data'
    input_path  = input_path.replace(data_path,docker_path)
    output_path  = output_path.replace(data_path,docker_path)

    call = [f'docker run -v {data_path}:/data nyudiffusionmri/designer2:v2.0.10 designer -denoise',
            f'{input_path}',
            f'{output_path}  -pf 0.75 -pe_dir i -algorithm veraart -extent 9,9,9 ']

    print(' '.join(call))
    os.system(' '.join(call))

def estim_SMI_designer(input_mif, mask_path, sigma_path, output_path, data_path):

    call = [f'docker run -v {data_path}:/data nyudiffusionmri/designer2:v2.0.10 tmi -SMI',
            f'{input_mif}',
            f'{output_path}',
            f'-sigma {sigma_path}',
            f'-mask {mask_path}']

    print(' '.join(call))
    os.system(' '.join(call))
    
    call = [f'docker run -v {data_path}:/data nyudiffusionmri/designer2:v2.0.10 chmod -R 777 {output_path}']
    print(' '.join(call))
    os.system(' '.join(call))
    
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

        call = [f'/home/localadmin/SOFTWARES/ants-2.5.3/bin/antsMotionCorr',
                f'-d {dim}',
                f'-a {input}',
                f'-o {output}']

        os.system(' '.join(call))


def calc_snr(input_path1, input_path2, output_path):
    binary_op(input_path1, input_path2, '-div', output_path)


def calc_noise_floor(input_path1, output_path):
    a = np.sqrt(np.pi/2)
    call = [f'fslmaths',
            f'{input_path1}',
            f'-sqrt',
            f' {output_path}']

    os.system(' '.join(call))

    binary_op(output_path, a, '-mul', output_path)


def brain_extract_RATS(input_path):

    RATS_path = '/home/localadmin/SOFTWARES/Rodent_Seg/distribution2/RATS_MM'
 
    call = [f'{RATS_path}',
            f'{input_path}',
            f'{input_path.replace(".nii.gz", "_brain_mask.nii.gz")}',
            f'-t 1500']
    

    print(' '.join(call))
    os.system(' '.join(call))
    
    binary_op(input_path,input_path.replace(".nii.gz", "_brain_mask.nii.gz"), '-mul', input_path.replace(".nii.gz", "_brain.nii.gz"))
    

def brain_extract_BREX(input_path):

    BREX_path = '/home/localadmin/SOFTWARES/atlasBREX-master/'
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


##### NIFTI HANDLE #####


def nifti_to_mif(nifti_path, bvecs_path, bvals_path, mif_path):

    call = [f'mrconvert',
            f'{nifti_path}',
            f'-fslgrad {bvecs_path} {bvals_path} ',
            f'{mif_path} -force']

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
    ''' The function concatenates either all volumes in the list of input niftis
        either an integer number of volumes for each of the niftis in the list '''

    if opt == 'all':
        combined_nifti = nib.concat_images(
            list_niftis, check_affines=True,axis=3)
        nib.save(combined_nifti, output_path)
    elif type(opt) == int:
        volumes = []
        for ii in range(len(list_niftis)):
            temp_nifti = nib.load(list_niftis[ii])
            temp_img = temp_nifti.get_fdata()
            if len(temp_img.shape)==4: # if there is more than one volume take the first opt volumes
                temp_nifti_vols = temp_img[:, :, :, :opt]
            elif len(temp_img.shape)==3: # if there is only one volume already take that one
                temp_nifti_vols = np.expand_dims(temp_img[:, :, :],axis=-1) 
            volumes.append(temp_nifti_vols)
        combined_nifti = np.concatenate(volumes, axis=3)

        array_to_nii(temp_nifti, combined_nifti, output_path)
    else:
        raise ValueError(
            'Please enter either all or an integer number of volumes!')


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
        call = [f'/home/localadmin/anaconda3/envs/Dicomifier/bin/dicomifier ',
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

def antsreg(fixed_path, moving_path, out_transform):

    out_im = out_transform + '.nii.gz'
    
    call = [f'/home/localadmin/SOFTWARES/ants-2.5.3/bin/antsRegistration -d 3 --interpolation Linear',
            f'--winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0 ',
            f'--initial-moving-transform [{fixed_path}, {moving_path},1]',
            f'--transform Rigid[0.1] --convergence [1000x500x250x0,1e-7,10] --shrink-factors 12x8x4x1 --smoothing-sigmas 5x4x3x1vox ',
            f'--metric MI[{fixed_path}, {moving_path},0.5,32,Regular,0.25]',
            #f'--metric CC[{fixed_path}, {moving_path},0.5,4]',
            f'--transform Affine[0.15] --convergence [1000x500x250x0,1e-7,10] --shrink-factors 12x8x4x1 --smoothing-sigmas 5x4x3x1vox ',
            f'--metric MI[{fixed_path}, {moving_path},1.25,32,Random,0.25]', \
            f'--metric CC[{fixed_path}, {moving_path},0.5,4]' ,\
            f'--transform SyN[0.1,4,0] --convergence [100x70x50x20,1e-7,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox ', \
            f'--metric MI[{fixed_path}, {moving_path},1.25,32,Random,0.25]' ,\
            f'--metric CC[{fixed_path}, {moving_path},1,4]', \
            f'-o [{out_transform},{out_im}] ']

    print(' '.join(call))
    os.system(' '.join(call))
    
def antsreg_simple(fixed_path, moving_path, out_transform):

    out_im = out_transform + '.nii.gz'
    
    call = [f'/home/localadmin/SOFTWARES/ants-2.5.3/bin/antsRegistration -d 3 --interpolation Linear',
            f'--winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0 ',
            f'--initial-moving-transform [{fixed_path}, {moving_path},1]',
            f'--transform Rigid[0.1] --convergence [1000x500x250x0,1e-7,10] --shrink-factors 12x8x4x1 --smoothing-sigmas 5x4x3x1vox ',
            f'--metric MI[{fixed_path}, {moving_path},0.5,32,Regular,0.25]',
            #f'--metric CC[{fixed_path}, {moving_path},0.5,4]',
            f'--transform Affine[0.15] --convergence [1000x500x250x0,1e-7,10] --shrink-factors 12x8x2x1 --smoothing-sigmas 5x4x3x0vox ',
            f'--metric MI[{fixed_path}, {moving_path},1.25,32,Random,0.25]', \
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
    call = [f'/home/localadmin/SOFTWARES/ants-2.5.3/bin/antsRegistration -d 3 --interpolation Linear',
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
    
# def antsreg_simple(moving_path, fixed_path, out_transform):

#     out_im = out_transform + '.nii.gz'


#     call = [f'/home/localadmin/SOFTWARES/ants-2.5.3/bin/antsRegistration -d 3 --interpolation Linear',
#             f'--winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0 ',
#             f'--initial-moving-transform [{fixed_path}, {moving_path},1]',
#             f'--transform Rigid[0.1] --convergence [1000x500x250x0,1e-7,10] --shrink-factors 12x8x4x1 --smoothing-sigmas 5x4x3x1vox --masks [NULL,NULL]',
#             f'--metric MI[{fixed_path}, {moving_path},0.5,32,Regular,0.25]',
#             f'--transform Affine[0.15] --convergence [1000x500x250x0,1e-7,10] --shrink-factors 12x8x4x1 --smoothing-sigmas 5x4x3x1vox --masks [NULL,NULL]',
#             f'--metric MI[{fixed_path}, {moving_path},1.25,32,Random,0.25]', \
#             f'-o [{out_transform},{out_im}] ']

#     print(' '.join(call))
#     os.system(' '.join(call))


def antsreg_syn(fixed_path, moving_path, output_prefix, transformation):

    call = [f'antsRegistrationSyN.sh -d 3',
            f'-f {fixed_path}',
            f'-m {moving_path}',
            f'-o {output_prefix}',
            f'-t {transformation}']

    print(' '.join(call))
    os.system(' '.join(call))


def ants_apply_transforms_simple(input_path, ref_path, output_path, transf_1):  # input_type

    for ii in range(len(input_path)):

        input_temp = input_path[ii]
        output_temp = output_path[ii]

        call = [f'antsApplyTransforms',
                f'-d 3', \
                # f'-e {input_type}', \
                f'-i {input_temp}', \
                f'-r {ref_path}', \
                f'-t {transf_1}', \
                f'-o {output_temp}']
        print(' '.join(call))
        os.system(' '.join(call))
        
def ants_apply_transforms(input_path, ref_path, output_path, transf_1, transf_2):  # input_type

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
        print(' '.join(call))
        os.system(' '.join(call))

##### MODELS #####

def get_param_names_model(model):
    
    if model=='Nexi':
        patterns = ["*t_ex*", "*di*","*de*","*f*"]
        lims = [(0, 100), (0, 4), (0, 2),  (0, 0.85)]
        
    elif model=='Sandi':
        patterns = ["*fs*", "*di*","*de*","*f*"]
        lims = [(0, 0.1), (0, 4), (0, 2),  (0, 0.85)]
        
    elif model=='SMI':
        patterns = ["*Da*", "*DePar*", "*DePerp*", "*f*", "*fw*", "*p2*", "*p4*"]
        lims = [(0, 4), (0, 4), (0, 4),  (0, 0.85), (0, 3), (0, 0.5), (0,0.5)]
    
    return patterns, lims
    