import numpy as np
import dipy.reconst.dti as dti
import dipy.reconst.dki as dki
from dipy.core.gradients import gradient_table
import nibabel as nib
import pandas as pd
import math
import os
import warnings
import plotly.graph_objs as go
from plotly.offline import plot
import glob
import matplotlib.pyplot as plt

giro_ratio = 2.6751525e8 # Gyromagnetic radio [rad/(s*T)]

def get_bvals(scheme_file_path):
    """
    Function that gets b-values from scheme file

    Args:
        scheme_file_path (str) : path of the scheme file

    Returns:
        (np.ndarray) : b-values [ms/um²]
    """

    b_values = []  # Initialize an empty list to store the values from the 4th column

    try:
        with open(scheme_file_path, 'r') as file:
            for line in file:
                # Split each line into space-separated values and extract the 4th column
                columns = line.strip().split()
                if len(columns) >= 5:
                    G     = float(columns[3]) * 1e-6 # [T/um]
                    giro  = giro_ratio * 1e-3 # Gyromagnetic radio [rad/(ms*T)]
                    delta = float(columns[5]) * 1e3 # [ms]
                    Delta = float(columns[4]) * 1e3 # [ms]
                    b     = pow(G * giro * delta, 2) * (Delta - delta/3) # [ms/um²]
                    b_values.append(b)  # Assuming columns are 0-based

    except FileNotFoundError:
        print(f"File not found: {scheme_file_path}")
    
    return np.array(b_values)

def get_bvecs(scheme_file_path):
    """
    Function that gets b-vectors from scheme file

    Args:
        scheme_file_path (str) : path of the scheme file

    Returns:
        (np.ndarray) : Unitary b_vectors, 3D
    """

    b_vectors = []  # Initialize an empty list to store the values from the first three columns

    try:
        with open(scheme_file_path, 'r') as file:
            for line in file:
                # Split each line into space-separated values and extract the first three columns
                b_vec = line.strip().split()[:3]  # Assuming columns are 0-based
                # Skip the header
                if len(b_vec) == 3:
                    b_value = float(line.strip().split()[3])
                    if b_value != 0:
                        b_vec = [float(val) for val in b_vec]  # Convert to float if needed
                        b_vectors.append(b_vec)
                    else :
                        b_vectors.append([0 ,0, 0])

    except FileNotFoundError:
        print(f"File not found: {scheme_file_path}")
    
    return np.array(b_vectors)

def get_dwi(dwi_path):
    """
    Function that gets the DWI values from the simulation output file
    
    Args:
        dwi_path (pathlib.PoxisPath) : path of the output DWI file
        
    Returns:
        (np.ndarray) : DWI values
    """

    if ".bfloat" in str(dwi_path):
        return np.fromfile(dwi_path, dtype="float32")
    elif ".txt" in str(dwi_path):
        signal = []
        with open(dwi_path) as f:
            [signal.append(float(line)) for line in f.readlines()]
        return np.array(signal)

def array_to_nifti(dwi_array):

    # Create an empty 4x4 affine matrix with ones on the diagonal
    affine = np.eye(4)

    img = nib.Nifti1Image(dwi_array, affine)
    return img

def calculate_DKI(scheme_file_path, dwi):
    """
    Function that calculates DTI / DKI metrics. 
    For DTI, 2 b-vals (<= 1) and 6 directions are needed. 
    For DKI, 3 b-vals (<= 2-3) and 21 directions are needed
    
    Args:
        scheme_file_path (str) : path of the scheme file
        dwi (np.ndarray)       : 4D DWI image

    Returns:
        FA (np.float64) : Fractional Anisotropy
        MD (np.float64) : Mean diffusivity
        AD (np.float64) : Axial diffusivity
        RD (np.float64) : Radial diffusivity
        MK (np.float64) : Mean Kurtosis
        AK (np.float64) : Axial Kurtosis
        RK (np.float64) : Radial Kurtosis

    """

    # Get b-values <= 1 [ms/um^2] (DTI has Gaussian assumption => small b needed)
    bvalues = get_bvals(scheme_file_path)      
    bvecs = get_bvecs(scheme_file_path)
   
    
    # DTI fit
    idx       = bvalues <= 1
    bvals_dti = bvalues[idx]  
    bvecs_dti = bvecs[idx]    
    gtab      = gradient_table(bvals_dti, bvecs_dti)
    # build model
    dkimodel  = dki.DiffusionKurtosisModel(gtab)
    # Create an empty 4x4 affine matrix with ones on the diagonal
    affine = np.eye(4)
    dwi_nii   = nib.Nifti1Image(dwi[idx], affine)
    dkifit    = dkimodel.fit(dwi_nii.get_fdata())
    # save maps
    FA = dkifit.fa
    MD = dkifit.md
    AD = dkifit.ad
    RD = dkifit.rd

    # DKI fit
    idx       = bvalues <= 3
    bvals_dki = bvalues[idx]  
    bvecs_dki = bvecs[idx]    
    gtab      = gradient_table(bvals_dki, bvecs_dki)
    # build model
    dkimodel  = dki.DiffusionKurtosisModel(gtab)
    dki_nii   = nib.Nifti1Image(dwi[idx], affine)
    dkifit    = dkimodel.fit(dki_nii.get_fdata())
    MK = dkifit.mk(0, 10)
    AK = dkifit.ak(0, 10)
    RK = dkifit.rk(0, 10)

    return FA, MD, AD, RD, MK, AK, RK

def create_conf_MCSim(N,T,dur,Di,De,scheme_name,out_path,substract,simulator_folder,vx_size):
    
    conf_filename   = os.path.join(simulator_folder,'instructions','conf','model_test.conf')
    scheme_filename = os.path.join(simulator_folder,'instructions','scheme', scheme_name)

    N = str(N)
    vx_size= vx_size*1e-3
    
    call = [f'N {N}',
            f'T {T} ',
            f'duration {dur}',
            f'diffusivity_intra {Di}',
            f'diffusivity_extra {De}',
            '',
            f'scheme_file {scheme_filename}',
            f'exp_prefix {out_path}',
            '',
            f'<obstacle>',
            f'<axons_list>',
            f'{substract}',
            f'permeability global 0',
            f'</axons_list>',
            f'</obstacle>',
            f'ini_walkers_pos intra',
            '',
            f'scale_from_stu 1',
            f'write_txt 0',
            f'write_bin 1',
            f'write_traj_file 1',
            f'num_process 10',
            '',
            f'<voxel>',
            f'0 0 0',
            f'{vx_size} {vx_size} {vx_size}',
            f'</voxel>',
            '',
            f'<END>']
    
    with open(conf_filename, 'w') as file:
        for content in call:
            file.write(content + '\n')  # Write each entry followed by a newline
            print(f"Written: {content}")  # Optional: Print confirmation for each line

    print('Finished writing conf file')
    
   
def run_sim(simulator_folder):
    
    filename        = os.path.join(simulator_folder,'run_mcds.sh')
    conf_filename   = os.path.join(simulator_folder,'instructions','conf','model_test.conf')

    with open(filename, 'r') as file:
        lines = file.readlines()

    # Look for the line containing "MC-DC_Simulator" and modify the content after it
    for i, line in enumerate(lines):
        if "MC-DC_Simulator" in line:
            # Replace the content of the line
            lines[i] = f"./MC-DC_Simulator {conf_filename}\n"
            break  
    
    # Write the updated lines back to the file
    with open(filename, 'w') as file:
        file.writelines(lines)
        
    os.chdir(simulator_folder)
    
    print('Running simulator ...')
    os.system("./run_mcds.sh")

def plot_traj(output_folder, sub_file):
    """
    Function that plots the trajectory and some example substracts 
    Args:
        output_folder (str)  : path of the output folder
        sub_file (str)      : path of the substract_file file
    """
    
    ## Trajectory
    traj_file = glob.glob(os.path.join(output_folder, "*.traj"))[0]
    
    lines = np.fromfile(traj_file, dtype="float32")
    xp = []
    yp = []
    zp = []
    for i in range(int(len(lines))):
        if i%3 == 0:
            xp.append(float(lines[i]))
        elif i%3 == 1:
            yp.append(float(lines[i]))
        elif i%3 == 2:
            zp.append(float(lines[i]))
    
    
    df_traj = pd.DataFrame(columns=["id_ax", "x", "y", "z", "r"])
    d = {'x':xp[::5000], 'y': yp[::5000], 'z': zp[::5000], 'r':0.1, 'traj': 1}
    df_traj = pd.concat([df_traj, pd.DataFrame(d)])

    ## Substract
    chosen_ax = 3
    df_subs = pd.DataFrame(columns=["id_ax", "x", "y", "z", "r"])
    
    for chosen_ax in range(3):
        with open(sub_file) as f:
            lines = f.readlines()
            lines = lines[2:]
        
            for i in range(len(lines)):
                coords = lines[i].split(' ')
                coords_n = [float(coords[i]) for i in [0, 4, 5, 6, 7]]
        
                if  int(coords[0]) == chosen_ax:
                    d = {"id_ax": coords_n[0], 
                            "x": coords_n[1]/1000, 
                            "y": coords_n[2]/1000, 
                            "z": coords_n[3]/1000, 
                            "r": coords_n[4],
                            "traj": 0}
        
                    df_avg_data = pd.DataFrame(d, index=[chosen_ax])
                    df_subs = pd.concat([df_subs, df_avg_data])
    
    
    df_all = pd.concat([df_traj, df_subs])
    print(df_all)
    
    ## Plot
    colors = [f"rgb({int(255*t)}, 0, {int(255*(1-t))})" for t in df_all["traj"]]
    fig = go.Figure()    
    fig.add_trace(go.Scatter3d(
                                x=df_all["x"],
                                y=df_all["y"],
                                z=df_all["z"],
                                type="scatter3d",
                                mode="markers",
                                marker=dict(
                                    sizemode="diameter",
                                    size=df_all["r"]*10,
                                    color=colors,
                                    line=dict(
                                        color="rgba(0, 0, 0, 0)",
                                        width=0
                                    )
                                )
                            )
                    )
    
    plot(fig, auto_open=True,filename=os.path.join(output_folder,'Traj_example'))

