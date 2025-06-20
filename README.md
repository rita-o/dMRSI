# dMRI-MRS Processing Toolbox

This package includes:
   - Codes to preprocess and analyse dMRI data
   - Codes to preprocess and analyse dMRS(I) data

> 🛠️ **This toolbox is a continuous work in progress.**  
> If you encounter issues or have questions, **please let us know** so we can improve it.
Contacts: ana.veiga-de-oliveira@chuv.ch

 <br> 

## USAGE

Download this toolbox to your computer and ensure all dependencies are installed as described above.

### dMRI Processing: 

There are three main scripts for diffusion MRI (dMRI) processing (in `processing_dwi`). Each script corresponds to a specific type of data:

- **Main_rat.py** — For rodent dMRI data  
- **Main_human.py** — For human dMRI data  
- **Main_organoid.py** — For organoid dMRI data

Each script contains the complete pipeline for dMRI preprocessing and analysis, organized into sequential steps.

Instructions:
1. Open the relevant script based on your dataset.
2. Review and customize the parameters at the beginning of the script to suit your experimental setup and processing needs.
3. Run the entire script (not advisable) for full pipeline execution,  
   or
   Run individual steps (e.g., `StepX`) if you want more control or are rerunning specific stages.

This pipeline is designed to process **multi-shell** diffusion data with **multiple diffusion times**, supporting both **Linear Tensor Encoding (LTE)** and **Spherical Tensor Encoding (STE)** for processing and analysis, along with an **anatomical** reference image (T1- or T2-weighted). Several images to control for quality are generated along the processing and saved under (`QA_X`).

### dMRS Processing: 

(coming soon (in `processing_dmrs`))

 <br>
 
## NOTES

1. Each dataset must include an Excel file with **cohort metadata** (e.g., subject ID, group, scan date).  
Some columns will be automatically filled during **Step1** of the pipeline.  
An example file is provided in the `common` folder.


2. The `common` folder contains essential resources shared across processing pipelines:

- Configuration files for **dMRS fitting**
- Basis dataset for **dMRS fitting**
- A toolbox to convert MRS data from **Bruker format to NIfTI**: `nifti_mrs_from_raw`
- Atlas files for anatomical segmentation and registration:
  - A **template file** (T1- or T2-weighted) → the filename must include `'template'`
  - A corresponding **atlas file** where each region has a number associated to it → the filename must include `'atlas'`
  - A **label file** containing the mapping between the numbers of the atlas and the region names → the filename must include `'labels'`

> Atlases used:  
> - **WHS_v4** for rodents  
> - **DKT** and **Juelich** for humans  
> Contact us if you want to have access to these atlas folders (too large for GitHub)


3. There are some commmon python scripts that support multiple components of the pipeline:

- **`auxiliar_modelling.py`**  
  Auxiliary functions for modeling routines (e.g., fitting options like `nexi`).

- **`bids_structure.py`**  
  Functions for organizing data in BIDS format.

- **`custom_functions.py`**  
  General-purpose helper functions used throughout the pipeline.

 <br> 
  
## PREREQUISITES (not provided here)

This package runs in Python and uses Conda to deal with multiple environments.

The following tools must be installed:

- [**RATS_MM**](https://iibi.uiowa.edu/rats-rodent-brain-mri) Add to your system's `PATH` after installation. Need only for brain extraction of rodent data. 

- [**ANTs (Advanced Normalization Tools)**](https://github.com/ANTsX/ANTs) Add to your system's `PATH` after installation. Needed throughout preprocessing.

- [**FSL (FMRIB Software Library)**](https://fsl.fmrib.ox.ac.uk/fsl/docs/#/) Add to your system's `PATH` after installation. Needed throughout preprocessing.

- [**MRtrix3**](https://www.mrtrix.org/) Needed throughout preprocessing.

- [**DESIGNER**](https://nyu-diffusionmri.github.io/DESIGNER-v2/) Needed for denoising and DTI/DKI fitting.

 <br> 

These tools should each be installed in their own dedicated Conda environments:

- [**Dicomifier**](https://github.com/lamyj/dicomifier) Environment name: `Dicomifier`; Purpose: Conversion of Bruker data to NIfTI. Only needed for dMRI data acquired with Bruker scanner - on rodents or organoids for example.
 
- [**dcm2niix**](https://github.com/rordenlab/dcm2niix) Environment name: `niix2bids`; Purpose: Conversion of Siemens data to NIfTI. Only needed for dMRI data acquired with human Siemens scanner.

- [**SwissKnife**](https://github.com/QuentinUhl/graymatter_swissknife) Environment name: `SwissKnife`; Purpose: Apply microstructural models to the dMRI data. Needed to apply NEXI, SANDI or SMEX.
 
- [**FSL MRS**](https://open.win.ox.ac.uk/pages/fsl/fsl_mrs/) Environment name: `fsl_mrs`; Purpose: analyze dMRS data

- **Matlab** (with [**MPPCA**](https://github.com/Neurophysics-CFIN/MP-PCA-Denoising) and [**tMPPCA**](https://github.com/Neurophysics-CFIN/Tensor-MP-PCA) toolboxes) is required if denosing with the matlab options; if you don't have matlab you can denoise with mrtrix or designer options (check `Step3.py`).

 <br> 
 
## AUTHORS

Authors: Rita Oliveira & Malte Brammerloh

Supervisor: Ileana Jelescu

Microstructure Mapping Lab (mic-map),
Department of Radiology,
Lausanne University Hospital and University of Lausanne (CHUV),
Rue Pépinet 3, 1003 Lausanne, Switzerland

Email: ana.veiga-de-oliveira@chuv.ch, malte.brammerloh@chuv.ch 

Last updated: June 2025

