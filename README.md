# dMRSI

This package includes:
   - Codes to preprocess and analyse dMRI data
   - Codes to preprocess and analyse dMRS(I) data
  

## USAGE:

Download this toolbox to your computer and ensure all dependencies are installed as described above.

### Cohort Metadata

For each dataset, you must provide an Excel file containing metadata about your cohort. This file serves as a central reference for subject information and processing.

- The file should include basic information about each subject (e.g., ID, group, scan date, etc.).
- Some of the columns will be automatically filled during **Step1** of the processing pipeline.
- An example file is included in the `common` folder

### dMRI Processing Scripts 
(in `processing_dwi`)

There are three main scripts for diffusion MRI (dMRI) processing. Each script corresponds to a specific type of data:

- **Main_rat.py** — For rodent dMRI data  
- **Main_human.py** — For human dMRI data  
- **Main_organoid.py** — For organoid dMRI data

Each script contains the complete pipeline for dMRI preprocessing and analysis, organized into sequential steps.

Instructions:
1. Open the relevant script based on your dataset.
2. Review and customize the parameters at the beginning of the script to suit your experimental setup and processing needs.
3. Run the entire script** for full pipeline execution,  
   or
   Run individual steps (e.g., `StepX`) if you want more control or are rerunning specific stages.

### dMRS Processing Scripts 
(in `processing_dmrs`)
(coming soon)

## PREREQUISITES (not provided here):

This package runs in Python with Conda.

The following tools must be installed and available in your system's `PATH`:

- [**RATS_MM**](https://iibi.uiowa.edu/rats-rodent-brain-mri)  (used for brain extraction of rat data)

- [**ANTs (Advanced Normalization Tools)**](https://github.com/ANTsX/ANTs)

- [**FSL (FMRIB Software Library)**](https://fsl.fmrib.ox.ac.uk/fsl/docs/#/)

⚠️ Make sure all the above programs are added to your system's `PATH` after installation.

- [**MRtrix3**](https://www.mrtrix.org/)

- [**DESIGNER**](https://nyu-diffusionmri.github.io/DESIGNER-v2/)
 

These tools should each be installed in their own dedicated Conda environments:

- [**Dicomifier**](https://github.com/lamyj/dicomifier)  
  - Environment name: `Dicomifier`  
  - *Purpose:* Conversion of Bruker data to NIfTI (not needed for dMRS)
 
- [**dcm2niix**](https://github.com/rordenlab/dcm2niix)  
  - Environment name: `niix2bids`  
  - *Purpose:* Conversion of Siemens data to NIfTI (not needed for dMRS)

- [**SwissKnife**](https://github.com/QuentinUhl/graymatter_swissknife)  
  - Environment name: `SwissKnife`
  - *Purpose:* Apply microstructural models to the dMRI data
 
- [**FSL MRS**](https://open.win.ox.ac.uk/pages/fsl/fsl_mrs/)  
  - Environment name: `fsl_mrs`
  - *Purpose:* analyze dMRS data


## AUTHORS:

Authors: Rita Oliveira & Malte Brammerloh

Supervisor: Ileana Jelescu

Microstructure Mapping Lab (mic-map),
Department of Radiology,
Lausanne University Hospital and University of Lausanne (CHUV),
Rue Pépinet 3, 1003 Lausanne, Switzerland

Email: ana.veiga-de-oliveira@chuv.ch, malte.brammerloh@chuv.ch 

Last updated: June 2025
