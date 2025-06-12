# dMRSI

This package includes:
   - Codes to preprocess and analyse dMRI data
   - Codes to preprocess and analyse dMRS(I) data
  

## USAGE:



## PREREQUISITES (not provided here):

This package runs in Python with Conda

The following tools must be installed and available in your system's `PATH`:

- [**RATS_MM**](https://iibi.uiowa.edu/rats-rodent-brain-mri)  

- **ANTs (Advanced Normalization Tools)**  
  - [GitHub Repository](https://github.com/ANTsX/ANTs)

- **FSL (FMRIB Software Library)**  
  - [Documentation](https://fsl.fmrib.ox.ac.uk/fsl/docs/#/)

> ⚠️ Make sure all the above programs are added to your system's `PATH` after installation.

- **MRtrix3**  
  - [Official Website](https://www.mrtrix.org/)

- **DESIGNER**  
  - [DESIGNER v2](https://nyu-diffusionmri.github.io/DESIGNER-v2/)
 

## Programs Installed in Conda Environments

These tools should each be installed in their own dedicated Conda environments:

- **Dicomifier**  
  - [GitHub Repository](https://github.com/lamyj/dicomifier)  
  - Environment name: `Dicomifier`  
  - **Purpose:** Conversion of Bruker data to NIfTI

- **SwissKnife**  
  - [GitHub Repository](https://github.com/QuentinUhl/graymatter_swissknife)  
  - Environment name: `SwissKnife`

- **dcm2niix**  
  - [GitHub Repository](https://github.com/rordenlab/dcm2niix)  
  - Environment name: `niix2bids`  
  - **Purpose:** Conversion of Siemens data to NIfTI


## AUTHORS:

Authors: Rita Oliveira & Malte Brammerloh
Supervisor: Ileana Jelescu

Microstructure Mapping Lab (mic-map),
Department of Radiology,
Lausanne University Hospital and University of Lausanne (CHUV),
Rue Pépinet 3, 1003 Lausanne, Switzerland

Email: ana.veiga-de-oliveira@chuv.ch, malte.brammerloh@chuv.ch 

Last updated: June 2025
