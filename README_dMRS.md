# 🟡 dMRS Processing overview

This pipeline is designed to process **multi–diffusion-time** and **multi–b-value** dMRS data acquired using the **SPECIAL** sequence.

## 🗒️ Description of analysis steps:

- **Step1_preproc**
Processes raw dMRS data using MATLAB scripts developed by the **EPFL group of Cristina Cudalbu** (by: Jessie Mosso, Toi Phan, and Éloïse Mougel), curated by Malte Brammerloh and integrated by Rita Oliveira, to:

  - Convert Bruker raw data to MATLAB (`.mat`) format  
  - Pre-process dMRS data  
  - Quantify metabolites using **LCModel**


- **Step2_fitting**
Fits the quantified dMRS data using diffusion models (developed by Malte Brammerloh), including:

  - `DTI`
  - `Stick`
  - `DKI`
  - `Cylinder`
  - `Cylinder–Sphere`
  - `Stick–Sphere`

This step requires the **SwissKnife** Conda environment.

