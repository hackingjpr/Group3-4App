# MRT Risk Calculator App

# Overview
This script "app.R" encodes a shiny app that upon uploading idat files, will give a relative risk score for patient samples. 
This risk score can be measured using ATRT, ECRT, or MRT (ATRT & ECRT) metagenes.  

# Background
Malignant rhabdoid tumours (MRT) are aggressive paediatric tumours with a poor prognosis despite aggressive multi-modal therapies. These tumours are named according to tumour location: Atypical Teratoid Rhabdoid Tumour (ATRT) if located within the central nervous system (CNS) and Extra-Cranial Rhabdoid Tumours (ECRT) when located outside the CNS.
This app provides a relative risk score to identify different levels of biological risk for patient samples generated using metagene signatures for either ATRT or ECRT specifically. There is an additional signature also provided for use across all rhabdoid tumours (MRT). 
The higher the risk score the higher the predicted biological risk for the patient. These scores require prospective validation in upcoming clinical trials and are currently solely for research purposes only.
Further information regarding the details of these signatures and the setting where they could have possible utility can be found in the following publication: paper link.


# System Requirements
### Hardware Requirements
Functions provided here are compatible with any standard computer with enough RAM to support the in-memory operations.
                
### Software Requirements
Functions provided import and depend on a number of R packages. Functionality has been tested on *R 4.0.4* with *Ubuntu 20.04.2 LTS*
                  
                  
                  
# Disclaimer : This app is designed exclusively for research purposes and is strictly not for diagnostic use.
                  