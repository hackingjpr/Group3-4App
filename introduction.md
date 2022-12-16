# Group 3-4 Continuum Score

## Overview

This shiny app, upon uploading idat files, will give a group 3/4 score for patient samples.
This app is derived from Williamson et als. paper: *Medulloblastoma group 3 and 4 tumors comprise a clinically and biologically significant expression continuum reflecting human cerebellar development*, which can be found [here](https://doi.org/10.1016/j.celrep.2022.111162), or navigate to the "Paper" tab above. The app and protocol in general is covered in the STAR Protocols paper entitled *"Deriving a continuum score for group 3 and 4 medulloblastoma tumour samples analysed via RNA-sequencing or DNA methylation microarray"*


  ![Graphical abstract](https://ars.els-cdn.com/content/image/1-s2.0-S2211124722009718-fx1.jpg)

## Background

The script behind this app codes for the specific in-silico steps for deriving a continuum score from group3/group4 medulloblastoma sample data analysed using RNA-sequencing or DNA methylation microarray (either 450k or EPIC platforms). Determination of sample subgroup is therefore essential and can be determined via the upload of DNA methylation microarray sample data (IDAT format) to the Molecular Neuropathology Classifier website (v11.4 and v12.5 subgroup calls are compatible – [link](https://www.molecularneuropathology.org/mnp/)). Prior knowledge of data processing for tumour sample data from RNA-sequencing and DNA methylation Microarray sources is also assumed: please refer to the original publication (Williamson et al., 2022) and Key Resources Table section for detail on platform-specific pre-processing strategies and software used.  

To derive a continuum score for group 3/group 4 tumours analysed via RNA-sequencing, an NMF (Non-negative Matrix Factorization) model trained using 331 medulloblastoma tumours is projected onto incoming sample data using a procedure broadly outline by Tamayo et al (Tamayo et al 2007). This procedure is somewhat resistant to noise, differences in platform or even species. Two of the metagenes which correspond to Group3 and Group4 patients are logistically transformed and a ratio between the two produced which is then scaled between 0 and 1 (G3/G4 score). Whilst knowledge of metagenes and NMF is not expressly required for this protocol, more information about NMF and how it is applied in this study can be found in this study’s publication (Williamson et al., 2022).  

By utilising paired-methylation profiles where RNA-based continuum score was known (n=192/331), a regression random forest model trained using 400 CpGs was created that is capable of reverse engineering a continuum score for group 3/group 4 sample data as a proxy from both the Illumina Methylation450k and MethylationEPIC platform with high accuracy (RMSE = 0.036). More detail on how the random forest model classifier was trained is provided in the original manuscript (Williamson et al., 2022). 

## System Requirements

### Hardware Requirements

Functions provided here are compatible with any standard computer with enough RAM to support the in-memory operations.

### Software Requirements

Functions provided import and depend on a number of R packages. Functionality has been tested on *R 4.0.4* with *Ubuntu 20.04.2 LTS*
