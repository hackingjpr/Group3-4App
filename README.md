# Deriving a continuum score for group 3 and 4 medulloblastoma tumour samples analysed via RNA-sequencing or DNA methylation microarray


- [Overview](#overview)
- [Tutorial](#tutorial)
  - [Step One - Select Expression or Methylation and Upload Data](#step-one---select-expression-or-methylation-and-upload-data)
  - [Step Two - Generate Group 3-4 Scores](#step-two---generate-group-3-4-scores)
  - [Step Three - Results](#step-three---results)
  - [Step Four - Export Data](#step-four---export-data)
  - [Step Five - Reset](#step-five---reset)
- [Run script without using Shiny App (with example data)](#run-script-without-using-shiny-app-with-example-data)
- [System Requirements](#system-requirements)
- [Troubleshooting](#troubleshooting)

# Overview
This script "app.R" encodes a shiny app that, upon uploading idat files, will give a group 3/4 score for patient samples.   
The repository can be cloned onto your RStudio or the base code can be run independently, this is shown underneath.
This app is derived from the paper of Williamson et al.: *Medulloblastoma group 3 and 4 tumors comprise a clinically and biologically significant expression continuum reflecting human cerebellar development*, which can be found [here](https://doi.org/10.1016/j.celrep.2022.111162).
And is covered in the STAR Protocols paper: *Deriving a continuum score for group 3 and 4 medulloblastoma tumour samples analysed via RNA-sequencing or DNA methylation microarray*. For those not wanting to run their samples on the app the script: [Grp34NoShinyScript.R](https://github.com/hackingjpr/Group3-4App/blob/main/Grp34NoShinyScript.R) can be used, this will provide graphs and outputs but is less interactive. The scripts can also be run in the terminal and not RStudio but the tutorial will be for RStudio. 

# Installation Instructions
To install and begin using the app on your own machine you will first need to install RStudio and then clone the repository. Cloning can either be done in the command line or directly in RStudio. The tutorial for how to do so using command line is [here](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository). 
And for using RStudio directly a tutorial is [here](https://resources.github.com/github-and-rstudio/).
Once you have cloned the repository you should see the contained files in the box in the bottom right quadrant. From these files double click on "app.R", this will open the script behind the shiny app. RStudio should recognise that this is a shiny app and display a "Run App" button in the top right corner of the upper left quadrant. If there is no button simply click anywhere in the script and press Ctrl+A to select all text and then Ctrl+Return to run what is selected. The first time you run the script it may take a while as it will install all of the needed packages, it may also ask you if you want to update your other packages, this choice is up to you. If it still doesn't work then check to see if there is any help in the [troubleshooting](#troubleshooting) section of this README. If running the app on Apple silicon you may need to install lgfortran library, the information about that can be found here: https://mac.r-project.org/tools/.

# Tutorial

## Step One - Select Expression or Methylation and Upload Data

<details>
  <summary>Select Expression or Methylation and Upload Data</summary>

Depending on whether you are uploading Expression or Methylation data select the appropriate option.

Upload your idat files including both red and green files for Methylation, or RDS/TXT/CSV files for Expression.

![upload.png](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/app_upload.png)

Increasing the number of samples will of course increase the length of time for the upcoming processes so we recommend ~10 sample batches. This will make looking through the results easier and will speed up the process.

If uploading Expression data you will be asked to give up to two further inputs:  
1. Selecting whether to scale your results against the data frame of Williamson et al. or against your own uploaded data.
2.  If you selected scaling against your own uploaded data you will be asked if you want to filter out any outliers. This is done via a sliding scale from one to four, for removing samples more than one to four standard deviations from the mean. 

</details>

## Step Two - Generate Group 3-4 Scores
<details>
  <summary>Generate Group 3-4 Scores</summary>
Click the "Generate Group 3/4 Score" button. This will start the process of generating Group 3/4 Continuum Scores and a loading bar should begin filling underneath the "Reset" button.

![Generate Scores](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/app_generate.png)
</details>

## Step Three - Results
<details>
  <summary>Results</summary>
Once the calculation has been completed you should be brought to the Results tab. This tab will show a data table at the top which displays your sample names on the left and their Group 3/4 Scores on the right.

![Results Table](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/app_score_table.png)

It will also show a number of graphs, an example is below.

### Group 3/4 Plot
<details>
  <summary>Plot</summary>


Places your sample data on a cumulative frequency plot based on data from the Cell Reports paper of Williamson et al. It tells you whether the patient is Group 3 or Group 4 and allows you to see where the patient ranks against this large dataset. 

 ![E1](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/E1.png)

</details>

### Survival Plot: No Risk Factors Considered

<details>
  <summary>Plot</summary>

Shows patients expected five year survival based on only their group 3/4 score and no other risk factors.

 ![E2](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/E2.png)

</details>

### Survival Plot: Age Considered

<details>
  <summary>Plot</summary>

Shows patients expected five year survival based on their group 3/4 score but also taking into account their age.

 ![E3](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/E3.png)

</details>

There will also be a box entitled "Selected Sample Information", this will inform you of the currently selected sample (which will also be highlighted on the graphs), and give some information about the sample and the expected survival.

![Selections](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/app_sample_selected.png)

</details>

## Step Four - Export Data

<details>
  <summary>Export Data</summary>

Once you have your results you can download your data as a CSV file (data table results), or as a PDF (data table and graphs).

![Download](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/app-download.png)

</details>

## Step Five - Reset

<details>
  <summary>Reset</summary>

Once you have looked over or downloaded your data you can reset the app to upload other samples. 

To do this click the "Reset" button in the bar on the left of the app.

![Reset](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/app_reset.png)

</details>


# Run script without using Shiny App (with example data)

 <details>
  <summary>Script</summary>
  
```
#####################################
############ METHYLATION ############
#####################################

# Install/Load required packages and their dependencies 

install.packages('mlbench', dependencies = TRUE) 
install.packages('caret', dependencies = TRUE) 
install.packages('randomForest', dependencies = TRUE) 

# For specific package versions, see Key Resource Table section. 

library(mlbench) 
library(caret) 
library(randomForest) 
# This loads each package into your working environment 

# CRITICAL: You MUST update ‘/your/directory/’ to the location which you cloned the GitHub repository in step 1 of Continuum score assignment (RNA-Sequencing). 

# Load in the prediction object 

load(file = "/your/directory/Group3-4App/StarProtocols_Guide/data/g3.g4.cont.rfe.Rdata") 
# This loads in the precalculated random forest model 

# Load in example methylation dataset. 

mvals.mat <- read.delim("/your/directory/Group3-4App/StarProtocols_Guide/data/mvals.mat.txt") 

# CRITICAL: The random forest model in this protocol requires that test data be provided as a matrix of M-values (logit-transformed beta values),
# where columns correspond to sample ID and rows correspond to probes. If your data is a matrix of beta values (object below named as “your.betas”), 
# you can easily convert these to M-values using the following: 

#mvals.mat <- log2(your.betas/(1-your.betas)) 
# logit-transformation 

# Subset M-Value matrix to probes used as predictors in model 

mvals.mat <- as.matrix(mvals.mat[predictors(g3.g4.cont.rfe),]) 
# Removes probes that are not used for prediction 

# Apply test set to model and get predicted continuum scores using predict() 

pred.cont.rand.for <- as.data.frame(predict(g3.g4.cont.rfe, t(mvals.mat))) 

write.csv(pred.cont.rand.for, file = '/your/directory/my_continuum_scores_Methylation.csv', row.names = TRUE) 
# Export as .csv 

# Expected outcome: A data.frame object where rows correspond to sample ID and column corresponds to each sample's respective continuum score value.  
  
source("./AppSourceFunctions1.13.R")

# Example graphs displaying the first 10 samples, graphs can get cluttered if too many samples are displayed
# currently displays samples 1-10 and highlights sample 1.
samples.to.display <- c(1:10) #change this to display different samples, currently 1-10.

### Generate Group3/4 score graph selecting the first sample to highlight
generate_figure_highlight_g3g4Expression(pred.cont.rand.for[samples.to.display,1]
                                         , 1)

### Generate Survival Plot selecting the first sample to highlight
survivalcurveplot(pred.cont.rand.for[samples.to.display,1]
                  ,1)

### Generate Survival Plot selecting the first sample to highlight
SurvivalAgePlot(pred.cont.rand.for[samples.to.display,1],
                1)  
 
####################################
############ EXPRESSION ############
####################################

# Install/load required R packages and their dependencies.  

install.packages("NMF", dependencies = TRUE) 
install.packages("MASS", dependencies = TRUE) 
BiocManager::install("biomaRt") 

# For specific package versions, see Key Resource Table section. When confronted with yes/no questions, answer yes to install dependency packages. 

library(NMF) 
library(MASS) 
# This loads the packages required into your working environment. 

# Load required data objects. 

# CRITICAL: You MUST update ‘/your/directory/’ to the location which you cloned the GitHub repository in step 1. 

nmf.res <- readRDS(file = "/your/directory/Group3-4App/StarProtocols_Guide/data/nmf.res.rds") 
# This loads in the precalculated NMF model. 

# Load the required custom functions. 

source(file = "/your/directory/Group3-4App/StarProtocols_Guide/R/Project_NMF.R") 
# Wrapper function used to project NMF model onto unseen group3/group4 sample data. A function breakdown is provided below (see figure 1.). 


# Load sample data as a matrix object. 
tpms.mat <- read.delim("/your/directory/Group3-4App/StarProtocols_Guide/data/tpms.mat.txt") 

# CRITICAL: If you wish to use your own RNA-sequencing sample data, you must ensure that it follows
# the same format as tpms.mat. This object is a matrix, where columns correspond to samples and rows correspond to genes,
# with expression counts presented in the transcripts per million (TPM) format or equivalent. 
# All genes (rows) must use HUGO gene nomenclature i.e., gene symbols.
# If your input dataset is not annotated correctly, please see Problem 1 in the Troubleshooting section.
# Note that a column-rank normalization procedure is employed, this coupled with the NMF projection and other
# normalisation procedures renders the results somewhat resistant to noise and compatible with representations of expression other than TPM.
# We have for example used Rlog, or variance stabilised transforms from DESeq or even other platforms such as Affymetrix microarray or nanostring data with success.
# Note when projecting onto platforms other than bulk RNA-seq appropriate filtering strategies to remove invariant genes/probes may be necessary. 

# Project NMF model onto sequencing data  

tpms.H <- project.NMF(input.array = as.matrix(tpms.mat), nmf.result = nmf.res) 
# Apply project.NMF function to input dataset.
  
# Extract Group 3 and Group 4 metagenes from data and transpose matrix. 

g3g4.tpms <- t(tpms.H[c(3,1),])  
# Rows 3 and 1 in tpms.H correspond to the metagenes for Groups 4 and 3 respectively. 

# Apply logistic transformation to metagenes. 

logistic.g3g4.tpms <- apply(g3g4.tpms,2,function(x){(1 / (1 + exp(-x)))}) 
# Apply a logistic transformation  

logistic.g3g4.tpms.score <- apply(logistic.g3g4.tpms,1,function(x){x[2]/(x[1]+x[2])}) 
# Calculate a ratio between logistically transformed Group3 and Group4 metagene  

# Scale values between 0 and 1. 

scaling.function <- function(x){(x-min(x)) / (max(x)-min(x))} 
# Create a function to scale values between 0 and 1 

logistic.g3g4.tpms.continuum.score <- scaling.function(logistic.g3g4.tpms.score) 
# Apply the function to the unscaled g3g4 scores  

# CRITICAL: If you are using a small dataset or one that does not represent the full spectrum of Group3/Group4 
# medulloblastomas you may want to omit this step and present unscaled G3/G4 ratios in which case the following command should be used. 

# Alternatively, you may wish to append to the precalculated G3/G4 ratios from Williamson et al 
# and then scale together with your new samples in which case the following alternative command should be used:  

scaling.function1 <- function(x){(x - 0.3953062) / (0.5964371 - 0.3953062)} 
# Create a function to scale values between 0 and 1 using Williamson et al. data) 

logistic.g3g4.tpms.continuum.score <- scaling.function1(logistic.g3g4.tpms.score) 
# Apply scaling 

# Note that theoretically this could lead to some sample returning values under 0 or over 1.
# The user would need to take a considered view on such samples. They could simply be producing
# values close to 1 or 0 but otherwise consistent with samples at the extreme limits of the G3/G4 continuum,
# in which case manually assigning them the maximum 1 or minimum 0 value may be a valid approach.
# Should they massively exceed previous limits they may simply be outliers or technical artefacts 
# that would be best noted but excluded from further analysis. 

# Present output as data.frame for export. 

logistic.g3g4.tpms.continuum.score <- as.data.frame(logistic.g3g4.tpms.continuum.score) 

colnames(logistic.g3g4.tpms.continuum.score) <- 'Continuum Score' 
# Renaming for easier interpretation 

write.csv(logistic.g3g4.tpms.continuum.score, file = '/your/directory/my_continuum_scores.csv ', row.names = TRUE) 
#Export as .csv table 

# Expected outcome: A data.frame object where rows correspond to sample ID and column corresponds to each samples respective continuum score value.
  
source("/your/directory/Group3-4App/AppSourceFunctions1.13.R")

# Example graphs displaying the first 10 samples, graphs can get cluttered if too many samples are displayed
# currently displays samples 1-10 and highlights sample 1.
samples.to.display <- c(1:10) #change this to display different samples, currently 1-10.

### Generate Group3/4 score graph selecting the first sample to highlight
generate_figure_highlight_g3g4Expression(logistic.g3g4.tpms.continuum.score[samples.to.display,1]
                                         , 1)

### Generate Survival Plot selecting the first sample to highlight
survivalcurveplot(logistic.g3g4.tpms.continuum.score[samples.to.display,1]
                  ,1)

### Generate Survival Plot selecting the first sample to highlight
SurvivalAgePlot(logistic.g3g4.tpms.continuum.score[samples.to.display,1],
                1)  
                
```

 </details>
 
## Example Script Results
<details>
  <summary>Methylation</summary>


 For methylation the expected results are:  
 
 ```
        predict(g3.g4.cont.rfe, t(mvals.mat))
NMB_109                            0.07934780
NMB_110                            0.17779087
NMB_111                            0.31949452
NMB_118                            0.41626223
NMB_119                            0.15780596
NMB_125                            0.05701327
NMB_130                            0.09525883
NMB_132                            0.29072448
NMB_134                            0.20310554
NMB_136                            0.37450195

```
  
 The following graphs will also be created:  
 ![Methylation1](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/Methylation1.png)
 ![Methylation2](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/Methylation2.png)
 ![Methylation3](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/Methylation3.png)
 
 </details>
 
 <details>
  <summary>Expression</summary>
  
  For expression the expected results are:
  
```
> logistic.g3g4.tpms.continuum.score
          Continuum Score
Sample1         0.5991949
Sample2         0.7221333
Sample3         0.8062744
Sample4         0.5188520
Sample5         0.7622158
Sample6         0.7608163
Sample7         0.5246359
Sample8         0.8032780
Sample9         0.7450637
Sample10        0.8553960...

```

 The following graphs will also be created:  
 ![Expression1](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/Expression1.png)
 ![Expression2](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/Expression2.png)
 ![Expression3](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/Expression3.png)

   </details>
 

# System Requirements
## Hardware Requirements
Functions provided here are compatible with any standard computer with enough RAM to support the in-memory operations.

## Software Requirements
Functions provided import and depend on a number of R packages. Functionality has been tested on *R 4.2.2* with *Ubuntu 20.04.2 LTS*

# Troubleshooting
## App won't install
1. Check if you have selected all of the script before running it.
2. Check to see if RStudio is asking you if you want to update your packages, this will appear in the bottom left quadrant and look like "Update (Y/N)". Click in the quadrant and type Y then press Return. If this doesn't sort it then next time press N and press return, RStudio can be temperamental.

## App won't open after running
Pop ups may be blocked, try running the app again and see if something appears saying "pop up blocked", click this and tell it to allow pop ups.

## App closes after uploading files and clicking "Generate Group 3/4 Scores"
If running Methylation idat files then make sure you have uploaded both red and green files for the sample.

## App won't open after closing itself due to an error
If the app closes itself down then you may need to press stop. The button for this is a red hexagon that can be found in the top right corner of the box in the lower left quadrant of the RStudio screen.




# *Disclaimer: This app is designed exclusively for research purposes and strictly not for diagnostic or clinical use.*
