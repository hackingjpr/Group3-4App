### Set working directory to wherever "source_functions.R" is
setwd("/your/directory/Group3-4App/")
source("./AppSourceFunctions1.13.R")


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