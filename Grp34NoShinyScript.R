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

# Load in the prediction object 

load(file = "/your/directory/Group3-4App/StarProtocols_Guide/data/g3.g4.cont.rfe.Rdata") 
# This loads in the precalculated random forest model 

# Load in example methylation dataset. 

mvals.mat <- read.delim("/your/directory/Group3-4App/StarProtocols_Guide/data/mvals.mat.txt") 

#If your data is a matrix of beta values (object below named as “your.betas”), 
# you can convert these to M-values using the following: 

#mvals.mat <- log2(your.betas/(1-your.betas)) 
# logit-transformation 

# Subset M-Value matrix to probes used as predictors in model 

mvals.mat <- as.matrix(mvals.mat[predictors(g3.g4.cont.rfe),]) 
# Removes probes that are not used for prediction 

# Apply test set to model and get predicted continuum scores using predict() 

pred.cont.rand.for <- as.data.frame(predict(g3.g4.cont.rfe, t(mvals.mat))) 

write.csv(pred.cont.rand.for, file = '/your/directory/my_continuum_scores_Methylation.csv', row.names = TRUE) 
# Export as .csv 


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


library(NMF) 
library(MASS) 
# This loads the packages required into your working environment. 

# Load required data objects. 


nmf.res <- readRDS(file = "/your/directory/Group3-4App/StarProtocols_Guide/data/nmf.res.rds") 
# This loads in the precalculated NMF model. 

# Load the required custom functions. 


# Load sample data as a matrix object. 
tpms.mat <- read.delim("/your/directory/Group3-4App/StarProtocols_Guide/data/tpms.mat.txt") 


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


# If you are using a small dataset you may want to omit this step and present unscaled G3/G4 ratios in which case the following command should be used. 
# Scale values between 0 and 1. 
scaling.function <- function(x){(x-min(x)) / (max(x)-min(x))} 
# Create a function to scale values between 0 and 1 

logistic.g3g4.tpms.continuum.score <- scaling.function(logistic.g3g4.tpms.score) 

# Alternatively, you may wish to append to the precalculated G3/G4 ratios from Williamson et al 
# and then scale together with your new samples. If so the following alternative command should be used:  

scaling.function1 <- function(x){(x - 0.3953062) / (0.5964371 - 0.3953062)} 
# Create a function to scale values between 0 and 1 using Williamson et al. data

logistic.g3g4.tpms.continuum.score <- scaling.function1(logistic.g3g4.tpms.score) 
# Apply scaling 

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