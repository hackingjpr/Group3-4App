#### Continuum Score Generation - RNA-Seq ####
#
# NOTE: Please change "/your/directory/" to the path where you have cloned this git repository.
#
#
# Install required packages and their dependencies:
install.packages("NMF", dependencies = TRUE)
install.packages("MASS", dependencies = TRUE)
BiocManager::install("biomaRt")
#
# Load packages
library(NMF)
library(biomaRt)
library(MASS)
#
# Load project.NMF function:
source(file = "/your/directory/Group3-4App-main/StarProtocols_Guide/R/Project_NMF.R")
#
#
# Load nmf.res
nmf.res <- readRDS(file = "/your/directory/Group3-4App-main/StarProtocols_Guide/data/nmf.res.rds")

# Load match selected tpms.mat
tpms.mat <- read.delim("/your/directory/Group3-4App-main/StarProtocols_Guide/data/tpms.mat.txt")

# Project NMF model onto tpms.mat
tpms.H <- project.NMF(input.array = as.matrix(tpms.mat),
                      nmf.result = nmf.res)

# Isolate Group 3 and Group 4 metagenes values
g3g4.tpms <- t(tpms.H[c(3,1),]) 


#Logistic transformation:               
logistic.g3g4.tpms <- apply(g3g4.tpms,2,function(x){(1 / (1 + exp(-x)))}) 
logistic.g3g4.tpms.score <- apply(logistic.g3g4.tpms,1,function(x){x[2]/(x[1]+x[2])})


#Scale between 0-1 and present output as a data.frame for export:
logistic.g3g4.tpms.continuum.score <- as.data.frame((logistic.g3g4.tpms.score-min(logistic.g3g4.tpms.score)) / (max(logistic.g3g4.tpms.score) - min(logistic.g3g4.tpms.score)))
colnames(logistic.g3g4.tpms.continuum.score) <- "Continuum Score"

# Preview output:
head(logistic.g3g4.tpms.continuum.score)

#Export as .csv
write.csv(logistic.g3g4.tpms.continuum.score, file = "/your/directory/my_continuum_scores.csv",
          row.names = TRUE)

