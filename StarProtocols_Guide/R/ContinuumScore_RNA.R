#### Continuum Score Generation - RNA-Seq - Easiest Reproducible Example ####
#
# ToDo : Update filepaths when testing gitclone
#
# Load packages
library(NMF)
library(MASS)
#
#
# Note: You MUST update '/your/directory/' to the directory that you clone this repository in to.
#
# Load nmf.res
nmf.res <- readRDS(file = "/your/directory/Group3-4App-main/StarProtocols_Guide/data/nmf.res.rds")

# tpms.mat
tpms.mat <- read.delim("/your/directory/Group3-4App-main/StarProtocols_Guide/data/tpms.mat.txt")

#Load project.NMF function
source(file = "/your/directory/Group3-4App-main/StarProtocols_Guide/R/Project_NMF.R")



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

head(logistic.g3g4.tpms.continuum.score)

