#### Continuum Score Generation - RNA-Seq - Easiest Reproducible Example ####
#
#
# Install/load packages
library(NMF)
library(MASS)
#
#
# Load nmf.res
nmf.res <- readRDS(file = "~/Documents/STARMETHODS/nmf.res.rds")

# tpms.mat
tpms.mat <- read.delim("~/Documents/STARMETHODS/tpms.mat.txt")

#Load project.NMF function
source(file = "~/Documents/Project_NMF.R")



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
