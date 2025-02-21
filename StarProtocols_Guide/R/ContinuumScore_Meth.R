### Continuum score generation - Methylation Array - Easy Reproducible Example ###
#
# Install required packages and their dependencies:
install.packages("mlbench", dependencies = TRUE)
install.packages("caret", dependencies = TRUE)
install.packages("randomForest", dependencies = TRUE)
#
# Load required packages (For installation, please see main text)
library(mlbench)
library(caret)
library(randomForest)
#
#
# NOTE: You MUST update '/your/directory/' to the path where you have cloned this repository.
#
# Load in the prediction object
load(file = "/your/directory/Group3-4App-main/StarProtocols_Guide/data/g3.g4.cont.rfe.Rdata")
#
# Load in example methylation dataset
mvals.mat <- read.delim("/your/directory/Group3-4App-main/StarProtocols_Guide/data/mvals.mat.txt")
#
#
# Subset matrix to probes used as predictors in model
mvals.mat <- as.matrix(mvals.mat[predictors(g3.g4.cont.rfe),])
#
# Apply test set to model and get predicted continuum scores using predict()
Pred.cont.rand.for <- as.data.frame(predict(g3.g4.cont.rfe, t(mvals.mat)))
#
# Format output for export
colnames(Pred.cont.rand.for) <- "Continuum Score"
#
# Preview of output
head(Pred.cont.rand.for)
#
# Export as .csv 
write.csv(Pred.cont.rand.for, file = "/your/directory/my_continuum_scores_Methylation.csv",
          row.names = TRUE)
