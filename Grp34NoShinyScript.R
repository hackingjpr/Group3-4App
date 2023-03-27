### Set working directory to wherever "source_functions.R" is
setwd("~/Group3-4App")
source("./AppSourceFunctions1.12.R")

#####################################
############ METHYLATION ############
#####################################



### load in the prediction object
load(file = "./AppExtraFiles/Inputs/g3.g4.cont.rfe.Rdata")

### Choose folder containing idats to be processed (make sure you set this to yours)
idats <- "~/your/idat/file/location"


### Get Basenames
temp.base <- get_basenames(idats)

### Process Idats
temp.processed <- process_idats(temp.base)


# Obtain MValues
beta2m(temp.processed$betas) -> M.values

if(ncol(M.values)==1){
  ### for single sample
  t(data.frame(t(M.values)[,predictors(g3.g4.cont.rfe)])) -> input.df
  colnames(M.values) -> rownames(input.df)
}else{
  t(M.values)[,predictors(g3.g4.cont.rfe)] -> input.df
}

### Round results to 3 figures
metagene <- round(predict(g3.g4.cont.rfe, input.df), digits = 3)
metagene.df <- data.frame('Group.3.4.Score' = metagene)


### This is your Group3/4 Scores
metagene.df

### Select Risk values column
# figure.input <- test.res$Risk_Value
figure.input <- metagene.df$Group.3.4.Score
names(figure.input) <- rownames(metagene.df)
print(figure.input)


### Name the rows
names(figure.input) <- rownames(metagene.df)
print(figure.input)

### Generate Group3/4 score graph selecting the first sample to highlight
generate_figure_highlight_g3g4(figure.input,
                              1)
### Generate Survival Plot selecting the first sample to highlight
survivalcurveplot(
  figure.input
  ,1)

### Generate Age Survival plot selecting the first sample to highlight
SurvivalAgePlot(figure.input,
                1)


####################################
############ EXPRESSION ############
####################################

## Load in your samples (currently is an example file available from the GitHub repository, if you cloned the GitHub you will already have this)
in.files <- "./AppExtraFiles/Inputs/subsetTpms.mat10.rds"

input.file <- in.files
if (file_ext(input.file) == "rds") {
  in.files <- readRDS(file = input.file)
} else if (file_ext(input.file) == "csv") {
  in.files <- read.csv(file = input.file, row.names = 1)
} else if (file_ext(input.file) == "txt") {
  in.files <- read.delim(file = input.file)
} else {
  message("file not right format!")
}

nmb.mat <- nmb.mat.prepped

# ## interset common genes / probes
tpms.mat <- match.select(nmb.mat, in.files)

## project using pseudo-inverse & post-projection normalise
# project back onto the same dataset
rnaseq.H <- project.NMF(input.array = nmb.mat,
                        nmf.result = nmf.res)

tpms.matrix <- as.matrix(tpms.mat)
if (ncol(tpms.matrix) == 1) {
  colnames(in.files) -> colnames(tpms.matrix)
}

# project onto fresh dataset
tpms.H <- project.NMF(input.array = tpms.matrix,
                      nmf.result = nmf.res)

### define new g3g4 score for projection back onto the original data
# t(rnaseq.H[c(3, 1), ]) -> g3g4.rnaseq

# apply(g3g4.rnaseq, 2, function(x) {
#   (1 / (1 + exp(-x)))
# }) -> logistic.g3g4.rnaseq
# 
# apply(logistic.g3g4.rnaseq, 1, function(x) {
#   x[2] / (x[1] + x[2])
# }) -> logistic.g3g4.rnaseq.score

## Scale to Williamson et al. dataset
# scaling.function3(logistic.g3g4.rnaseq.score) -> logistic.g3g4.rnaseq.score
ScalingChoice <- "ours"

## Scale to uploaded dataset 
#scaling.function(logistic.g3g4.rnaseq.score) -> logistic.g3g4.rnaseq.score
#ScalingChoice <- "yours"

## Outlier removal
outlier <- 0


t(tpms.H[c(3,1),]) -> g3g4.tpms

apply(g3g4.tpms, 2, function(x) {
  (1 / (1 + exp(-x)))
}) -> logistic.g3g4.tpms

if(is.null(dim(logistic.g3g4.tpms))){
  logistic.g3g4.tpms[2] / (logistic.g3g4.tpms[1] + logistic.g3g4.tpms[2]) -> logistic.g3g4.tpms.score
  apply(logistic.g3g4.tpms, 1, function(x) {
    x[2] / (x[1] + x[2])
  }) ->  logistic.g3g4.tpms.score
  message("is.null(dim(logistic.g3g4.tpms))")
  scaling.function3(logistic.g3g4.tpms.score) -> logistic.g3g4.tpms.score
}else{
  apply(logistic.g3g4.tpms, 1, function(x) {
    x[2] / (x[1] + x[2])
  }) -> logistic.g3g4.tpms.score
  
  mean(logistic.g3g4.tpms.score) -> mean.logistic.g3g4.tpms.score
  message(mean.logistic.g3g4.tpms.score)
  sd(logistic.g3g4.tpms.score) -> sd.logistic.g3g4.tpms.score
  message(sd.logistic.g3g4.tpms.score)
  
  if (outlier == 0){
    upper.limit <- 1
    lower.limit <- 0
  }
  else{
    upper.limit <-
      ((outlier) * sd.logistic.g3g4.tpms.score) +  mean.logistic.g3g4.tpms.score
    lower.limit <-
      mean.logistic.g3g4.tpms.score - ((outlier) * sd.logistic.g3g4.tpms.score)
  }
  
  outlier.idx <-
    which(logistic.g3g4.tpms.score > upper.limit |
            logistic.g3g4.tpms.score < lower.limit)
  
  
  if (length(outlier.idx) != 0 & ScalingChoice == "yours") {
    apply(logistic.g3g4.tpms, 1, function(x) {
      x[2] / (x[1] + x[2])
    }) ->  logistic.g3g4.tpms.score
    
    removed <- logistic.g3g4.tpms.score[-outlier.idx]
    
    # scaling.function(logistic.g3g4.tpms.score
    #                  [-outlier.idx]) -> logistic.g3g4.tpms.score
    scaling.function(removed) -> logistic.g3g4.tpms.score
    
  }else if (length(outlier.idx) == 0 & ScalingChoice == "yours") {
    apply(logistic.g3g4.tpms, 1, function(x) {
      x[2] / (x[1] + x[2])
    }) ->  logistic.g3g4.tpms.score
    
    scaling.function(logistic.g3g4.tpms.score) -> logistic.g3g4.tpms.score
    
  } else{
    apply(logistic.g3g4.tpms, 1, function(x) {
      x[2] / (x[1] + x[2])
    }) ->  logistic.g3g4.tpms.score
    scaling.function3(logistic.g3g4.tpms.score) -> logistic.g3g4.tpms.score

  }
}

round(logistic.g3g4.tpms.score, digits = 3) -> logistic.g3g4.tpms.score

## This is your Group 3/4 Values
data.frame('Group.3.4.Score' = logistic.g3g4.tpms.score) -> logistic.g3g4.tpms.score.df

### G3/4 Graph selecting the first sample to highlight
generate_figure_highlight_g3g4Expression(logistic.g3g4.tpms.score
                                         , 1)

## Survival plot selecting the first sample to highlight
survivalcurveplot(
  logistic.g3g4.tpms.score
  ,1)

## Age Survival Plot selecting the first sample to highlight
SurvivalAgePlot(logistic.g3g4.tpms.score,
                1)
