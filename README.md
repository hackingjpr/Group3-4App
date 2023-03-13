# Deriving a continuum score for group 3 and 4 medulloblastoma tumour samples analysed via RNA-sequencing or DNA methylation microarray


- [Overview](#overview)
- [Tutorial](#tutorial)
  - [Step One - Select Expression or Methylation and Upload Data](#step-one---select-expression-or-methylation-and-upload-data)
  - [Step Two - Generate Group 3-4 Scores](#step-two---generate-group-3-4-scores)
  - [Step Three - Results](#step-three---results)
  - [Step Four - Export Data](#step-four---export-data)
  - [Step Five - Reset](#step-five---reset)
- [Run script without using Shiny App](#run-script-without-using-shiny-app)
- [Example script to run](#example-script-to-run)
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
Once you have cloned the repository you should see the contained files in the box in the bottom right quadrant. From these files double click on "app.R", this will open the script behind the shiny app. RStudio should recognise that this is a shiny app and display a "Run App" button in the top right corner of the upper left quadrant. If there is no button simply click anywhere in the script and press Ctrl+A to select all text and then Ctrl+Return to run what is selected. The first time you run the script it may take a while as it will install all of the needed packages, it may also ask you if you want to update your other packages, this choice is up to you. If it still doesn't work then check to see if there is any help in the [troubleshooting](#troubleshooting) section of this README. 

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


# Run script without using Shiny App

 <details>
  <summary>Script</summary>
  
```
### Set working directory to wherever "source_functions.R" is
setwd("~/Group3-4App")
source("./AppSourceFunctions1.11.R")

#####################################
############ METHYLATION ############
#####################################



### load in the prediction object
load(file = "./AppExtraFiles/Inputs/g3.g4.cont.rfe.Rdata")

### Choose folder containing idats to be processed
idats <- "~/your/idat/folder"

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

## Load in your samples
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
t(rnaseq.H[c(3, 1), ]) -> g3g4.rnaseq

apply(g3g4.rnaseq, 2, function(x) {
  (1 / (1 + exp(-x)))
}) -> logistic.g3g4.rnaseq

apply(logistic.g3g4.rnaseq, 1, function(x) {
  x[2] / (x[1] + x[2])
}) -> logistic.g3g4.rnaseq.score

## Scale to Williamson et al. dataset
scaling.function3(logistic.g3g4.rnaseq.score) -> logistic.g3g4.rnaseq.score
ScalingChoice <- "ours"

## Scale to uploaded dataset 
#scaling.function(logistic.g3g4.rnaseq.score) -> logistic.g3g4.rnaseq.score
#ScalingChoice <- "yours"

## Outlier removal
outlier <- 1


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
                
```
</details>

# Example script to run
 <details>
  <summary>Example Script</summary>

```

### Set working directory to wherever "source_functions.R" is
setwd("~/Group3-4App")
source("./AppSourceFunctions1.11.R")

#####################################
############ METHYLATION ############
#####################################

# Obtain MValues (these have been obtained from idat files, for full workflow read above)
M.values <- read.delim("~/Group3-4App/StarProtocols_Guide/data/mvals.mat.txt")

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
t(rnaseq.H[c(3, 1), ]) -> g3g4.rnaseq

apply(g3g4.rnaseq, 2, function(x) {
  (1 / (1 + exp(-x)))
}) -> logistic.g3g4.rnaseq

apply(logistic.g3g4.rnaseq, 1, function(x) {
  x[2] / (x[1] + x[2])
}) -> logistic.g3g4.rnaseq.score

## Scale to Williamson et al. dataset
scaling.function3(logistic.g3g4.rnaseq.score) -> logistic.g3g4.rnaseq.score
ScalingChoice <- "ours"

## Scale to uploaded dataset 
#scaling.function(logistic.g3g4.rnaseq.score) -> logistic.g3g4.rnaseq.score
#ScalingChoice <- "yours"

## Outlier removal
outlier <- 1


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


```
 </details>
 
## Example Script Results
<details>
  <summary>Methylation</summary>


 For methylation the expected results are:  
 
 ```
> metagene.df
        Group.3.4.Score
NMB_109           0.079
NMB_110           0.178
NMB_111           0.319
NMB_118           0.416
NMB_119           0.158
NMB_125           0.057
NMB_130           0.095
NMB_132           0.291
NMB_134           0.203
NMB_136           0.375
```
  
 The following graphs will also be created:  
 ![M1](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/M1.png)
 ![M2](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/M2.png)
 ![M3](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/M3.png)
 
 </details>
 
 <details>
  <summary>Expression</summary>
  
```
> logistic.g3g4.tpms.score.df
                                 Group.3.4.Score
GSM3905407_BCH825_BCH825.P02.A01           0.597
GSM3905407_BCH825_BCH825.P02.A03           0.722
GSM3905407_BCH825_BCH825.P02.A05           0.806
GSM3905407_BCH825_BCH825.P02.A07           0.518
GSM3905407_BCH825_BCH825.P02.A08           0.763
GSM3905407_BCH825_BCH825.P02.A09           0.761
GSM3905407_BCH825_BCH825.P02.A10           0.521
GSM3905407_BCH825_BCH825.P02.A11           0.802
GSM3905407_BCH825_BCH825.P02.B01           0.744
GSM3905407_BCH825_BCH825.P02.B02           0.856

```

 The following graphs will also be created:  
 ![E1](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/E1.png)
 ![E2](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/E2.png)
 ![E3](https://github.com/hackingjpr/Group3-4App/blob/main/AppExtraFiles/Tutorial/E3.png)

   </details>
 

# System Requirements
## Hardware Requirements
Functions provided here are compatible with any standard computer with enough RAM to support the in-memory operations.

## Software Requirements
Functions provided import and depend on a number of R packages. Functionality has been tested on *R 4.0.4* with *Ubuntu 20.04.2 LTS*

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




# *Disclaimer: This app is designed exclusively for research purposes and strictly not for diagnostic use.*
