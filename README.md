# WHAT IS THE TITLE, WILL CHANGE IF THIS IS JUST FOR APP OR ALSO FOR DEAN'S STUFF


- [Overview](#overview)
- [Tutorial](#tutorial)
  - [Step One - Select Expression or Methylation](#step-one---Select-Expression-or-Methylation)
  - [Step Two - Uploading Data](#step-two---Uploading-Data)
  - [Step Three - Generate Group 3/4 Scores](#step-three---Generate-Group-3/4-Scores)
  - [Step Four - Results](#step-four---results)
  - [Step Five - Export Data and Reset](#step-five---export-data-and-reset)
- [System Requirements](#system-requirements)

# Overview
This script "app.R" encodes a shiny app that upon uploading idat files, will give a group 3/4 score for patient samples.   
The repository can be cloned onto your RStudio or the base code can be run independently, this is shown underneath.

# Tutorial
## Step One - Select Expression or Methylation

Depending on whether you are uploading Expression or Methylation data select the appropriate option.

## Step Two - Uploading Data

Upload your idat files (for now unzipped idat files only) including both red and green files for Methylation, or RDS file for Expression.

![upload](upload.png)

Upload a minimum of two samples if running idats. Increasing the number of samples will of course increase the length of time for the upcoming processes so we recommend \~10 sample batches. This will make looking through the results easier and will speed up the process.

# Step Three - Generate Group 3/4 Scores

This will start the process of generating Group 3/4 Continuum Scores.

# Step Four - Results

Once the calculation has been completed you should be brought to the Results tab. This tab will show a data table at the top which displays your sample names on the left and their Group 3/4 Scores on the right.

![Risk Values](risk_values.png)

It will also show a graph (as shown below), the details of this graph and what it shows can be found in the corresponding paper (link)\
(insert graph example)\
There will also be a green box in the bottom right which will inform you of the currently selected sample (which will be highlighted in orange on the graph), and the selected metagene set.

![Selections](selections.png)

# Step Five - Export Data and Reset

Once you have your results you can either reset the app by pressing the reset button:

![Reset](Reset.png)

Or download your data as a CSV file (data table results), or as a PDF (data table and graph).\
**Download functionality is not currently supported but will be by the time the paper is published.**


# Run script without using Shiny App
```
### Set working directory to wherever "source_functions.R" is
setwd("/example/Idat-Shiny")
source("./source_functions.R")

### Choose folder containing idats to be processed
idats <- "/example/idatfolder/"

### Get Basenames
temp.base <- get_basenames(idats)

### Process Idats
temp.processed <- process_idats(temp.base)

### Select metagene set
metagene <- ALL
#metagene <- ATRT
#metagene <- ECRT

### Extract Metagenes (This will be your risk values result)
test.res <- extract.metagene(
  as.character(metagene[[1]]$genes),
  as.numeric(metagene[[1]]$weights),
  beta2m(temp.processed$betas),
  as.numeric(metagene[[2]])
)

### Round results to 3 figures
round(test.res, digits = 3)
```

# Example script to run
```
### Set working directory to wherever "source_functions.R" is
setwd("/example/Idat-Shiny")
source("./source_functions.R")

### Choose folder containing idats to be processed

baseDir <- system.file("extdata", package = "minfiData")
list.files(baseDir)

idats <- (file.path(baseDir, "5723646052"))


### Get Basenames
temp.base <- get_basenames(idats)

### Process Idats
temp.processed <- process_idats(temp.base)

### Select metagene set
# For MRT
metagene <- ALL
#For ATRT
#metagene <- ATRT
#For ECRT
#metagene <- ECRT

### Extract Metagenes (This will be your risk values result)
test.res <- extract.metagene(
  as.character(metagene[[1]]$genes),
  as.numeric(metagene[[1]]$weights),
  beta2m(temp.processed$betas),
  as.numeric(metagene[[2]])
)

### Round results to 3 figures
round(test.res, digits = 3)

#For ALL should give : -0.325, -0.416, 0.222
#For ATRT should give : -1.423, -1.373, -1.678
#For ECRT should give : -0.071, -0.033, -0.122

### Select Risk values column
figure.input <- test.res$Risk_Value

### Name the rows
names(figure.input) <- rownames(test.res)
print(figure.input)

### Pick the generate_figure_highlight that is required
generate_figure_highlight_mrt(figure.input,
                              NA)

# generate_figure_highlight_atrt(figure.input,
#                                NA)

# generate_figure_highlight_ecrt(figure.input,
#                                NA)

```
## Example Script Results
For ALL(MRT) should give : -0.325, -0.416, 0.222  
For ATRT should give : -1.423, -1.373, -1.678  
For ECRT should give : -0.071, -0.033, -0.122  

For MRT the plot should look like this:

![MRT Example Plot](https://github.com/hackingjpr/Idat-Shiny/blob/main/Tutorial/MRTExamplePlot.png?raw=true)  

# System Requirements
## Hardware Requirements
Functions provided here are compatible with any standard computer with enough RAM to support the in-memory operations.

## Software Requirements
Functions provided import and depend on a number of R packages. Functionality has been tested on *R 4.0.4* with *Ubuntu 20.04.2 LTS*



# *Disclaimer : This app is designed exclusively for research purposes and strictly not for diagnostic use.*
