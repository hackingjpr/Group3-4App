### Set working directory to wherever "source_functions.R" is
setwd("~/Group3-4App")
source("./AppSourceFunctions1.3.R")

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