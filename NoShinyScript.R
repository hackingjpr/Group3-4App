### Set working directory to where "source_functions.R" is
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

### Extract Metagenes
test.res <- extract.metagene(
  as.character(metagene[[1]]$genes),
  as.numeric(metagene[[1]]$weights),
  beta2m(temp.processed$betas),
  as.numeric(metagene[[2]])
)

### Round results to 3 figures
round(test.res, digits = 3)

### Select Risk values
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