### Set working directory to wherever "source_functions.R" is
setwd("/example/Idat-Shiny")
source("./source_functions.R")


### load in the prediction object
load(file = "./g3.g4.cont.rfe.Rdata")
pred.cont.rand.for <- predict(g3.g4.cont.rfe, t(M.values)[,predictors(g3.g4.cont.rfe)])

### I have attached some know values you can read here
pred.cont.rand.for <- readRDS(file = "./pred.cont.rand.for.rds")

### Choose folder containing idats to be processed
idats <- "./Mix"

### Get Basenames
temp.base <- get_basenames(idats)

### Process Idats
temp.processed <- process_idats(temp.base)

# Obtain MValues
beta2m(temp.processed$betas) -> M.values

### Round results to 3 figures
round(test.res, digits = 3)

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