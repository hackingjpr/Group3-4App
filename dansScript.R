source("./source_projection_scripts.R")

# load in original metagene H values
avg.h.val <- readRDS("./avg.h.val.input.rds")

# load in original g3g4 metagene values
g3g4 <- readRDS("./g3g4.input.rds")

#load in filtered expression values from NMB MB samples used to generate H Values
nmb.mat <- readRDS("./nmb.mat.input.rds")

# load in a second object to project on to
tpms.mat <- readRDS(file = "./tpms.mat.rds")

# convert g3 and g4 metagenes into g3g4 score
apply(g3g4,2,function(x){(1 / (1 + exp(-x)))}) -> logistic.g3g4
apply(logistic.g3g4,1,function(x){x[2]/(x[1]+x[2])}) -> logistic.g3g4.score
scaling.function(logistic.g3g4.score) -> logistic.g3g4.score

# annotate nmb.mat to change into hgnc symbols
annotate.HTseq.IDs(rownames(nmb.mat)) -> annotation
nmb.mat[-which(annotation$hgnc_symbol==""|is.na(annotation$hgnc_symbol)),] -> nmb.mat
annotation$hgnc_symbol[-which(annotation$hgnc_symbol==""|is.na(annotation$hgnc_symbol))] -> rownames(nmb.mat)

## pre-filter datasets
idx <- param.filter(2^nmb.mat)
nmb.mat <- 2^nmb.mat[idx, ]

## pre-projection normalisation
nmb.mat <- prep.data(nmb.mat)

## interset common genes / probes
tpms.mat <- match.select(nmb.mat, tpms.mat)

## train NMF model
library(NMF)
## NEED PACKAGE NMF
init <- NMF::nmfModel(4, 
                 nmb.mat, 
                 W = 0.5, 
                 H = t(avg.h.val))

## generate NMF seeded with model
nmf.res <- NMF::nmf(nmb.mat, 
               4, 
               seed = init, 
               nrun = 8, 
               .pbackend = 20)

## project using pseudo-inverse & post-projection normalise
# project back onto the same dataset
rnaseq.H <- project.NMF(input.array = nmb.mat, 
                        nmf.result = nmf.res)

test <- as.matrix(tpms.mat)
# project onto fresh dataset
tpms.H <- project.NMF(input.array = test, 
                      nmf.result = nmf.res)
# tpms.H <- project.NMF(input.array = tpms.mat, 
#                       nmf.result = nmf.res)

### define new g3g4 score for projection back onto the original data
t(rnaseq.H[c(3,1),]) -> g3g4.rnaseq
apply(g3g4.rnaseq,2,function(x){(1 / (1 + exp(-x)))}) -> logistic.g3g4.rnaseq
apply(logistic.g3g4.rnaseq,1,function(x){x[2]/(x[1]+x[2])}) -> logistic.g3g4.rnaseq.score
scaling.function(logistic.g3g4.rnaseq.score) -> logistic.g3g4.rnaseq.score

## join and plot the two together (original g3g4 score and g3g4 score projected back onto the same data) kind of a control that it is working
## NEED TIDYFT
library(tidyft)
df <- inner_join(data.frame(logistic.g3g4.score,ids = names(logistic.g3g4.score)),
                 data.frame(logistic.g3g4.rnaseq.score, ids = names(logistic.g3g4.rnaseq.score)))
plot(df$logistic.g3g4.score,df$logistic.g3g4.rnaseq.score)

# extract only the g3 and g4 metagenes and convert to g3g4 score
t(tpms.H[c(3,1),]) -> g3g4.tpms
apply(g3g4.tpms,2,function(x){(1 / (1 + exp(-x)))}) -> logistic.g3g4.tpms
apply(logistic.g3g4.tpms,1,function(x){x[2]/(x[1]+x[2])}) -> logistic.g3g4.tpms.score

# as.data.frame(logistic.g3g4.tpms.score)

# some times helpful to remove outliers prior to scaling
outlier.idx <- c(head(order(logistic.g3g4.tpms.score), round(5*(length(logistic.g3g4.tpms.score)/100))),
                 tail(order(logistic.g3g4.tpms.score), round(5*(length(logistic.g3g4.tpms.score)/100)))
)

# scale to create final score
scaling.function(logistic.g3g4.tpms.score[-outlier.idx]) -> logistic.g3g4.tpms.score

round(logistic.g3g4.tpms.score, digits = 3)


as.data.frame(logistic.g3g4.tpms.score)
