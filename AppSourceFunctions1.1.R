########## Methylation ################


print(getwd())

library(BiocManager)
options(repos = BiocManager::repositories())
# if(!require("devtools"))
#   install.packages("devtools")
# library(rsconnect)
# if(!require(devtools)) {
#   BiocManager::install("devtools")
#   library(devtools)
# }
# 
# if(!require(rsconnect)) {
#   BiocManager::install("rsconnect")
#   library(rsconnect)
# }


if(!require(minfiData)) {
  BiocManager::install("minfiData")
  library(minfiData)
}
message("minfidata done")
if(!require(sva)){
  BiocManager::install("sva")
  library(sva)
}
message("sva done")
if(!require(bumphunter)){
  BiocManager::install("bumphunter")
  library(bumphunter)
}
message("bumphunter done")
# if(!require(lumi)){
#   BiocManager::install("lumi")
#   library(lumi)
# }
# message("lumi done")
# 
# if(!require(minfi)){
#   BiocManager::install("minfi")
#   library(minfi)
# }
# 
# message("minfi done")
if(!require(IlluminaHumanMethylationEPICmanifest)){
  BiocManager::install("IlluminaHumanMethylationEPICmanifest")
  library(IlluminaHumanMethylationEPICmanifest)
}
# message("illumina1 done")
if(!require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)){
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
}
#message("manifest loaded")
if(!require(shiny)){
  install.packages("shiny")
  library(shiny)
}
message("shiny done")

if(!require(DT)){
  install.packages("DT")
  library(DT)
}

message("DT done")

if(!require(shinyWidgets)){
  install.packages("shinyWidgets")
  library(shinyWidgets)
}
message("shinyWidgets done")


if(!require(shinydashboard)){
  install.packages("shinydashboard")
  library(shinydashboard)
}
message("shinydashboard done")


if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
message("ggplot2 done")


if(!require(ggpubr)){
  install.packages("ggpubr")
  library(ggpubr)
}
message("ggpubr done")

if(!require(ggnewscale)){
  install.packages("ggnewscale")
  library(ggnewscale)
}
message("ggnewscale done")

if(!require(foreach)){
  install.packages("foreach")
  library(foreach)
}
message("foreach done")

if(!require(ggrepel)){
  install.packages("ggrepel")
  library(ggrepel)
}
message("ggrepel done")

if(!require(pheatmap)){
  install.packages("pheatmap")
  library(pheatmap)
}
message("pheatmap done")

if(!require(shinycssloaders)){
  install.packages("shinycssloaders")
  library(shinycssloaders)
}
message("shinycssloaders done")

if(!require(shinybusy)){
  install.packages("shinybusy")
  library(shinybusy)
}
message("shinybusy done")

if(!require(waiter)){
  install.packages("waiter")
  library(waiter)
}
message("waiter done")

if(!require(tinytex)){
  install.packages('tinytex')
  library(tinytex)
}

if(!require(mlbench)){
  install.packages('mlbench')
  library(mlbench)
}

if(!require(caret)){
  install.packages('caret')
  library(caret)
}

if(!require(NMF)){
  install.packages('NMF')
  library(NMF)
}

if(!require(tidyft)){
  install.packages('tidyft')
  library(tidyft)
}

if(!require(plotly)){
  install.packages('plotly')
  library(plotly)
}

if(!require(randomForest)){
  packageurl <- "https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-14.tar.gz"
  install.packages(packageurl, repos=NULL, type="source")
  library(randomForest)
}


library(gridExtra)
library(grid)


load(file = "./g3.g4.cont.rfe.Rdata")



beta2m <- function (beta) {
  m <- log2(beta/(1 - beta))
  return(m)
}



process_idats <- function(basenames){
  
  #check file sizes, epic will be over 1e7
  # Split into mixed or single test e.g. 450k & Epic or one or the other
  file.sizes <- file.info(paste0(basenames, "_Red.idat"))$size
  
  length(which(file.sizes> 1e7))>0 & length(which(file.sizes< 1e7))>0 -> mixed.test
  length(which(file.sizes> 1e7))>0 | length(which(file.sizes< 1e7))>0 -> single.test
  
  if(mixed.test) {
    #  as before where you process them seperately then combine
    idx.epic <-
      file.info(paste0(basenames, "_Red.idat"))$size > 1e7
    m450k <- read.metharray(basenames[!idx.epic],
                            verbose = TRUE,
                            force = TRUE) 
    mEpic <- read.metharray(basenames[idx.epic],
                            verbose = TRUE,
                            force = TRUE) 
    comb.rgSet <- combineArrays(m450k, 
                                mEpic,
                                outType = "IlluminaHumanMethylation450k")
    comb.mSet <- preprocessNoob(comb.rgSet, 
                                dyeMethod = "single")
    comb.detP <- detectionP(comb.rgSet)
    comb.gmrSet <- mapToGenome(ratioConvert(comb.mSet))
    comb.beta <- getBeta(comb.gmrSet)
    comb.MRT.pd <- pData(comb.gmrSet)
  }
  
  if((!mixed.test)&single.test){
    # only need the one
    single <- read.metharray(basenames,
                             verbose = TRUE,
                             force = TRUE)
    comb.mSet <- preprocessNoob(single, dyeMethod = "single")
    comb.detP <- detectionP(single)
    comb.gmrSet <- mapToGenome(ratioConvert(comb.mSet))
    comb.beta <- getBeta(comb.gmrSet)
    comb.MRT.pd <- pData(comb.gmrSet)
  }
  
  return(list(betas = comb.beta, pvals = comb.detP))
}

get_basenames <- function(dir){
  list.files(path = dir, pattern = "_Red.idat", full.names = T) -> temp.files
  return(gsub("_Red.idat","", temp.files))
}

extract.metagene <- function(index, weights, exp.matrix, scaling) {
  exp.matrix[index, ] -> temp.exp
  apply(temp.exp, 2, function(x) {
    mean(x * weights)
  }) -> raw.metagene
  round(raw.metagene, digits = 3)
  as.numeric(scale(raw.metagene, center = scaling[1], scale = scaling[2])) -> Risk_Value
  round(Risk_Value, digits = 3)
  names(Risk_Value) <- colnames(exp.matrix)
  df <- as.data.frame(Risk_Value)
  return(df)
}

digits = 0:9
createRandString<- function() {
  v = c(sample(LETTERS, 5, replace = TRUE),
        sample(digits, 4, replace = TRUE),
        sample(LETTERS, 1, replace = TRUE))
  return(paste0(v,collapse = ""))
}



############# Expression ########################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require(biomaRt)){
  BiocManager::install("biomaRt")
  library(biomaRt)
}

# load in original metagene H values
avg.h.val <- readRDS("./avg.h.val.input.rds")

# load in original g3g4 metagene values
g3g4 <- readRDS("./g3g4.input.rds")

#load in filtered expression values from NMB MB samples used to generate H Values
#nmb.mat <- readRDS("./nmb.mat.input.rds")


nmb.mat.prepped <- readRDS(file = "./nmb.mat.prepped.rds")

#tpms.mat <- readRDS(file = "./tpms.mat.rds")

### annotate using library(biomaRt)
annotate.HTseq.IDs<-function(HTseq.IDs){
  ENSEMBL_DB_HOST = "uswest.ensembl.org" # Set back to default, once they are up and running again
  ENSEMBL_VERSION = "Ensembl Genes 101"  # Try to fix https://support.bioconductor.org/p/104454/
  
  #mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = ENSEMBL_DB_HOST, version = ENSEMBL_VERSION)
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  ensemblID <- gsub("\\..*", '', HTseq.IDs)
  #fetch the annotation information from the Ensembl database
  #ensembl <- useMart('ensembl', dataset='hsapiens_gene_ensembl')
  symbols <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), filters='ensembl_gene_id', ensemblID, mart=mart)
  #combine the annotation with the RNA-Seq data
  annotatedix <- match(ensemblID,symbols$ensembl_gene_id)
  symbols[annotatedix,] -> annotatedGenes
  return(cbind(HTseq.IDs,annotatedGenes))
}


### loading in all functions needed for projection
match.select <- function(input1.eset, input2.eset){
  
  m1 <- input1.eset
  gs.names1 <- row.names(m1)
  gs.descs1 <- row.names(m1)
  sample.names1 <- names(m1)
  
  m2 <- input2.eset
  gs.names2 <- row.names(m2)
  gs.descs2 <- row.names(m2)
  sample.names2 <- names(m2)
  
  gs.names3 <- intersect(gs.names1, gs.names2)
  
  locations2 <- match(gs.names3, gs.names2, nomatch=0)
  gs.names2 <- gs.names2[locations2]
  gs.descs2 <- gs.descs2[locations2]
  m2 <- m2[locations2, ]
  
  return(m2)
}

prep.data <- function(input.array){
  
  m <- input.array
  col.names <- colnames(input.array)
  gs.names <- row.names(input.array)
  gs.descs <- row.names(input.array)
  sample.names <- colnames(input.array)
  
  cols <- length(m[1,])
  for (j in 1:cols) {  # column rank normalization from 0 to N - 1
    m[,j] <- rank(m[,j], ties.method = "average") - 1
  }
  m <- 10000*m/(length(m[,1]) - 1)
  
  v <- m
  row.names(v) <- gs.names
  colnames(v) <- col.names
  
  return(v)
}

param.filter <- function(m = NULL,
                         thres = 20,
                         ceil = 10000,
                         shift = 1,
                         
                         fold = 3,
                         delta = 300){
  if (thres != "NULL") {
    m[m < thres] <- thres
  }
  if (ceil != "NULL") {
    m[m > ceil] <- ceil
  }
  if (shift != "NULL") {
    m <- m + shift
  }
  
  cols <- ncol(m)
  rows <- nrow(m)
  row.max <- apply(m, MARGIN = 1, FUN = max)
  plot(sort(row.max, decreasing = TRUE), type = "l", main = "row max", ylab = "Value", xlab = "Index")
  
  row.min <- apply(m, MARGIN = 1, FUN = min)
  
  flag <- array(dim = rows)
  flag <- (row.max / row.min >= fold) & (row.max - row.min >= delta)
  cat(length(which(flag)), "genes / probes retained", "\n")
  flag
}

scaling.function <- function(x){(x-min(x))/(max(x)-min(x))} 



param.filter <- function(m = NULL,
                         thres = 20,
                         ceil = 10000,
                         shift = 1,
                         fold = 3,
                         delta = 300){
  if (thres != "NULL") {
    m[m < thres] <- thres
  }
  if (ceil != "NULL") {
    m[m > ceil] <- ceil
  }
  if (shift != "NULL") {
    m <- m + shift
  }
  
  cols <- ncol(m)
  rows <- nrow(m)
  row.max <- apply(m, MARGIN = 1, FUN = max)
  plot(sort(row.max, decreasing = TRUE), type = "l", main = "row max", ylab = "Value", xlab = "Index")
  
  row.min <- apply(m, MARGIN = 1, FUN = min)
  
  flag <- array(dim = rows)
  flag <- (row.max / row.min >= fold) & (row.max - row.min >= delta)
  cat(length(which(flag)), "genes / probes retained", "\n")
  flag
}

prep.data <- function(input.array){
  
  m <- input.array
  col.names <- colnames(input.array)
  gs.names <- row.names(input.array)
  gs.descs <- row.names(input.array)
  sample.names <- colnames(input.array)
  
  cols <- length(m[1,])
  for (j in 1:cols) {  # column rank normalization from 0 to N - 1
    m[,j] <- rank(m[,j], ties.method = "average") - 1
  }
  m <- 10000*m/(length(m[,1]) - 1)
  
  v <- m
  row.names(v) <- gs.names
  colnames(v) <- col.names
  
  return(v)
}
match.select <- function(input1.eset, input2.eset){
  
  m1 <- input1.eset
  gs.names1 <- row.names(m1)
  gs.descs1 <- row.names(m1)
  sample.names1 <- names(m1)
  
  m2 <- input2.eset
  gs.names2 <- row.names(m2)
  gs.descs2 <- row.names(m2)
  sample.names2 <- names(m2)
  
  gs.names3 <- intersect(gs.names1, gs.names2)
  
  locations2 <- match(gs.names3, gs.names2, nomatch=0)
  gs.names2 <- gs.names2[locations2]
  gs.descs2 <- gs.descs2[locations2]
  m2 <- m2[locations2, ]
  
  return(m2)
}
project.NMF <- function(input.array, nmf.result){
  require(MASS)
  m <- input.array
  gs.names <- rownames(input.array)
  gs.descs <- rownames(input.array)
  sample.names <- colnames(input.array)
  
  W <- as.data.frame(basis(nmf.result))
  W.row.names <- row.names(W)
  W.row.descs <- row.names(W)
  W.names <- names(W)
  
  overlap <- intersect(gs.names, W.row.names)
  
  cat("Size of Input dataset:", length(gs.names), "genes\n")
  cat("Size of W matrix (rows):", length(W.row.names), "genes\n")
  cat("Size of overlap:", length(overlap), "genes\n")
  
  locations.m <- match(overlap, gs.names, nomatch=0)
  m2 <- m[locations.m, ]
  locations.W <- match(overlap, W.row.names, nomatch=0)
  W2 <- W[locations.W, ]
  W2 <- as.matrix(W2)
  
  H <- ginv(W2) %*% m2
  
  n.col <- length(H[1,])
  for (i in 1:n.col) {
    S.2 <- sqrt(sum(H[,i]*H[,i]))
    #        S.2 <- sum(H[,i])
    H[,i] <- H[,i]/S.2
  }
  
  V <- data.frame(H)
  names(V) <- sample.names
  row.names(V) <- W.names
  return(V)
}


annotate.HTseq.IDs<-function(HTseq.IDs){
  ENSEMBL_DB_HOST = "uswest.ensembl.org" # Set back to default, once they are up and running again
  ENSEMBL_VERSION = "Ensembl Genes 101"  # Try to fix https://support.bioconductor.org/p/104454/
  
  #mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = ENSEMBL_DB_HOST, version = ENSEMBL_VERSION)
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  ensemblID <- gsub("\\..*", '', HTseq.IDs)
  #fetch the annotation information from the Ensembl database
  #ensembl <- useMart('ensembl', dataset='hsapiens_gene_ensembl')
  symbols <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), filters='ensembl_gene_id', ensemblID, mart=mart)
  #combine the annotation with the RNA-Seq data
  annotatedix <- match(ensemblID,symbols$ensembl_gene_id)
  symbols[annotatedix,] -> annotatedGenes
  return(cbind(HTseq.IDs,annotatedGenes))
}

################ graph test##################
generate_figure_highlight_ecrt <- function(new.sample.meta.score, indexRow) {
  if(is.null(indexRow)){indexRow=1}
  #temp.df <- readRDS(file = "https://github.com/hackingjpr/Idat-Shiny/blob/main/temp.df.rds")
  temp.df <- readRDS(file = "./ecrt20.dist.rds")
  #df.cat.atrt <- readRDS(file = "https://github.com/hackingjpr/Idat-Shiny/blob/main/df.cat.atrt.rds")
  # df.cat.atrt <- readRDS(file = "./df.cat.atrt.rds")
  #comb.SDb.atrt <- readRDS(file = "https://github.com/hackingjpr/Idat-Shiny/blob/main/comb.SDb.atrt.rds")
  # comb.SDb.atrt <- readRDS(file = "./comb.SDb.atrt.rds")
  ggplot(aes(x = 1:nrow(temp.df), y = ecrt20), data = temp.df) +
    geom_line() +
    scale_shape_manual(values = c(1, 4,  3)) +
    scale_color_manual(values = c('#E69F00', '#999999', "white")) +
    ylab("ATRT-8") +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    ) -> b
  
  df.lines.hor <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = 0,
                xend = max(which(
                  temp.df$ecrt20 < new.sample.meta.score[i]
                )),
                y = new.sample.meta.score[i],
                yend = new.sample.meta.score[i]
              )
            }
  df.lines.hor$labels <- names(new.sample.meta.score)
  
  df.lines.ver <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = max(which(
                  temp.df$ecrt20 < new.sample.meta.score[i]
                )),
                xend = max(which(
                  temp.df$ecrt20 < new.sample.meta.score[i]
                )),
                y = new.sample.meta.score[i],
                yend = min(temp.df$ecrt20)
              )
            }
  df.lines.ver$perc <-
    paste0(round(df.lines.ver$xend / length(temp.df$ecrt20) * 100), "th")
  
  df.lines.ver$colour <- factor(ifelse(1:nrow(df.lines.ver)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  df.lines.hor$colour <- factor(ifelse(1:nrow(df.lines.hor)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  
  b <- b +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = as.character(colour)
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.hor
    ) +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = colour
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.ver
    ) +
    #geom_text_repel(aes(x = x+10, y = y+0.1, label = labels), direction = "y", data = df.lines.hor)
    geom_text(aes(
      x = x + 10,
      y = y + 0.1,
      label = labels,
      colour = colour
    ), data = df.lines.hor) 
  #scale_color_discrete("red", "lightgrey")
  
  df <- data.frame()
  c <- ggplot() + theme_void()
  
  
  #ggarrange(a,a2,ggarrange(
  # c,b, ncol = 2, nrow = 1, widths = c(0.015,1)),ncol=1,nrow=3)
  ggarrange(c,
            b,
            ncol = 2,
            nrow = 1,
            widths = c(0.015, 1))
  #return(data.frame(perc = df.lines.ver[,5], row.names=rownames(df.lines.ver)))
}

df1 <- data.frame(green = c(0,0.68),
                  yellow = c(0.68,1),
                  colour = c("G", "Y"))


generate_figure_highlight_g3g4 <- function(new.sample.meta.score, indexRow){
  if(is.null(indexRow)){indexRow=1}
  y <- readRDS(file = "./y.vals.rds")
  df <- as.data.frame(y)
  ecdf(y) -> model
  model(y) -> y2
  b <- ggplot(df, aes(x=y2, y=y)) +
    geom_rect(data = df1, aes(NULL,NULL, xmin=green, xmax=yellow, fill=colour),
              ymin=0,ymax=1, colour="white", size=0.5, alpha=0.2) +
    scale_fill_manual(values=c("G"= "green", "Y" = "yellow")) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_hline(yintercept=0.5, linetype="dashed") +
    geom_hline(yintercept=1, linetype="dashed") +
    geom_vline(xintercept=0.68, linetype="dashed") +
    geom_line() +
    theme(legend.position = "none")
  
  
  
  
  
  df.lines.hor <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                # x = 0,
                x = 0,
                xend = y2[max(which(
                  y < new.sample.meta.score[i]
                ))],
                y = y[max(which(
                  y < new.sample.meta.score[i]
                ))],
                yend = y[max(which(
                  y < new.sample.meta.score[i]
                ))]
              )
            }
  df.lines.hor$labels <- names(new.sample.meta.score)
  
  df.lines.ver <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = y2[max(which(
                  y < new.sample.meta.score[i]
                ))],
                xend = y2[max(which(
                  y < new.sample.meta.score[i]
                ))],
                y =  y[max(which(
                  y < new.sample.meta.score[i]
                ))],
                yend = 0
              )
            }
  df.lines.ver$perc <-
    paste0(round(df.lines.ver$xend / length(y2) * 100), "th")
  
  df.lines.ver$colour <- factor(ifelse(1:nrow(df.lines.ver)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  df.lines.hor$colour <- factor(ifelse(1:nrow(df.lines.hor)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  
  
  b <- b +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = as.character(colour)
      ),
      #colour = "red",
      linetype = "dashed",
      size = 1,
      data = df.lines.hor
    ) +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = colour
      ),
      #colour = "red",
      linetype = "dashed",
      size = 1,
      data = df.lines.ver
      
    ) +
    xlab("") + 
    ylab("G3/G4 Score")
  #geom_text_repel(aes(x = x+10, y = y+0.1, label = labels), direction = "y", data = df.lines.hor)
  # geom_text(aes(
  #   x = x + 10,
  #   y = y + 0.1,
  #   label = labels,
  #   colour = colour
  # ), data = df.lines.hor)
  #scale_color_discrete("red", "lightgrey")
  
  df <- data.frame()
  c <- ggplot() + theme_void()
  
  
  #ggarrange(a,a2,ggarrange(
  # c,b, ncol = 2, nrow = 1, widths = c(0.015,1)),ncol=1,nrow=3)
  d<- ggarrange(c,
                b,
                ncol = 2,
                nrow = 1,
                widths = c(0.015, 1))
  
  # d <- subplot(c,
  #              b,
  #              ncol=2,
  #              nrows = 1,
  #              widths = c(0.015,1))
  
  # d <- ggplotly(d)
  d
  #return(data.frame(perc = df.lines.ver[,5], row.names=rownames(df.lines.ver)))
}



generate_figure_highlight_g3g4PERC <- function(new.sample.meta.score, indexRow){
  if(is.null(indexRow)){indexRow=1}
  y <- readRDS(file = "./y.vals.rds")
  df <- as.data.frame(y)
  ecdf(y) -> model
  model(y) -> y2
  b <- ggplot(df, aes(x=y2, y=y)) +
    geom_rect(data = df1, aes(NULL,NULL, xmin=green, xmax=yellow, fill=colour),
              ymin=0,ymax=1, colour="white", size=0.5, alpha=0.2) +
    scale_fill_manual(values=c("G"= "green", "Y" = "yellow")) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_hline(yintercept=0.5, linetype="dashed") +
    geom_hline(yintercept=1, linetype="dashed") +
    geom_vline(xintercept=0.68, linetype="dashed") +
    geom_line() +
    theme(legend.position = "none")
  
  
  
  
  
  df.lines.hor <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                # x = 0,
                x = 0,
                xend = y2[max(which(
                  y < new.sample.meta.score[i]
                ))],
                y = y[max(which(
                  y < new.sample.meta.score[i]
                ))],
                yend = y[max(which(
                  y < new.sample.meta.score[i]
                ))]
              )
            }
  df.lines.hor$labels <- names(new.sample.meta.score)
  
  df.lines.ver <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = y2[max(which(
                  y < new.sample.meta.score[i]
                ))],
                xend = y2[max(which(
                  y < new.sample.meta.score[i]
                ))],
                y =  y[max(which(
                  y < new.sample.meta.score[i]
                ))],
                yend = 0
              )
            }
  df.lines.ver$perc <-
    # paste0(round(df.lines.hor$xend / length(y2) * 100), "th")
    paste0(round(df.lines.hor$xend / 1 * 100), "th")
  
  return(df.lines.ver$perc[indexRow])
}

time.comb <- readRDS(file = "./time.comb.rds")
age.comb <- readRDS(file = "./age.comb.rds")
status.comb <- readRDS(file = "./status.comb.rds")
comb.cont <-readRDS(file = "./comb.cont.rds")

library(survival)


# new.sample.meta.score <- c(0.5,0.3)
# indexRow <- 0.3

coxph(Surv(time.comb, status.comb) ~ comb.cont) -> train.fit
message("Train fit complete")
summary(survfit(train.fit, data.frame(g3g4.values=comb.cont)), time = 5) -> x
df2 <- data.frame(pred = comb.cont,
                  surv = as.numeric(x$surv),
                  up = as.numeric(x$upper),
                  lo = as.numeric(x$lower))

survivalcurveplot <- function(new.sample.meta.score,indexRow){
  if(is.null(indexRow)){indexRow=1}
  df2$pred -> pred
  df2$surv -> surv
  b<- ggplot(df2, aes(x=pred, y=surv)) +
    geom_line() +
    geom_point(alpha = 1/20) +
    geom_line(aes(x=pred, y=lo),linetype="dotted") +
    geom_line(aes(x=pred, y=up),linetype="dotted") +
    theme_classic() + xlab("G3/G4 Score") + ylab("Survival") +
    labs(title = "New plot title", subtitle = "A subtitle") +
    ylim(0,1) +
    theme(legend.position = "none")
  
  ggplotly(b)
  
  df.lines.hor <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              surv[which(
                pred < new.sample.meta.score[i]
              )] -> temp.surv
              data.frame(
                x = 0,
                xend = new.sample.meta.score[i],
                y = temp.surv[which.min(temp.surv)],
                yend = temp.surv[which.min(temp.surv)]
              )
            }
  
  df.lines.hor$labels <- names(new.sample.meta.score)

  
  df.lines.ver <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              surv[which(
                pred < new.sample.meta.score[i]
              )] -> temp.surv
              data.frame(
                x = new.sample.meta.score[i],
                xend = new.sample.meta.score[i],
                y = 0,
                yend = temp.surv[which.min(temp.surv)]
              )
            }
  message(df.lines.ver)
  
  
  df.lines.ver$perc <-
    paste0(round(df.lines.ver$xend / length(pred) * 100), "th")
  
  df.lines.ver$colour <- factor(ifelse(1:nrow(df.lines.ver)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  df.lines.hor$colour <- factor(ifelse(1:nrow(df.lines.hor)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  
  
  b <- b +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = as.character(colour)
      ),
      #colour = "red",
      linetype = "dashed",
      size = 1,
      data = df.lines.hor
    ) +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = colour
      ),
      #colour = "red",
      linetype = "dashed",
      size = 1,
      data = df.lines.ver
      
    ) 
  
  
  df <- data.frame()
  c <- ggplot() + theme_void()
  
  
  d<- ggarrange(c,
                b,
                ncol = 2,
                nrow = 1,
                widths = c(0.015, 1))
  
  
  d
}
# 
# new.sample.meta.score <- c(0.5,0.3)
# indexRow <- 0.3
# 
#survivalcurveplot(new.sample.meta.score,indexRow)


