################# Packages ################# 
print(getwd())

if (!require(bumphunter)) {
  BiocManager::install("bumphunter")
  library(bumphunter)
}

if(!require(minfi)){
  BiocManager::install("minfi")
  library(minfi)
}

if (!require(shiny)) {
  install.packages("shiny")
  library(shiny)
}

if (!require(DT)) {
  install.packages("DT")
  library(DT)
}


if (!require(shinyWidgets)) {
  install.packages("shinyWidgets")
  library(shinyWidgets)
}

if (!require(shinydashboard)) {
  install.packages("shinydashboard")
  library(shinydashboard)
}

if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require(ggpubr)) {
  install.packages("ggpubr")
  library(ggpubr)
}

if (!require(foreach)) {
  install.packages("foreach")
  library(foreach)
}

if (!require(waiter)) {
  install.packages("waiter")
  library(waiter)
}

if (!require(caret)) {
  install.packages('caret')
  library(caret)
}
if (!require(NMF)) {
  install.packages('NMF')
  library(NMF)
}

if (!require(randomForest)) {
  install.packages("https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-14.tar.gz", repos = NULL, type = "source")
  library(randomForest)
}
if (!require(biomaRt)) {
  BiocManager::install("biomaRt")
  library(biomaRt)
}
if (!require(gridExtra)) {
  install.packages('gridExtra')
  library(gridExtra)
}

if (!require(survival)) {
  install.packages('survival')
  library(survival)
}

if (!require(tools)) {
  install.packages('tools')
  library(tools)
}

if (!require(markdown)) {
  install.packages('markdown')
  library(markdown)
}
################# Methylation ################# 


beta2m <- function (beta) {
  m <- log2(beta / (1 - beta))
  return(m)
}

process_idats <- function(basenames) {
  #check file sizes, epic will be over 1e7
  # Split into mixed or single test e.g. 450k & Epic or one or the other
  file.sizes <- file.info(paste0(basenames, "_Red.idat"))$size
  
  length(which(file.sizes > 1e7)) > 0 &
    length(which(file.sizes < 1e7)) > 0 -> mixed.test
  length(which(file.sizes > 1e7)) > 0 |
    length(which(file.sizes < 1e7)) > 0 -> single.test
  
  if (mixed.test) {
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
  
  if ((!mixed.test) & single.test) {
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

get_basenames <- function(dir) {
  list.files(path = dir,
             pattern = "_Red.idat",
             full.names = T) -> temp.files
  return(gsub("_Red.idat", "", temp.files))
}

extract.metagene <- function(index, weights, exp.matrix, scaling) {
  exp.matrix[index,] -> temp.exp
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
createRandString <- function() {
  v = c(
    sample(LETTERS, 5, replace = TRUE),
    sample(digits, 4, replace = TRUE),
    sample(LETTERS, 1, replace = TRUE)
  )
  return(paste0(v, collapse = ""))
}

load("./AppExtraFiles/Inputs/g3.g4.cont.rfe.Rdata")

################# Expression ################# 


# load in original metagene H values
avg.h.val <- readRDS("./AppExtraFiles/Inputs/avg.h.val.input.rds")

# load in original g3g4 metagene values
g3g4 <- readRDS("./AppExtraFiles/Inputs/g3g4.input.rds")

#load in filtered expression values from NMB MB samples used to generate H Values
#nmb.mat <- readRDS("./nmb.mat.input.rds")


nmb.mat.prepped <- readRDS(file = "./AppExtraFiles/Inputs/nmb.mat.prepped.rds")
nmf.res <- readRDS("./AppExtraFiles/Inputs/nmf.res.rds")

#tpms.mat <- readRDS(file = "./tpms.mat.rds")

### annotate using library(biomaRt)
annotate.HTseq.IDs <- function(HTseq.IDs) {
  ENSEMBL_DB_HOST = "uswest.ensembl.org" # Set back to default, once they are up and running again
  ENSEMBL_VERSION = "Ensembl Genes 101"  # Try to fix https://support.bioconductor.org/p/104454/
  
  #mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = ENSEMBL_DB_HOST, version = ENSEMBL_VERSION)
  mart <-
    useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  ensemblID <- gsub("\\..*", '', HTseq.IDs)
  #fetch the annotation information from the Ensembl database
  #ensembl <- useMart('ensembl', dataset='hsapiens_gene_ensembl')
  symbols <-
    getBM(
      attributes = c('ensembl_gene_id', 'hgnc_symbol'),
      filters = 'ensembl_gene_id',
      ensemblID,
      mart = mart
    )
  #combine the annotation with the RNA-Seq data
  annotatedix <- match(ensemblID, symbols$ensembl_gene_id)
  symbols[annotatedix, ] -> annotatedGenes
  return(cbind(HTseq.IDs, annotatedGenes))
}


### loading in all functions needed for projection
match.select <- function(input1.eset, input2.eset) {
  m1 <- input1.eset
  gs.names1 <- row.names(m1)
  gs.descs1 <- row.names(m1)
  sample.names1 <- names(m1)
  
  m2 <- input2.eset
  gs.names2 <- row.names(m2)
  gs.descs2 <- row.names(m2)
  sample.names2 <- names(m2)
  
  gs.names3 <- intersect(gs.names1, gs.names2)
  
  locations2 <- match(gs.names3, gs.names2, nomatch = 0)
  gs.names2 <- gs.names2[locations2]
  gs.descs2 <- gs.descs2[locations2]
  if (ncol(m2) == 1) {
    m2.out <- data.frame(m2[locations2,])
    rownames(m2.out) <- gs.descs2
    colnames(m2.out) <- sample.names2
  } else{
    m2.out <- m2[locations2,]
  }
  return(m2.out)
}

prep.data <- function(input.array) {
  m <- input.array
  col.names <- colnames(input.array)
  gs.names <- row.names(input.array)
  gs.descs <- row.names(input.array)
  sample.names <- colnames(input.array)
  
  cols <- length(m[1, ])
  for (j in 1:cols) {
    # column rank normalization from 0 to N - 1
    m[, j] <- rank(m[, j], ties.method = "average") - 1
  }
  m <- 10000 * m / (length(m[, 1]) - 1)
  
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
                         delta = 300) {
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
  plot(
    sort(row.max, decreasing = TRUE),
    type = "l",
    main = "row max",
    ylab = "Value",
    xlab = "Index"
  )
  
  row.min <- apply(m, MARGIN = 1, FUN = min)
  
  flag <- array(dim = rows)
  flag <- (row.max / row.min >= fold) & (row.max - row.min >= delta)
  cat(length(which(flag)), "genes / probes retained", "\n")
  flag
}

scaling.function <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
# scaling.function2 <- function(x){(x-0.3608444)/(0.6065306-0.3608444)}
scaling.function3 <-
  function(x) {
    (x - 0.3953062) / (0.5964371 - 0.3953062)
  }




param.filter <- function(m = NULL,
                         thres = 20,
                         ceil = 10000,
                         shift = 1,
                         fold = 3,
                         delta = 300) {
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
  plot(
    sort(row.max, decreasing = TRUE),
    type = "l",
    main = "row max",
    ylab = "Value",
    xlab = "Index"
  )
  
  row.min <- apply(m, MARGIN = 1, FUN = min)
  
  flag <- array(dim = rows)
  flag <- (row.max / row.min >= fold) & (row.max - row.min >= delta)
  cat(length(which(flag)), "genes / probes retained", "\n")
  flag
}

prep.data <- function(input.array) {
  m <- input.array
  col.names <- colnames(input.array)
  gs.names <- row.names(input.array)
  gs.descs <- row.names(input.array)
  sample.names <- colnames(input.array)
  
  cols <- length(m[1, ])
  for (j in 1:cols) {
    # column rank normalization from 0 to N - 1
    m[, j] <- rank(m[, j], ties.method = "average") - 1
  }
  m <- 10000 * m / (length(m[, 1]) - 1)
  
  v <- m
  row.names(v) <- gs.names
  colnames(v) <- col.names
  
  return(v)
}

project.NMF <- function(input.array, nmf.result) {
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
  
  locations.m <- match(overlap, gs.names, nomatch = 0)
  m2 <- m[locations.m,]
  locations.W <- match(overlap, W.row.names, nomatch = 0)
  W2 <- W[locations.W,]
  W2 <- as.matrix(W2)
  
  H <- ginv(W2) %*% m2
  
  n.col <- length(H[1, ])
  for (i in 1:n.col) {
    S.2 <- sqrt(sum(H[, i] * H[, i]))
    #        S.2 <- sum(H[,i])
    H[, i] <- H[, i] / S.2
  }
  
  V <- data.frame(H)
  names(V) <- sample.names
  row.names(V) <- W.names
  return(V)
}


annotate.HTseq.IDs <- function(HTseq.IDs) {
  ENSEMBL_DB_HOST = "uswest.ensembl.org" # Set back to default, once they are up and running again
  ENSEMBL_VERSION = "Ensembl Genes 101"  # Try to fix https://support.bioconductor.org/p/104454/
  
  #mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = ENSEMBL_DB_HOST, version = ENSEMBL_VERSION)
  mart <-
    useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  ensemblID <- gsub("\\..*", '', HTseq.IDs)
  #fetch the annotation information from the Ensembl database
  #ensembl <- useMart('ensembl', dataset='hsapiens_gene_ensembl')
  symbols <-
    getBM(
      attributes = c('ensembl_gene_id', 'hgnc_symbol'),
      filters = 'ensembl_gene_id',
      ensemblID,
      mart = mart
    )
  #combine the annotation with the RNA-Seq data
  annotatedix <- match(ensemblID, symbols$ensembl_gene_id)
  symbols[annotatedix, ] -> annotatedGenes
  return(cbind(HTseq.IDs, annotatedGenes))
}

################# Graphs ################# 

df1 <- data.frame(
  green = c(0, 0.68),
  yellow = c(0.68, 1),
  colour = c("G", "Y")
)


generate_figure_highlight_g3g4 <-
  function(new.sample.meta.score, indexRow) {
    if (is.null(indexRow)) {
      indexRow = 1
    }
    y <- readRDS(file = "./AppExtraFiles/Inputs/y.vals.rds")
    df <- as.data.frame(y)
    ecdf(y) -> model
    model(y) -> y2
    b <- ggplot(df, aes(x = y2, y = y)) +
      geom_rect(
        data = df1,
        aes(
          NULL,
          NULL,
          xmin = green,
          xmax = yellow,
          fill = colour
        ),
        ymin = 0,
        ymax = 1,
        colour = "white",
        size = 0.5,
        alpha = 0.2
      ) +
      scale_fill_manual(values = c("G" = "green", "Y" = "yellow")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_hline(yintercept = 0.5, linetype = "dashed") +
      geom_hline(yintercept = 1, linetype = "dashed") +
      geom_vline(xintercept = 0.68, linetype = "dashed") +
      geom_line() +
      theme(legend.position = "none") +
      theme(text = element_text(size = 15))
    
    
    
    
    
    df.lines.hor <-
      foreach(i = 1:length(new.sample.meta.score),
              .combine = rbind) %do% {
                data.frame(
                  # x = 0,
                  x = 0,
                  xend = y2[max(which(y < new.sample.meta.score[i]))],
                  y = y[max(which(y < new.sample.meta.score[i]))],
                  yend = y[max(which(y < new.sample.meta.score[i]))]
                )
              }
    df.lines.hor$labels <- names(new.sample.meta.score)
    
    df.lines.ver <-
      foreach(i = 1:length(new.sample.meta.score),
              .combine = rbind) %do% {
                data.frame(
                  x = y2[max(which(y < new.sample.meta.score[i]))],
                  xend = y2[max(which(y < new.sample.meta.score[i]))],
                  y =  y[max(which(y < new.sample.meta.score[i]))],
                  yend = 0
                )
              }
    df.lines.ver$perc <-
      paste0(round(df.lines.ver$xend / length(y2) * 100), "th")
    
    df.lines.ver$colour <-
      factor(
        ifelse(
          1:nrow(df.lines.ver) == indexRow,
          "highlight",
          "no.highlight"
        ),
        levels = c("highlight", "no.highlight")
      )
    df.lines.hor$colour <-
      factor(
        ifelse(
          1:nrow(df.lines.hor) == indexRow,
          "highlight",
          "no.highlight"
        ),
        levels = c("highlight", "no.highlight")
      )
    
    
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
        size = 0.5,
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
        size = 0.5,
        data = df.lines.ver
        
      ) +
      xlab("") +
      ylab("Group 3/4 Score")
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
    d <- ggarrange(c,
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

generate_figure_highlight_g3g4Expression <-
  function(new.sample.meta.score, indexRow) {
    if (is.null(indexRow)) {
      indexRow = 1
    }
    y <- readRDS(file = "./AppExtraFiles/Inputs/y.vals.rds")
    df <- as.data.frame(y)
    ecdf(y) -> model
    model(y) -> y2
    b <- ggplot(df, aes(x = y2, y = y)) +
      geom_rect(
        data = df1,
        aes(
          NULL,
          NULL,
          xmin = green,
          xmax = yellow,
          fill = colour
        ),
        ymin = 0,
        ymax = 1,
        colour = "white",
        size = 0.5,
        alpha = 0.2
      ) +
      scale_fill_manual(values = c("G" = "green", "Y" = "yellow")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_hline(yintercept = 0.5, linetype = "dashed") +
      geom_hline(yintercept = 1, linetype = "dashed") +
      geom_vline(xintercept = 0.68, linetype = "dashed") +
      geom_line() +
      theme(legend.position = "none") +
      theme(text = element_text(size = 15))
    
    df.lines.hor <-
      foreach(i = 1:length(new.sample.meta.score),
              .combine = rbind) %do% {
                data.frame(
                  # x = 0,
                  x = 0,
                  xend = y2[max(which(y < new.sample.meta.score[i]))],
                  y = y[max(which(y < new.sample.meta.score[i]))],
                  yend = y[max(which(y < new.sample.meta.score[i]))]
                )
              }
    df.lines.hor$labels <- names(new.sample.meta.score)
    
    df.lines.ver <-
      foreach(i = 1:length(new.sample.meta.score),
              .combine = rbind) %do% {
                data.frame(
                  x = y2[max(which(y < new.sample.meta.score[i]))],
                  xend = y2[max(which(y < new.sample.meta.score[i]))],
                  y =  y[max(which(y < new.sample.meta.score[i]))],
                  yend = 0
                )
              }
    df.lines.ver$perc <-
      paste0(round(df.lines.ver$xend / length(y2) * 100), "th")
    
    df.lines.ver$colour <-
      factor(
        ifelse(
          1:nrow(df.lines.ver) == indexRow,
          "highlight",
          "no.highlight"
        ),
        levels = c("highlight", "no.highlight")
      )
    df.lines.hor$colour <-
      factor(
        ifelse(
          1:nrow(df.lines.hor) == indexRow,
          "highlight",
          "no.highlight"
        ),
        levels = c("highlight", "no.highlight")
      )
    
    
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
        size = 0.5,
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
        size = 0.5,
        data = df.lines.ver
        
      ) +
      xlab("") +
      ylab("Group 3/4 Score")
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
    d <- ggarrange(c,
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


generate_figure_highlight_g3g4PERC <-
  function(new.sample.meta.score, indexRow) {
    if (is.null(indexRow)) {
      indexRow = 1
    }
    y <- readRDS(file = "./AppExtraFiles/Inputs/y.vals.rds")
    df <- as.data.frame(y)
    ecdf(y) -> model
    model(y) -> y2
    b <- ggplot(df, aes(x = y2, y = y)) +
      geom_rect(
        data = df1,
        aes(
          NULL,
          NULL,
          xmin = green,
          xmax = yellow,
          fill = colour
        ),
        ymin = 0,
        ymax = 1,
        colour = "white",
        size = 0.5,
        alpha = 0.2
      ) +
      scale_fill_manual(values = c("G" = "green", "Y" = "yellow")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_hline(yintercept = 0.5, linetype = "dashed") +
      geom_hline(yintercept = 1, linetype = "dashed") +
      geom_vline(xintercept = 0.68, linetype = "dashed") +
      geom_line() +
      theme(legend.position = "none") +
      theme(text = element_text(size = 15))
    
    
    
    
    
    df.lines.hor <-
      foreach(i = 1:length(new.sample.meta.score),
              .combine = rbind) %do% {
                data.frame(
                  # x = 0,
                  x = 0,
                  xend = y2[max(which(y < new.sample.meta.score[i]))],
                  y = y[max(which(y < new.sample.meta.score[i]))],
                  yend = y[max(which(y < new.sample.meta.score[i]))]
                )
              }
    df.lines.hor$labels <- names(new.sample.meta.score)
    
    df.lines.ver <-
      foreach(i = 1:length(new.sample.meta.score),
              .combine = rbind) %do% {
                data.frame(
                  x = y2[max(which(y < new.sample.meta.score[i]))],
                  xend = y2[max(which(y < new.sample.meta.score[i]))],
                  y =  y[max(which(y < new.sample.meta.score[i]))],
                  yend = 0
                )
              }
    df.lines.ver$perc <-
      # paste0(round(df.lines.hor$xend / length(y2) * 100), "th")
      paste0(round(df.lines.hor$xend / 1 * 100), "%")
    
    return(df.lines.ver$perc[indexRow])
  }


################# FOR AGE BASED SURVIVAL PLOT #################
#
time.comb <- readRDS(file = "./AppExtraFiles/Inputs/time.comb.rds")
age.comb <- readRDS(file = "./AppExtraFiles/Inputs/age.comb.rds")
status.comb <- readRDS(file = "./AppExtraFiles/Inputs/status.comb.rds")
comb.cont <- readRDS(file = "./AppExtraFiles/Inputs/comb.cont.rds")


# new.sample.meta.score <- c(0.5,0.3)
# indexRow <- 0.3

coxph(Surv(time.comb, status.comb) ~ comb.cont) -> train.fit
# message("Train fit complete")

summary(survfit(train.fit, data.frame(g3g4.values = as.numeric(comb.cont))), time = 5) -> x

df2 <- data.frame(
  pred = comb.cont,
  surv = as.numeric(x$surv),
  up = as.numeric(x$upper),
  lo = as.numeric(x$lower)
)

# df2$pred -> x
# df2$surv  -> y
# data = data.frame(x, y)
# ggplot(aes(x = x, y = y), data = data) +
#   geom_line()
# 

fit.extrap <- lm(surv~poly(pred,2,raw=TRUE), data=df2)
fit.extrap.up <- lm(up~poly(pred,2,raw=TRUE), data=df2)
fit.extrap.lo <- lm(lo~poly(pred,2,raw=TRUE), data=df2)
extrap.surv <-  predict(fit.extrap, data.frame(pred = c(0,0.01,0.02,0.03,0.98,0.99,1)))
extrap.pred <- c(0,0.01,0.02,0.03,0.98,0.99,1)
extrap.hi <-  predict(fit.extrap.up, data.frame(pred = c(0,0.01,0.02,0.03,0.98,0.99,1)))
extrap.lo <- predict(fit.extrap.lo, data.frame(pred = c(0,0.01,0.02,0.03,0.98,0.99,1)))

data.frame(pred = extrap.pred,
           surv = extrap.surv,
           up = extrap.hi,
           lo = extrap.lo) -> df2.extrap
# #define x-axis values
# x_axis <- seq(0, 1, length=15)
# plot(data$x, data$y, pch=19, xlab='x', ylab='y')
# #add curve of each model to plot
# lines(x_axis, predict(fit2, data.frame(x=x_axis)), col='green')


# summary(fit)
#Formula: y ~ SSlogis(x, Asym, xmid, scal)
#
#Parameters:
#      Estimate Std. Error t value Pr(>|t|)
#Asym 1.473e+04  2.309e+04   0.638    0.551
#xmid 4.094e+00  2.739e+00   1.495    0.195
#scal 9.487e-01  5.851e-01   1.622    0.166
#
#Residual standard error: 941.9 on 5 degrees of freedom
#
#Number of iterations to convergence: 0 
#Achieved convergence tolerance: 4.928e-06

# lines(seq(0.5, 4, length.out = 100), 
#       predict(fit, newdata = data.frame(x = seq(0.5, 4, length.out = 100))))


coxph(Surv(time.comb, c(status.comb)) ~ comb.cont + age.comb) -> train.fit.age

summary(survfit(
  train.fit.age,
  data.frame(comb.cont = comb.cont,
             age.comb = age.comb)
), time = 5) -> x

summary(survfit(train.fit.age, data.frame(
  comb.cont = c(seq(0, 1, 0.01), seq(0, 1, 0.01))
  ,
  age.comb = c(rep("TRUE", 101), rep("FALSE", 101))
)), time = 5) -> y


# fit.extrap3 <- lm(surv~poly(pred,2,raw=TRUE), data=df3)
# fit.extrap.up3 <- lm(up~poly(pred,2,raw=TRUE), data=df3)
# fit.extrap.lo3 <- lm(lo~poly(pred,2,raw=TRUE), data=df3)
# extrap.surv3 <-  predict(fit.extrap3, data.frame(pred = c(0,0.01,0.02,0.03,0.98,0.99,1)))
# extrap.pred3 <- c(0,0.01,0.02,0.03,0.98,0.99,1)
# extrap.hi3 <-  predict(fit.extrap.up3, data.frame(pred = c(0,0.01,0.02,0.03,0.98,0.99,1)))
# extrap.lo3 <- predict(fit.extrap.lo3, data.frame(pred = c(0,0.01,0.02,0.03,0.98,0.99,1)))
# 
# data.frame(pred = extrap.pred3,
#            surv = extrap.surv3,
#            up = extrap.hi3,
#            lo = extrap.lo3,
#            age = age.comb[row.names(x$table)]) -> df3.extrap

df3 <- data.frame(
  pred = comb.cont[row.names(x$table)],
  surv = as.numeric(x$surv),
  up = as.numeric(x$upper),
  lo = as.numeric(x$lower),
  age = age.comb[row.names(x$table)]
)


df3.y <- data.frame(
  pred = c(seq(0, 1, 0.01), seq(0, 1, 0.01)),
  surv = as.numeric(y$surv),
  up = as.numeric(y$upper),
  lo = as.numeric(y$lower),
  age = c(rep("TRUE", 101), rep("FALSE", 101))
)

# df3$point <- rep("yes", nrow(df3))
# df3.y$point <- rep("no", nrow(df3.y))



################# SURVIVAL PLOT AND PERCENTAGE ################# 



survivalcurveplot <- function(new.sample.meta.score, indexRow) {
  if (is.null(indexRow)) {
    indexRow = 1
  }
  
  df2 <- rbind(df2, df2.extrap)
  df2$pred -> pred
  df2$surv -> surv
  b <- ggplot(df2, aes(x = pred, y = surv)) +
    geom_line() +
    geom_point(alpha = 1 / 20) +
    geom_line(aes(x = pred, y = lo), linetype = "dotted") +
    geom_line(aes(x = pred, y = up), linetype = "dotted") +
    theme_classic() + xlab("Group 3/4 Score") + ylab("Survival") +
    # labs(title = "New plot title", subtitle = "A subtitle") +
    ylim(0, 1) +
    xlim(0, 1) +
    theme(legend.position = "none") +
    theme(text = element_text(size = 15))
  
  df.lines.hor <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              if(length(which(pred < new.sample.meta.score[i]))==0){
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = 0.8,
                  yend = 0.8
                )  
              }else{
                surv[which(pred < new.sample.meta.score[i])] -> temp.surv
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = temp.surv[which.min(temp.surv)],
                  yend = temp.surv[which.min(temp.surv)]
                )
              }
            }
  
  df.lines.hor$labels <- names(new.sample.meta.score)
  
  
  df.lines.ver <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              if(length(which(pred < new.sample.meta.score[i]))==0){
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = 0,
                  yend = 0.8
                )  
              }else{
                surv[which(pred < new.sample.meta.score[i])] -> temp.surv
                data.frame(
                  x = new.sample.meta.score[i],
                  xend = new.sample.meta.score[i],
                  y = 0,
                  yend = temp.surv[which.min(temp.surv)]
                )
              }
            }
  message(df.lines.ver)
  
  
  df.lines.ver$perc <-
    paste0(round(df.lines.ver$xend / length(pred) * 100), "th")
  
  df.lines.ver$colour <-
    factor(
      ifelse(
        1:nrow(df.lines.ver) == indexRow,
        "highlight",
        "no.highlight"
      ),
      levels = c("highlight", "no.highlight")
    )
  df.lines.hor$colour <-
    factor(
      ifelse(
        1:nrow(df.lines.hor) == indexRow,
        "highlight",
        "no.highlight"
      ),
      levels = c("highlight", "no.highlight")
    )
  
  
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
      size = 0.5,
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
      size = 0.5,
      data = df.lines.ver
      
    )
  
  
  df <- data.frame()
  c <- ggplot() + theme_void()
  
  
  d <- ggarrange(c,
                 b,
                 ncol = 2,
                 nrow = 1,
                 widths = c(0.015, 1))
  
  
  d
}





survivalcurveplotPERC <- function(new.sample.meta.score, indexRow) {
  
  
  if (is.null(indexRow)) {
    indexRow = 1
  }
  df2 <- rbind(df2, df2.extrap)
  df2$pred -> pred
  df2$surv -> surv
  b <- ggplot(df2, aes(x = pred, y = surv)) +
    geom_line() +
    geom_point(alpha = 1 / 20) +
    geom_line(aes(x = pred, y = lo), linetype = "dotted") +
    geom_line(aes(x = pred, y = up), linetype = "dotted") +
    # theme_classic() +
    xlab("Group 3/4 Score") + ylab("Survival") +
    # labs(title = "New plot title", subtitle = "A subtitle") +
    ylim(0, 1) +
    theme(legend.position = "none") +
    theme(text = element_text(size = 15))
  
  
  df.lines.hor <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              if(length(which(pred < new.sample.meta.score[i]))==0){
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = 0.8,
                  yend = 0.8
                )  
              }else{
                surv[which(pred < new.sample.meta.score[i])] -> temp.surv
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = temp.surv[which.min(temp.surv)],
                  yend = temp.surv[which.min(temp.surv)]
                )
              }
            }
  
  df.lines.hor$labels <- names(new.sample.meta.score)
  
  
  df.lines.ver <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              if(length(which(pred < new.sample.meta.score[i]))==0){
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = 0,
                  yend = 0.8
                )  
              }else{
                surv[which(pred < new.sample.meta.score[i])] -> temp.surv
                data.frame(
                  x = new.sample.meta.score[i],
                  xend = new.sample.meta.score[i],
                  y = 0,
                  yend = temp.surv[which.min(temp.surv)]
                )
              }
            }
  message(df.lines.ver)
  
  df.lines.hor$perc <-
    paste0(round(df.lines.ver$yend* 100), "%")
  
  return(df.lines.hor$perc[indexRow])
  
  
}
#
# new.sample.meta.score <- c(0.5,0.3)
# indexRow <- 0.3
#

################# AGE PLOT AND PERCENTAGES ################# 
new.sample.meta.score <- c(0.01, 0, 0.98, 1)
indexRow <- 2



SurvivalAgePlot <- function (new.sample.meta.score, indexRow) {
  if (is.null(indexRow)) {
    indexRow = 1
  }
  
  # df3 <- rbind(df3, df3.extrap)
  
  
  age <- df3$age
  age <- ('True' = 'Over Five')
  df3$pred -> pred
  df3$surv -> surv
  df4 <- subset(df3.y, age == TRUE)
  df5 <- subset(df3.y, age != TRUE)
  df4$pred -> pred4
  df4$surv -> surv4
  df5$pred -> pred5
  df5$surv -> surv5
  b <- ggplot(df3, aes(
    x = pred,
    y = surv,
    group = age,
    color = age
  )) +
    #geom_line() +
    geom_point(alpha = 1 / 10) +
    geom_line(data = df3.y, aes(
      x = pred,
      y = surv,
      group = age,
      color = age
    )) +
    geom_line(data = df3.y,
              aes(x = pred, y = lo, group = age),
              linetype = "dotted") +
    geom_line(data = df3.y,
              aes(x = pred, y = up, group = age),
              linetype = "dotted") +
    theme_classic() + xlab("Group 3/4 Score") + ylab("Survival") +
    scale_color_manual(
      name = 'Age',
      labels = c("Over Three", "Under Three", "no.highlight", "highlight"),
      values = c('red', 'dodgerblue', "grey", "orange")
    ) +
    # scale_color_manual(name='Age',
    #                    labels=c("Over Five", "Under Five", "no.highlight", "highlight"),
    #                    values=c('red','dodgerblue', "grey", "orange")) +
    #
    theme(legend.position = "none") +
    # labs(title = "New plot title", subtitle = "A subtitle") +
    ylim(0, 1) +
    theme(text = element_text(size = 15))
  
  df.lines.hor4 <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              surv4[which(pred4 < new.sample.meta.score[i])] -> temp.surv
              if(length(which(pred4 < new.sample.meta.score[i]))==0){
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = 0.8,
                  yend = 0.8
                )  
              }else{
                surv4[which(pred4 < new.sample.meta.score[i])] -> temp.surv
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = temp.surv[which.min(temp.surv)],
                  yend = temp.surv[which.min(temp.surv)]
                )
              }
            }
  
  df.lines.hor4$labels <- names(new.sample.meta.score)
  
  
  # df.lines.ver4 <-
  #   foreach(i = 1:length(new.sample.meta.score),
  #           .combine = rbind) %do% {
  #             surv4[which(
  #               pred4 < new.sample.meta.score[i]
  #             )] -> temp.surv
  #             data.frame(
  #               x = new.sample.meta.score[i],
  #               xend = new.sample.meta.score[i],
  #               y = 0,
  #               yend = temp.surv[which.min(temp.surv)]
  #             )
  #           }
  # message(df.lines.ver4)
  
  df.lines.hor5 <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              surv4[which(pred5 < new.sample.meta.score[i])] -> temp.surv
              if(length(which(pred4 < new.sample.meta.score[i]))==0){
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = 0.8,
                  yend = 0.8
                )  
              }else{
                surv5[which(pred5 < new.sample.meta.score[i])] -> temp.surv
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = temp.surv[which.min(temp.surv)],
                  yend = temp.surv[which.min(temp.surv)]
                )
              }
            }
  
  df.lines.hor5$labels <- names(new.sample.meta.score)
  
  
  df.lines.ver5 <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              if(length(which(pred5 < new.sample.meta.score[i]))==0){
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = 0,
                  yend = 0.8
                )  
              }else{
                surv5[which(pred5 < new.sample.meta.score[i])] -> temp.surv
                data.frame(
                  x = new.sample.meta.score[i],
                  xend = new.sample.meta.score[i],
                  y = 0,
                  yend = temp.surv[which.min(temp.surv)]
                )
              }
            }
  
  
  # df.lines.ver4$perc <-
  #   paste0(round(df.lines.ver4$xend / length(pred4) * 100), "th")
  
  df.lines.ver5$perc <-
    paste0(round(df.lines.ver5$xend / length(pred5) * 100), "th")
  
  # df.lines.ver4$cols <- factor(ifelse(1:nrow(df.lines.ver4)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  df.lines.hor4$cols <-
    factor(
      ifelse(
        1:nrow(df.lines.hor4) == indexRow,
        "highlight",
        "no.highlight"
      ),
      levels = c("highlight", "no.highlight")
    )
  df.lines.ver5$cols <-
    factor(
      ifelse(
        1:nrow(df.lines.ver5) == indexRow,
        "highlight",
        "no.highlight"
      ),
      levels = c("highlight", "no.highlight")
    )
  df.lines.hor5$cols <-
    factor(
      ifelse(
        1:nrow(df.lines.hor5) == indexRow,
        "highlight",
        "no.highlight"
      ),
      levels = c("highlight", "no.highlight")
    )
  
  
  b <- b +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = as.character(cols)
      ),
      #colour = "red",
      linetype = "dashed",
      size = 0.5,
      data = df.lines.hor4
    ) +
    # geom_segment(
    #   aes(
    #     x = x,
    #     y = y,
    #     xend = xend,
    #     yend = yend,
    #     colour = as.character(cols)
    #   ),
    #   #colour = "red",
    #   linetype = "dashed",
    #   size = 1,
  #   data = df.lines.ver4
  # ) +
  geom_segment(
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      colour = as.character(cols)
    ),
    #colour = "red",
    linetype = "dashed",
    size = 0.5,
    data = df.lines.hor5
  ) +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = as.character(cols)
      ),
      #colour = "red",
      linetype = "dashed",
      size = 0.5,
      data = df.lines.ver5
    )
  
  
  df <- data.frame()
  c <- ggplot() + theme_void()
  
  
  d <- ggarrange(c,
                 b,
                 ncol = 2,
                 nrow = 1,
                 widths = c(0.015, 1))
  
  
  d
  # b
}



SurvivalAgePlotPerc5 <- function (new.sample.meta.score, indexRow) {
  if (is.null(indexRow)) {
    indexRow = 1
  }
  
  # df3 <- rbind(df3, df3.extrap)
  
  
  age <- df3$age
  age <- ('True' = 'Over Five')
  df3$pred -> pred
  df3$surv -> surv
  df4 <- subset(df3.y, age == TRUE)
  df5 <- subset(df3.y, age != TRUE)
  df4$pred -> pred4
  df4$surv -> surv4
  df5$pred -> pred5
  df5$surv -> surv5
  b <- ggplot(df3, aes(
    x = pred,
    y = surv,
    group = age,
    color = age
  )) +
    #geom_line() +
    geom_point(alpha = 1 / 10) +
    geom_line(data = df3.y, aes(
      x = pred,
      y = surv,
      group = age,
      color = age
    )) +
    geom_line(data = df3.y,
              aes(x = pred, y = lo, group = age),
              linetype = "dotted") +
    geom_line(data = df3.y,
              aes(x = pred, y = up, group = age),
              linetype = "dotted") +
    theme_classic() + xlab("Group 3/4 Score") + ylab("Survival") +
    scale_color_manual(
      name = 'Age',
      labels = c("Over Three", "Under Three", "no.highlight", "highlight"),
      values = c('red', 'dodgerblue', "grey", "orange")
    ) +
    # scale_color_manual(name='Age',
    #                    labels=c("Over Five", "Under Five", "no.highlight", "highlight"),
    #                    values=c('red','dodgerblue', "grey", "orange")) +
    #
    # theme(legend.position = "none") +
    # labs(title = "New plot title", subtitle = "A subtitle") +
    ylim(0, 1) +
    theme(text = element_text(size = 15))
  
  df.lines.hor4 <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              surv4[which(pred4 < new.sample.meta.score[i])] -> temp.surv
              if(length(which(pred4 < new.sample.meta.score[i]))==0){
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = 0.8,
                  yend = 0.8
                )  
              }else{
                surv4[which(pred4 < new.sample.meta.score[i])] -> temp.surv
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = temp.surv[which.min(temp.surv)],
                  yend = temp.surv[which.min(temp.surv)]
                )
              }
            }
  
  df.lines.hor4$labels <- names(new.sample.meta.score)
  
  
  # df.lines.ver4 <-
  #   foreach(i = 1:length(new.sample.meta.score),
  #           .combine = rbind) %do% {
  #             surv4[which(
  #               pred4 < new.sample.meta.score[i]
  #             )] -> temp.surv
  #             data.frame(
  #               x = new.sample.meta.score[i],
  #               xend = new.sample.meta.score[i],
  #               y = 0,
  #               yend = temp.surv[which.min(temp.surv)]
  #             )
  #           }
  # message(df.lines.ver4)
  
  df.lines.hor5 <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              surv4[which(pred5 < new.sample.meta.score[i])] -> temp.surv
              if(length(which(pred4 < new.sample.meta.score[i]))==0){
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = 0.8,
                  yend = 0.8
                )  
              }else{
                surv5[which(pred5 < new.sample.meta.score[i])] -> temp.surv
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = temp.surv[which.min(temp.surv)],
                  yend = temp.surv[which.min(temp.surv)]
                )
              }
            }
  
  df.lines.hor5$labels <- names(new.sample.meta.score)
  
  
  df.lines.ver5 <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              if(length(which(pred5 < new.sample.meta.score[i]))==0){
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = 0,
                  yend = 0.8
                )  
              }else{
                surv5[which(pred5 < new.sample.meta.score[i])] -> temp.surv
                data.frame(
                  x = new.sample.meta.score[i],
                  xend = new.sample.meta.score[i],
                  y = 0,
                  yend = temp.surv[which.min(temp.surv)]
                )
              }
            }
  
  df.lines.hor5$perc <-
    paste0(round(df.lines.hor5$y* 100), "%")
  
  return(df.lines.hor5$perc[indexRow])
}



SurvivalAgePlotPerc4 <- function (new.sample.meta.score, indexRow) {
  if (is.null(indexRow)) {
    indexRow = 1
  }
  
  # df3 <- rbind(df3, df3.extrap)
  
  
  age <- df3$age
  age <- ('True' = 'Over Five')
  df3$pred -> pred
  df3$surv -> surv
  df4 <- subset(df3.y, age == TRUE)
  df5 <- subset(df3.y, age != TRUE)
  df4$pred -> pred4
  df4$surv -> surv4
  df5$pred -> pred5
  df5$surv -> surv5
  b <- ggplot(df3, aes(
    x = pred,
    y = surv,
    group = age,
    color = age
  )) +
    #geom_line() +
    geom_point(alpha = 1 / 10) +
    geom_line(data = df3.y, aes(
      x = pred,
      y = surv,
      group = age,
      color = age
    )) +
    geom_line(data = df3.y,
              aes(x = pred, y = lo, group = age),
              linetype = "dotted") +
    geom_line(data = df3.y,
              aes(x = pred, y = up, group = age),
              linetype = "dotted") +
    theme_classic() + xlab("Group 3/4 Score") + ylab("Survival") +
    # scale_color_manual(
    #   name = 'Age',
    #   labels = c("Over Three", "Under Three", "no.highlight", "highlight"),
    #   values = c('red', 'dodgerblue', "grey", "orange")
    # ) +
    # scale_color_manual(name='Age',
    #                    labels=c("Over Five", "Under Five", "no.highlight", "highlight"),
    #                    values=c('red','dodgerblue', "grey", "orange")) +
    #
    theme(legend.position = "none") +
    # labs(title = "New plot title", subtitle = "A subtitle") +
    ylim(0, 1) +
    theme(text = element_text(size = 15))
  
  df.lines.hor4 <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              surv4[which(pred4 < new.sample.meta.score[i])] -> temp.surv
              if(length(which(pred4 < new.sample.meta.score[i]))==0){
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = 0.8,
                  yend = 0.8
                )  
              }else{
                surv4[which(pred4 < new.sample.meta.score[i])] -> temp.surv
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = temp.surv[which.min(temp.surv)],
                  yend = temp.surv[which.min(temp.surv)]
                )
              }
            }
  
  df.lines.hor4$labels <- names(new.sample.meta.score)
  
  df.lines.hor5 <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              surv4[which(pred5 < new.sample.meta.score[i])] -> temp.surv
              if(length(which(pred4 < new.sample.meta.score[i]))==0){
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = 0.8,
                  yend = 0.8
                )  
              }else{
                surv5[which(pred5 < new.sample.meta.score[i])] -> temp.surv
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = temp.surv[which.min(temp.surv)],
                  yend = temp.surv[which.min(temp.surv)]
                )
              }
            }
  
  df.lines.hor5$labels <- names(new.sample.meta.score)
  
  
  df.lines.ver5 <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              if(length(which(pred5 < new.sample.meta.score[i]))==0){
                data.frame(
                  x = 0,
                  xend = new.sample.meta.score[i],
                  y = 0,
                  yend = 0.8
                )  
              }else{
                surv5[which(pred5 < new.sample.meta.score[i])] -> temp.surv
                data.frame(
                  x = new.sample.meta.score[i],
                  xend = new.sample.meta.score[i],
                  y = 0,
                  yend = temp.surv[which.min(temp.surv)]
                )
              }
            }
  
  
  df.lines.hor4$perc <-
    paste0(round(df.lines.hor4$y* 100), "%")
  
  return(df.lines.hor4$perc[indexRow])
}