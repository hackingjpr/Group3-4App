if(!require(minfiData)) {
  BiocManager::install("minfiData")
  library(minfiData)
}

if(!require(sva)){
  BiocManager::install("sva")
  library(sva)
}

if(!require(devtools)){
  BiocManager::install("devtools")
  library(devtools)
}

if(!require(bumphunter)){
  BiocManager::install("bumphunter")
  library(bumphunter)
}

if(!require(lumi)){
  BiocManager::install("lumi")
  library(lumi)
}

if(!require(minfi)){
  BiocManager::install("minfi")
  library(minfi)
}

if(!require(IlluminaHumanMethylationEPICmanifest)){
  BiocManager::install("IlluminaHumanMethylationEPICmanifest")
  library(IlluminaHumanMethylationEPICmanifest)
}

if(!require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)){
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
}

if(!require(shiny)){
  BiocManager::install("shiny")
  library(shiny)
}

if(!require(DT)){
  BiocManager::install("DT")
  library(DT)
}

if(!require(shinyWidgets)){
  BiocManager::install("shinyWidgets")
  library(shinyWidgets)
}

if(!require(shinydashboard)){
  BiocManager::install("shinydashboard")
  library(shinydashboard)
}

if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}

if(!require(ggpubr)){
  install.packages("ggpubr")
  library(ggpubr)
}

if(!require(ggnewscale)){
  install.packages("ggnewscale")
  library(ggnewscale)
}

if(!require(foreach)){
  install.packages("foreach")
  library(foreach)
}

if(!require(ggrepel)){
  install.packages("ggrepel")
  library(ggrepel)
}

if(!require(pheatmap)){
  install.packages("pheatmap")
  library(pheatmap)
}



load(file = "./ATRT.v3.abs.chun.Rdata")
atrt.meth.os.meta.n8.extract -> ATRT
load(file = "./ECRT.v3.abs.chun.Rdata")
ecrt.meth.os.meta.n20.extract -> ECRT
load(file = "./ALL.v3.abs.chun.Rdata")
all.meth.os.meta.n54.extract -> ALL

process_idats <- function(basenames){
  
  #check file sizes, epic will be over 1e7
  # idx.epic <-
  #   file.info(paste0(basenames, "_Red.idat"))$size > 1e7
  
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
  as.numeric(scale(raw.metagene, center = scaling[1], scale = scaling[2])) -> scaled.metagene
  round(scaled.metagene, digits = 3)
  names(scaled.metagene) <- colnames(exp.matrix)
  df <- as.data.frame(scaled.metagene)
  return(df)
}

digits = 0:9
createRandString<- function() {
  v = c(sample(LETTERS, 5, replace = TRUE),
        sample(digits, 4, replace = TRUE),
        sample(LETTERS, 1, replace = TRUE))
  return(paste0(v,collapse = ""))
}


generate_figure <- function(new.sample.meta.score){
  temp.df <- readRDS(file = "~/Idat-Shiny/temp.df.rds")
  df.cat.atrt <- readRDS(file = "~/Idat-Shiny/df.cat.atrt.rds")
  comb.SDb.atrt <- readRDS(file = "~/Idat-Shiny/comb.SDb.atrt.rds")
  ggplot(aes(x=1:nrow(temp.df),y=atrt8), data = temp.df) +
    geom_line() +
    scale_shape_manual(values=c(1, 4,  3)) +
    scale_color_manual(values=c('#E69F00','#999999',"white")) +
    ylab("ATRT-8") +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) -> b
  
  df.lines.hor <-foreach(i = 1:length(new.sample.meta.score), .combine = rbind)%do%{
    data.frame(x=0,
               xend=max(which(temp.df$atrt8 <new.sample.meta.score[i])),
               y=new.sample.meta.score[i],
               yend=new.sample.meta.score[i])}
  df.lines.hor$labels <- names(new.sample.meta.score)
  
  df.lines.ver <-foreach(i = 1:length(new.sample.meta.score), .combine = rbind)%do%{
    data.frame(x=max(which(temp.df$atrt8 < new.sample.meta.score[i])),
               xend=max(which(temp.df$atrt8 < new.sample.meta.score[i])),
               y=new.sample.meta.score[i],
               yend=min(temp.df$atrt8))}
  df.lines.ver$perc <- paste0(round(df.lines.ver$xend/ length(temp.df$atrt8)*100),"th")
  b<- b + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), colour = "red", linetype = "dashed", data = df.lines.hor) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), colour = "red", linetype = "dashed", data = df.lines.ver) +
    #geom_text_repel(aes(x = x+10, y = y+0.1, label = labels), direction = "y", data = df.lines.hor)
    geom_text(aes(x = x+10, y = y+0.1, label = labels), data = df.lines.hor) 
  
  df <- data.frame()
  c <- ggplot() + theme_void()
  
  
  #ggarrange(a,a2,ggarrange(
   # c,b, ncol = 2, nrow = 1, widths = c(0.015,1)),ncol=1,nrow=3)
  ggarrange(
    c,b, ncol = 2, nrow = 1, widths = c(0.015,1))
  #return(data.frame(perc = df.lines.ver[,5], row.names=rownames(df.lines.ver)))
}




