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

load(file = "/home/njh264/Idat-Shiny/atrt.meth.os.meta.n8.extract.v3.abs.chun.Rdata")
load(file = "/home/njh264/Idat-Shiny/ecrt.meth.os.meta.n20.extract.v3.abs.chun.Rdata")
load(file = "/home/njh264/Idat-Shiny/all.meth.os.meta.n54.extract.v3.abs.chun.Rdata")

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
  as.numeric(scale(raw.metagene, center = scaling[1], scale = scaling[2])) -> scaled.metagene
  names(scaled.metagene) <- colnames(exp.matrix)
  return(scaled.metagene)
}


