


BiocManager::install("minfiData")

BiocManager::install("sva")

BiocManager::install("devtools")

BiocManager::install("bumphunter")

BiocManager::install("lumi")



library(bumphunter)
library(minfi)
library(minfiData)
library(sva)
library(lumi)



library(devtools)

install_github("hansenlab/tutorial.450k")





baseDir <- system.file("extdata", package="minfiData")

samp.df <- read.metharray.sheet(baseDir)



# turn idats into beta values



require(minfi)



## load samplesheet for NCL samples

samp.df <- read.delim(baseDir,
                      
                      stringsAsFactors = FALSE,
                      
                      sep = ",",
                      
                      head = TRUE,
                      
                      row.names = 1)



idx.epic <- file.info(paste0(samp.df$Basename, "_Red.idat"))$size > 1e7



m450k.rgSet <- read.metharray(samp.df$Basename[!idx.epic],
                              
                              verbose = TRUE,
                              
                              force = TRUE)

epic.rgSet <- read.metharray(samp.df$Basename[idx.epic],
                             
                             verbose = TRUE,
                             
                             force = TRUE)



comb.rgSet <- combineArrays(m450k.rgSet, epic.rgSet,
                            
                            outType = "IlluminaHumanMethylation450k")



comb.mSet <- preprocessNoob(m450k.rgSet, dyeMethod = "single")



# dir.create("/home/yuri/Rhabdoid_Analysis/metagene_survival/raw/")

# saveRDS(comb.rgSet,

#         "/home/yuri/Rhabdoid_Analysis/metagene_survival/raw/combined_450K_rgSet.rds")



comb.detP <- detectionP(m450k.rgSet)

comb.gmrSet <- mapToGenome(ratioConvert(comb.mSet))



comb.beta <- getBeta(comb.gmrSet)

comb.MRT.pd <- pData(comb.gmrSet)



#######



#premade beta values



comb.SDb.mrt <- read.csv("./data/comb.SDb.mrt.v3.chun.csv", row.names = 1)

comb.SDb.ecrt <- read.csv("./data/comb.SDb.ecrt.v3.chun.csv", row.names = 1)

comb.SDb.atrt <- read.csv("./data/comb.SDb.atrt.v3.chun.csv", row.names = 1)





comb.SDb.mrt



extract.metagene <- function(index, weights, exp.matrix, scaling){
  
  exp.matrix[index,] -> temp.exp
  
  apply(temp.exp,2,function(x){mean(x*weights)}) -> raw.metagene
  
  as.numeric(scale(raw.metagene, center = scaling[1], scale = scaling[2])) -> scaled.metagene
  
  names(scaled.metagene) <- colnames(exp.matrix)
  
  return(scaled.metagene)
  
}



load(file = "/home/cancer/mrt_survival/coxboost_extract/atrt.meth.os.meta.n8.extract.v3.abs.chun.Rdata")



extract.metagene(as.character(atrt.meth.os.meta.n8.extract[[6]]$genes),
                 
                 as.numeric(atrt.meth.os.meta.n8.extract[[6]]$weights),
                 
                 beta2m(comb.SDb.atrt),
                 
                 as.numeric(atrt.meth.os.meta.n8.extract[[7]])) -> atrt.meth.os.meta.n8.extract.meta





load(file = "/home/cancer/mrt_survival/coxboost_extract/ecrt.meth.os.meta.n20.extract.v3.abs.chun.Rdata")



extract.metagene(as.character(ecrt.meth.os.meta.n20.extract[[6]]$genes),
                 
                 as.numeric(ecrt.meth.os.meta.n20.extract[[6]]$weights),
                 
                 beta2m(comb.SDb.ecrt),
                 
                 as.numeric(ecrt.meth.os.meta.n20.extract[[7]])) -> ecrt.meth.os.meta.n20.extract.meta





load(file = "/home/cancer/mrt_survival/coxboost_extract/all.meth.os.meta.n54.extract.v3.abs.chun.Rdata")



extract.metagene(as.character(all.meth.os.meta.n54.extract[[6]]$genes),
                 
                 as.numeric(all.meth.os.meta.n54.extract[[6]]$weights),
                 
                 beta2m(comb.SDb.mrt),
                 
                 as.numeric(all.meth.os.meta.n54.extract[[7]])) -> all.meth.os.meta.n54.extract.meta