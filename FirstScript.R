
# samp.df <- read.delim("/home/njh264/Idats/James_Methylation_Samplesheet_Apr2022.csv", 
#                       stringsAsFactors = FALSE, 
#                       sep = ",", 
#                       head = TRUE,
#                       row.names = 1)
# 
# 
# samp.df <- read.delim("/home/njh264/Idats/James_Methylation_Samplesheet_450k_Apr2022.csv", 
#                       stringsAsFactors = FALSE, 
#                       sep = ",", 
#                       head = TRUE,
#                       row.names = 1)
# 
# 
# samp.df <- read.delim("/home/njh264/Idats/James_Methylation_Samplesheet_Epic_Apr2022.csv", 
#                       stringsAsFactors = FALSE, 
#                       sep = ",", 
#                       head = TRUE,
#                       row.names = 1)

source("./source_functions.R")
temp.base <- get_basenames("/home/njh264/Idats/Epic")

temp.processed <- process_idats(temp.base)


final.metageneATRT <- extract.metagene(
  as.character(atrt.meth.os.meta.n8.extract[[6]]$genes),
  as.numeric(atrt.meth.os.meta.n8.extract[[6]]$weights),
  beta2m(temp.processed$betas),
  as.numeric(atrt.meth.os.meta.n8.extract[[7]])
) 

final.metageneECRT <- extract.metagene(
  as.character(ecrt.meth.os.meta.n20.extract[[6]]$genes),
  as.numeric(ecrt.meth.os.meta.n20.extract[[6]]$weights),
  beta2m(temp.processed$betas),
  as.numeric(ecrt.meth.os.meta.n20.extract[[7]])
) 

final.metageneALL <- extract.metagene(
  as.character(all.meth.os.meta.n54.extract[[6]]$genes),
  as.numeric(all.meth.os.meta.n54.extract[[6]]$weights),
  beta2m(temp.processed$betas),
  as.numeric(all.meth.os.meta.n54.extract[[7]])
) 
