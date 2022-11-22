###################################################################
# Get Ligands that bind to multiple GPCRs of different families
###################################################################
library(foreach)
library(doParallel)
library(dplyr)

# set directory
setwd("/Users/kwabena/Research/GPCR/Entire_work_organized/GPCR_Ligand_Data")

Final_Data <- read.table(file = "Final_Data.tsv", sep = "\t", fill = F, stringsAsFactors = F,
                       header = T)
# Find unique ligands
ligand_uniq <- as.vector(unique(Final_Data$Ligand.InChI.Key))
l.na <- which(is.na(ligand_uniq), arr.ind = T)
if (sum(l.na) > 0){
  ligand_uniq <- ligand_uniq[-c(l.na)]
}

# Parallel approach to finding ligands that binds to GPCRs of multiple families
# We save the data  
cores=detectCores()
cl <- makeForkCluster(cores[1]-4) #not to overload the computer
registerDoParallel(cl)

ligand_mult_GPCR_fam <- foreach(i=1:length(ligand_uniq), .combine = rbind) %dopar% {
  cond <- Final_Data$Ligand.InChI.Key == ligand_uniq[i]
  if (sum(cond, na.rm = T) > 1){
    subs <- Final_Data[cond,c(1,2,4,5)]
    subs <- distinct(subs)
    len <- length(unique(subs$IUPHAR.Class.of.Target))
    if (len > 1){
      subs
    }
  }  
}
stopCluster(cl)

# write data to tsv file
write.table(ligand_mult_GPCR_fam, sep = "\t",
            file = "/Users/kwabena/Research/GPCR/Entire_work_organized/ligand_mult_GPCR_fam_Data.tsv",
            row.names=FALSE)

cond <- ligand_mult_GPCR_fam$Target.Common.name == "Human"
sum(cond, na.rm = T)
ligand_mult_GPCR_fam <- ligand_mult_GPCR_fam[cond,]
ligand_uniq <- as.vector(unique(ligand_mult_GPCR_fam$Ligand.InChI.Key))
rm(cond)

cores=detectCores()
cl <- makeForkCluster(cores[1]-4) #not to overload the computer
registerDoParallel(cl)

ligand_mult_GPCR_fam_Human <- foreach(i=1:length(ligand_uniq), .combine = rbind) %dopar% {
  cond <- ligand_mult_GPCR_fam$Ligand.InChI.Key == ligand_uniq[i]
  if (sum(cond, na.rm = T) > 1){
    subs <- ligand_mult_GPCR_fam[cond,]
    subs <- distinct(subs)
    len <- length(unique(subs$IUPHAR.Class.of.Target))
    if (len > 1){
      subs
    }
  }  
}

stopCluster(cl)

# write data to tsv file
write.table(ligand_mult_GPCR_fam_Human, sep = "\t",
            file = "/Users/kwabena/Research/GPCR/Entire_work_organized/ligand_mult_GPCR_fam_Human_Data.tsv",
            row.names=FALSE)

rm(list = ls())
