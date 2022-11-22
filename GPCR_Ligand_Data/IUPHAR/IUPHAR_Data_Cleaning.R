# This code cleans the data downloaded from IUPHAR
# We only excludes some of the columns 
# The cleaning is done before it is combined with data
# from GLASS and BindingDB

# set directory
setwd("/Users/kwabena/Research/GPCR/Entire_work_organized/GPCR_Ligand_Data/IUPHAR")

interactions <- read.csv(file = "interactions.tsv", sep = '\t', fill = F, stringsAsFactors = F,
                         header = T, skip = 1)
colnames(interactions)
if(colnames(interactions)[16] == "X.Ligand.ID"){
  names(interactions)[16] <- "Ligand.ID"
}

ligands <- read.csv(file = "ligands.tsv", sep = '\t', 
                    fill = F, stringsAsFactors = F, header = T, skip = 1)
colnames(ligands)

iuphar.all <- merge(interactions, ligands, by.x = "Ligand.ID", by.y = "Ligand.id")

colnames(iuphar.all)

iuphar.all <- iuphar.all[, c(6,29,30,31,32,55,56)]

write.table(iuphar.all, sep = "\t",
          file = "/Users/kwabena/Research/GPCR/Entire_work_organized/GPCR_Ligand_Data/IUPHAR/iuphar.all.tsv",
          row.names = F)

remove(list = ls(all.names = T))