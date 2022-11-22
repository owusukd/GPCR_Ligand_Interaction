# This code cleans the data downloaded from GLASS
# We only excludes some of the columns 
# The cleaning is done before it is combined with data
# from BindingDB and IUPHAR

# set directory
setwd("/Users/kwabena/Research/GPCR/Entire_work_organized/GPCR_Ligand_Data/GLASS")

glass.all <- read.csv(file = "interactions_total.tsv", sep = '\t',
                          fill = T, stringsAsFactors = F, header = T)

colnames(glass.all) 
# Columns to retain
glass.all <- glass.all[, c(1:6)]
cond <- which(glass.all$UniProt.ID == "")
if(length(cond) != 0){
  glass.all <- glass.all[-cond,]
}

names(glass.all) <- c("UniProt.ID.of.Target","Ligand.InChI.Key","Parameter",
                      "Value","Unit","Database.Source")

write.table(glass.all, sep = "\t",
          file = "/Users/kwabena/Research/GPCR/Entire_work_organized/GPCR_Ligand_Data/GLASS/glass.all.tsv",
          row.names = FALSE)

remove(list = ls(all.names = T))
