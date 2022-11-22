# This code cleans the data downloaded from BindingDB
# We only excludes some of the columns 
# The cleaning is done before it is combined with data
# from GLASS and IUPHAR

# set directory
setwd("/Users/kwabena/Research/GPCR/Entire_work_organized/GPCR_Ligand_Data/BindingDB")

bindingDb.all <- read.csv(file = "BindingDB_All.tsv", sep = '\t',
                          fill = T, stringsAsFactors = F, header = T)

colnames(bindingDb.all) 
# Columns to retain
names(bindingDb.all)[c(42)] <- c("UniProt.ID.of.Target")
bindingDb.all <- bindingDb.all[, c(2,4,9:12,27,28,33,39,42)]
cond <- which(bindingDb.all$UniProt.ID.of.Target == "")
if(length(cond) != 0){
  bindingDb.all <- bindingDb.all[-cond,] 
}

names(bindingDb.all) <- c("Ligand.SMILES","Ligand.InChI.Key",
                          "Ki.nM","IC50.nM","Kd.nM","EC50.nM","Ligand.HET.ID.in.PDB",
                          "PDB.ID.for.Ligand.Target.Complex","DrugBank.ID.of.Ligand","PDB.ID.of.Target.Chain",
                          "UniProt.ID.of.Target")

write.table(bindingDb.all, sep = "\t",
          file = "/Users/kwabena/Research/GPCR/Entire_work_organized/GPCR_Ligand_Data/BindingDB/bindingDb.all.tsv",
          row.names = FALSE)

remove(list = ls(all.names = T))
