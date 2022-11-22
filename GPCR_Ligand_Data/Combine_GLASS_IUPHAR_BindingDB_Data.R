####################
# Combine data sets
####################

# set directory
setwd("/Users/kwabena/Research/GPCR/Entire_work_organized/GPCR_Ligand_Data")

# get data
GPCR_Pen <-read.csv("GPCR_PEnDB.tsv",header = T, fill=F, sep = "\t", skip = 3,
                    stringsAsFactors = F)
GPCR_Pen1 <- GPCR_Pen[,c(1,11)]
names(GPCR_Pen1) <- c("UniProt.ID.of.Target","Target.Common.name")

GPCR_Pen <-as.data.frame(GPCR_Pen)
colnames(GPCR_Pen)
GPCR_Pen <- GPCR_Pen[,c(1,12)]
names(GPCR_Pen) <- c("UniProt.ID.of.Target","IUPHAR.Class.of.Target")

# set directory
setwd("/Users/kwabena/Research/GPCR/Entire_work_organized/GPCR_Ligand_Data/BindingDB")

bindingDb.all <- read.csv(file = "bindingDb.all.tsv", sep = '\t', 
                          fill = T, stringsAsFactors = F, header = T)

# set directory
setwd("/Users/kwabena/Research/GPCR/Entire_work_organized/GPCR_Ligand_Data/IUPHAR")

iuphar.all <- read.csv(file = "iuphar.all.tsv", sep = '\t', fill = F, stringsAsFactors = F,
                       header = T)

# set directory
setwd("/Users/kwabena/Research/GPCR/Entire_work_organized/GPCR_Ligand_Data/GLASS")

glassData <- read.csv(file = "glass.all.tsv", sep = '\t', fill = F, stringsAsFactors = F,
                      header = T)

# set directory
setwd("/Users/kwabena/Research/GPCR/Entire_work_organized/GPCR_Ligand_Data")

colnames(bindingDb.all)
colnames(iuphar.all)
colnames(glass.all)


# check for the GPCRs
cond <- bindingDb.all$UniProt.ID.of.Target %in% GPCR_Pen$UniProt.ID.of.Target
sum(cond)
bindingDb.all <- bindingDb.all[cond, ]

cond <- glass.all$UniProt.ID.of.Target %in% GPCR_Pen$UniProt.ID.of.Target
sum(cond)
glass.all <- glass.all[cond, ]

cond <- iuphar.all$Target.UniProt.ID %in% GPCR_Pen$UniProt.ID.of.Target
sum(cond)
iuphar.all <- iuphar.all[cond, ]
names(iuphar.all)[c(1,6,7)] <- c("UniProt.ID.of.Target", "Ligand.SMILES", "Ligand.InChI.Key")

colnames(bindingDb.all)
colnames(iuphar.all)
colnames(glass.all)

# get them to have the same columns
bindingDbNewNames <- c("Affinity.Units","Affinity.High","Affinity.Median","Affinity.Low",
                       "Parameter","Value","Unit")
bindingDb.all[, bindingDbNewNames] <- NA
bindingDb.all[["Database.Source"]] <- "BindingDB"

glassNewNames <- c("Ligand.SMILES","Ki.nM","IC50.nM","Kd.nM","EC50.nM","Ligand.HET.ID.in.PDB",
                   "PDB.ID.for.Ligand.Target.Complex","DrugBank.ID.of.Ligand","PDB.ID.of.Target.Chain",
                   "Affinity.Units","Affinity.High","Affinity.Median","Affinity.Low")
glass.all[, glassNewNames] <- NA
glass.all[["Database.Source"]] <- "GLASS"

iupharNewNames <- c("Ki.nM","IC50.nM","Kd.nM","EC50.nM","Ligand.HET.ID.in.PDB",
                    "PDB.ID.for.Ligand.Target.Complex","DrugBank.ID.of.Ligand","PDB.ID.of.Target.Chain",
                    "Parameter","Value","Unit")
iuphar.all[, iupharNewNames] <- NA
iuphar.all[["Database.Source"]] <- "IUPHAR"


# check if they have same columns
colnames(bindingDb.all)
colnames(iuphar.all)
colnames(glass.all)

# get the columns in the same order
ord_glass <- c(1,2,7,3,4,5,8,9,10,11,17,18,19,16,12,14,13,15,6)
ord_iuphar <- c(1,7,6,16,17,18,8,9,10,11,3,4,5,2,12,14,13,15,19)
ord_bindingDb <- c(11,2,1,16,17,18,3,4,5,6,13,14,15,12,7,9,8,10,19)
bindingDb.all <- bindingDb.all[, ord_bindingDb ]
iuphar.all <- iuphar.all[, ord_iuphar]
glass.all <- glass.all[, ord_glass]


# append them together
Final_Data <- rbind(iuphar.all, bindingDb.all)
Final_Data <- rbind(Final_Data, glass.all)

# check and keep only GPCR data
cond <- Final_Data$UniProt.ID.of.Target == ""
sum(cond)
Final_Data <- Final_Data[!cond, ]

cond <- Final_Data$UniProt.ID.of.Target %in% GPCR_Pen$UniProt.ID.of.Target
sum(cond)
Final_Data <- Final_Data[cond, ]

# Checking for duplicates
suppressPackageStartupMessages(library(dplyr))
Final_Data <- distinct(Final_Data)

# ligands with no InChIKey
cond <- Final_Data$Ligand.InChI.Key == ""
sum(cond)
Final_Data <- Final_Data[-cond, ]

Final_Data <- merge(Final_Data, GPCR_Pen, by.x = "UniProt.ID.of.Target", 
                     by.y = "UniProt.ID.of.Target")
cond <- Final_Data$IUPHAR.Class.of.Target ==""
sum(cond)


#######################################################################################
# Restructure Final_Data set to have all affinity parameters on a ligand-GPCR on one row
#######################################################################################

# keep only Set 1: Top 7 highest (Ki, IC50, Potency, EC50, Inhibition, Activity, Kd, pKB,
# pKi, pIC50, pKd, pEC50) parameter
table(Final_Data$Parameter)
table(Final_Data$Affinity.Units)
table(Final_Data$Unit)

## combine some the parameters
# p[A50] = p[A]50 = pEC50 
cond <- (Final_Data$Parameter == "p[A50]")|(Final_Data$Parameter == "p[A]50")|
  (Final_Data$Parameter == "pEC50")
sum(cond, na.rm = T)
Final_Data$Parameter[cond] <- "pEC50"
# -Log KD = pKD
cond <- (Final_Data$Parameter == "pKD")|(Final_Data$Parameter == "-Log KD")
sum(cond, na.rm = T)
Final_Data$Parameter[cond] <- "pKd"
# pKi(uM) = pKi
cond <- (Final_Data$Parameter == "pKi(uM)")|(Final_Data$Parameter == "pKi")
sum(cond, na.rm = T)
Final_Data$Parameter[cond] <- "pKi"
# pKB = PKb = PkB = pkB = pKb = -Log K B = -Log KB 
cond <- (Final_Data$Parameter=="pKb")|(Final_Data$Parameter=="pKB")|(Final_Data$Parameter=="PKb")|
  (Final_Data$Parameter=="PkB")|(Final_Data$Parameter=="pkB")|(Final_Data$Parameter == "-Log KB")|
  (Final_Data$Parameter == "-Log K B")
sum(cond, na.rm = T)
Final_Data$Parameter[cond] <- "pKB"
# Ki = KI_MICROM
cond <- (Final_Data$Parameter == "KI_MICROM")|(Final_Data$Parameter == "Ki")
sum(cond, na.rm = T)
Final_Data$Parameter[cond] <- "Ki"

# subset Final_Data for data from glass, iuphar, bindingdb
cond <- Final_Data$Database.Source == "GLASS"
sub.glass <- Final_Data[cond,]
cond <- Final_Data$Database.Source == "IUPHAR"
sub.iuphar <- Final_Data[cond,]
cond <- Final_Data$Database.Source == "BindingDB"
sub.bindingdb <- Final_Data[cond,]

# remove duplicates in glass data, iuphar, and bindingdb
sub.glass <- distinct(sub.glass)
sub.iuphar <- distinct(sub.iuphar)
sub.bindingdb <- distinct(sub.bindingdb)

cond <- (sub.glass$Parameter == "Ki")|(sub.glass$Parameter == "IC50")|
  (sub.glass$Parameter == "Potency")|(sub.glass$Parameter == "EC50")|
  (sub.glass$Parameter == "Inhibition")|(sub.glass$Parameter == "Activity")|
  (sub.glass$Parameter == "Kd")|(sub.glass$Parameter == "pKB")|
  (sub.glass$Parameter == "pKi")|(sub.glass$Parameter == "pIC50")|
  (sub.glass$Parameter == "pKd")|(sub.glass$Parameter == "pEC50")
sum(cond, na.rm = T)
sub.glass <- sub.glass[cond,]

# keep units nM(Ki, Kd, IC50, EC50, Potency), %(Activity, Inhibition), -(pKB) 
table(sub.glass$Unit, sub.glass$Parameter)

cond <- (sub.glass$Unit == "nM")|(sub.glass$Unit == "%")|(sub.glass$Unit == "-")
sum(cond, na.rm = T)
sub.glass <- sub.glass[cond,]

# keep only pKi, pIC50, pKd, pEC50, pKB of Affinity.Units
cond <- (sub.iuphar$Affinity.Units=="pKi")|(sub.iuphar$Affinity.Units=="pIC50")|
  (sub.iuphar$Affinity.Units=="pKd")|(sub.iuphar$Affinity.Units=="pEC50")|
  (sub.iuphar$Affinity.Units=="pKB")
sum(cond, na.rm = T)
sub.iuphar <- sub.iuphar[cond,]

#### Glass
colnames(Final_Data)
table(sub.glass$Parameter)
table(sub.iuphar$Affinity.Units)
table(sub.bindingdb)
# create the variables Activity, Inhibition, Potency, pKB, pKi, pIC50, pEC50, pKd
var.names <- c("Activity.%","Inhibition.%","Potency.nM","pKB","pKi","pIC50","pEC50","pKd")
sub.glass[,var.names] <- NA
colnames(sub.glass)

var_indx <- c(21,22,23,7,8,10,9,24,27,25,26,28)
parameters <- c("Activity", "Inhibition", "Potency", "Ki", "IC50", "EC50", "Kd", "pKB",
                "pEC50", "pKi", "pIC50", "pKd")
for (i in 1:length(parameters)){
  cond <- which((sub.glass$Parameter == parameters[i]), arr.ind = T)
  sub.glass[cond,var_indx[i]] <- sub.glass$Value[cond]
}

#### IUPHAR
sub.iuphar[,var.names] <- NA
parameters <- c("pKB","pEC50","pIC50","pKd","pKi")
colnames(sub.glass)

var_indx <- c(24,27,26,28,25)
for (i in 1:length(parameters)){
  cond <- which((sub.iuphar$Affinity.Units == parameters[i]), arr.ind = T)
  sub.iuphar[cond,var_indx[i]] <- paste(sub.iuphar$Affinity.Low[cond],
                                    sub.iuphar$Affinity.High[cond], sep = " - ")
}

#### Bindingdb
sub.bindingdb[,var.names] <- NA

# Rearrange column
Final_Data <- rbind(sub.bindingdb,sub.iuphar,sub.glass)
Final_Data <- Final_Data[,-c(4,5,6,11,12,13,14)]
colnames(Final_Data)
var.nm <- c(1,11,13,2,3,9,8,4:7,14:20,10,12)
Final_Data <- Final_Data[, var.nm]
colnames(Final_Data)

######### Separate operators from affinity values
library(stringr)
operator <- c("<" , ">", "=")
Final_Data$Affinity_relation <- NA

### Ki.nM do for "<" , ">", "="
for (op in operator) {
  symb <- str_which(Final_Data$Ki.nM, op)
  if (length(symb) >= 1){
    ki <- Final_Data$Ki.nM[symb]
    ki <- unlist(ki)
    ki_List <- strsplit(ki,split= op, fixed=T)
    
    ki.val <- data.frame(stringsAsFactors = F)
    for (q in 1:length(ki_List)) {
      ki.val.list <- data.frame(val = ki_List[[q]][2], stringsAsFactors = F)
      ki.val <- rbind(ki.val, ki.val.list, stringsAsFactors = F)
    }
    ki.val <- as.vector(unlist(ki.val))
    
    Final_Data$Ki.nM[symb] <- ki.val
    Final_Data$Affinity_relation[symb] <- op
    
    rm(ki_List, ki.val, ki.val.list,symb,ki,q)
  }
}

### Kd.nM do for "<" , ">", "="
for (op in operator) {
  symb <- str_which(Final_Data$Kd.nM, op)
  if (length(symb) >= 1){
    kd <- Final_Data$Kd.nM[symb]
    kd <- unlist(kd)
    kd_List <- strsplit(kd,split=op, fixed=T)
    
    kd.val <- data.frame(stringsAsFactors = F)
    for (q in 1:length(kd_List)) {
      kd.val.list <- data.frame(val = kd_List[[q]][2], stringsAsFactors = F)
      kd.val <- rbind(kd.val, kd.val.list, stringsAsFactors = F)
    }
    kd.val <- as.vector(unlist(kd.val))
    
    Final_Data$Kd.nM[symb] <- kd.val
    Final_Data$Affinity_relation[symb] <- op
    
    rm(kd_List, kd.val, kd.val.list,symb,kd,q)
  }
}

### IC50.nM do for "<" , ">"
for (op in operator) {
  symb <- str_which(Final_Data$IC50.nM, op)
  if (length(symb) >= 1){
    IC50 <- Final_Data$IC50.nM[symb]
    IC50 <- unlist(IC50)
    IC50_List <- strsplit(IC50,split= op, fixed=T)
    
    IC50.val <- data.frame(stringsAsFactors = F)
    for (q in 1:length(IC50_List)) {
      IC50.val.list <- data.frame(val = IC50_List[[q]][2], stringsAsFactors = F)
      IC50.val <- rbind(IC50.val, IC50.val.list, stringsAsFactors = F)
    }
    IC50.val <- as.vector(unlist(IC50.val))
    
    Final_Data$IC50.nM[symb] <- IC50.val
    Final_Data$Affinity_relation[symb] <- op
    
    rm(IC50_List, IC50.val, IC50.val.list,symb,IC50,q)
  }
}

### EC50.nM do for "<" , ">", "="
for (op in operator) {
  symb <- str_which(Final_Data$EC50.nM, op)
  if (length(symb) >= 1){
    EC50 <- Final_Data$EC50.nM[symb]
    EC50 <- unlist(EC50)
    EC50_List <- strsplit(EC50,split=op, fixed=T)
    
    EC50.val <- data.frame(stringsAsFactors = F)
    for (q in 1:length(EC50_List)) {
      EC50.val.list <- data.frame(val = EC50_List[[q]][2], stringsAsFactors = F)
      EC50.val <- rbind(EC50.val, EC50.val.list, stringsAsFactors = F)
    }
    EC50.val <- as.vector(unlist(EC50.val))
    
    Final_Data$EC50.nM[symb] <- EC50.val
    Final_Data$Affinity_relation[symb] <- op
    
    rm(EC50_List, EC50.val, EC50.val.list,symb,EC50,q)
  }
}

# make some variables numeric
Final_Data$Ki.nM <- as.numeric(Final_Data$Ki.nM)
Final_Data$Kd.nM <- as.numeric(Final_Data$Kd.nM)
Final_Data$IC50.nM <- as.numeric(Final_Data$IC50.nM)
Final_Data$EC50.nM <- as.numeric(Final_Data$EC50.nM)

cond <- is.na(Final_Data$Affinity_relation)
sum(cond, na.rm = T)
Final_Data$Affinity_relation[cond] <- "="

# Rearrange columns 
colnames(Final_Data)
Final_Data <- Final_Data[,c(1:7,21,8:20)]

Final_Data <- merge(Final_Data, GPCR_Pen1, by.x = "UniProt.ID.of.Target", 
                    by.y = "UniProt.ID.of.Target", suffixes = c("",".y"))
Final_Data <- Final_Data[,c(1,22,2:21)]
cond <- Final_Data$Ligand.InChI.Key == ""
sum(cond, na.rm = T)
Final_Data <- Final_Data[!cond,]

# write data to tsv file
write.table(Final_Data, sep = "\t",
            file = "/Users/kwabena/Research/GPCR/Entire_work_organized/GPCR_Ligand_Data/Final_Data.tsv",
            row.names=FALSE)

rm(list = ls())
