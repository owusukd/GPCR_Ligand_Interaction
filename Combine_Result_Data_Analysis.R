# Adding 3D similar Overlap Pocket and AA of Pocket
# Read data and merges them  in a special 
# way that assigns the right value to the 
# right row and column
# 
library(xlsx)
library(readxl)
library(tibble)

# Read data
# data files have the same naming format: 
# 1. for docking it is "ligang_name"_AutoDock_vina_Results.tsv
# e.g. AJLF_AutoDock_vina_Results.tsv
# 2. for pocket it is "ligang_name"_APoc_Binding_pocket_comparison.tsv
# e.g. AJLF_APoc_Binding_pocket_comparison.tsv

ligand <- c("AJLF","CLQV","DTZD","IKSH","NKOP","USZP","XLWJ","YKMS")

data_comb_final <- data.frame()

for (lig in ligand) {
  # Read docking results data
  path_dock <- c("/Users/kwabena/Research/GPCR/AutoDock_ligands_GPCR") 
  path_dock <- paste(path_dock, lig, sep = '/')
  file_dock <- paste(lig, c("AutoDock_vina_Results.tsv"), sep = '_')
  setwd(path_dock)
  data_dock <- read.table(file = file_dock, header = T, sep = "\t", 
                          stringsAsFactors = F)
  
  # Read pocket comparison results data
  path_pock <- c("/Users/kwabena/Research/GPCR/Binding_pocket_comparison")
  path_pock <- paste(path_pock, lig, c("APoc"), sep = '/')
  file_pock <- paste(lig, c("APoc_Binding_pocket_comparison.tsv"), sep = '_')
  setwd(path_pock)
  data_pock <- read.table(file = file_pock, header = T, sep = "\t", 
                          stringsAsFactors = F)
  
  # add columns to data_pock for each ligand
  data_pock <- add_column(data_pock, ligand_ID = lig, .before = "pocket_1" )
  data_pock <- add_column(data_pock, pocket_1_Affinity_kcal_mol = NA, 
                          .after = "pocket_1" )
  data_pock <- add_column(data_pock, pocket_2_Affinity_kcal_mol = NA, 
                          .after = "pocket_2" )
  
  # extract the affinities for each ligand and add it to data_pock
  pocket_1_pdbID <- substr(data_pock[,2], 10, 13)
  pocket_1_num <- as.numeric(substr(data_pock[,2], 22, 22))
  pocket_2_pdbID <- substr(data_pock[,4], 10, 13)
  pocket_2_num <- as.numeric(substr(data_pock[,4], 22, 22))
  affinity_1 <- c()
  affinity_2 <- c()
  for (row in 1:NROW(data_pock)) {
    cond1 <- (data_dock$protein_PDB_ID == pocket_1_pdbID[row] & 
                data_dock$pocket_number == pocket_1_num[row])
    affinity_1[row] <- data_dock[cond1,4]
    
    cond2 <- (data_dock$protein_PDB_ID == pocket_2_pdbID[row] & 
                data_dock$pocket_number == pocket_2_num[row])
    affinity_2[row] <- data_dock[cond2,4]
  }
  # find the adsolute difference of the affinities of the two pockets
  data_pock$abs_diff_aff_1_aff2 <- abs(affinity_1-affinity_2)
  data_pock$pocket_1_Affinity_kcal_mol <- affinity_1
  data_pock$pocket_2_Affinity_kcal_mol <- affinity_2
  
  # combine data for all ligands
  data_comb_final <- rbind(data_comb_final, data_pock)
}

setwd("/Users/kwabena/Research/GPCR/AutoDock_ligands_GPCR")
write.table(data_comb_final, file = "Combine_docking_pocket_comp_results_data.tsv",
            sep = "\t", row.names = F)

setwd("/Users/kwabena/Research/GPCR/AutoDock_ligands_GPCR")
data_comb_final <- read.table("Combine_docking_pocket_comp_results_data.tsv",sep = "\t",
                             stringsAsFactors = F)

setwd("/Users/kwabena/Research/GPCR/3D_Similar_Coincide_Pkt")
Similar_Coincide_Pkt <- read.table(file = "3D_Similar_Overlap_Pkt.csv", header = T, 
                                   sep = ",", stringsAsFactors = F)
AA_Pkt <- read.table(file = "AA_Pkt.csv", header = T, 
                     sep = "\t", stringsAsFactors = F)

setwd("/Users/kwabena/Research/GPCR/AutoDock_ligands_GPCR")
Ligs_Align_Pkt <- read.table(file = "Ligs_Align_Pkt.csv", header = T, 
                             sep = "\t", stringsAsFactors = F)

######## 3D similar Overlap Pocket 
for (lig in ligand) {
  cond_data_comb <- which(data_comb_final$ligand_ID == lig)
  cond_similar_coincide <- which(Similar_Coincide_Pkt$ligand_ID == lig)
  for (i in cond_similar_coincide) {
    for (j in cond_data_comb) {
      if (substr(data_comb_final[j,2], 10, 13)== substr(Similar_Coincide_Pkt[i,2], 1, 4)){
        if (substr(data_comb_final[j,4], 10, 13)== substr(Similar_Coincide_Pkt[i,2], 11, 14)){
          if (substr(data_comb_final[j,2], 10, 13)== substr(Similar_Coincide_Pkt[i,3], 1, 4)){
            if (substr(data_comb_final[j,2], 22, 22)== substr(Similar_Coincide_Pkt[i,3], 10, 10)){
              data_comb_final[j,15] <- Similar_Coincide_Pkt[i,4]
              data_comb_final[j,16] <- Similar_Coincide_Pkt[i,5]
            }
          }
        }
      }
      
      if (substr(data_comb_final[j,2], 10, 13)== substr(Similar_Coincide_Pkt[i,2], 11, 14)){
        if (substr(data_comb_final[j,4], 10, 13)== substr(Similar_Coincide_Pkt[i,2], 1, 4)){
          if (substr(data_comb_final[j,4], 10, 13)== substr(Similar_Coincide_Pkt[i,3], 1, 4)){
            if (substr(data_comb_final[j,4], 22, 22)== substr(Similar_Coincide_Pkt[i,3], 10, 10)){
              data_comb_final[j,17] <- Similar_Coincide_Pkt[i,4]
              data_comb_final[j,18] <- Similar_Coincide_Pkt[i,5]
            }
          }
        }
      }
    }
  }
}

######## Alignment of ligands for each pockets compared
for (lig in ligand) {
  cond_data_comb <- which(data_comb_final$ligand_ID == lig)
  cond_ligs_align <- which(Ligs_Align_Pkt$ligand_ID == lig)
  for (i in cond_data_comb) {
    for (j in cond_ligs_align) {
      if (substr(data_comb_final[i,2], 10, 13)== substr(Ligs_Align_Pkt[j,2], 1, 4)){
        if (substr(data_comb_final[i,2], 22, 22)== substr(Ligs_Align_Pkt[j,2], 6, 6)){
          if (substr(data_comb_final[i,4], 10, 13)== substr(Ligs_Align_Pkt[j,2], 8, 11)){
            if (substr(data_comb_final[i,4], 22, 22)== substr(Ligs_Align_Pkt[j,2], 13, 13)){
              data_comb_final$Pkt_ligs_align_rmsd[i] <- Ligs_Align_Pkt$RMSD[j]
              break      
            }
          }
        }
      }
    }
  }
}

######## AA of Pocket
for (ligs in ligand) {
  cond_data_comb_2 <- which(data_comb_final$ligand_ID == ligs)
  cond_AA_Pkt <- which(AA_Pkt$ligand_ID == ligs)
  for (s in cond_AA_Pkt) {
    for (n in cond_data_comb_2) {
      
      if (substr(data_comb_final[n,2], 10, 13) == AA_Pkt[s,2]){
        if (substr(data_comb_final[n,2], 22, 22) == AA_Pkt[s,3]){
          data_comb_final$GPCR_Resi_pkt_1[n] <- AA_Pkt[s,4]
        }
      }
      if (substr(data_comb_final[n,4], 10, 13) == AA_Pkt[s,2]){
        if (substr(data_comb_final[n,4], 22, 22) == AA_Pkt[s,3]){
          data_comb_final$GPCR_Resi_pkt_2[n] <- AA_Pkt[s,4]
        }
      }
      
    }
  }
}

# Correlation analysis
corr_table <- data.frame(row.names = F)
for (i in 1:length(ligand)) {
  cond <- which(data_comb_final$ligand_ID == ligand[i])
  data <- data_comb_final[cond,]
  corr <- cor.test(data$PS_score, data$abs_diff_aff_1_aff2,
                   alternative = c("two.sided"), method = c("pearson"), 
                   conf.level = 0.95)
  corr_list <- list(Ligand = ligand[i], 
                    Correlation_PS_score_abs_diff_aff1_aff2_pockets=as.vector(corr$estimate),
                    P_value =corr$p.value)
  corr_table <- rbind(corr_table, corr_list, stringsAsFactors = F)
}



setwd("/Users/kwabena/Research/GPCR/Manuscript/MDPI biomolecules/New_Data_Control")
data_control <- read_excel("Combined_Data.xlsx")

colnames(data_control)

######## Convert AA 3 letter code to 1 letter code
AA <- c("ALA"="A","ARG"="R","ASN"="N","ASP"="D","CYS"="C","GLU"="E","GLN"="Q","GLY"="G","HIS"="H",
        "ILE"="I","LEU"="L","LYS"="K","MET"="M","PHE"="F","PRO"="P","SER"="S","THR"="T","TRP"="W",
        "TYR"="Y","VAL"="V")
data_control$Protein_AA_Interact_pkt_1 <- as.character(data_control$Protein_AA_Interact_pkt_1)
data_control$Protein_AA_Interact_pkt_2 <- as.character(data_control$Protein_AA_Interact_pkt_2)
data_control$Resi_pkt_1 <- as.character(data_control$Resi_pkt_1)
data_control$Resi_pkt_2 <- as.character(data_control$Resi_pkt_2)

for (k in 1:NROW(data_control)) {
  ####### Begin: Residues that are involved in interaction 
  if (data_control$Protein_AA_Interact_pkt_1[k] != "-"){
    pkt_1 <- unlist(strsplit(data_control$Protein_AA_Interact_pkt_1[k], split = ","))
    pkt_1_seq = ""
    for (l in 1:length(pkt_1)) {
      pkt_1_seq <- paste0(pkt_1_seq, AA[toupper(pkt_1[l])])
    }
    if (nchar(pkt_1_seq) == length(pkt_1)){
      data_control$AA_Interact_pkt_1[k] <- pkt_1_seq 
    }
  }else{
    data_control$AA_Interact_pkt_1[k] <- "-"
  }
  
  if (data_control$Protein_AA_Interact_pkt_2[k] != "-"){
    pkt_2 <- unlist(strsplit(data_control$Protein_AA_Interact_pkt_2[k], split = ","))
    pkt_2_seq = ""
    for (m in 1:length(pkt_2)) {
      pkt_2_seq <- paste0(pkt_2_seq, AA[toupper(pkt_2[m])])
    }
    if (nchar(pkt_2_seq) == length(pkt_2)){
      data_control$AA_Interact_pkt_2[k] <- pkt_2_seq
    }
  }else{
    data_control$AA_Interact_pkt_2[k] <- "-"
  }
  ####### End: Residues that are involved in interaction 
  
  ####### Begin: Residues in the pocket
  pkts_1 <- unlist(strsplit(data_control$Resi_pkt_1[k], split = ","))
  pkts_1_seq = ""
  for (t in 1:length(pkts_1)) {
    pkts_1_seq <- paste0(pkts_1_seq, AA[toupper(pkts_1[t])])
  }
  if (nchar(pkts_1_seq) == length(pkts_1)){
    data_control$AA_pkt_1[k] <- pkts_1_seq 
    
  }
  
  pkts_2 <- unlist(strsplit(data_control$Resi_pkt_2[k], split = ","))
  pkts_2_seq = ""
  for (u in 1:length(pkts_2)) {
    pkts_2_seq <- paste0(pkts_2_seq, AA[toupper(pkts_2[u])])
  }
  if (nchar(pkts_2_seq) == length(pkts_2)){
    data_control$AA_pkt_2[k] <- pkts_2_seq
  }
  ####### End: Residues in the pocket 
}

######## Create more features
library(Peptides)
for (v in 1:NROW(data_control)) {
  ############################################### Pkt 1
  seqs1 <- data_control$AA_pkt_1[v]
  
  # MS-WHIM scores of a protein sequence: MS-WHIM scores were derived from 36 
  # electrostatic potential properties derived from the three- dimensional 
  # structure of the 20 natural amino acids
  c = mswhimScores(seqs1)[[1]][c(1,2,3)]
  data_control$Pkt_1_MSWHIM1[v] <- c[1]; data_control$Pkt_1_MSWHIM2[v] <- c[2]; 
  data_control$Pkt_1_MSWHIM3[v] <- c[3]; 
  
  ############################################### Pkt 2
  seqs2 <- data_control$AA_pkt_2[v]
  
  # MS-WHIM scores of a protein sequence: MS-WHIM scores were derived from 36 
  # electrostatic potential properties derived from the three- dimensional 
  # structure of the 20 natural amino acids
  cc = mswhimScores(seqs2)[[1]][c(1,2,3)]
  data_control$Pkt_2_MSWHIM1[v] <- cc[1]; data_control$Pkt_2_MSWHIM2[v] <- cc[2]; 
  data_control$Pkt_2_MSWHIM3[v] <- cc[3]; 
  
}

# MSWHIM distance

lf_dist_MSWHIM <- function(x,y){
  sum_abs_diff_mswhim <- max(abs(x[1]-y[1]), abs(x[2]-y[2]), abs(x[3]-y[3]))
  sum_abs_diff_mswhim
}

colnames(data_control)
for (h in 1:NROW(data_control)) {
  x <- unlist(as.vector(data_control[h,c(21,22,23)]))
  y <- unlist(as.vector(data_control[h,c(24,25,26)]))
  data_control$Lf_Dist_MSWHIM[h] <- lf_dist_MSWHIM(x,y)
}

# Affinity
l1_dist_Affinity <- function(x,y){
  sum_abs_diff_Affinity <- abs(x-y)
  sum_abs_diff_Affinity
}

colnames(data_control)
for (h in 1:NROW(data_control)) {
  x <- unlist(as.vector(data_control[h,c(15)]))
  y <- unlist(as.vector(data_control[h,c(16)]))
  data_control$L1_Dist_Affinity[h] <- l1_dist_Affinity(x,y)
}


# Num_Same_Resi_Interact
for (s in 1:NROW(data_control)) {
  string <- strsplit(c(data_control$Protein_AA_Interact_pkt_1[s],
                       data_control$Protein_AA_Interact_pkt_2[s]), ",")
  count_aa <- 0
  if(length(string[[1]]) <= length(string[[2]])){
    for (t in string[[1]]) {
      indx <- which(string[[2]] == t)
      if (length(indx) >= 1){
        count_aa = count_aa + 1
        string[[2]] <- string[[2]][-c(indx[1])]
      }
    }
  }else{
    for (t in string[[2]]) {
      indx <- which(string[[1]] == t)
      if (length(indx) >= 1){
        count_aa = count_aa + 1
        string[[1]] <- string[[1]][-c(indx[1])]
      }
    }
  }
  data_control$Num_Same_Resi_Interact[s] <- count_aa
}

# P3D_similar_Coincide_pkt
data_control$Average_P3D_similar_Coincide_pkt_1 <- 0.5*(data_control$P3D_similar_Coincide_pkt_1_flex + 
                                                         data_control$P3D_similar_Coincide_pkt_1_rigid)
data_control$Average_P3D_similar_Coincide_pkt_2 <- 0.5*(data_control$P3D_similar_Coincide_pkt_2_flex + 
                                                          data_control$P3D_similar_Coincide_pkt_2_rigid)

data_control$Sum_P3D_similar_Coincide_pkts_pair <- (data_control$Average_P3D_similar_Coincide_pkt_1 + 
                                                      data_control$Average_P3D_similar_Coincide_pkt_2)


######## Write data to excel
library(xlsx)
setwd("/Users/kwabena/Research/GPCR/Manuscript/MDPI biomolecules/New_Data_Control")
write.table(data_control[,c(1:5,8,15,16,28)], file = "Control_Combined_Data.tsv",sep = "\t",
            col.names = TRUE, row.names = FALSE, append = FALSE)

######### Plot

# Actual_pkt_ligs_aligned_RMSD vs Docked_pkt_ligs_aligned_RMSD
plot(data_control$Actual_Pkt_Ligs_Aligned_RMSD, type="l",col="red", 
     ylab = "RMSD of Aligned Ligands in Pairs of Pockets", xlab = "Index of Pairs of Pockets")
lines(data_control$Docked_Pkt_Ligs_Aligned_RMSD, col="blue", type = "l")
legend(1, 3, legend=c("RMSD Acutal", "RMSD Docked"),
       col=c("red", "blue"), lty=c(1,1), cex=0.8)

######### Correlation Analysi
# PS_score vs:
# Actual_pkt_ligs_aligned_RMSD
# Lf_Dist_MSWHIM
# Num_Same_Resi_Interact
# L1_Dist_Affinity
# Sum_P3D_similar_Coincide_pkts_pair

# PS_score
mean(data_control$PS_score)
sd(data_control$PS_score)
min(data_control$PS_score)
max(data_control$PS_score)


# Actual_pkt_ligs_aligned_RMSD --- Negative and Significant
corr <- cor.test(data_control$PS_score, data_control$Actual_Pkt_Ligs_Aligned_RMSD,
                 conf.level = 0.95, alternative = c("two.sided"), method = c("pearson"))
corr

# Actual_pkt_ligs_aligned_RMSD vs Docked_pkt_ligs_aligned_RMSD --- Significant
var.test(x=data_control$Actual_Pkt_Ligs_Aligned_RMSD, y=data_control$Docked_Pkt_Ligs_Aligned_RMSD,
         alternative = "two.sided")
t.test(x=data_control$Actual_Pkt_Ligs_Aligned_RMSD, y=data_control$Docked_Pkt_Ligs_Aligned_RMSD, 
       alternative = "two.sided", var.equal = F, paired = T)

mean(data_control$Actual_Pkt_Ligs_Aligned_RMSD)
sd(data_control$Actual_Pkt_Ligs_Aligned_RMSD)

mean(data_control$Docked_Pkt_Ligs_Aligned_RMSD)
sd(data_control$Docked_Pkt_Ligs_Aligned_RMSD)

# L1_Dist_Affinity --- Negative and Significant
corr <- cor.test(data_control$PS_score, data_control$L1_Dist_Affinity,
                 alternative = c("two.sided"), method = c("pearson"), conf.level = 0.95)
corr

# Ignore
# Sum_P3D_similar_Coincide_pkts_pair --- Positive but Not Significant
corr <- cor.test(data_control$PS_score, data_control$Sum_P3D_similar_Coincide_pkts_pair,
                 alternative = c("two.sided"), method = c("pearson"), conf.level = 0.95)
corr

# Ignore
# Lf_Dist_MSWHIM --- Negative and Significant
corr <- cor.test(data_control$PS_score, data_control$Lf_Dist_MSWHIM,
                 alternative = c("two.sided"), method = c("pearson"), conf.level = 0.95)
corr

# Ignore
# Num_Same_Resi_Interact --- Positive but Not Significant
corr <- cor.test(data_control$PS_score, data_control$Num_Same_Resi_Interact,
                 alternative = c("two.sided"), method = c("pearson"), conf.level = 0.95)
corr

