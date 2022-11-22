# cut sequences by adding 5 amino acid at end of the 1st extracellular loop
# and adds 5 amino acid at beginning and the end of the remaining extracellular loops
# and writes all the extracellular loops as one file for each ligand
# set directory
setwd("/home/owusu/Desktop/GPCR/final_data/fasta_MSA/")
library(seqinr)

ligID <- c("FQUAFMNPXPXOJE-UHFFFAOYSA-N","MLQFOEOUNIRULR-UHFFFAOYSA-N","YKMSTUDOGGAEJH-UHFFFAOYSA-N",
           "USZPQRMQYJIDII-UHFFFAOYSA-N","BYBLEWFAAKGYCD-UHFFFAOYSA-N","NKOPNLUYOHOGFZ-UHFFFAOYSA-N",
           "AJLFQFYMLRXVHV-UHFFFAOYSA-N","CLQVVBPDAXJGBV-UHFFFAOYSA-N","DTZDSNQYNPNCPK-UHFFFAOYSA-N",
           "IKSHHOBCJKJKOG-UHFFFAOYSA-N","XLWJPQQFJNGUPA-UHFFFAOYSA-N")

num.ligID <- length(ligID)

for (k in 1:num.ligID) {
  # read the positions of the extracellular loops, helices, etc
  filName <- paste(ligID[k],"txt", sep = ".")
  positn <- read.table(file = filName, header = T, sep = "", stringsAsFactors = F)
  # read the fasta file of the gpcrs
  fastaName <- paste(ligID[k], "fasta", sep = ".")
  fasta <- read.fasta(file = fastaName, seqtype = "AA", as.string = T, whole.header = F, 
                      strip.desc = F)
  # get the number of gpcr in fasta file
  num.seq <- (NROW(positn))/2
  
  Eis <- c(3, 7, 11, 15) # column numbers of the Ei's in the positions of the extracellular loops, helices, etc
  seq.fasta <- NULL
  seq.ID <- c() 
   
  for (j in Eis) {
    
    # row begin = i+(i-1), row end = i+((i+1)-1), i = 1:num.seq, rows = 1:NROW(positn) 
    for (i in 1:num.seq) {
      beg <- i+(i-1)
      end <- i+((i+1)-1)
      strt <- positn[beg, j]
      ends <- positn[end, j]
    
      if (j == 3){
        Positn_strt <- strt
      } else {
          Positn_strt <- strt - 5
      }
      Positn_ends <- ends + 5
      seq.sub <- substring(fasta[[i]][1], Positn_strt, Positn_ends)
      seq.InID <- paste(attr(fasta[[i]], "name"), names(positn)[j], sep = "")
      seq.ID <- c(seq.ID, seq.InID)
      seq.list <- list(seq = seq.sub)
      seq.fasta <- rbind(seq.fasta, seq.list)
    }
  }
      
  # write fasta file
  file.out <- paste(ligID[k], "E", sep = "_")
  file.out <- paste(file.out, "txt", sep = ".")
  write.fasta(sequences = seq.fasta, names = seq.ID, file.out = file.out, open = "w", 
              nbchar = 60, as.string = T)
}
rm(list = ls())
