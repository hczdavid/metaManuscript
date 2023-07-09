#---------------------------------------#
#  
#  Remove no Primer reads in Brown Study
#  
#  Author: David Huang
#  Update: 6/7/2021
#---------------------------------------#

# This code is used to remove the no primer reads for Ba and Bb study
# The first 20 samples have no primer so keep the original file

library(ShortRead)

# get the primer seq from the paper
FWD <- "GAGTTTGATCNTGGCTCAG" 
REV <- "TGCTGCCTCCCGTAGGAGT" 

path_ba <- "/home/chuang19/Ben_data/Brown2018a/Brown2018a_fastq" 
path_bb <- "/home/chuang19/Ben_data/Brown2018b/Brown2018b_fastq" 

fnFs_ba <- sort(list.files(path_ba, pattern="1.fastq", full.names = TRUE))
fnRs_ba <- sort(list.files(path_ba, pattern="2.fastq", full.names = TRUE))
sname_ba <- sapply(strsplit(basename(fnFs_ba), "_"), `[`, 1)

fnFs_bb <- sort(list.files(path_bb, pattern="1.fastq", full.names = TRUE))
fnRs_bb <- sort(list.files(path_bb, pattern="2.fastq", full.names = TRUE))
sname_bb <- sapply(strsplit(basename(fnFs_bb), "_"), `[`, 1)


prop_ba <- c()
for(i in 1:length(sname_ba)){
  hit1 <- vcountPattern(FWD, sread(readFastq(fnFs_ba[i])), fixed = FALSE)
  hit2 <- vcountPattern(REV, sread(readFastq(fnRs_ba[i])), fixed = FALSE)
  prop_ba[i] <- sum(hit1 == 1 & hit2 == 1)/length(hit1)
}

prop_bb <- c()
for(i in 1:length(sname_bb)){
  hit1 <- vcountPattern(FWD, sread(readFastq(fnFs_bb[i])), fixed = FALSE)
  hit2 <- vcountPattern(REV, sread(readFastq(fnRs_bb[i])), fixed = FALSE)
  prop_bb[i] <- sum(hit1 == 1 & hit2 == 1)/length(hit1)
}
save(prop_ba, prop_bb, file = "keep_reads.Rdata")

fun <- function(x, primer, file2, primer2) {
  allhit1 <- vcountPattern(primer, sread(x), fixed = FALSE)
  allhit2 <- vcountPattern(primer2, sread(readFastq(file2)), fixed = FALSE)
  if(sum(allhit1) < 10){
    return(x)
  }else{
    return(x[allhit1==1 & allhit2==1])
  }
}

outpath_ba <- "/home/chuang19/Ben_data/Brown2018a/fastq_rm_noprimer/"
outpath_bb <- "/home/chuang19/Ben_data/Brown2018b/fastq_rm_noprimer/"

outfile_ba_1 <- paste0(outpath_ba, sname_ba, "_1_f.fastq")
outfile_ba_2 <- paste0(outpath_ba, sname_ba, "_2_f.fastq")
outfile_bb_1 <- paste0(outpath_bb, sname_bb, "_1_f.fastq")
outfile_bb_2 <- paste0(outpath_bb, sname_bb, "_2_f.fastq")


for(i in 1:length(sname_ba)){
  print(paste0(i, " Ba"))
  filterFastq(fnFs_ba[i], outfile_ba_1[i], file2=fnRs_ba[i], primer =FWD, primer2 =REV, filter=fun, compress=F)
  filterFastq(fnRs_ba[i], outfile_ba_2[i], file2=fnFs_ba[i], primer =REV, primer2 =FWD, filter=fun, compress=F)
}

for(i in 1:length(sname_bb)){
  print(paste0(i, " Bb"))
  filterFastq(fnFs_bb[i], outfile_bb_1[i], file2=fnRs_bb[i], primer =FWD, primer2 =REV, filter=fun, compress=F)
  filterFastq(fnRs_bb[i], outfile_bb_2[i], file2=fnFs_bb[i], primer =REV, primer2 =FWD, filter=fun, compress=F)
}

