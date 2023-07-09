

# Check Primer for Ba and Bb
library(ShortRead)

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# get the primer seq from the paper
FWD <- "GAGTTTGATCNTGGCTCAG" 
REV <- "TGCTGCCTCCCGTAGGAGT" 

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# path_ba <- "/home/chuang19/Ben_data/Brown2018a/Brown2018a_fastq" 
# path_bb <- "/home/chuang19/Ben_data/Brown2018b/Brown2018b_fastq" 

path_ba <- "/home/chuang19/Ben_data/Brown2018a/fastq_rm_noprimer_nomatch" 
path_bb <- "/home/chuang19/Ben_data/Brown2018b/fastq_rm_noprimer_nomatch" 


# path_ba <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/Example data/Ba/" 
# path_bb <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/Example data/Bb/" 

fnFs_ba <- sort(list.files(path_ba, pattern="1_f.fastq", full.names = TRUE))
fnRs_ba <- sort(list.files(path_ba, pattern="2_f.fastq", full.names = TRUE))
sname_ba <- sapply(strsplit(basename(fnFs_ba), "_"), `[`, 1)

fnFs_bb <- sort(list.files(path_bb, pattern="1_f.fastq", full.names = TRUE))
fnRs_bb <- sort(list.files(path_bb, pattern="2_f.fastq", full.names = TRUE))
sname_bb <- sapply(strsplit(basename(fnFs_bb), "_"), `[`, 1)

allfile    <- list(fnFs_ba, fnRs_ba,fnFs_bb, fnRs_bb)

# all_hits <- lapply(allfile, function(y){
#   sapply(y, function(x){sapply(FWD.orients, primerHits, fn = x)})
# })

all_hits <- list()
for (i in 1:4) {
  if(i == 1 | i ==3){orient <- FWD.orients}else{orient <- REV.orients}
  all_hits[[i]] <- sapply(allfile[[i]], function(x){sapply(orient, primerHits, fn = x)})
}

# check the read length
allfile    <- list(fnFs_ba, fnRs_ba,fnFs_bb, fnRs_bb)
all_fastq  <- lapply(allfile, function(x){lapply(x, readFastq)})
all_rl     <- lapply(all_fastq, function(x){
  lapply(x, function(y){sread(y)@ranges@width})})
  

all_ltable <- lapply(all_rl, function(y){
  
  sapply(y, function(x){
    
    myt <- c(sum(x <= 246),
      sum(x <= 247),
      sum(x <= 248),
      sum(x <= 249), 
             sum(x == 250),
             sum(x == 251),
             sum(x > 251 & x <= 300),
             sum(x == 301))
    names(myt) <- c("<=246","<=247","<=248","<=249", "250", "251", "252-300", "301")
    myt
  })
  
})


alltable <-  mapply(function(x,y,z){xx <- rbind(x,y); colnames(xx) <- z;t(xx)}, 
                    all_hits, all_ltable, list(sname_ba, sname_ba,sname_bb,sname_bb))
#alltable1 <- lapply(alltable, function(x){prop = round(x[,1]/rowSums(x[,5:9]),2); cbind(x, prop)})

names(alltable) <- c("Ba_FWR", "Ba_REV","Bb_FWR", "Bb_REV")
save(alltable, file = "check_primer_length_Brown_rmNOprimer.Rdata")
