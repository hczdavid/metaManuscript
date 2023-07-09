#-----------------------------------#
# 
# Su DADA2 pipeline
# 
# Author: David
# Update: 06/12/2021
# 
#-----------------------------------#


library(ShortRead)
library(dada2)

# ---------- Check Primer --------------- #
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


FWD <- "GTGCCAGCMGCCGCGGTAA"  
REV <- "GGACTACHVGGGTWTCTAAT"  

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# path    <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/Example data/Su//"
# example <- list.files(path, pattern = ".fastq",full.names = T)
# 
# # NO primer included; also double-check by eye
# lapply(example, function(x){
#   rbind(FWDreads = sapply(FWD.orients, primerHits, fn = x),
#         REVreads = sapply(REV.orients, primerHits, fn = x))
# })

# test <- sread(readFastq(example[3]))
# length(test@ranges@width)
# allhit <- which(vcountPattern("ATTAGAWACCCBDGTAGTCC", test, fixed = FALSE)==1)
# allhit_seq <- as.character(test[allhit])
# sapply(allhit_seq, function(x)grepRaw("ATTAGA", x))


# -------------- DADA2 ----------------#

# only use the forward reads since the reverse reads is not stable and forward is good for truncate window.
path <- "/home/chuang19/Ben_data/Subramaniam2016/fastq" 

# Forward fastq filenames have format: 1.fastq
fnFs <- sort(list.files(path, pattern="1.fastq", full.names = TRUE))
#fnRs <- sort(list.files(path, pattern="2.fastq", full.names = TRUE))

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <-  sample.names <- sapply(strsplit(basename(fnFs), "[.]"), `[`, 1)

#plotQualityProfile(fnFs[1:2])

#Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
#filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#filter and trimming: all reads are 251 length
out <- filterAndTrim(fnFs, filtFs, truncLen=c(250),
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
out
errF <- learnErrors(filtFs, multithread=TRUE)
#errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

#get the ESV table
seqtab <- makeSequenceTable(dadaFs)
#table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out,sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input","filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)

ASV_Su   <- seqtab.nochim
track_Su <- track

outpath <- "/home/chuang19/Ben_data/Subramaniam2016/" 
save(ASV_Su, track_Su, file = paste0(outpath,"Su_dada2.RData"))
