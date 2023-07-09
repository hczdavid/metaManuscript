#-----------------------------------#
# 
# Fettweis2019 DADA2 pipeline
# 
# Author: David
# Update: 05/31/2021
# 
#-----------------------------------#

library(ShortRead)
library(dada2)

# ---------- Check Primer --------------- #
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# 
# FWD1 <- "AGAGTTYGATYMTGGCTYAG"
# FWD2 <- "AGARTTTGATCYTGGTTCAG"
# REV  <- "ATTACCGCGGCTGCTGG"  
# 
# FWD1.orients <- allOrients(FWD1)
# FWD2.orients <- allOrients(FWD2)
# REV.orients <- allOrients(REV)
# 
# path    <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/Example data/Fe/"
# example <- list.files(path, full.names = T)
# 
# # Primer is included; also double-check by eye to make sure the primer is in the begaining.
# lapply(example, function(x){
#   rbind(FWD1reads = sapply(FWD1.orients, primerHits, fn = x),
#         FWD2reads = sapply(FWD2.orients, primerHits, fn = x),
#         REVreads = sapply(REV.orients, primerHits, fn = x))
# })



# -------------- DADA2 ----------------#
# The full data is stored at BRC clustere "/home/chuang19/ncbi/public/fastq_fettwis"
# Here we use the example data to as example.
# We only use the forward reads for two reasons:
#  (1) Fe study have V1-V3 and the forward reads should already cover V1-V2 with good QC
#  (2) The QS of reverse reads is not stable.

# path <- "/home/chuang19/ncbi/public/fastq_fettwis" 
# plotQualityProfile(example[c(3,5)])

# Forward fastq filenames have format: 1.fastq
fnFs <- sort(list.files(path, pattern="1.fastq", full.names = TRUE))

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))

#filter and trimming; use trimleft = 20 to remove the primers
out <- filterAndTrim(fnFs, filtFs, truncLen=c(260),
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,trimLeft = 20,
                     compress=TRUE, multithread=FALSE) 

errF <- learnErrors(filtFs, multithread=TRUE)
# plotErrors(errF, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

#get the ESV table
seqtab <- makeSequenceTable(dadaFs)

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)


ASV_Fe   <- seqtab.nochim
track_Fe <- track

outpath <- "/home/chuang19/ncbi/public/" 
save(ASV_Fe, track_Fe, file = paste0(outpath,"Fe_dada2.RData"))

