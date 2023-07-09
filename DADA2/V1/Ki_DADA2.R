#-----------------------------------#
# 
# Kindinger2017 DADA2 pipeline
# 
# Author: David
# Update: 05/31/2021
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


FWD <- "GAGTTTGATCNTGGCTCA"  
REV <- "GTNTTACNGCGGCKGCTG"  

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

path    <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/Example data/Ki/"
example <- list.files(path, full.names = T)

# NO primer included; also double-check by eye
lapply(example, function(x){
  rbind(FWDreads = sapply(FWD.orients, primerHits, fn = x),
        REVreads = sapply(REV.orients, primerHits, fn = x))
})


# -------------- DADA2 ----------------#

# The full data is stored at BRC clustere "/home/chuang19/Ben_data/Kindinger2017/fastq"
# Here we use the example data to as example. 
# We only download the forward reads for two reasons:
#  (1) Ki study have V1-V3 and the forward reads should already cover V1-V2 with good QC
#  (2) The QS of reverse reads is not stable.

# path <- "/home/chuang19/Ben_data/Kindinger2017/fastq" 



# Forward fastq filenames have format: 1.fastq
fnFs <- sort(list.files(path, pattern="1.fastq", full.names = TRUE))

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# plotQualityProfile(example[1:2]) #trimming after 250bp

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
#filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#filter and trimming
out <- filterAndTrim(fnFs, filtFs, truncLen=c(250),
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

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
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)

ASV_Ki   <- seqtab.nochim
track_Ki <- track

outpath <- "/home/chuang19/Ben_data/Kindinger2017/" 
save(ASV_Ki, track_Ki, file = paste0(outpath,"Ki_dada2.RData"))

