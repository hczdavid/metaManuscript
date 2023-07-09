#-----------------------------------#
# 
# St DADA2 pipeline
# 
# Author: David
# Update: 06/11/2021
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



FWD <- "AGAGTTTGATYMTGGCTCAG"  
REV <- "GWATTACCGCGGCKGCTG"  

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

path <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/Example data/St/fastq/"
example <- list.files(path, pattern = ".fastq",full.names = T)


# there seems primers in the file but it's not.
lapply(example, function(x){
  rbind(FWDreads = sapply(FWD.orients, primerHits, fn = x),
        REVreads = sapply(REV.orients, primerHits, fn = x))
})

# test <- sread(readFastq(example[1]))
# length(test@ranges@width)
# which(vcountPattern(FWD.orients[1], test, fixed = FALSE) == 1)
# myseq <- as.character(test[38])
# grep(myseq, pattern = "AGAGTTTGAT")

# path1 <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/Example data/St/fastq/filtered/"
# example1 <- list.files(path1, pattern = ".fastq",full.names = T)
# lapply(example1, function(x){
#   rbind(FWDreads = sapply(FWD.orients, primerHits, fn = x),
#         REVreads = sapply(REV.orients, primerHits, fn = x))
# })

# -------------- DADA2 ----------------#

# We can directly use fastq file instead the 

# list all forward reads file
path <- "/home/chuang19/Ben_data/Stafford2017/Stafford2017_fastq" 

# Forward fastq filenames have format: 1.fastq
fnFs <- sort(list.files(path, pattern=".fastq", full.names = TRUE))

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <-  sample.names <- sapply(strsplit(basename(fnFs), "[.]"), `[`, 1)

# plotQualityProfile(fnFs[3:4])
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))

# Place filtered files in filtered/ subdirectory
out <- filterAndTrim(fnFs, filtFs, truncLen=350, maxLen=650,
                     maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                     verbose=TRUE, compress=TRUE)

err <- learnErrors(filtFs, BAND_SIZE=16, HOMOPOLYMER_GAP_PENALTY=-1, 
                   multithread=TRUE, verbose=TRUE)
# plotErrors(err, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(filtFs, err=err, multithread=TRUE)

#get the ESV table
seqtab <- makeSequenceTable(dadaFs)

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out,sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input","filtered", "denoised", "nonchim")
rownames(track) <- sample.names
head(track)

ASV_St   <- seqtab.nochim
track_St <- track

outpath <- "/home/chuang19/Ben_data/Stafford2017/" 
save(ASV_St, track_St, file = paste0(outpath,"St_dada2.RData"))

