#-----------------------------------#
# 
# Brown2018a DADA2 pipeline
# 
# Author: David
# Update: 05/31/2021
# 
#-----------------------------------#


library(ShortRead)
library(dada2)

# ---------- Check Primer --------------- #

# An script "check_primer_Brown.R" is used to check the primers for Ba and Ba study (same primer for
# both studies)
# We found that the first 40 samples (20 forward and 20 reverse) don't inlucde primers
# and all other samples include primers with length 19


# -------------- DADA2 ----------------#

# The full data is stored at BRC clustere "/home/chuang19/Ben_data/Brown2018a_fastq" 
# Here we use the example data to as example. 

path <- "/home/chuang19/Ben_data/Brown2018a/fastq_rm_noprimer_nomatch" 
# path    <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/Example data/Ba/"
# example <- list.files(path, full.names = T)

# Check the Quality Profile for forward and reverse reads
# We can consider use both forwrd and reverse reads based on the quality score
# Trim forward at 240 bp and Reverse at 150 bp
# plotQualityProfile(example[1:2]) 

# Forward fastq filenames have format: 1.fastq
# Reverse fastq filenames have format: 2.fastq
fnFs <- sort(list.files(path, pattern="1_f.fastq", full.names = TRUE))
#fnRs <- sort(list.files(path, pattern="2_f.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
#filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# Filter and trimming
# no trimleft for the first 20 samples
out1 <- filterAndTrim(fnFs[1:20], filtFs[1:20], 
                      truncLen=c(245),
                      maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE) 
# trimleft = 19 for all other samples
out2 <- filterAndTrim(fnFs[-c(1:20)], filtFs[-c(1:20)],  
                      truncLen=c(245), trimLeft = 19,
                      maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=FALSE) 

out <- rbind(out1, out2)

errF <- learnErrors(filtFs, multithread=TRUE)
#errR <- learnErrors(filtRs, multithread=TRUE)

# plotErrors(errF, nominalQ=TRUE)

# Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
#dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#merger
#mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Get the ESV table
seqtab <- makeSequenceTable(dadaFs)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)


getN <- function(x) sum(getUniques(x))
#track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
#colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)

ASV_Ba   <- seqtab.nochim
track_Ba <- track

outpath <- "/home/chuang19/Ben_data/Brown2018a/" 
save(ASV_Ba, track_Ba, file = paste0(outpath,"Ba_dada2.RData"))

