#-----------------------------------#
# 
# Peru DADA2 pipeline
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

FWD <- "GTGYCAGCMGCCGCGGTAA"  
REV <- "GGACTACNVGGGTWTCTAAT"  

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)


# path    <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/Example data/Pe/Orient35/FWD/"
# path    <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/Example data/Pe/Orient35/REV/"
# path    <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/Example data/Pe/Orient53/FWD/"
# path    <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/Example data/Pe/Orient53/REV/"

# example <- list.files(path, full.names = T)
# # primers in the fastq file: use trimleft to remove
# lapply(example, function(x){
#   rbind(FWDreads = sapply(FWD.orients, primerHits, fn = x),
#         REVreads = sapply(REV.orients, primerHits, fn = x))
# })



# -------------- DADA2 ----------------#

# Every sample is in both orient53 and orient35 though, 
# and it should be possible to process both (independently) 
# and then merge the tables together, summing reads for each sample.

# The samples are from 3 different runs, samples 1-48, 49-96, and 97-137 
# are each from a different run, and thus those groups of samples should be denoised together.
# The primers are still on the reads. They are using the 515F/806R EMP primers, 
# you can look at demux_peru.R to see the exact sequences. They should be removed with `trimLeft`.


#-------------- process orient 53 -----------------#

pathf_35 <- "/home/chuang19/Ben_data/Peru2020/demux/orient35/FWD" 
pathr_35 <- "/home/chuang19/Ben_data/Peru2020/demux/orient35/REV" 
pathf_53 <- "/home/chuang19/Ben_data/Peru2020/demux/orient53/FWD" 
pathr_53 <- "/home/chuang19/Ben_data/Peru2020/demux/orient53/REV" 

# Forward fastq filenames have format: 1.fastq
fnFs_35 <- sort(list.files(pathf_35, pattern=".fastq", full.names = TRUE))
fnRs_35 <- sort(list.files(pathr_35, pattern=".fastq", full.names = TRUE))
fnFs_53 <- sort(list.files(pathf_53, pattern=".fastq", full.names = TRUE))
fnRs_53 <- sort(list.files(pathr_53, pattern=".fastq", full.names = TRUE))

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs_35), "[_]"), `[`, 1)

# plotQualityProfile(example[c(1:2)])
# Place filtered files in filtered/ subdirectory

filtFs_35 <- file.path("/home/chuang19/Ben_data/Peru2020/filtered/FWD35", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs_35 <- file.path("/home/chuang19/Ben_data/Peru2020/filtered/REV35", paste0(sample.names, "_R_filt.fastq.gz"))
filtFs_53 <- file.path("/home/chuang19/Ben_data/Peru2020/filtered/FWD53", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs_53 <- file.path("/home/chuang19/Ben_data/Peru2020/filtered/REV53", paste0(sample.names, "_R_filt.fastq.gz"))

# filter and trimming
out_35 <- filterAndTrim(fnFs_35, filtFs_35, fnRs_35, filtRs_35, trimLeft = c(19,20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

out_53 <- filterAndTrim(fnFs_53, filtFs_53, fnRs_53, filtRs_53, trimLeft = c(19,20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
trackall_35 <- c()
trackall_53 <- c()

#three runs 1-48, 49-96, and 97-137 
grp <- list(1:48,49:96,97:137)

seqtab.nochim35 <- c()
seqtab.nochim53 <- c()

for(i in 1:3){
  
  nfiltFs_35 <- filtFs_35[grp[[i]]]
  nfiltRs_35 <- filtRs_35[grp[[i]]]
  nfiltFs_53 <- filtFs_53[grp[[i]]]
  nfiltRs_53 <- filtRs_53[grp[[i]]]

  errF_35 <- learnErrors(nfiltFs_35, multithread=TRUE)
  errR_35 <- learnErrors(nfiltRs_35, multithread=TRUE)
  errF_53 <- learnErrors(nfiltFs_53, multithread=TRUE)
  errR_53 <- learnErrors(nfiltRs_53, multithread=TRUE)
  
  # plotErrors(errF_53, nominalQ=TRUE)
  
  # Sample Inference
  dadaFs_35 <- dada(nfiltFs_35, err=errF_35, multithread=TRUE)
  dadaRs_35 <- dada(nfiltRs_35, err=errR_35, multithread=TRUE)
  dadaFs_53 <- dada(nfiltFs_53, err=errF_53, multithread=TRUE)
  dadaRs_53 <- dada(nfiltRs_53, err=errR_53, multithread=TRUE)
  
  mergers_35 <- mergePairs(dadaFs_35, nfiltFs_35, dadaRs_35, nfiltRs_35, verbose=TRUE)
  mergers_53 <- mergePairs(dadaFs_53, nfiltFs_53, dadaRs_53, nfiltRs_53, verbose=TRUE)
  
  seqtab_35 <- makeSequenceTable(mergers_35)
  seqtab_53 <- makeSequenceTable(mergers_53)
  
  seqtab.nochim_35 <- removeBimeraDenovo(seqtab_35, method="consensus", multithread=TRUE, verbose=TRUE)
  seqtab.nochim_53 <- removeBimeraDenovo(seqtab_53, method="consensus", multithread=TRUE, verbose=TRUE)
  
  track_35 <- cbind(out_35[grp[[i]],], sapply(dadaFs_35, getN), sapply(dadaRs_35, getN), sapply(mergers_35, getN), rowSums(seqtab.nochim_35))# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
  track_53 <- cbind(out_53[grp[[i]],], sapply(dadaFs_53, getN), sapply(dadaRs_53, getN), sapply(mergers_53, getN), rowSums(seqtab.nochim_53))# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
  
  colnames(track_35) <- c("input","filtered", "denoisedF", "denoisedR", "merge", "nonchim")
  rownames(track_35) <- sample.names[grp[[i]],]
  colnames(track_53) <- c("input","filtered", "denoisedF", "denoisedR", "merge", "nonchim")
  rownames(track_53) <- sample.names[grp[[i]],]
  
  trackall_35 <- rbind(trackall_35,track_35)
  trackall_53 <- rbind(trackall_53,track_53)
  
  seqtab.nochim35[[i]] <- seqtab.nochim_35
  seqtab.nochim53[[i]] <- seqtab.nochim_53
  
}

outpath <- "/home/chuang19/Ben_data/Peru2020/" 
save(seqtab.nochim35, seqtab.nochim53, trackall_35, trackall_53, file = paste0(outpath,"Pe_dada2.RData"))





