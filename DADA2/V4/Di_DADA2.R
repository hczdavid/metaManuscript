#-----------------------------------#
# 
# Di DADA2 pipeline
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

FWD <- "ACTCCTACGGGAGGCAGCA"  
REV <- "CCGTCAATTCCTTTGAGTTT" 

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)


# path    <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/Example data/Di/"
# example <- list.files(path, pattern = ".fastq",full.names = T)
# 
# 
# # NO primer included; also double-check by eye
# lapply(example, function(x){
#   rbind(FWDreads = sapply(FWD.orients, primerHits, fn = x),
#         REVreads = sapply(REV.orients, primerHits, fn = x))
# })



# -------------- DADA2 ----------------#

path <- "/home/chuang19/Ben_data/Digiulio2015/fastq_vaginal"

# Forward fastq filenames have format: 1.fastq
fnFs <- sort(list.files(path, pattern=".fastq", full.names = TRUE))

#Extract sample names
bfnFs <- basename(fnFs)
sample.names <- c()
for(i in 1:length(bfnFs)){
  mn <- unlist(strsplit(bfnFs[i],"[.]"))
  if(length(mn) == 3){
    sample.names[i] <- mn[1]
  }else{
    sample.names[i] <- paste0(mn[1],"_rs")
  }
}

# plotQualityProfile(example[1:2]) #trimming after 250bp

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))

#filter and trimming
out <- filterAndTrim(fnFs, filtFs, truncLen=c(400),
                     maxN=0, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE


errF <- learnErrors(filtFs, BAND_SIZE=16, HOMOPOLYMER_GAP_PENALTY=-1, 
                    verbose=TRUE, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(filtFs, err=errF, BAND_SIZE=16, HOMOPOLYMER_GAP_PENALTY=-1, multithread=TRUE)

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

ASV_Di   <- seqtab.nochim
track_Di <- track

outpath <- "/home/chuang19/Ben_data/Digiulio2015/" 
save(ASV_Di, track_Di, file = paste0(outpath,"Di_dada2.RData"))
