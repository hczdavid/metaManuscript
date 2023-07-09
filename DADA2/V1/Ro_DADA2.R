#-----------------------------------#
# 
# Ro DADA2 pipeline
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



FWD1 <- "AGAGTTTGATYMTGGCTCAG"  
FWD2 <- "AGGGTTCGATTCTGGCTCAG"  
FWD3 <- "AGAGTTTGATCCTGGCTTAG"  
FWD4 <- "AGAATTTGATCTTGGTTCAG"  
REV <- "TTACCGCGGCTGCTGGCA"  

FWD1.orients <- allOrients(FWD1)
FWD2.orients <- allOrients(FWD2)
FWD3.orients <- allOrients(FWD3)
FWD4.orients <- allOrients(FWD4)
REV.orients  <- allOrients(REV)

path <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/Example data/Ro/fastq/" 
example <- list.files(path, full.names = T)

# there are rev primer in the file and we need to remove it use cutadapt.
lapply(example, function(x){
  rbind(FWDreads1 = sapply(FWD1.orients, primerHits, fn = x),
        FWDreads2 = sapply(FWD2.orients, primerHits, fn = x),
        FWDreads3 = sapply(FWD3.orients, primerHits, fn = x),
        FWDreads4 = sapply(FWD4.orients, primerHits, fn = x),
        REVreads  = sapply(REV.orients, primerHits, fn = x))
})


# remove the primer
#--------------------------------------------------#
# ! /bin/bash
# outdir="/home/chuang19/Ben_data/Romero2014_trim/"
# inpdir="/home/chuang19/Ben_data/Romero2014_fastq/"
# names=$(</home/chuang19/Ben_data/run1.txt)
# fastq=".fastq"
# fp="TTACCGCGGCTGCTGGCA"
# 
# for name in $names
# do
# cutadapt -g  $fp -o $outdir$name$fastq  $inpdir$name$fastq
# done
# echo All done
#--------------------------------------------------#

# file after remove the primer; sanity check
path1 <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/Example data/Ro/fastq_noprimer/" 
example1 <- list.files(path1, full.names = T)
lapply(example1, function(x){
  rbind(FWDreads1 = sapply(FWD1.orients, primerHits, fn = x),
        FWDreads2 = sapply(FWD2.orients, primerHits, fn = x),
        FWDreads3 = sapply(FWD3.orients, primerHits, fn = x),
        FWDreads4 = sapply(FWD4.orients, primerHits, fn = x),
        REVreads  = sapply(REV.orients, primerHits, fn = x))
})

# test <- sread(readFastq(example[1]))
# test1 <- sread(readFastq(example1[1]))
# length(test@ranges@width)
# length(test1@ranges@width)
# seq_length <- cbind(test@ranges@width, test1@ranges@width)
# as.character(test1[1])
# as.character(test[1])
# grep(as.character(test[1]), pattern = "TTACCGCGGCTGCTGGCA")



# -------------- DADA2 ----------------#

# We can use fastq file after trimming primer 
path <- "/home/chuang19/Ben_data/Romero2014/Romero2014_trim" 

# Forward fastq filenames have format: 1.fastq
fnFs <- sort(list.files(path, pattern=".fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <-  sample.names <- sapply(strsplit(basename(fnFs), "[.]"), `[`, 1)

# plotQualityProfile(fnFs[1:2]) #
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))

# filter and trimming, we need to reverse the ASV later; we tried 450 and 350, after explore the truncate window
# we decide to use 450 cutoff.
out <- filterAndTrim(fnFs, filtFs, truncLen=c(450),
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

errF <- learnErrors(filtFs, multithread=TRUE, BAND_SIZE=16, HOMOPOLYMER_GAP_PENALTY=-1,verbose=TRUE)
# plotErrors(errF, nominalQ=TRUE)

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

#-------------------------#
# Reverse the ESV sequence
#-------------------------#

ESVname <- colnames(seqtab.nochim)
newEVSname <- c()
for(i in 1:length(ESVname)){newEVSname[i] <- allOrients(ESVname[i])[4]}
colnames(seqtab.nochim) <- newEVSname

ASV_Ro   <- seqtab.nochim
track_Ro <- track

outpath <- "/home/chuang19/Ben_data/Romero2014/" 
save(ASV_Ro, track_Ro, file = paste0(outpath,"Ro_dada2.RData"))


