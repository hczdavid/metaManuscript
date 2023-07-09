#-----------------------------------#
# 
# El DADA2 pipeline
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


FWD <- "ACTCCTRCGGGAGGCAGCAG"  
REV <- "GGACTACHVGGGTWTCTAAT"  

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# path    <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/Example data/El/fastq//"
# example <- list.files(path, full.names = T)
# 
# # the files include primers
# lapply(example, function(x){
#   rbind(FWDreads = sapply(FWD.orients, primerHits, fn = x),
#         REVreads = sapply(REV.orients, primerHits, fn = x))
# })

# use cutadapt the remove the primers
# ! /bin/bash
# outdir="/home/chuang19/ncbi/dbGaP-24369/Elovitz2019/trimfastq/"
# inpdir="/home/chuang19/ncbi/dbGaP-24369/Elovitz2019/fastq/"
# names=$(</home/chuang19/ncbi/dbGaP-24369/Elovitz2019/SRR_Acc_List.txt)
# fastq1="_1.fastq"
# fastq2="_2.fastq"
# fp="ACTCCTRCGGGAGGCAGCAG"
# rp="GGACTACHVGGGTWTCTAAT"
# 
# for name in $names
# do
# cutadapt -g  $fp -o $outdir$name$fastq1  $inpdir$name$fastq1
# cutadapt -g  $rp -o $outdir$name$fastq2  $inpdir$name$fastq2
# done
# echo All done


# path1    <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/Example data/El/fastq_noprimer/"
# example1 <- list.files(path1, full.names = T)
# 
# # NO primer included; also double-check by eye
# lapply(example1, function(x){
#   rbind(FWDreads = sapply(FWD.orients, primerHits, fn = x),
#         REVreads = sapply(REV.orients, primerHits, fn = x))
# })


# -------------- DADA2 ----------------#

path <- "/home/chuang19/ncbi/dbGaP-24369/Elovitz2019/trimfastq/" 



# Forward fastq filenames have format: 1.fastq
fnFs <- sort(list.files(path, pattern="1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="2.fastq", full.names = TRUE))

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "[_]"), `[`, 1)

# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])

#Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#filter and trimming
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,270),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
# plotErrors(errF, nominalQ=TRUE)

# Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#merger
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample

#get the ESV table
seqtab <- makeSequenceTable(mergers)

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input","filtered", "denoisedF", "denoisedR", "merge", "nonchim")
rownames(track) <- sample.names
head(track)

ASV_El   <- seqtab.nochim
track_El <- track

outpath <- "/home/chuang19/ncbi/dbGaP-24369/Elovitz2019/" 
save(ASV_El, track_El, file = paste0(outpath,"El_dada2.RData"))



