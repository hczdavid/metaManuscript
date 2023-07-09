#-----------------#
# merge peru data
#-----------------#

library(tidyverse)

load("/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/output/dada2/Pe/Pe_dada2_notmerge.RData")

# expane m col to v col
# m is OTU table
# v is the colname name of new OTU table
exp_matrix <- function(m,v){
  n_name <- v[!v%in%colnames(m)]
  m0 <- matrix(0,nrow = nrow(m),ncol = length(n_name),dimnames = list(rownames(m),n_name))
  allm <- cbind(m,m0)
  allm <- allm[,v]
}


ASVname53 <- lapply(seqtab.nochim53, colnames) %>% unlist() %>% unique()
ASV_53    <- Reduce(rbind, lapply(seqtab.nochim53, exp_matrix, ASVname53))

ASVname35 <- lapply(seqtab.nochim35, colnames) %>% unlist() %>% unique()
ASV_35    <- Reduce(rbind, lapply(seqtab.nochim35, exp_matrix, ASVname35))

# merge ASV_53 and ASV_35
ASVname <- unique(c(ASVname53,ASVname35))
ASV_Pe  <- exp_matrix(ASV_35,ASVname) + exp_matrix(ASV_53,ASVname) 

# merge track data
track_Pe <- trackall_35 + trackall_53

outpath <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/output/dada2/"
save(ASV_Pe, track_Pe, file = paste0(outpath, "Pe_dada2.RData"))


