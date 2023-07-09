#------------------#
# Track ASV
#------------------#

# function used to load object to a variable
load_obj <- function(f,n){
  env <- new.env()
  nm <- load(f, env)[n]
  env[[nm]]
}

path <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/output/dada2/"
sname <- sapply(strsplit(list.files(path,pattern = ".RData"), "_"), `[`,1)

ASV_all    <- sapply(list.files(path, pattern = ".RData",full.names = T), load_obj,1);names(ASV_all)   <- sname
track_all  <- sapply(list.files(path, pattern = ".RData",full.names = T), load_obj,2);names(track_all) <- sname

# For El and Pe they both have 6 columns, we remove column 3 and 4 (denoisedF and denoisedR), and show merge
track_all$El <- track_all$El[,-c(3:4)]
track_all$Pe <- track_all$Pe[,-c(3:4)]

#"input"     "filtered"  "denoisedF" "nonchim"  
prop_nonchim <- sapply(track_all, function(x){x[,4]/x[,1]}) %>% unlist()
prop_denoise <- sapply(track_all, function(x){x[,3]/x[,1]}) %>% unlist()
prop_filter  <- sapply(track_all, function(x){x[,2]/x[,1]}) %>% unlist()

plotdata <- data.frame(nonchim = prop_nonchim, 
                       denoise_merge = prop_denoise,
                       filter  = prop_filter,
                       study = sapply(strsplit(names(prop_nonchim), "[.]"), `[`,1)) %>% 
  gather(key = "step", value = "prop",1:3) %>% 
  mutate(step = factor(step, levels = c("filter", "denoise_merge", "nonchim")))
  

readbox <- ggplot(plotdata, aes(x=study, y=prop, col = step)) + 
  geom_boxplot() + 
  theme_bw()+
  geom_hline(yintercept = 1, linetype = "dashed", col = "red", cex = 1) + 
  ylab("Proportion of Input reads") + xlab("Study") + 
  theme(legend.position = "top")

outpath <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/output/"
ggsave(filename = paste0(outpath, "trackReads.png"), readbox, width = 8, height = 5, dpi = 600)

# check the small proportion sample
# test <- track_all[["Fe"]]
# test[(test[,2]/test[,1]) < 0.5,]

# change rownames: remove _F...
sapply(ASV_all, function(x){rownames(x)[1]})
ASV_all_rename <- sapply(ASV_all, function(x){rownames(x) <- stringr::str_sub(rownames(x), end = -17);x})
sapply(ASV_all_rename, function(x){rownames(x)[1]})
rownames(ASV_all_rename$Su) <- stringr::str_sub(rownames(ASV_all_rename$Su), end = -3)
sapply(ASV_all_rename, function(x){rownames(x)[1]})

# add UAB and STC to the ASV table list
load("/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/output/dada2/Ca/processed.rda")

# check if name is matched
sum(rownames(st) != rownames(df))
UAB <- st[df$Cohort == "UAB", ]
STC <- st[df$Cohort == "Stanford", ]

ASV_all_rename <- append(ASV_all_rename, list(UAB=UAB,STC=STC))

v1study <- c("Ba", "Bb", "Fe", "Ki", "Ro", "St")
v4study <- c("Di", "El", "Pe", "STC", "Su", "Ta", "UAB")

ASV_all_rename <- ASV_all_rename[c(v1study, v4study)]
names(ASV_all_rename) <- c(paste0(v1study, "_v1"), paste0(v4study, "_v4"))

ASV_all <- ASV_all_rename
outpath <- "/Users/huangc9/Documents/Ben research/Publishable/1_DADA2/output/ASV/"
save(ASV_all, file = paste0(outpath, "ASV_all.RData"))

