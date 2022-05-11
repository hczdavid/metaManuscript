#--------------------------------------------------#
# 
# Meta-analysis of Vaginal microbiome and PTB
#
# Analysis 3: Data Transformation and Classifiers
# 
# Author: David Huang
# Update: 06/28/2021
#
#--------------------------------------------------#

# Goal: 
# 1. Compare different definitions of PTB


# Data setup:
# 1. Average for longitudinal data
# 2. Use the CLR transformed data

# Classifier setup
# 1. Use random forest classifier
# 2. Do intra-, cross- and loso-analysis
# 3. 5-fold CV and 20 interations for intra-analysis

# Input: genus_data.Rdata
# Output: a3_intra.Rdata; a3_cross.Rdata; a3_loso.Rdata

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(randomForest)
library(caret)
library(pROC)
library(glmnet)
library(foreach)
library(doParallel)
library(tidyverse)


load("/Users/huangc9/Documents/Ben research/Publishable/4_prepare_data/output/genus_data.Rdata")

# get the basic variables
metaSuball  <- genus_data$meta
genus_clr   <- genus_data$trans$clr
genus_clr   <- genus_data$trans$log


snum  <-  length(genus_clr)
sname <-  names(genus_clr)

registerDoParallel(detectCores()-1)


# --------------------------------- #
#   Obtain the data based on GAAD 
# --------------------------------- #

genus_clr_l32   <- list()
genus_clr_l34   <- list()
genus_clr_34_37 <- list()
genus_clr_g39   <- list()

for(i in 1:snum){
  genus_clr_l32[[i]]   <- genus_clr[[i]][metaSuball[[i]][6] <  32 & !is.na(metaSuball[[i]][6]),,drop=F]
  genus_clr_l34[[i]]   <- genus_clr[[i]][metaSuball[[i]][6] <  34 & !is.na(metaSuball[[i]][6]),,drop=F]
  genus_clr_34_37[[i]] <- genus_clr[[i]][metaSuball[[i]][6] >= 34 & metaSuball[[i]][6] < 37 & !is.na(metaSuball[[i]][6]),,drop=F]
  genus_clr_g39[[i]]   <- genus_clr[[i]][metaSuball[[i]][6] >= 39 & !is.na(metaSuball[[i]][6]),,drop=F]
}

samplenumber <- rbind(sapply(genus_clr_l32, nrow), sapply(genus_clr_l34, nrow), sapply(genus_clr_34_37, nrow))
colnames(samplenumber) <- sname
# write.csv(samplenumber, file = "mysample.csv")

# exclude Ro, Di, Pe, STC, Su
sindex          <- c(1:3, 5, 7, 11, 12)
fsamplenumber   <- samplenumber[,sindex]
fstudy          <- colnames(fsamplenumber)
minsample       <- apply(fsamplenumber, 2, min)
minsample_which <- apply(fsamplenumber, 2, which.min)

genus_clr_l32   <- genus_clr_l32[sindex]
genus_clr_l34   <- genus_clr_l34[sindex]
genus_clr_34_37 <- genus_clr_34_37[sindex]
genus_clr_g39   <- genus_clr_g39[sindex]



#x1 <- rowMeans(genus_clr_g39[[1]])
# genus_clr_g39[[1]][,"Bifidobacterium"]
# # plot the CLR data
# 
# allclr <- mapply(function(x, y){
#   rbind(x, y) %>% as.data.frame() %>%  mutate(birth = c(rep("PTB", nrow(x)), rep("TB", nrow(y))))
# }, genus_clr_l32, genus_clr_g39, SIMPLIFY = F)
# 
# 
# 
# allclr_plot <- lapply(allclr, function(x){
#   x %>%  gather(key = "genus", value = "value", 1:25) %>% mutate(genus = factor(genus, levels = unique(genus)))
# })
# 
# box1 <- ggplot(allclr_plot[[1]], aes(x = birth, y = value)) +
#   geom_boxplot(width = 0.5) + 
#   geom_jitter(width = 0.3, size = 1, alpha = 0.5) +
#   facet_wrap(~genus, nrow = 5, scales = "free_y") + theme_bw()
# 
# ggsave(box1,filename = "box1_clr4.png",width = 12, height = 8, units = "in", dpi = 600)

# allrank <- list()
# for(jj in 1:10){ # repeat 10 times
#   set.seed(jj)
#   decrea_accu <- matrix(NA, nrow = 7, ncol = 25, dimnames = list(fstudy, genusname))
#   
#   for(i in 1:7){
#     
#     cdata <- rbind(genus_clr_l32[[i]], genus_clr_g39[[i]])
#     pheno <- c(rep(1, nrow(genus_clr_l32[[i]])), rep(0, nrow(genus_clr_g39[[i]])))
#     
#     myrf <- randomForest(cdata,y=as.factor(pheno) ,ntree = 1000, importance = TRUE)
#     decrea_accu[i,] <- importance(myrf, type = 1, scale = F)
#   }
#   decrea_accu_rank <- apply(decrea_accu*-1, 1, rank ) %>% t
#   allrank[[jj]] <- decrea_accu_rank
# }
# rank_averg <- Reduce("+", allrank)/10
# 
# 
# allrank <- rank_averg[,rev(genusorder)] %>% as.data.frame() %>% mutate(study = fstudy) %>% 
#   gather(key = "genus", value = "rank", 1:25) %>% mutate(genus = factor(genus, levels = c(genusorder))) %>% 
#   mutate(study = factor(study, levels = fstudy))
# 
# 
# 
# plotrank <- ggplot(allrank, aes(x=study , y= genus, fill = rank)) + 
#   geom_tile(color = "gray")+ # can adjust height and width
#   scale_fill_gradientn(colours = rev(heat.colors(40)[-c(1:5)]), limits = c(1,25),guide = guide_colourbar(order = 1)) +
#   geom_text(aes(label = round(rank)), size = 1.5, col="blue") + 
#   theme_bw() + labs(fill = "Rank")+
#   ylab("Dataset") + xlab("Genus") 
# 
# ggsave(plotrank,filename = "a5_rank37prop.png",width = 6, height = 4, units = "in", dpi = 600)
# 



# --------------------------------- #
#            Intra Study 
# --------------------------------- #

control      <- trainControl(method="repeatedcv", number=5, repeats=5,search = "random", allowParallel = FALSE)
iter         <- 20

a3_intra     <- c()

for(kk in 1:3){ # loop for < 32, < 34 and 34 - 37 #

  myintra <- list()

  for(iii in 1:10){ # 10 repeats loop #

  print(paste("Iteration", iii))

  genus_clr_l32_new   <- list()
  genus_clr_l34_new   <- list()
  genus_clr_34_37_new <- list()

  # get the equal sample for different PTB group
  for(i in 1:length(fstudy)){ 
    set.seed(iii*100)
    if(minsample_which[i] == 1){
      genus_clr_l32_new[[i]] <- genus_clr_l32[[i]]
      genus_clr_l34_new[[i]] <- genus_clr_l34[[i]][sample(fsamplenumber[2, i], fsamplenumber[1, i]),]
      genus_clr_34_37_new[[i]] <- genus_clr_34_37[[i]][sample(fsamplenumber[3, i], fsamplenumber[1, i]),]
    }else{
      genus_clr_l32_new[[i]] <- genus_clr_l32[[i]][sample(fsamplenumber[1, i], fsamplenumber[3, i]),]
      genus_clr_l34_new[[i]] <- genus_clr_l34[[i]][sample(fsamplenumber[2, i], fsamplenumber[3, i]),]
      genus_clr_34_37_new[[i]] <- genus_clr_34_37[[i]]
    }
  }

  # define the matrix to save AUC for each study and iteration
  auctable_rf <- matrix(NA,nrow = iter,ncol = length(fstudy),dimnames = list(NULL,fstudy))

  for(jj in 1:length(fstudy)){ # loop for 7 studies #

    print(paste("Study", jj))

    if(kk == 1){
      genus_clr_new <- genus_clr_l32_new
    }else if(kk == 2){
      genus_clr_new <- genus_clr_l34_new
    }else{
      genus_clr_new <- genus_clr_34_37_new
    }

    cdata <- rbind(genus_clr_new[[jj]], genus_clr_g39[[jj]])
    pheno <- c(rep(1, nrow(genus_clr_new[[jj]])), rep(0, nrow(genus_clr_g39[[jj]])))

    y.sample.matrix <- matrix(nrow=length(pheno),ncol=iter)
    y.test.matrix   <- matrix(nrow=length(pheno),ncol=iter)

    rf.pred.matrix  <- matrix(nrow=length(pheno),ncol=iter)
    set.seed(100)

    hyper.model <- train(cdata,as.factor(pheno),method="rf",ntree = 1000,tuneLength=15,trControl=control)
    bestmtry <- as.matrix(hyper.model$bestTune)[1]

    result <- foreach(n = 1:iter, .inorder = FALSE)%dopar%{ # parallel computing for 20 iterations #

      flds    <- createFolds(pheno, k = 5, list = TRUE, returnTrain = FALSE)
      y.test  <- NULL
      rf.pred <- NULL

      for (i in 1:5) {
    
        y.test <- c(y.test, pheno[flds[[i]]])
        y.train <- pheno[-flds[[i]]]
        x.train <- cdata[-flds[[i]],]
        if(sum(colSums(x.train)==0)!=0){x.train[,which(colSums(x.train)==0)] <- 1e-6}
        x.test <- cdata[flds[[i]],]

        rf.train <- randomForest(x.train,y=as.factor(y.train) ,ntree = 1000, mtry = bestmtry)
        rf.pred <- c(rf.pred,predict(rf.train,x.test,type="prob")[,2])
      }

      y.sample <- unlist(flds)
      myresult <- list(y.sample,y.test,rf.pred)

      return(myresult)
    }

    for(n in 1:iter){
      y.sample.matrix[,n]       <- result[[n]][[1]]
      y.test.matrix[,n]         <- result[[n]][[2]]
      rf.pred.matrix[,n]        <- result[[n]][[3]]
    }

    for (n in 1:iter) {
      y.test.temp <- cbind(y.sample.matrix[,n],y.test.matrix[,n])
      y.test.temp <- y.test.temp[order(y.test.temp[,1]),]
      y.test.matrix[,n] <- y.test.temp[,2]

      rf.pred.temp <- cbind(y.sample.matrix[,n],rf.pred.matrix[,n])
      rf.pred.temp <- rf.pred.temp[order(rf.pred.temp[,1]),]
      rf.pred.matrix[,n] <- rf.pred.temp[,2]
    }
    
    auctable_rf[,jj]  <- apply(rf.pred.matrix,2,FUN = function(x){roc(pheno,x, quiet = T)$auc})
    
  }
  
  myintra[[iii]]  <- list(rf = auctable_rf)
  
  }
  
  a3_intra[[kk]] <- myintra
}

save(a3_intra, file = "a3_intra.Rdata")



# --------------------------------- #
#            Cross Study 
# --------------------------------- #

a3_cross <- c()

for(kk in 1:3){
  
  result <- foreach(iii = 1:10, .inorder = FALSE)%dopar%{ # loop for < 32, < 34 and 34 - 37 #
    
    genus_clr_l32_new   <- list()
    genus_clr_l34_new   <- list()
    genus_clr_34_37_new <- list()
    
    for(i in 1:length(fstudy)){
      set.seed(iii*100)
      if(minsample_which[i] == 1){
        genus_clr_l32_new[[i]] <- genus_clr_l32[[i]]
        genus_clr_l34_new[[i]] <- genus_clr_l34[[i]][sample(fsamplenumber[2, i], fsamplenumber[1, i]),]
        genus_clr_34_37_new[[i]] <- genus_clr_34_37[[i]][sample(fsamplenumber[3, i], fsamplenumber[1, i]),]
      }else{
        genus_clr_l32_new[[i]] <- genus_clr_l32[[i]][sample(fsamplenumber[1, i], fsamplenumber[3, i]),]
        genus_clr_l34_new[[i]] <- genus_clr_l34[[i]][sample(fsamplenumber[2, i], fsamplenumber[3, i]),]
        genus_clr_34_37_new[[i]] <- genus_clr_34_37[[i]]
      }
    }
    
    # define the matrix to save AUC for each study 
    AUC_rf       <- matrix(NA,nrow = length(fstudy),ncol = length(fstudy),dimnames = list(fstudy,fstudy))
    
    for(jj in 1:length(fstudy)){ # loop for training study #
      
      print(paste("training study ", jj))
      
      if(kk == 1){
        genus_clr_new <- genus_clr_l32_new
      }else if(kk == 2){
        genus_clr_new <- genus_clr_l34_new
      }else{
        genus_clr_new <- genus_clr_34_37_new
      }
      
      set1       <- rbind(genus_clr_new[[jj]], genus_clr_g39[[jj]])
      set1_pheno <- c(rep(1, nrow(genus_clr_new[[jj]])), rep(0, nrow(genus_clr_g39[[jj]])))
      
      rf.train  <- train(set1,as.factor(set1_pheno),method="rf",importance=T,ntree = 1000,tuneLength=15,trControl=control)
      
      for(mm in 1:length(fstudy)){  # loop for testing study #
        
        set2 <- rbind(genus_clr_new[[mm]], genus_clr_g39[[mm]])
        set2_pheno <- c(rep(1, nrow(genus_clr_new[[mm]])), rep(0, nrow(genus_clr_g39[[mm]])))
        
        rf.pred             <- predict(rf.train,set2,type="prob")
        AUC_rf[jj,mm]       <- roc(set2_pheno,rf.pred[,2], quiet = T)$auc
      }
      
    }
    
    return(list(AUC_rf))
    
  }
  
  a3_cross[[kk]] <- result
  
}

save(a3_cross,file = "a3_cross.Rdata")


# --------------------------------- #
#            LOSO Study 
# --------------------------------- #

a3_loso <- c()

for(kk in 1:3){ 
  
  result <- foreach(iii = 1:10, .inorder = FALSE)%dopar%{ # parallel computing for 10 repeats #
    
    genus_clr_l32_new   <- list()
    genus_clr_l34_new   <- list()
    genus_clr_34_37_new <- list()
    
    for(i in 1:length(fstudy)){
      set.seed(iii*100)
      if(minsample_which[i] == 1){
        genus_clr_l32_new[[i]] <- genus_clr_l32[[i]]
        genus_clr_l34_new[[i]] <- genus_clr_l34[[i]][sample(fsamplenumber[2, i], fsamplenumber[1, i]),]
        genus_clr_34_37_new[[i]] <- genus_clr_34_37[[i]][sample(fsamplenumber[3, i], fsamplenumber[1, i]),]
      }else{
        genus_clr_l32_new[[i]] <- genus_clr_l32[[i]][sample(fsamplenumber[1, i], fsamplenumber[3, i]),]
        genus_clr_l34_new[[i]] <- genus_clr_l34[[i]][sample(fsamplenumber[2, i], fsamplenumber[3, i]),]
        genus_clr_34_37_new[[i]] <- genus_clr_34_37[[i]]
      }
    }
    
    AUC_rf <- rep(NA,length(fstudy)); names(AUC_rf) <- fstudy
    
    for(jj in 1:length(fstudy)){ # loop for excluded study #
      
      print(paste("Exclued Study ", jj))
      
      if(kk == 1){
        genus_clr_new <- genus_clr_l32_new
      }else if(kk == 2){
        genus_clr_new <- genus_clr_l34_new
      }else{
        genus_clr_new <- genus_clr_34_37_new
      }
      
      # testing set
      set1       <- rbind(genus_clr_new[[jj]], genus_clr_g39[[jj]])
      set1_pheno <- c(rep(1, nrow(genus_clr_new[[jj]])), rep(0, nrow(genus_clr_g39[[jj]])))
      
      # training set
      set2_ptb     <- Reduce(rbind, genus_clr_new[-jj])
      set2_tb      <- Reduce(rbind, genus_clr_g39[-jj])
      set2         <- rbind(set2_ptb, set2_tb)
      set2_pheno   <- c(rep(1, nrow(set2_ptb)), rep(0, nrow(set2_tb)))
      
      
      rf.train        <- train(set2,as.factor(set2_pheno),method="rf",ntree = 1000,tuneLength=15,trControl=control)
      rf.pred         <- predict(rf.train,set1,type="prob")
      AUC_rf[jj]      <- roc(set1_pheno,rf.pred[,2], quiet = T)$auc
      
    }
    
    return(list(AUC_rf))
    
  }
  
  a3_loso[[kk]] <- result
}

save(a3_loso,file = "a3_loso3.Rdata")
stopImplicitCluster() # end parallel
