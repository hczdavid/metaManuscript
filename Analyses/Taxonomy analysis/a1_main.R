#-------------------------------------------#
# 
# Meta-analysis of Vaginal microbiome and PTB
#
# Analysis 1: Taxonomic Resolution Analysis
# 
# Author: David Huang
# Update: 06/28/2021
#
#-------------------------------------------#


# Goal: 
# 1. Compare different taxonomy levels (ASV, Genus with Lac.sp, Family, Order, Class and Phylum).
# 2. Compare the common table and top 0.1% table.

# Data setup:
# 1. Average for longitudinal data
# 2. Use the proportial data
# 3. Define PTB based on < 37 GAAD

# Classifier setup
# 1. Use random forest classifier
# 2. Do intra-, cross- and loso-analysis
# 3. 5-fold CV and 20 interations for intra-analysis

# Input: taxa_data.Rdata
# Output: a1_intra.Rdata; a1_cross.Rdata; a1_loso.Rdata

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(randomForest)
library(caret)
library(pROC)
library(glmnet)
library(foreach)
library(doParallel)
library(tidyverse)

load("/Users/huangc9/Documents/Ben research/Publishable/4_prepare_data/output/taxa_data.Rdata")


# get the basic variables
alldata <- taxa_data$prop

metaSub_v1  <- taxa_data$meta[1:5]
metaSub_v4  <- taxa_data$meta[6:12]
sname_v1    <- names(metaSub_v1)
sname_v4    <- names(metaSub_v4)
snum_v1     <- length(sname_v1)
snum_v4     <- length(sname_v4)
taxa_level  <- names(alldata$com_V1$Br)
ln          <- length(taxa_level)

registerDoParallel(detectCores()-1)


# --------------------------------- #
#            Intra Study 
# --------------------------------- #

iter     <- 20
control  <- trainControl(method="repeatedcv", number=5, repeats=5,search = "random")
a1_intra <- list()

for(ww in 1:length(alldata)){ # loop for com_v1/top_v1/com_v4/top_v4 #
  
  print(paste0("Current running the cohort ** ", ww, " **" ))
  
  taxadata<- alldata[[ww]]
  
  # different parameters for v1 and v4
  if(ww==1 | ww==2){
    metaSub <- metaSub_v1; snum <- snum_v1; sname <- sname_v1
  }else{metaSub <- metaSub_v4;snum <- snum_v4; sname <- sname_v4}
  
  
  auctable <- list()
  
  for(jj in 1:snum){ # loop for studies #
    
    print(paste0("Current running the study ** ", sname[jj], " **" ))
    
    # define matrix to save AUC for all iteration and all taxa level
    auctable_rf <- matrix(NA,nrow = iter,ncol = ln,dimnames = list(NULL,taxa_level))
    
    for(kk in 1:ln){ # loop for taxa level #
      
      print(paste0("Current running the taxa ** ", taxa_level[kk], " **" ))
      
      cdata <- taxadata[[jj]][[kk]]
      pheno <- metaSub[[jj]][,15]
      
      y.sample.matrix <- matrix(nrow=length(pheno),ncol=iter)
      y.test.matrix   <- matrix(nrow=length(pheno),ncol=iter)
      rf.pred.matrix  <- matrix(nrow=length(pheno),ncol=iter)
      
      set.seed(100)
      
      # tunning mtry variable
      hyper.model <- train(cdata,as.factor(pheno),method="rf",ntree = 1000,tuneLength=15,trControl=control)
      bestmtry <- as.matrix(hyper.model$bestTune)[1]
      
      result <- foreach(n = 1:iter, .inorder = FALSE)%dopar%{ # start paralell computing for iteration #
        
        flds    <- createFolds(pheno, k = 5, list = TRUE, returnTrain = FALSE)
        y.test  <- NULL
        rf.pred <- NULL
        
        for(i in 1:5){ # start 5-fold CV
          
          y.test <- c(y.test, pheno[flds[[i]]])
          y.train <- pheno[-flds[[i]]]
          x.train <- cdata[-flds[[i]],]
          
          #remove OTUs that have all zeroes in training set
          x.train <- subset(x.train,select=apply(x.train,2,var)>0)
          x.test <- cdata[flds[[i]],]
          
          rf.train <- randomForest(x.train,y=as.factor(y.train), ntree = 1000,mtry = bestmtry)
          rf.pred <- c(rf.pred,predict(rf.train,x.test,type="prob")[,2])
        }
        
        y.sample <- unlist(flds)
        myresult <- list(y.sample,y.test,rf.pred)
      }

      for(n in 1:iter){
        y.sample.matrix[,n]       <- result[[n]][[1]]
        y.test.matrix[,n]         <- result[[n]][[2]]
        rf.pred.matrix[,n]        <- result[[n]][[3]]
      }
      
      #calculate average predicted y by sample
      for (n in 1:iter) {
        y.test.temp <- cbind(y.sample.matrix[,n],y.test.matrix[,n])
        y.test.temp <- y.test.temp[order(y.test.temp[,1]),]
        y.test.matrix[,n] <- y.test.temp[,2]
        
        rf.pred.temp <- cbind(y.sample.matrix[,n],rf.pred.matrix[,n])
        rf.pred.temp <- rf.pred.temp[order(rf.pred.temp[,1]),]
        rf.pred.matrix[,n] <- rf.pred.temp[,2]
      }
      
      auctable_rf[,kk]  <- apply(rf.pred.matrix,2, function(x){roc(pheno,x,quiet=T)$auc})
    }
    
    auctable[[jj]] <- auctable_rf
  }
  a1_intra[[ww]] <- auctable
}
save(a1_intra, file = "a1_intra.Rdata")

# --------------------------------- #
#            Cross Study 
# --------------------------------- #

a1_cross <- foreach(ww = 1:length(alldata), .inorder = FALSE)%dopar%{ # start parallel computing #
  
  taxadata <- alldata[[ww]]
  
  # different parameters for v1 and v4
  if(ww==1 | ww==2){
    metaSub <- metaSub_v1; snum <- snum_v1; sname <- sname_v1
  }else{metaSub <- metaSub_v4;snum <- snum_v4; sname <- sname_v4}

  # list to save all taxa results 
  taxa_cross_result <- list()
  
  for(kk in 1:ln){ # loop for taxa level #
    
    print(paste0("Current running the taxa ** ", taxa_level[kk], " **" ))
    
    # define matrix to save AUC for all cross study
    AUC_rf <- matrix(NA,nrow = snum,ncol = snum,dimnames = list(sname,sname))
    
    for(jj in 1:snum){ # start training loop #
      
      set1       <- taxadata[[jj]][[kk]]
      set1_pheno <- metaSub[[jj]][,15]
      
      # tuning mtry and get the best model for testing
      rf.train   <- train(set1,as.factor(set1_pheno),method="rf",importance=T,ntree = 1000,tuneLength=15,trControl=control)
      
      for(mm in 1:snum){ # start testing loop #
        
        set2       <- taxadata[[mm]][[kk]]
        set2_pheno <- metaSub[[mm]][,15]
        
        rf.pred         <- predict(rf.train,set2,type="prob")
        AUC_rf[jj,mm]   <- roc(set2_pheno,rf.pred[,2],quiet=T)$auc
      }
    }
    taxa_cross_result[[kk]] <- AUC_rf
  }
  return(taxa_cross_result)
}
names(a1_cross) <- names(alldata)
save(a1_cross, file = "a1_cross.Rdata")


# --------------------------------- #
#            LOSO Study 
# --------------------------------- #


a1_loso <- foreach(ww = 1:length(alldata), .inorder = FALSE)%dopar%{ # start parallel computing #
  
    taxadata<- alldata[[ww]]
    
    # different parameters for v1 and v4
    if(ww==1 | ww==2){
      metaSub <- metaSub_v1; snum <- snum_v1; sname <- sname_v1
    }else{metaSub <- metaSub_v4;snum <- snum_v4; sname <- sname_v4}
    
   
    # define matrix to save AUC for each study and taxa
    losoAUC <- matrix(NA,nrow = ln,ncol = snum,dimnames = list(taxa_level,sname))
    
    for(kk in 1:ln){ # loop for taxa level #
      
      print(paste0("Current running the taxa ** ", taxa_level[kk], " **" ))
      
      for(jj in 1:snum){ # loop for the study which is excluded in the training set
        
        
          # get the testing set
          set1       <- taxadata[[jj]][[kk]]
          set1_pheno <- metaSub[[jj]][,15]
          
          # combine other studies and get the training set
          set2       <- Reduce(rbind, lapply(taxadata[-jj], `[[`, kk))
          set2_pheno <- unlist(sapply(metaSub[-jj], `[`, 15))
          
          # traing and tuning mtry
          rf.train <- train(set2,as.factor(set2_pheno),method="rf",ntree = 1000,tuneLength=15,trControl=control)
          
          # testing and get AUC
          rf.pred  <- predict(rf.train,set1,type="prob")
          losoAUC[kk,jj] <- roc(set1_pheno,rf.pred[,2],quiet=T)$auc
          
        }
    }
    
    return(losoAUC)
}
names(a1_loso) <- names(alldata)
save(a1_loso, file = "a1_loso.Rdata")

stopImplicitCluster()



