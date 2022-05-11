#--------------------------------------------------#
# 
# Meta-analysis of Vaginal micribome and PTB
#
# Analysis 2: Data Transformation and Classifiers
# 
# Author: David Huang
# Update: 06/28/2021
#
#--------------------------------------------------#

# Goal: 
# 1. Compare classifiers: Random Forest, Logistic regression, LASSO, Ridge and Elastic.
# 2. Compare data transformation methods: proportion, CLR, ALR, Log, Rank.
# 3. Compare if standardize the each genus to mean 0 and variance 1


# Data setup:
# 1. Avarage for longitudinal data
# 2. Define PTB based on < 37 GAAD
# 3. Have both scaled and unscaled data for each transformation

# Classifier setup
# 1. Tuning hyperparameter of RF (mtry) and LASSO (lambda)
# 2. Do intra-, cross- and loso-analysis
# 3. 5-fold CV and 20 interations for intra-analysis
# 4. For glmnet, set standardize = FALSE

# Input: genus_data.Rdata
# Output: a2_intra.Rdata; a2_cross.Rdata; a2_loso.Rdata

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
metaSuball <- genus_data$meta
comgenus   <- genus_data$trans

ntrans    <-  length(comgenus)
nametrans <-  names(comgenus)
snum      <-  length(comgenus$prop)
sname     <-  names(comgenus$prop)

registerDoParallel(detectCores()-1)

# --------------------------------- #
#            Intra Study 
# --------------------------------- #

control  <- trainControl(method="repeatedcv", number=5, repeats=5,search = "random")
iter     <- 20

a2_intra <- c()
for(ww in 1:ntrans){ # transformation loop #
  
  # define matrices to store AUC for five classifiers
  auctable_rf    <- matrix(NA,nrow = iter, ncol = snum, dimnames = list(NULL,sname))
  auctable_lasso <- auctable_logistic <- auctable_ridge <- auctable_enet <- auctable_rf
 
  for(jj in 1:snum){ # study loop #
    
    print(paste("Transformation ww = ", ww, "Study jj = ", jj)) 
    
    cdata <- comgenus[[ww]][[jj]]   
    pheno <- metaSuball[[jj]][,15]
    
    y.sample.matrix <- matrix(nrow=length(pheno),ncol=iter) 
    y.test.matrix   <- matrix(nrow=length(pheno),ncol=iter) 
    
    rf.pred.matrix    <- matrix(nrow=length(pheno),ncol=iter)
    lasso.pred.matrix <- logistic.pred.matrix <- ridge.pred.matrix <- enet.pred.matrix <- rf.pred.matrix

    # tuning mtry
    set.seed(100)
    hyper.model <- train(cdata,as.factor(pheno),method="rf",ntree = 1000,tuneLength=15,trControl=control)
    bestmtry <- as.matrix(hyper.model$bestTune)[1]
    
   
    result <- foreach(n = 1:iter, .inorder = FALSE)%dopar%{  # start parallel computing loop #
      
      set.seed(n)
      flds          <- createFolds(pheno, k = 5, list = TRUE, returnTrain = FALSE)
      y.test        <- NULL; rf.pred    <- NULL; lasso.pred <- NULL
      logistic.pred <- NULL; ridge.pred <- NULL; enet.pred  <- NULL
      
      for(i in 1:5){ # start 5-fold CV #
        y.test  <- c(y.test, pheno[flds[[i]]])
        y.train <- pheno[-flds[[i]]]
        x.train <- cdata[-flds[[i]],]
        if(sum(colSums(x.train)==0)!=0){x.train[,which(colSums(x.train)==0)] <- 1e-6}
        x.test  <- cdata[flds[[i]],]
        
        rf.train <- randomForest(x.train,y=as.factor(y.train) ,ntree = 1000,mtry = bestmtry)
        rf.pred  <- c(rf.pred,predict(rf.train,x.test,type="prob")[,2])
        
        cv.lasso.train <- cv.glmnet(x.train,as.factor(y.train), family = "binomial", alpha = 1,standardize = FALSE)
        lasso.pred     <- c(lasso.pred, predict(cv.lasso.train,x.test,s = "lambda.min",type = "response"))
        
        logistic.train <- train(x.train,as.factor(y.train),method="glm")
        logistic.pred  <- c(logistic.pred,predict(logistic.train,x.test,type="prob")[,2])
        
        cv.ridge.train <- cv.glmnet(x.train,as.factor(y.train), family = "binomial", alpha = 0,standardize = FALSE)
        ridge.pred     <- c(ridge.pred, predict(cv.ridge.train,x.test,s = "lambda.min",type = "response"))
        
        enet.train <- train(x.train,as.factor(y.train),method="glmnet",tuneLength=15,trControl=control,standardize = FALSE)
        enet.pred  <- c(enet.pred,predict(enet.train,x.test,type="prob")[,2])
      }
      
      y.sample <- unlist(flds)
      myresult <- list(y.sample,y.test,rf.pred,lasso.pred,logistic.pred,ridge.pred,enet.pred)
      
      return(myresult)
    }

    for(n in 1:iter){
      y.sample.matrix[,n]       <- result[[n]][[1]]
      y.test.matrix[,n]         <- result[[n]][[2]]
      rf.pred.matrix[,n]        <- result[[n]][[3]]
      lasso.pred.matrix[,n]     <- result[[n]][[4]]
      logistic.pred.matrix[,n]  <- result[[n]][[5]]
      ridge.pred.matrix[,n]     <- result[[n]][[6]]
      enet.pred.matrix[,n]      <- result[[n]][[7]]
    }
    
    for (n in 1:iter) {
      
      y.test.temp <- cbind(y.sample.matrix[,n],y.test.matrix[,n])
      y.test.temp <- y.test.temp[order(y.test.temp[,1]),]
      y.test.matrix[,n] <- y.test.temp[,2]
      
      rf.pred.temp <- cbind(y.sample.matrix[,n],rf.pred.matrix[,n])
      rf.pred.temp <- rf.pred.temp[order(rf.pred.temp[,1]),]
      rf.pred.matrix[,n] <- rf.pred.temp[,2]
      
      lasso.pred.temp <- cbind(y.sample.matrix[,n],lasso.pred.matrix[,n])
      lasso.pred.temp <- lasso.pred.temp[order(lasso.pred.temp[,1]),]
      lasso.pred.matrix[,n] <- lasso.pred.temp[,2]
      
      logistic.pred.temp <- cbind(y.sample.matrix[,n],logistic.pred.matrix[,n])
      logistic.pred.temp <- logistic.pred.temp[order(logistic.pred.temp[,1]),]
      logistic.pred.matrix[,n] <- logistic.pred.temp[,2]
      
      ridge.pred.temp <- cbind(y.sample.matrix[,n],ridge.pred.matrix[,n])
      ridge.pred.temp <- ridge.pred.temp[order(ridge.pred.temp[,1]),]
      ridge.pred.matrix[,n] <- ridge.pred.temp[,2]
      
      enet.pred.temp <- cbind(y.sample.matrix[,n],enet.pred.matrix[,n])
      enet.pred.temp <- enet.pred.temp[order(enet.pred.temp[,1]),]
      enet.pred.matrix[,n] <- enet.pred.temp[,2]
    }
    
    auctable_rf[,jj]       <- apply(rf.pred.matrix,2, function(x){roc(pheno,x)$auc})
    auctable_lasso[,jj]    <- apply(lasso.pred.matrix,2, function(x){roc(pheno,x)$auc})
    auctable_logistic[,jj] <- apply(logistic.pred.matrix,2, function(x){roc(pheno,x)$auc})
    auctable_ridge[,jj]    <- apply(ridge.pred.matrix,2, function(x){roc(pheno,x)$auc})
    auctable_enet[,jj]     <- apply(enet.pred.matrix,2, function(x){roc(pheno,x)$auc})
    
  }
  
  a2_intra[[ww]]       <- list(rf        = auctable_rf,
                               lasso     = auctable_lasso,
                               logistic  = auctable_logistic,
                               ridge     = auctable_ridge,
                               enet      = auctable_enet)
}
save(a2_intra, file = "a2_intra.Rdata")

# --------------------------------- #
#            Cross Study 
# --------------------------------- #

a2_cross <- foreach(ww = 1:ntrans, .inorder = FALSE)%dopar%{ # data transformation parallel computing loop #

     AUC_rf       <- matrix(NA,nrow = snum,ncol = snum,dimnames = list(sname,sname))
     AUC_lasso    <- AUC_logistic <- AUC_ridge <- AUC_enet <- AUC_rf

    for(jj in 1:snum){ # study loop for training #
      
      print(paste("Transformation ww = ", ww, "Study jj = ", jj)) 
      
      # training set
      set1       <- comgenus[[ww]][[jj]]
      set1_pheno <- metaSuball[[jj]][,15]
      
      rf.train        <- train(set1,as.factor(set1_pheno),method="rf",importance=T,ntree = 1000,tuneLength=15,trControl=control)
      cv.lasso.train  <- cv.glmnet(set1,as.factor(set1_pheno), family = "binomial", alpha = 1,standardize = FALSE)
      cv.ridge.train  <- cv.glmnet(set1,as.factor(set1_pheno), family = "binomial", alpha = 0,standardize = FALSE)
      enet.train      <- train(set1,as.factor(set1_pheno),method="glmnet",tuneLength=15,trControl=control,standardize = FALSE)
      logistic.train  <- train(set1,as.factor(set1_pheno),method="glm")
      
      for(mm in 1:snum){ # study loop for testing #
        
        # testing set
        set2       <- comgenus[[ww]][[mm]]
        set2_pheno <- metaSuball[[mm]][,15]

        rf.pred         <- predict(rf.train,set2,type="prob")
        AUC_rf[jj,mm]   <- roc(set2_pheno,rf.pred[,2])$auc
        
        lasso.pred        <- predict(cv.lasso.train,set2,s = "lambda.min",type = "response")
        AUC_lasso[jj,mm]  <- roc(set2_pheno,as.vector(lasso.pred))$auc

        logistic.pred         <- predict(logistic.train,set2,type="prob")
        AUC_logistic[jj,mm]   <- roc(set2_pheno,logistic.pred[,2])$auc

        ridge.pred          <- predict(cv.ridge.train,set2,s = "lambda.min",type = "response")
        AUC_ridge[jj,mm]    <- roc(set2_pheno,as.vector(ridge.pred))$auc

        enet.pred             <- predict(enet.train,set2,type="prob")
        AUC_enet[jj,mm]       <- roc(set2_pheno,enet.pred[,2])$auc
      }
      
    }
    return(list(rf        = AUC_rf,
                lasso     = AUC_lasso,
                logistic  = AUC_logistic,
                ridge     = AUC_ridge,
                enet      = AUC_enet))
}
save(a2_cross,file = "a2_cross.Rdata")

# --------------------------------- #
#            LOSO Study 
# --------------------------------- #

a2_loso <- foreach(ww = 1:ntrans, .inorder = FALSE)%dopar%{
    
    AUC_rf <- rep(NA,snum); names(AUC_rf) <- sname
    AUC_lasso <- AUC_logistic <- AUC_ridge <- AUC_enet <-  AUC_rf 

      for(jj in 1:snum){ #study
        
        
          print(paste("jj = ",jj, "ww =", ww))
        
          # testing set
          set1 <- comgenus[[ww]][[jj]]
          set1_pheno <- metaSuball[[jj]][,15]
          
          # training set
          set2       <- Reduce(rbind, comgenus[[ww]][-jj])
          set2_pheno <- unlist(sapply(metaSuball[-jj], `[`, 15))
          
          rf.train        <- train(set2,as.factor(set2_pheno),method="rf",ntree = 1000,tuneLength=15,trControl=control)
          rf.pred         <- predict(rf.train,set1,type="prob")
          AUC_rf[jj]      <- roc(set1_pheno,rf.pred[,2])$auc
          
          cv.lasso.train  <- cv.glmnet(set2,as.factor(set2_pheno), family = "binomial", alpha = 1,standardize = FALSE)
          lasso.pred      <- predict(cv.lasso.train,set1,s = "lambda.min")
          AUC_lasso[jj]   <- roc(set1_pheno,lasso.pred)$auc
          
          cv.ridge.train  <- cv.glmnet(set2,as.factor(set2_pheno), family = "binomial", alpha = 0,standardize = FALSE)
          ridge.pred      <- predict(cv.ridge.train,set1,s = "lambda.min")
          AUC_ridge[jj]   <- roc(set1_pheno,ridge.pred)$auc
          
          enet.train      <- train(set2,as.factor(set2_pheno),method="glmnet",tuneLength=15,trControl=control,standardize = FALSE)
          enet.pred       <- predict(enet.train,set1,type="prob")
          AUC_enet[jj]    <- roc(set1_pheno,enet.pred[,2])$auc
          
          
          logistic.train    <- train(set2,as.factor(set2_pheno),method="glm")
          logistic.pred     <- predict(logistic.train,set1,type="prob")
          AUC_logistic[jj]  <- roc(set1_pheno,logistic.pred[,2])$auc
          
        }
      
    return(list(rf        = AUC_rf,
                lasso     = AUC_lasso,
                logistic  = AUC_logistic,
                ridge     = AUC_ridge,
                enet      = AUC_enet))
}
save(a2_loso,file = "a2_loso.Rdata")

stopImplicitCluster() # end parallel
