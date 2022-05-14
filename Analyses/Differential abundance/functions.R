# functions for Bayesian Approach




# Contengency table

#         PTB TB
# Present
# Absent 

#' get_ctable: get contengency table based on p0, oddr, n and ratio
#' p0: probability of present given PTB
#' oddr: the present odds ratio of PTB to TB
#' n: sample size of PTB
#' ratio: ratio of TB sample size to PTB sample size
get_ctable <- function(p0, oddr, n, ratio){
  
  u  <- p0/(1-p0) 
  p1 <- u/oddr/(1+u/oddr)
  
  ctable <- matrix(NA,2,2)
  ctable[1,1] <- round(n*p0) # round can be used 
  ctable[2,1] <- n - ctable[1,1] 
  ctable[1,2] <- round(n*ratio*p1) # round can be used 
  ctable[2,2] <- n*ratio - ctable[1,2] 
  ctable
}

#' getoddr: get odds ratio given a contengency table
getoddr <- function(ctable)(ctable[1,1]*ctable[2,2]/ctable[2,1]/ctable[1,2])

#' get_ctable_dr: get the contengency table with detection ratio based on binomial distribution
get_ctable_dr <- function(ctable, dr, seed=1){
  
  set.seed(seed)
  newtable <- ctable
  newtable[1,1] <- rbinom(1, ctable[1,1], dr)
  newtable[1,2] <- rbinom(1, ctable[1,2], dr)
  newtable[2,] <- colSums(ctable) - newtable[1,]
  newtable
}

#' get_p1: get the probability of present given PB based on p0 and oddr
get_p1 <- function(p0,oddr){
  u  <- p0/(1-p0) 
  p1 <- u/oddr/(1+u/oddr)
  p1
}

#' get_loglikelihood_p0: get the likihood vector given a contengency table, value p0 and a r vector
get_loglikelihood_p0 <- function(contable, p0, r){
  
  if (p0==1) {stop("p0 is 1 and the odds is infinity")}
  if (p0==0) {stop("p0 is 0 and the odds is 0")}
  
  u  <- p0/(1-p0) 
  
  logout <- c()
  for(i in 1:length(r)){
    p1 <- u/r[i]/(1+u/r[i])
    logout[i] <- log(p0)*contable[1,1] + log(1-p0)*contable[2,1]+ log(p1)*contable[1,2] + log(1-p1)*contable[2,2]
  }
  logout
} 

#' get_loglikelihood_u: get the likihood matrix given a contengency table, value p0 and a r vector
get_loglikelihood_u <- function(contable, u, r){
  
  logout <- matrix(NA, nrow = length(r), ncol = length(u))
  
  for(i in 1:length(r)){
    for(j in 1:length(u)){
      p0 <- u[j]/(1+u[j])
      p1 <- u[j]/r[i]/(1+u[j]/r[i])
      logout[i,j] <- log(p0)*contable[1,1] + log(1-p0)*contable[2,1]+ log(p1)*contable[1,2] + log(1-p1)*contable[2,2]
    }
  }
  logout
} 

#' get a uniform prior given a r vector
get_prior_uniform <- function(x){data.frame(value = x, prior = 1/length(x))}

#' get the posterior distribution for r, fold-change of odds
#' the function only need the prior of r and loglikelihood matrix. 
#' we can ignor the prior for u beacuase we assume it's uniform distribution for any of the study
#' and mathmatically the prior of u can be eliminated.


get_posteior_r_u   <- function(prior_r, loglikelihood){ 
  
  # For a given r value, we first sum all possible u value
  # thus, we first need to calcuate log(u_1 + u_2 + ... u_n) given a r value.
  # to avoid underflow, we can first calcuate log(u_1 + u_2), then log(u_1 + u_2  + u_3), so on and so forth.
  # log(u_1 + u_2)  = log(u_1) + log(1 + exp(log(u_2) - log(u_1)))
  prior_r[prior_r[,2]==0,2] <- 3e-323
  
  if(dim(loglikelihood)[2]==1){
    log_sum_r <- loglikelihood[,1]
  }else{
    log_sum_r <-  apply(loglikelihood, 1, function(x){
      
      x <- x[order(x)]
      tempsum <- x[1]
      for(i in 1:(length(x)-1)){
        tempsum <- log( 1 + exp(tempsum - x[i+1]) ) + x[i+1]
      }
      tempsum
    })
    }
  
  # since it is in the log scale, we add the prior of r.
  log_sum_prior <- log_sum_r + log(prior_r[,2])
  
  log_sum_prior_1 <- log_sum_prior[order(log_sum_prior)]
  # use the same technic as above
  log_sum <- log_sum_prior_1[1]
  for(i in 1:(length(log_sum_prior_1)-1)){
    log_sum <- log( 1 + exp( log_sum - log_sum_prior_1[i+1]) ) + log_sum_prior_1[i+1]
  }
  
  data.frame(value = prior_r$value, posteiror = exp(log_sum_prior - log_sum))
  
  # The following code can be used to confirm the above caculate is correct when there is no underflow issue.
  # likihood <- exp(loglikelihood)
  # allr     <- rowSums(likihood)
  # data.frame(value = prior_r$value, posteiror = allr*prior_r[,2]/sum(allr*prior_r[,2])
  
}

#' for likelihood is a vector
get_posteior_r_p0  <- function(prior_r, loglikelihood){ 
  
  # For a given r value, we first sum all possible u value
  # thus, we first need to calcuate log(u_1 + u_2 + ... u_n) given a r value.
  # to avoid underflow, we can first calcuate log(u_1 + u_2), then log(u_1 + u_2  + u_3), so on and so forth.
  # log(u_1 + u_2)  = log(u_1) + log(1 + exp(log(u_2) - log(u_1)))
  
  
  prior_r[prior_r[,2]==0,2] <- 3e-323
  # since it is in the log scale, we add the prior of r.
  log_sum_prior <- loglikelihood + log(prior_r[,2])
  
  log_sum_prior_1 <- log_sum_prior[order(log_sum_prior)]
  # use the same technic as above
  log_sum <- log_sum_prior_1[1]
  for(i in 1:(length(log_sum_prior_1)-1)){
    log_sum <- log( 1 + exp( log_sum - log_sum_prior_1[i+1]) ) + log_sum_prior_1[i+1]
  }
  
  data.frame(value = prior_r$value, posteiror = exp(log_sum_prior - log_sum))
  
  # The following code can be used to confirm the above caculate is correct when there is no underflow issue.
  # likihood <- exp(loglikelihood)
  # allr     <- rowSums(likihood)
  # data.frame(value = prior_r$value, posteiror = allr*prior_r[,2]/sum(allr*prior_r[,2])
  
}

#' get_observe_p_log: get the likelihood use observed present count and true present counts and detection rate
get_observe_p_log <- function(obs_table, true_table, dr){
  log(dr)*(sum(obs_table[1,])) + log(1-dr)*(sum(true_table[1,]) - sum(obs_table[1,])) + 
    log(choose(obs_table[2,1], true_table[1,1] - obs_table[1,1])) + 
    log(choose(obs_table[2,2], true_table[1,2] - obs_table[1,2])) 
}


get_loglikelihood_p0_drate <- function(contable, p0, r, drate){
  
  u <- p0/(1-p0)
  logout <- matrix(NA, nrow = length(r), ncol = length(drate))
  
  for(i in 1:length(r)){
    for(j in 1:length(drate)){
      p0 <- u/(1+u)*drate[j]
      p1 <- u/r[i]/(1+u/r[i])*drate[j]
      logout[i,j] <- log(p0)*contable[1,1] + log(1-p0)*contable[2,1]+ log(p1)*contable[1,2] + log(1-p1)*contable[2,2]
    }
  }
  logout
} 
