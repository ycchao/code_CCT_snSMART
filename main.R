##==================================================================================================##
## This is the code for the simulation study of the manuscript "Dynamic Treatment Regimens in small ##
## n, Sequential, Multiple Assignment, Randomized Trials: an application in focal segmental         ##
## glomerulosclerosis".                                                                             ##
##                                                                                                  ##
## The code was written by Yan-Cheng Chao from the Department of Biostatistics at the University    ##
## of Michigan. The code was last run on Sept. 17, 2019.                                            ##
##==================================================================================================##

library(geepack)
library(rjags)
library(nCDunnett)

##### Estimation #####

n.sim <- 1000
NUM_ARMS <- 3
MCMC_SAMPLE <- 6000
BURN.IN <- 1000
n_MCMC_chain <- 1

# prior for BJSM
pi_prior.a <- c(0.4,0.4,0.4)
pi_prior.b <- c(1.6,1.6,1.6) #mean of pi of treatment = 0.4/(0.4+1.6) = 0.2
beta0_prior.a <- 1.6
beta0_prior.b <- 0.4         #mean of beta0 = 1.6/(1.6+0.4) = 0.8
beta1_prior.r <- 2         #mean of beta1 = r/mu = 2/2 = 1, beta1 ~ gamma(r,mu)
beta1_prior.mu <- 2        #beta1 ~ gamma

#scenario 1
pi_1A <- 0.4
pi_1B <- 0.4        
pi_1C <- 0.2
#scenario 2
pi_1A <- 0.45           
pi_1B <- 0.45        
pi_1C <- 0.2
#scenario 3
pi_1A <- 0.45           
pi_1B <- 0.2         
pi_1C <- 0.2
#scenario 4
pi_1A <- 0.45           
pi_1B <- 0.3        
pi_1C <- 0.2

#scenario a
discount_y <- c(1,1,1)
discount_n1 <- c(0.8,0.6,0.4)
discount_n2 <- c(0.8,0.6,0.4)
#scenario b
discount_y <- c(1.5,1,0.5)
discount_n1 <- c(0.8,0.6,0.4)
discount_n2 <- c(0.8,0.6,0.4)
#scenario c
discount_y <- c(1.5,1,0.5)
discount_n1 <- c(0.65,0.7,0.75)
discount_n2 <- c(0.75,0.6,0.45)

data_generation <- function(pi_1A,pi_1B,pi_1C,discount_y,discount_n1,discount_n2,n_A,n_B,n_C){
  pi_2A_y <<- pi_1A * discount_y[1]        # Second stage response rate of responders to A 
  pi_2B.A_n <<- pi_1A * discount_n1[2]     # Second stage response rate of non-responders to B who receive A in the second stage
  pi_2C.A_n <<- pi_1A * discount_n1[3]     # Second stage response rate of non-responders to C who receive A in the second stage
  pi_2B_y <<- pi_1B * discount_y[2]        # Second stage response rate of responders to B
  pi_2A.B_n <<- pi_1B * discount_n1[1]     # Second stage response rate of non-responders to A who receive B in the second stage
  pi_2C.B_n <<- pi_1B * discount_n2[3]     # Second stage response rate of non-responders to C who receive B in the second stage
  pi_2C_y <<- pi_1C * discount_y[3]        # Second stage response rate of responders to C
  pi_2A.C_n <<- pi_1C * discount_n2[1]     # Second stage response rate of non-responders to A who receive C in the second stage
  pi_2B.C_n <<- pi_1C * discount_n2[2]     # Second stage response rate of non-responders to B who receive C in the second stage
  ### stage 1 ###
  # trt A
  n_resp_A <- rbinom(n=1,size=n_A,prob=pi_1A)
  n_nonresp_A <- n_A - n_resp_A
  n_nonresp_A_stage2_B <<- ifelse(n_nonresp_A%%2==1,(n_nonresp_A+1)/2,n_nonresp_A/2)
  n_nonresp_A_stage2_C <<- n_nonresp_A - n_nonresp_A_stage2_B
  
  # trt B
  n_resp_B <- rbinom(n=1,size=n_B,prob=pi_1B)
  n_nonresp_B <- n_B - n_resp_B
  n_nonresp_B_stage2_A <<- ifelse(n_nonresp_B%%2==1,(n_nonresp_B+1)/2,n_nonresp_B/2)
  n_nonresp_B_stage2_C <<- n_nonresp_B - n_nonresp_B_stage2_A
  
  # trt C
  n_resp_C <- rbinom(n=1,size=n_C,prob=pi_1C)
  n_nonresp_C <- n_C - n_resp_C
  n_nonresp_C_stage2_A <<- ifelse(n_nonresp_C%%2==1,(n_nonresp_C+1)/2,n_nonresp_C/2)
  n_nonresp_C_stage2_B <<- n_nonresp_C - n_nonresp_C_stage2_A
  
  ### stage 2 ###
  # stage I = trt A
  n_resp_A.A.Y <- rbinom(n=1,size=n_resp_A,prob=pi_2A_y)
  data_A.A.Y <- data.frame(treatment_stageI = rep(1,n_resp_A),
                           response_stageI = rep(1,n_resp_A),
                           treatment_stageII = rep(1,n_resp_A),
                           response_stageII = c(rep(1,n_resp_A.A.Y),rep(0,n_resp_A-n_resp_A.A.Y)))
  
  n_resp_A.B.Y <- rbinom(n=1,size=n_nonresp_A_stage2_B,prob=pi_2A.B_n)
  data_A.B.Y <- data.frame(treatment_stageI = rep(1,n_nonresp_A_stage2_B),
                           response_stageI = rep(0,n_nonresp_A_stage2_B),
                           treatment_stageII = rep(2,n_nonresp_A_stage2_B),
                           response_stageII = c(rep(1,n_resp_A.B.Y),rep(0,n_nonresp_A_stage2_B-n_resp_A.B.Y)))
  
  n_resp_A.C.Y <- rbinom(n=1,size=n_nonresp_A_stage2_C,prob=pi_2A.C_n)
  data_A.C.Y <- data.frame(treatment_stageI = rep(1,n_nonresp_A_stage2_C),
                           response_stageI = rep(0,n_nonresp_A_stage2_C),
                           treatment_stageII = rep(3,n_nonresp_A_stage2_C),
                           response_stageII = c(rep(1,n_resp_A.C.Y),rep(0,n_nonresp_A_stage2_C-n_resp_A.C.Y)))
  
  # stage I = trt B
  n_resp_B.B.Y <- rbinom(n=1,size=n_resp_B,prob=pi_2B_y)
  data_B.B.Y <- data.frame(treatment_stageI = rep(2,n_resp_B),
                           response_stageI = rep(1,n_resp_B),
                           treatment_stageII = rep(2,n_resp_B),
                           response_stageII = c(rep(1,n_resp_B.B.Y),rep(0,n_resp_B-n_resp_B.B.Y)))
  
  n_resp_B.A.Y <- rbinom(n=1,size=n_nonresp_B_stage2_A,prob=pi_2B.A_n)
  data_B.A.Y <- data.frame(treatment_stageI = rep(2,n_nonresp_B_stage2_A),
                           response_stageI = rep(0,n_nonresp_B_stage2_A),
                           treatment_stageII = rep(1,n_nonresp_B_stage2_A),
                           response_stageII = c(rep(1,n_resp_B.A.Y),rep(0,n_nonresp_B_stage2_A-n_resp_B.A.Y)))
  
  n_resp_B.C.Y <- rbinom(n=1,size=n_nonresp_B_stage2_C,prob=pi_2B.C_n)
  data_B.C.Y <- data.frame(treatment_stageI = rep(2,n_nonresp_B_stage2_C),
                           response_stageI = rep(0,n_nonresp_B_stage2_C),
                           treatment_stageII = rep(3,n_nonresp_B_stage2_C),
                           response_stageII = c(rep(1,n_resp_B.C.Y),rep(0,n_nonresp_B_stage2_C-n_resp_B.C.Y)))
  
  # stage I = trt C
  n_resp_C.C.Y <- rbinom(n=1,size=n_resp_C,prob=pi_2C_y)
  data_C.C.Y <- data.frame(treatment_stageI = rep(3,n_resp_C),
                           response_stageI = rep(1,n_resp_C),
                           treatment_stageII = rep(3,n_resp_C),
                           response_stageII = c(rep(1,n_resp_C.C.Y),rep(0,n_resp_C-n_resp_C.C.Y)))
  
  n_resp_C.A.Y <- rbinom(n=1,size=n_nonresp_C_stage2_A,prob=pi_2C.A_n)
  data_C.A.Y <- data.frame(treatment_stageI = rep(3,n_nonresp_C_stage2_A),
                           response_stageI = rep(0,n_nonresp_C_stage2_A),
                           treatment_stageII = rep(1,n_nonresp_C_stage2_A),
                           response_stageII = c(rep(1,n_resp_C.A.Y),rep(0,n_nonresp_C_stage2_A-n_resp_C.A.Y)))
  
  n_resp_C.B.Y <- rbinom(n=1,size=n_nonresp_C_stage2_B,prob=pi_2C.B_n)
  data_C.B.Y <- data.frame(treatment_stageI = rep(3,n_nonresp_C_stage2_B),
                           response_stageI = rep(0,n_nonresp_C_stage2_B),
                           treatment_stageII = rep(2,n_nonresp_C_stage2_B),
                           response_stageII = c(rep(1,n_resp_C.B.Y),rep(0,n_nonresp_C_stage2_B-n_resp_C.B.Y)))
  
  data_stageI.II <- rbind(data_A.A.Y,data_A.B.Y,data_A.C.Y,
                          data_B.B.Y,data_B.A.Y,data_B.C.Y,
                          data_C.C.Y,data_C.A.Y,data_C.B.Y)
  return(data_stageI.II)
} 
Dunnett_test_SMART_6beta <- function(pi_1A,pi_1B,pi_1C,discount_y,discount_n1,discount_n2,n_A,n_B,n_C,n.sim){
  p_value <- beta <- c()
  error <- 0
  for (i in 1:n.sim){
    set.seed(i)
    mydata <- data_generation(pi_1A,pi_1B,pi_1C,discount_y,discount_n1,discount_n2,n_A,n_B,n_C)
    Y <- c(mydata$response_stageI, mydata$response_stageII)
    XA1 <- as.numeric(mydata$treatment_stageI==1)
    XB1 <- as.numeric(mydata$treatment_stageI==2)
    XC1 <- as.numeric(mydata$treatment_stageI==3)
    XA2 <- as.numeric(mydata$treatment_stageII==1)
    XB2 <- as.numeric(mydata$treatment_stageII==2)
    XC2 <- as.numeric(mydata$treatment_stageII==3)
    XA <- c(XA1, XA2)
    XB <- c(XB1, XB2)
    XC <- c(XC1, XC2)
    Z1A <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==1,mydata$response_stageI,0))
    Z2A <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==1,1-mydata$response_stageI,0))
    Z1B <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==2,mydata$response_stageI,0))
    Z2B <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==2,1-mydata$response_stageI,0))
    Z1C <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==3,mydata$response_stageI,0))
    Z2C <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==3,1-mydata$response_stageI,0))
    ptid <- rep(1:nrow(mydata), 2)
    geedata <- data.frame(ptid, XA, XB, XC, Z1A, Z2A, Z1B, Z2B, Z1C, Z2C, Y)
    geedata <- geedata[order(geedata$ptid),]
    rm(ptid, Y, Z1A, Z2A, Z1B, Z2B, Z1C, Z2C, XA, XB, XC)
    error_temp <- try({
      mod1 <- geeglm(Y~XA+XB+XC+Z1A+Z2A+Z1B+Z2B+Z1C+Z2C-1, family=poisson(link="log"), data=geedata, id=ptid,corstr = "independence")
      beta_hat <- mod1$coefficients;
      cov_hat <- mod1$geese$vbeta;
      L <- matrix(c(1,0,-1,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0),nrow=2,ncol=9,byrow=T);
      test_stat <- abs(L %*% beta_hat)/sqrt(diag(L %*% cov_hat %*% t(L)));
      R <- cov2cor(L %*% cov_hat %*% t(L));
      lambda <- fa(R,rotate="none",fm="pa")$loadings[1:2];
      p_value <- rbind(p_value,1-pNCDun(test_stat,nu=Inf,rho=lambda^2,delta=c(0,0)))
      beta <- rbind(beta,beta_hat)
    }
    )
    if(inherits(error_temp,"try-error"))
      error <- error + 1
  }
  return(list(p_value,error,beta))
}

n_A <- n_B <- n_C <- 40

pi_hat <- sd_pi_hat <- pi_hat_fsmle <- sd_pi_hat_fsmle <- pi_hat_bjsm <- c()
error_count <- 0
pi_DTR_hat <- pi_DTR_se <- NULL
pi_DTR_est <- c()
DTR_pi_mean_WRRM <- DTR_pi_se_WRRM <- c()
for (i in 1:n.sim){
  # JSRM
  set.seed(i+10000)
  mydata <- data_generation(pi_1A,pi_1B,pi_1C,discount_y,discount_n1,discount_n2,n_A,n_B,n_C)
  mydata$disc <- 2 * mydata$treatment_stageI - (mydata$response_stageI == 0)
  Y <- c(mydata$response_stageI, mydata$response_stageII)
  XA1 <- as.numeric(mydata$treatment_stageI==1)
  XB1 <- as.numeric(mydata$treatment_stageI==2)
  XC1 <- as.numeric(mydata$treatment_stageI==3)
  XA2 <- as.numeric(mydata$treatment_stageII==1)
  XB2 <- as.numeric(mydata$treatment_stageII==2)
  XC2 <- as.numeric(mydata$treatment_stageII==3)
  XA <- c(XA1, XA2)
  XB <- c(XB1, XB2)
  XC <- c(XC1, XC2)
  Z1A <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==1,mydata$response_stageI,0))
  Z2A <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==1,1-mydata$response_stageI,0))
  Z1B <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==2,mydata$response_stageI,0))
  Z2B <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==2,1-mydata$response_stageI,0))
  Z1C <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==3,mydata$response_stageI,0))
  Z2C <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==3,1-mydata$response_stageI,0))
  ptid <- rep(1:nrow(mydata), 2)
  geedata <- data.frame(ptid, XA, XB, XC, Z1A, Z2A, Z1B, Z2B, Z1C, Z2C, Y)
  geedata <- geedata[order(geedata$ptid),]
  rm(ptid, Y, Z1A, Z2A, Z1B, Z2B, Z1C, Z2C, XA, XB, XC)
  try({
    mod1 <- geeglm(Y~XA+XB+XC+Z1A+Z2A+Z1B+Z2B+Z1C+Z2C-1, family=poisson(link="log"), data=geedata, id=ptid,corstr = "independence")
    beta_hat <- mod1$coefficients[1:3];
    sd_beta_hat <- summary(mod1)$coef[1:3,2];
    pi_hat <- rbind(pi_hat,exp(beta_hat));
    sd_pi_hat <- rbind(sd_pi_hat,exp(beta_hat)*sd_beta_hat);
    b_hat <- coef(mod1);
    grad <- exp(rbind(c(2*b_hat[1]+b_hat[4], b_hat[2]+b_hat[5], b_hat[1]+b_hat[2]+b_hat[5]),
                      c(2*b_hat[1]+b_hat[4], b_hat[3]+b_hat[5], b_hat[1]+b_hat[3]+b_hat[5]),
                      c(2*b_hat[2]+b_hat[6], b_hat[1]+b_hat[7], b_hat[2]+b_hat[1]+b_hat[7]),
                      c(2*b_hat[2]+b_hat[6], b_hat[3]+b_hat[7], b_hat[2]+b_hat[3]+b_hat[7]),
                      c(2*b_hat[3]+b_hat[8], b_hat[1]+b_hat[9], b_hat[3]+b_hat[1]+b_hat[9]),
                      c(2*b_hat[3]+b_hat[8], b_hat[2]+b_hat[9], b_hat[3]+b_hat[2]+b_hat[9])));
    pi_DTR_hat <- rbind(pi_DTR_hat, c(1, 1, -1) %*% t(grad))  
    grad[,3] <- -grad[,3]
    sigma_b <- mod1$geese$vbeta
    L1 <- matrix(c(2, 0, 0, 1, 0, 0, 0, 0, 0,
                   0, 1, 0, 0, 1, 0, 0, 0, 0,
                   1, 1, 0, 0, 1, 0, 0, 0, 0), nrow=3, ncol=9, byrow=T)
    sigma_g <- L1 %*% sigma_b %*% t(L1)
    seAB <- sqrt(grad[1,] %*% sigma_g %*% grad[1,])
    
    L2 <- matrix(c(2, 0, 0, 1, 0, 0, 0, 0, 0,
                   0, 0, 1, 0, 1, 0, 0, 0, 0,
                   1, 0, 1, 0, 1, 0, 0, 0, 0), nrow=3, ncol=9, byrow=T)
    sigma_g <- L2 %*% sigma_b %*% t(L2)
    seAC <- sqrt(grad[2,] %*% sigma_g %*% grad[2,])
    
    L3 <- matrix(c(0, 2, 0, 0, 0, 1, 0, 0, 0,
                   1, 0, 0, 0, 0, 0, 1, 0, 0,
                   1, 1, 0, 0, 0, 0, 1, 0, 0), nrow=3, ncol=9, byrow=T)
    sigma_g <- L3 %*% sigma_b %*% t(L3)
    seBA <- sqrt(grad[3,] %*% sigma_g %*% grad[3,])
    
    L4 <- matrix(c(0, 2, 0, 0, 0, 1, 0, 0, 0,
                   0, 0, 1, 0, 0, 0, 1, 0, 0,
                   0, 1, 1, 0, 0, 0, 1, 0, 0), nrow=3, ncol=9, byrow=T)
    sigma_g <- L4 %*% sigma_b %*% t(L4)
    seBC <- sqrt(grad[4,] %*% sigma_g %*% grad[4,])
    
    L5 <- matrix(c(0, 0, 2, 0, 0, 0, 0, 1, 0,
                   1, 0, 0, 0, 0, 0, 0, 0, 1,
                   1, 0, 1, 0, 0, 0, 0, 0, 1), nrow=3, ncol=9, byrow=T)
    sigma_g <- L5 %*% sigma_b %*% t(L5)
    seCA <- sqrt(grad[5,] %*% sigma_g %*% grad[5,])
    
    L6 <- matrix(c(0, 0, 2, 0, 0, 0, 0, 1, 0,
                   0, 1, 0, 0, 0, 0, 0, 0, 1,
                   0, 1, 1, 0, 0, 0, 0, 0, 1), nrow=3, ncol=9, byrow=T)
    sigma_g <- L6 %*% sigma_b %*% t(L6)
    seCB <- sqrt(grad[6,] %*% sigma_g %*% grad[6,])
    pi_DTR_se <- rbind(pi_DTR_se, c(seAB, seAC, seBA, seBC, seCA, seCB))
  })
  #BJSM
  
  jags_path <- "Your Path to BUG file"
  jag.model.name <- "BJSM_6betas_gamma_beta1.bug"   # beta1 ~ gamma
  tryCatch({
    jag <- jags.model(file.path(jags_path,jag.model.name),
                      data=list(n = nrow(mydata),
                                num_arms = NUM_ARMS,
                                Y1 = mydata$response_stageI,
                                Y2 = mydata$response_stageII,
                                treatment_stageI = mydata$treatment_stageI,
                                treatment_stageII = mydata$treatment_stageII,
                                response_stageI_disc = mydata$disc,
                                #prior
                                pi_prior.a = pi_prior.a,
                                pi_prior.b = pi_prior.b,
                                beta0_prior.a = beta0_prior.a,
                                beta0_prior.b = beta0_prior.b,
                                # beta1_prior.a = beta1_prior.a,    # pareto
                                # beta1_prior.c = beta1_prior.c     # pareto
                                beta1_prior.a = beta1_prior.r,   # gamma
                                beta1_prior.c = beta1_prior.mu  # gamma
                                # beta1_prior.a = beta1_prior.mu,  # lognormal
                                # beta1_prior.c = beta1_prior.tau # lognormal
                      ),
                      n.chains=n_MCMC_chain,n.adapt = BURN.IN)   
    posterior_sample <- coda.samples(jag,
                                     c('pi','beta'),
                                     MCMC_SAMPLE)
  },
  warning = function(war){
    warning_count <- warning_count + 1
    err_war_message <- rbind(paste("The warning ", warning_count, " is: ", war))
  },
  error = function(err){
    error_count <- error_count + 1
    err_war_message <- rbind(paste("The error ", error_count, " is: ", err))
    error_ind <- 1
  },
  finally = {
    print(i)     # show the number of iterations run 
  }
  )
  out_post <- posterior_sample[[1]]
  pi_hat_bjsm <- rbind(pi_hat_bjsm,apply(out_post[,7:9],2,mean))
  pi_AB_tt <- out_post[,7]^2*out_post[,2]+(1-out_post[,7])*out_post[,8]*out_post[,1]    
  pi_AC_tt <- out_post[,7]^2*out_post[,2]+(1-out_post[,7])*out_post[,9]*out_post[,1]
  pi_BA_tt <- out_post[,8]^2*out_post[,4]+(1-out_post[,8])*out_post[,7]*out_post[,3]
  pi_BC_tt <- out_post[,8]^2*out_post[,4]+(1-out_post[,8])*out_post[,9]*out_post[,3]
  pi_CA_tt <- out_post[,9]^2*out_post[,6]+(1-out_post[,9])*out_post[,7]*out_post[,5]
  pi_CB_tt <- out_post[,9]^2*out_post[,6]+(1-out_post[,9])*out_post[,8]*out_post[,5]
  pi_DTR_est <- rbind(pi_DTR_est,c(mean(pi_AB_tt),mean(pi_AC_tt),mean(pi_BA_tt),mean(pi_BC_tt),mean(pi_CA_tt),mean(pi_CB_tt)))
  #WRRM
  mydata$id <- seq(1,nrow(mydata),by=1)
  data_stageI.II_rep <- mydata[rep(row.names(mydata),mydata$response_stageI+1),]
  if(length(data_stageI.II_rep[data_stageI.II_rep$response_stageI==1 & data_stageI.II_rep$treatment_stageI==1,3])!=0){
    data_stageI.II_rep[data_stageI.II_rep$response_stageI==1 & data_stageI.II_rep$treatment_stageI==1,3] <- rep(c(2,3))
  }
  if(length(data_stageI.II_rep[data_stageI.II_rep$response_stageI==1 & data_stageI.II_rep$treatment_stageI==2,3])!=0){
    data_stageI.II_rep[data_stageI.II_rep$response_stageI==1 & data_stageI.II_rep$treatment_stageI==2,3] <- rep(c(1,3))
  }
  if(length(data_stageI.II_rep[data_stageI.II_rep$response_stageI==1 & data_stageI.II_rep$treatment_stageI==3,3])!=0){
    data_stageI.II_rep[data_stageI.II_rep$response_stageI==1 & data_stageI.II_rep$treatment_stageI==3,3] <- rep(c(1,2))
  }
  
  data_stageI.II_rep$weight <- ifelse(data_stageI.II_rep$response_stageI==1,3,6)  ## second weights need to be changed if the second-stage randomization is unbalanced
  data_stageI.II_rep$trt_stageIA <- ifelse(data_stageI.II_rep$treatment_stageI==1,1,0)
  data_stageI.II_rep$trt_stageIB <- ifelse(data_stageI.II_rep$treatment_stageI==2,1,0)
  data_stageI.II_rep$trt_stageIC <- ifelse(data_stageI.II_rep$treatment_stageI==3,1,0)
  data_stageI.II_rep$trt_stageIIB <- ifelse(data_stageI.II_rep$treatment_stageII==2,1,0)
  data_stageI.II_rep$trt_stageIIC <- ifelse(data_stageI.II_rep$treatment_stageII==3,1,0)
  fitted <- geeglm(response_stageII~trt_stageIB+trt_stageIC+trt_stageIA:trt_stageIIC+
                     trt_stageIIC:trt_stageIB+trt_stageIIB:trt_stageIC,id=id,weights=weight,
                   family=poisson("log"),corstr = "independence",data=data_stageI.II_rep)
  contrast_matrix <- matrix(c(1,0,0,0,0,0,
                              1,0,0,1,0,0,
                              1,1,0,0,0,0,
                              1,1,0,0,1,0,
                              1,0,1,0,0,0,
                              1,0,1,0,0,1),nrow=6,ncol=6,byrow=T)
  beta_est <- fitted$coefficients
  vcov_beta <- fitted$geese$vbeta
  DTR_est <- rep(NA,6)
  se_DTR_est <- rep(NA,6)
  lower_CI <- rep(NA,6)
  upper_CI <- rep(NA,6)
  for (k in 1:6){
    cont <- contrast_matrix[k,]
    DTR_est[k] <- exp(cont%*%beta_est)
    se_DTR_est[k] <- DTR_est[k]*sqrt(cont%*%vcov_beta%*%cont)   #delta method is used
  }
  DTR_pi_mean_WRRM <- rbind(DTR_pi_mean_WRRM,DTR_est)
  DTR_pi_se_WRRM <- rbind(DTR_pi_se_WRRM,se_DTR_est)
  
  pi_hat_fsmle <- rbind(pi_hat_fsmle,aggregate(response_stageI~treatment_stageI,data=mydata,function(x) mean(x==1))[,2])
  sd_pi_hat_fsmle <- rbind(sd_pi_hat_fsmle,sqrt(pi_hat_fsmle * (1 - pi_hat_fsmle) / c(n_A,n_B,n_C)))
}

trt_effect_output <- data.frame(true_pi = c(pi_1A,pi_1B,pi_1C),
                             pi_BJSM = apply(pi_hat_bjsm,2,mean),
                             sd_BJSM = apply(pi_hat_bjsm,2,sd),
                             pi_JSRM = apply(pi_hat,2,mean),
                             sd_JSRM = apply(pi_hat,2,sd),
                             emp_sd_JSRM = apply(sd_pi_hat,2,mean),
                             pi_fsmle = apply(pi_hat_fsmle,2,mean),
                             sd_fsmle = apply(pi_hat_fsmle,2,sd),
                             emp_sd_fsmle = apply(sd_pi_hat_fsmle,2,mean))
trt_effect_output$bias_BJSM <- trt_effect_output$pi_BJSM - trt_effect_output$true_pi
trt_effect_output$bias_JSRM <- trt_effect_output$pi_JSRM - trt_effect_output$true_pi
trt_effect_output$bias_fsmle <- trt_effect_output$pi_fsmle - trt_effect_output$true_pi
trt_effect_output$rMSE_BJSM <- sqrt(trt_effect_output$sd_BJSM^2 + trt_effect_output$bias_BJSM^2)
trt_effect_output$rMSE_JSRM <- sqrt(trt_effect_output$emp_sd_JSRM^2 + trt_effect_output$bias_JSRM^2)
trt_effect_output$rMSE_fsmle <- sqrt(trt_effect_output$emp_sd_fsmle^2 + trt_effect_output$bias_fsmle^2)

expected <- c()  # expected DTR response rates
expected[1] <- pi_1A * pi_2A_y + (1 - pi_1A) * pi_2A.B_n
expected[2] <- pi_1A * pi_2A_y + (1 - pi_1A) * pi_2A.C_n
expected[3] <- pi_1B * pi_2B_y + (1 - pi_1B) * pi_2B.A_n
expected[4] <- pi_1B * pi_2B_y + (1 - pi_1B) * pi_2B.C_n
expected[5] <- pi_1C * pi_2C_y + (1 - pi_1C) * pi_2C.A_n
expected[6] <- pi_1C * pi_2C_y + (1 - pi_1C) * pi_2C.B_n
dtr_effect_output <- data.frame(true_pi_dtr = expected,
                                 pi_DTR_BJSM = apply(pi_DTR_est,2,function(x) mean(x,na.rm=T)),
                                 sd_DTR_BJSM = apply(pi_DTR_est,2,sd),
                                 pi_DTR_JSRM = apply(pi_DTR_hat,2,function(x) mean(x,na.rm=T)),
                                 sd_DTR_JSRM = apply(pi_DTR_hat,2,sd),
                                 emp_sd_DTR_JSRM = apply(pi_DTR_se,2,function(x) mean(x,na.rm=T)),
                                 pi_DTR_WRRM = apply(DTR_pi_mean_WRRM,2,function(x) mean(x,na.rm=T)),
                                 sd_DTR_WRRM = apply(DTR_pi_mean_WRRM,2,sd),
                                 emp_sd_DTR_WRRM = apply(DTR_pi_se_WRRM,2,function(x) mean(x,na.rm=T)))
dtr_effect_output$bias_BJSM <- dtr_effect_output$pi_DTR_BJSM - dtr_effect_output$true_pi_dtr
dtr_effect_output$bias_JSRM <- dtr_effect_output$pi_DTR_JSRM - dtr_effect_output$true_pi_dtr
dtr_effect_output$bias_WRRM <- dtr_effect_output$pi_DTR_WRRM - dtr_effect_output$true_pi_dtr
dtr_effect_output$rMSE_BJSM <- sqrt(dtr_effect_output$sd_DTR_BJSM^2 + dtr_effect_output$bias_BJSM^2)
dtr_effect_output$rMSE_JSRM <- sqrt(dtr_effect_output$emp_sd_DTR_JSRM^2 + dtr_effect_output$bias_JSRM^2)
dtr_effect_output$rMSE_WRRM <- sqrt(dtr_effect_output$emp_sd_DTR_WRRM^2 + dtr_effect_output$bias_WRRM^2)

##### Sample Size Calculation #####
n.sim <- 1000
## Scenario 1a: (0.4,0.4,0.2)
pi_1A_H0 <- 0.2
pi_1B_H0 <- 0.2
pi_1C_H0 <- 0.2
discount_y_H0 <- c(1.5,1,0.5)    # The pre-specified values of linkage parameters for responders to treatments A, B, C
discount_n1_H0 <- c(0.8,0.6,0.4)   # The pre-specified values of linkage parameters for non-responders who receive first second-stage treatment. 
discount_n2_H0 <- c(0.8,0.6,0.4)   # The pre-specified values of linkage parameters for non-responders who receive second second-stage treatment. 

pi_1A <- 0.45   # First-stage response rate for trt A
pi_1B <- 0.3     # First-stage response rate for trt B
pi_1C <- 0.2    # First-stage response rate for trt C
discount_y <- c(1.5,1,0.5)    # For the first-stage responders, their second-stage trt response rates for A, B and C, respectively, are (0.35*1,0.35*1,0.1*1)=(0.35,0.35,0.1)
discount_n1 <- c(0.8,0.6,0.4)   # For the first-stage non-responders, the second-stage trt response rates for A-->B,B-->A and C-->A are (0.35*5/7,0.35*5/7,0.35*5/7)=(0.25,0.25,0.25)
discount_n2 <- c(0.8,0.6,0.4)   # For the first-stage non-responders, the second-stage trt response rates for A-->C,B-->C and C-->B are (0.10*1/2,0.10*1/2,0.35*5/7)=(0.05,0.05,0.25)

# The sample sizes for each arm in our design
n_A <- n_B <- n_C <- c(25,30,40,50,75,100)

alpha <- 0.1
type_1_error_rate <- c()
power <- c()
error_count_H0 <- c()
error_count_experiment <- c()
beta_H0 <- beta_experiment <- c()
# power2 <- c()
for (j in 1:length(n_A)){
  temp_H0 <- Dunnett_test_SMART_6beta(pi_1A_H0,pi_1B_H0,pi_1C_H0,discount_y_H0,
                                      discount_n1_H0,discount_n2_H0,n_A[j],n_B[j],n_C[j],n.sim)
  p_value_H0 <- data.frame(temp_H0[[1]])
  error_count_H0 <- c(error_count_H0,temp_H0[[2]])
  beta_H0 <- rbind(beta_H0,apply(temp_H0[[3]],2,mean))
  p_value_H0$sig <- 1 - (p_value_H0$X1 > alpha) * (p_value_H0$X2 > alpha)
  type_1_error_rate <- c(type_1_error_rate,mean(p_value_H0$sig))
  temp_experiment <- Dunnett_test_SMART_6beta(pi_1A,pi_1B,pi_1C,discount_y,
                                              discount_n1,discount_n2,n_A[j],n_B[j],n_C[j],n.sim)
  p_value_experiment <- data.frame(temp_experiment[[1]])
  error_count_experiment <- c(error_count_experiment,temp_experiment[[2]])
  beta_experiment <- rbind(beta_experiment,apply(temp_experiment[[3]],2,mean))
  p_value_experiment$sig1 <- 1 - (p_value_experiment$X1 > alpha) * (p_value_experiment$X2 > alpha)
  power <- c(power,mean(p_value_experiment$sig1))
  print(j)
  # p_value_experiment$sig2 <- (p_value_experiment$X1 <= 0.05) * (p_value_experiment$X2 <= 0.05)
  # power2 <- c(power2,mean(p_value_experiment$sig2))
}
error_count_H0
error_count_experiment
result.table <- data.frame(n=3*n_A,typeIerror=type_1_error_rate,power=power)
