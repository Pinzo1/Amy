#### Testing armidillo
####
rm(list=ls())
require("Rcpp")
require("RcppArmadillo")
require("deSolve")

setwd("H:/My Documents/Model_Pneumo_simple")

Age_groups=c("1m","2m","3m","4m","5m","6m","7m","8m","9m","10m","11m","12m",
             "1-2y","2-3y","3-4y","4-5y","5-10","10-20","20-30","30-40","40-50","50+")    ## define the age groups

N_agegroups = length(Age_groups) ## how many age groups are we considering

N_comms = 4 ## how many communes are we looking at here


age.bins <<- c(1:11, seq(12,59,12), seq(60,119,60), seq(120,720,120)) #in months



generation.time <<- 0.5 #time step is 2 weeks (half a month)

# calculate aging rate

age.lows <- c(0,age.bins[2:length(age.bins)-1])

ac.sz <- age.bins - age.lows

aging.rate <<- generation.time/ac.sz

age.out <- aging.rate

clearVT = (c(.062,.062,.062,.062,.062,.062,.062,.062,.062,.062,.062,.062,
             .062,.12,.12,.12,.34,.34,.34,.34,.34,.34)*2)

clearVT =matrix(clearVT, nrow = N_comms, ncol=N_agegroups, byrow=TRUE) #weekly conversion to daily ## rate of clearing VTs


clearNVT = (c(.086,.086,.086,.086,.086,.086,.086,.086,.086,.086,.086,.086,.086,
               .15,.15,.15,.34,.34,.34,.34,.34,.34)*2)

clearNVT =matrix(clearNVT, nrow=N_comms, ncol=N_agegroups, byrow=TRUE) #weekly conversion to daily ## rate of clearing VTs


ageout <- matrix(age.out, nrow=N_comms, ncol=N_agegroups, byrow=TRUE)


mixing_matrix <- read.csv("AvNumberOfContacts_NT_AP_with1s.csv")

mixing_matrix <- mixing_matrix[,2:23]


sourceCpp("Fit_by_arm.cpp")

#state_prop=c(rep(c(2200,2100,2000,1900,1800,1700,1600,1500,1400,1300,1200,1100,1000,900,800,700,600,500,400,300,200,100), Ncomms),rep(1,N_agegroups*N_comms),rep(1,N_agegroups*N_comms),rep(0,N_agegroups*N_comms)) ## put a certian number of people in each state 


#state_prop=c((c(2200,2100,2000,1900,1800,1700,1600,1500,1400,1300,1200,1100,1000,900,800,700,600,500,400,300,200,100,4200,4100,4000,3900,3800,3700,3600,3500,3400,3300,3200,3100,3000,2900,2800,2700,2600,2500,2400,2300,2200,2100,5200,5100,5000,4900,4800,4700,4600,4500,4400,4300,4200,4100,4000,3900,3800,3700,3600,3500,3400,3300,3200,3100)),
#            rep(c(22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1), N_comms),
#            rep(c(22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1), N_comms),rep(0,N_agegroups*N_comms)) ## put a certian number of people in each state 

##############################
## Source the function to
## simulate

source("Vietnam_prev_function.R")

SIS_ode_cpp(1, state_prop, param.list)



run_mod_test = as.data.frame(lsoda(y=state_prop, times=times, func=SIS_ode_cpp, parms=param.list))

run_mod_test = run_mod_test[,-1]

# solve model steady state
res = runsteady(y=state_prop,
                fun=SIS_ode_cpp,
                parms=param.list,
                times=c(0,1e5))



###############################################
### read in the true data

setwd("H:/My Documents/Carriage_survey_baseline")


VT_all <- read.csv("VT_prev_bycommandage.csv")

VT_age1to2 <- rowSums(VT_all[,14:25], na.rm=TRUE)

VT_all_new <- cbind(VT_all[,5:13], VT_age1to2)

VT_all_new_2 <- cbind(VT_all_new, VT_all[,30:34])


VT_all_new_2[is.na(VT_all_new_2)] <- 0


NVT_all <- read.csv("NVT_prev_bycommandage.csv")

NVT_age1to2 <- rowSums(NVT_all[,14:25], na.rm=TRUE)

NVT_all_new <- cbind(NVT_all[,5:13], NVT_age1to2)

NVT_all_new_2 <- cbind(NVT_all_new, NVT_all[,30:34])

NVT_all_new_2[is.na(NVT_all_new_2)] <- 0



N_samples_tot <- read.csv("N_samples_by_ageandcomm_all.csv")

N_samples_tot <- N_samples_tot[,2:16]


VT_prev_dat <- VT_all_new_2/N_samples_tot

VT_prev_dat[is.na(VT_prev_dat)] <- 0


NVT_prev_dat <- NVT_all_new_2/N_samples_tot 

NVT_prev_dat[is.na(NVT_prev_dat)] <- 0

age1_VT <- rowSums(VT_prev_dat[,1:9])/9
age2_VT <- (VT_prev_dat[,10])/1
age3_VT <- rowSums(VT_prev_dat[,11:15])/5

VT_prev_dat <- cbind(age1_VT,age2_VT,age3_VT)


#########################
## ARM 1 = 0p+1 5,10,16,22,25
## ARM 2 = 1p+1 2,4,14,17,21,27
## ARM 3 = 2p+1 3,8,15,23,24,26
## ARM 4 = 3p+0 6,11,12,13,18,20

VT_prev_dat_arm1 <- rbind(VT_prev_dat[5,],VT_prev_dat[10,],VT_prev_dat[16,],VT_prev_dat[19,],VT_prev_dat[22,],VT_prev_dat[25,])
VT_prev_dat_arm2 <- rbind(VT_prev_dat[2,],VT_prev_dat[4,],VT_prev_dat[14,],VT_prev_dat[17,],VT_prev_dat[21,],VT_prev_dat[27,])
VT_prev_dat_arm3 <- rbind(VT_prev_dat[3,],VT_prev_dat[8,],VT_prev_dat[15,],VT_prev_dat[23,],VT_prev_dat[24,],VT_prev_dat[26,])
VT_prev_dat_arm4 <- rbind(VT_prev_dat[6,],VT_prev_dat[11,],VT_prev_dat[12,],VT_prev_dat[13,],VT_prev_dat[18,],VT_prev_dat[20,])


age1_NVT <- rowSums(NVT_prev_dat[,1:9])/9
age2_NVT <- (NVT_prev_dat[,10])/1
age3_NVT <- rowSums(NVT_prev_dat[,11:15])/5

NVT_prev_dat <- cbind(age1_NVT,age2_NVT,age3_NVT)

#########################
## ARM 1 = 0p+1 5,10,16,22,25
## ARM 2 = 1p+1 2,4,14,17,21,27
## ARM 3 = 2p+1 3,8,15,23,24,26
## ARM 4 = 3p+0 6,11,12,13,18,20

NVT_prev_dat_arm1 <- rbind(NVT_prev_dat[5,],NVT_prev_dat[10,],NVT_prev_dat[16,],NVT_prev_dat[19,],NVT_prev_dat[22,],NVT_prev_dat[25,])
NVT_prev_dat_arm2 <- rbind(NVT_prev_dat[2,],NVT_prev_dat[4,],NVT_prev_dat[14,],NVT_prev_dat[17,],NVT_prev_dat[21,],NVT_prev_dat[27,])
NVT_prev_dat_arm3 <- rbind(NVT_prev_dat[3,],NVT_prev_dat[8,],NVT_prev_dat[15,],NVT_prev_dat[23,],NVT_prev_dat[24,],NVT_prev_dat[26,])
NVT_prev_dat_arm4 <- rbind(NVT_prev_dat[6,],NVT_prev_dat[11,],NVT_prev_dat[12,],NVT_prev_dat[13,],NVT_prev_dat[18,],NVT_prev_dat[20,])


#########################
## ARM 1 = 0p+1 5,10,16,22,25
## ARM 2 = 1p+1 2,4,14,17,21,27
## ARM 3 = 2p+1 3,8,15,23,24,26
## ARM 4 = 3p+0 6,11,12,13,18,20


age1_N <- rowSums(N_samples_tot[,1:9])
age2_N <- (N_samples_tot[,10])
age3_N <- rowSums(N_samples_tot[,11:15])

N_samples_tot <- cbind(age1_N,age2_N,age3_N)

N_samples_tot_arm1 <- rbind(N_samples_tot[5,],N_samples_tot[10,],N_samples_tot[16,],N_samples_tot[19,],N_samples_tot[22,],N_samples_tot[25,])
N_samples_tot_arm2 <- rbind(N_samples_tot[2,],N_samples_tot[4,],N_samples_tot[14,],N_samples_tot[17,],N_samples_tot[21,],N_samples_tot[27,])
N_samples_tot_arm3 <- rbind(N_samples_tot[3,],N_samples_tot[8,],N_samples_tot[15,],N_samples_tot[23,],N_samples_tot[24,],N_samples_tot[26,])
N_samples_tot_arm4 <- rbind(N_samples_tot[6,],N_samples_tot[11,],N_samples_tot[12,],N_samples_tot[13,],N_samples_tot[18,],N_samples_tot[20,])

###############################
## calculate likelihood for prev



calc_Multinomial_log_LL = function(observed_count_of_carriage, modeled_prevalence_of_carriers){
  LL = sum(observed_count_of_carriage * log(modeled_prevalence_of_carriers)) 
  if(any(modeled_prevalence_of_carriers<=0)) LL=-1000001
  return(LL)
}


########################
## median pars
## 
Pneumo_pars <- c( 0.2638921,0.1602888,0.005841182 ,0.3308179, 0.1970636,0.005352119)

Pneumo_pars <- c(0.2558836,0.1476148,0.004430936,0.3205566,0.1840608, 0.003248876) # lower

Pneumo_pars <- c(0.2704093 ,0.1695973,0.007176028, 0.3390404 ,  0.2149428, 0.007949989)
  
##############################
## Evaluate the model likelihood
##

mod_like = function(Pneumo_pars){
  
  
  mod_prev <- prev.function(Pneumo_pars ) ## need to make the vector outputs of this into a matrix
  
  
  S_tot <- matrix (mod_prev[1:88], nrow = N_comms, ncol = N_agegroups, byrow=TRUE)
  VT_tot <- matrix (mod_prev[89:176], nrow = N_comms, ncol = N_agegroups, byrow=TRUE)
  NVT_tot <- matrix (mod_prev[177:264], nrow = N_comms, ncol = N_agegroups, byrow=TRUE)
  B_tot <- matrix (mod_prev[265:352], nrow = N_comms, ncol = N_agegroups, byrow=TRUE)
  
  N_by_ageandcomm <- S_tot + VT_tot + NVT_tot + B_tot
  
  VT_prev_all <- (cbind(VT_tot[,4:13],VT_tot[,18:22]) + cbind(B_tot[,4:13],B_tot[,18:22]))
  
  NVT_prev_all <- (cbind(NVT_tot[,4:13],NVT_tot[,18:22]) + cbind(B_tot[,4:13],B_tot[,18:22]))
  
  N_all <- cbind(N_by_ageandcomm[,4:13],N_by_ageandcomm[,18:22])
  
  VT_sum_prev_mod = cbind(rowSums(VT_tot[,4:12]), (VT_tot[,13]) ,rowSums(VT_tot[,18:22])) + cbind(rowSums(B_tot[,4:12]), cbind(B_tot[,13]), rowSums(B_tot[,18:22]))
  
  NVT_sum_prev_mod = cbind(rowSums(NVT_tot[,4:12]), (NVT_tot[,13]) ,rowSums(NVT_tot[,18:22])) + cbind(rowSums(B_tot[,4:12]), cbind(B_tot[,13]), rowSums(B_tot[,18:22]))
  
  N_all_mod = cbind(rowSums(N_by_ageandcomm[,4:12]) , (N_by_ageandcomm[,13]) , rowSums(N_by_ageandcomm[,18:22])) 
  
  #######################################
  ## mod prev of infection in each age group
  
  VT_prev_model = (VT_sum_prev_mod)/N_all_mod
  
  NVT_prev_model = (NVT_sum_prev_mod)/N_all_mod
  
  N = S_tot+VT_tot+NVT_tot + B_tot
  
  VT_prev_model[  VT_prev_model < 0] <- 0
  NVT_prev_model[  NVT_prev_model < 0] <- 0
  
  
  #### this complete calcualtion of the lieklihood can be done hre we can then just call this mod like function in the MCMC evalution 
  LL_prop=0
  
  comm_log_like_arm1 <- matrix(NA, nrow=6, ncol= 3 )
  
  for(comms in 1:6)
  {
    
    for (age in 1:3) 
    { 
      comm_log_like_arm1[comms,age] = LL_prop +  calc_Multinomial_log_LL(c(VT_prev_dat_arm1[comms,age], NVT_prev_dat_arm1[comms,age], 1-NVT_prev_dat_arm1[comms,age]-VT_prev_dat_arm1[comms,age])*N_samples_tot_arm1[comms,age],
                                                                    c(VT_prev_model[1,age],NVT_prev_model[1,age],1-VT_prev_model[1,age]-NVT_prev_model[1,age]))
    }
  }
  
  
  ######################## 
  ## LL calc for arm 2
  
  LL_prop=0
  
  comm_log_like_arm2 <- matrix(NA, nrow=6, ncol= 3 )
  
  for(comms in 1:6)
  {
    
    for (age in 1:3) 
    { 
      comm_log_like_arm2[comms,age] = LL_prop +  calc_Multinomial_log_LL(c(VT_prev_dat_arm2[comms,age], NVT_prev_dat_arm2[comms,age], 1-NVT_prev_dat_arm2[comms,age]-VT_prev_dat_arm2[comms,age])*N_samples_tot_arm2[comms,age],
                                                                         c(VT_prev_model[2,age],NVT_prev_model[2,age],1-VT_prev_model[2,age]-NVT_prev_model[2,age]))
    }
  }
  
  ######################## 
  ## LL calc for arm 3
  
  LL_prop=0
  
  comm_log_like_arm3 <- matrix(NA, nrow=6, ncol= 3 )
  
  for(comms in 1:6)
  {
    
    for (age in 1:3) 
    { 
      comm_log_like_arm3[comms,age] = LL_prop +  calc_Multinomial_log_LL(c(VT_prev_dat_arm3[comms,age], NVT_prev_dat_arm3[comms,age], 1-NVT_prev_dat_arm3[comms,age]-VT_prev_dat_arm3[comms,age])*N_samples_tot_arm3[comms,age],
                                                                         c(VT_prev_model[3,age],NVT_prev_model[3,age],1-VT_prev_model[3,age]-NVT_prev_model[3,age]))
    }
  }
  
  ######################## 
  ## LL calc for arm 4 
  
  LL_prop=0
  
  comm_log_like_arm4 <- matrix(NA, nrow=6, ncol= 3 )
  

  for(comms in 1:6)
  {
    
    for (age in 1:3) 
    { 
      comm_log_like_arm4[comms,age] = LL_prop +  calc_Multinomial_log_LL(c(VT_prev_dat_arm4[comms,age], NVT_prev_dat_arm4[comms,age], 1-NVT_prev_dat_arm4[comms,age]-VT_prev_dat_arm4[comms,age])*N_samples_tot_arm4[comms,age],
                                                                         c(VT_prev_model[4,age],NVT_prev_model[4,age],1-VT_prev_model[4,age]-NVT_prev_model[4,age]))
    }
  }
  
  
  ###############################
  ## calculate total log-like
  
  total_log_like <- sum(comm_log_like_arm1) + sum(comm_log_like_arm2) + sum(comm_log_like_arm3) + sum(comm_log_like_arm4)
  
  return(total_log_like)
  
  
  
}


library(MASS)

###################################################
## 1.3 PRIOR

prior_M2 <- function( par ){
  
  beta_VT_1    <- par[1]
  beta_VT_2    <- par[2]
  beta_VT_3    <- par[3]

  
  beta_NVT_1    <- par[4]
  beta_NVT_2    <- par[5]
  beta_NVT_3    <- par[6]

  
  ######################################
  ## Uniform prior on prior_beta_1 ~ U(0,1)
  
  if( beta_VT_1>0 && beta_VT_1<1 )
  {
    prior_beta_VT_1 <- 1/1 ### check this value is ok with other code
    
  }else{
    prior_beta_VT_1 <- -1e6
  }
  
  ######################################
  ## Uniform prior on prior_beta_2 ~ U(0,1)
  
  if( beta_VT_2>0 && beta_VT_2<1 )
  {
    prior_beta_VT_2 <- 1/1 ### check this value is ok with other code
    
  }else{
    prior_beta_VT_2 <- -1e6
  }
  
  ######################################
  ## Uniform prior on prior_beta_3 ~ U(0,1)
  
  if( beta_VT_3>0 && beta_VT_3<3 )
  {
    prior_beta_VT_3 <- 1/1 ### check this value is ok with other code
    
  }else{
    prior_beta_VT_3 <- -1e6
  }

  
  ######################################
  ## Uniform prior on prior_beta_1 ~ U(0,1)
  
  if( beta_NVT_1>0 && beta_NVT_1<1 )
  {
    prior_beta_NVT_1 <- 1/1 ### check this value is ok with other code
    
  }else{
    prior_beta_NVT_1 <- -1e6
  }
  
  ######################################
  ## Uniform prior on prior_beta_2 ~ U(0,1)
  
  if( beta_NVT_2>0 && beta_NVT_2<1 )
  {
    prior_beta_NVT_2 <- 1/1 ### check this value is ok with other code
    
  }else{
    prior_beta_NVT_2 <- -1e6
  }
  
  ######################################
  ## Uniform prior on prior_beta_3 ~ U(0,1)
  
  if( beta_NVT_3>0 && beta_NVT_3<3 )
  {
    prior_beta_NVT_3 <- 1/1 ### check this value is ok with other code
    
  }else{
    prior_beta_NVT_3 <- -1e6
  }
  
  
  
  prior <- prior_beta_VT_1 +  prior_beta_VT_2 +  prior_beta_VT_3 +
    prior_beta_NVT_1 +  prior_beta_NVT_2 +  prior_beta_NVT_3 
  
  prior
}




#################################################
#################################################
##          ##                                 ##
##   ####   ##  #     #  ####  #     #  ####   ##
##  ##  ##  ##  ##   ## ##  ## ##   ## ##  ##  ##
##     ##   ##  ####### ##     ####### ##      ##
##    ##    ##  ## # ## ##  ## ## # ## ##  ##  ##
##   #####  ##  ##   ##  ####  ##   ##  ####   ##
##          ##                                 ##
#################################################
#################################################


N_mcmc       <- 10000      ## Number of MCMC iterations
N_tune_start <- 300
N_tune_end   <- 1000
N_adapt      <- 2000

#################################################
## 2.1 Robbins-munro step scaler


step_scale  <- 1
MCMC_accept <- 0


rm_scale <- function(step_scale, mc, log_prob){
  
  dd <- exp(log_prob)
  if( dd < -30 ){ dd <- 0 }
  dd <- min( dd, 1 )
  
  rm_temp <- ( dd - 0.23 )/( (mc+1)/(0.01*N_adapt+1) )
  
  out <- step_scale*exp(rm_temp)
  
  out <- max( out, 0.05 )
  out <- min( out, 5)
  out
}


#################################################
## 2.2 Prepare object for MCMC fitting

MCMC_par           <- matrix(NA, nrow=N_mcmc, ncol=7)
colnames(MCMC_par) <- c("beta_VT_1", "beta_VT_2", "beta_VT_3", "beta_NVT_1", "beta_NVT_2", "beta_NVT_3", "loglike" )


#########################################################
## 2.3 Implement MCMC iterations


par_MC <- c(0.2590076, 0.1652772, 0.006122574,  0.3225996,  0.2149428, 0.005619053)  ## (initial guesses for the difference values of beta)

Sigma_MC <- diag( c(0.008,0.009,0.005,0.008,0.009,0.005)^2 )

min_det <- 0.01*abs(det( Sigma_MC ))


## Ideally fill in your existing best guess
## based on the estimated posterior from a previous burn-in.
## The most important thing is to make sure that things are
## approximately the correct order of magnitude.



loglike_MC <- mod_like( par_MC ) + prior_M2( par_MC )



for(mc in 1:N_mcmc)
{
  par_MCp1 <- mvrnorm(n=1, mu=par_MC, Sigma=step_scale*Sigma_MC)
  
  
  if( par_MCp1[1] > 0 &&
      par_MCp1[2] > 0 &&
      par_MCp1[3] > 0 &&
      par_MCp1[4] > 0 &&
      par_MCp1[5] > 0 &&
      par_MCp1[6] > 0 
        ){
    
    loglike_MCp1 <- mod_like( par_MCp1 ) + prior_M2( par_MCp1 )
    
    log_prob = min( loglike_MCp1-loglike_MC, 0 )           
    
    if( log(runif(1)) < log_prob ) 
    {
      par_MC <- par_MCp1
      
      loglike_MC  <- loglike_MCp1
      MCMC_accept <- MCMC_accept + 1                       
    }
    
    if( mc < N_adapt ){
      step_scale <- rm_scale( step_scale, mc, log_prob)
    }
    
    if( (mc > N_tune_start) && (mc < N_tune_end) )
    {
      cov_MC <- cov( MCMC_par[1:(mc-1),1:6] )
      
      if( abs(det( cov_MC )) > min_det )
      { 
        Sigma_MC <- cov_MC
      }
    }
    
    if( mc%%500 == 0 )
    {
      write.table( MCMC_par, "MCMC_NT_fit.txt")
    }
  }
  
  MCMC_par[mc,1:6] <- par_MC
  MCMC_par[mc,7]   <- loglike_MC
}



#########################################################
## 2.4 Examine MCMC chains


par(mfrow=c(3,3))



#####################################
## PANEL 1: alpha_0 MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,1], 
     pch=19, col="grey", cex=0.25,
     xlab="MCMC iteration", ylab="beta_1_VT", 
     main="beta_1_VT")



#####################################
## PANEL 2: gamma MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,2], 
     ylim=c(0,1),
     pch=19, col="grey", cex=0.25,
     xlab="MCMC iteration", ylab="beta_2_VT", 
     main="beta_2_VT")



#####################################
## PANEL 3 rr MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,3], 
     pch=19, col="grey", cex=0.25,
     xlab="MCMC iteration", ylab="beta_3_VT", 
     main="beta_3_VT")



#####################################
## PANEL 4: time_c MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,4], 
     pch=19, col="grey", cex=0.25,
     xlab="MCMC iteration", ylab="beta_1_NVT", 
     main="beta_1_NVT" )

#####################################
## PANEL 5: time_c MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,5], 
     pch=19, col="grey", cex=0.25,
     xlab="MCMC iteration", ylab="beta_2_NVT", 
     main="beta_2_NVT" )



plot(x=1:N_mcmc, y=MCMC_par[,6], 
     pch=19, col="grey", cex=0.25,
     xlab="MCMC iteration", ylab="beta_3_NVT", 
     main="beta_3_NVT" )

plot(x=1:40000, y=MCMC_burn[,7], 
     pch=19, col="grey", cex=0.25,
     xlab="MCMC iteration", ylab="likelihood", 
     main="likelihood" )





#########################################################
## 2.5 Examine posterior distribution

MCMC_burn = MCMC_par[floor(0.2*nrow(MCMC_par)):(nrow(MCMC_par)-1),]


par(mfrow=c(2,3))


#####################################
## PANEL 1: MCMC posterior

for(k in 1:6)
{
  DEN = density( MCMC_burn[,k] )
  
  QUANT = quantile( MCMC_burn[,k], prob=c(0.025, 0.5, 0.975) )
  
  plot(x=DEN$x, y=DEN$y, type='l',
       xlim=c(0, max(DEN$x)),
       xlab="alpha_0", ylab="", 
       main=colnames(MCMC_par)[k] )
  
  
  low_index  = which(DEN$x<QUANT[1])
  mid_index  = intersect( which(DEN$x>=QUANT[1]), which(DEN$x<=QUANT[3]) )
  high_index = which(DEN$x>QUANT[3])
  
  polygon( x=c( DEN$x[low_index], rev(DEN$x[low_index]) ),
           y=c( rep(0,length(low_index)), rev(DEN$y[low_index]) ), 
           col="pink")
  
  polygon( x=c( DEN$x[mid_index], rev(DEN$x[mid_index]) ),
           y=c( rep(0,length(mid_index)), rev(DEN$y[mid_index]) ), 
           col="grey")
  
  polygon( x=c( DEN$x[high_index], rev(DEN$x[high_index]) ),
           y=c( rep(0,length(high_index)), rev(DEN$y[high_index]) ), 
           col="pink")
  
  points(x=rep(QUANT[2],2), y=c(0,max(DEN$y)), type='l', lty="dashed", lwd=2)
}




par(mfrow=c(3,3))



#####################################
## PANEL 1: alpha_0 MCMC chain

plot(x=1:40000, y=MCMC_burn[,1], 
     pch=19, col="grey", cex=0.25,
     xlab="MCMC iteration", ylab="beta_1_VT", 
     main="beta_1_VT")



#####################################
## PANEL 2: gamma MCMC chain

plot(x=1:40000, y=MCMC_burn[,2], 
     
     pch=19, col="grey", cex=0.25,
     xlab="MCMC iteration", ylab="beta_2_VT", 
     main="beta_2_VT")



#####################################
## PANEL 3 rr MCMC chain

plot(x=1:40000, y=MCMC_burn[,3], 
     pch=19, col="grey", cex=0.25,
     xlab="MCMC iteration", ylab="beta_3_VT", 
     main="beta_3_VT")



#####################################
## PANEL 4: time_c MCMC chain

plot(x=1:40000, y=MCMC_burn[,4], 
     pch=19, col="grey", cex=0.25,
     xlab="MCMC iteration", ylab="beta_1_NVT", 
     main="beta_1_NVT" )

#####################################
## PANEL 5: time_c MCMC chain

plot(x=1:40000, y=MCMC_burn[,5], 
     pch=19, col="grey", cex=0.25,
     xlab="MCMC iteration", ylab="beta_2_NVT", 
     main="beta_2_NVT" )



plot(x=1:40000, y=MCMC_burn[,6], 
     pch=19, col="grey", cex=0.25,
     xlab="MCMC iteration", ylab="beta_3_NVT", 
     main="beta_3_NVT" )













S_inf_comm1_age1 <- (run_mod_test[,(1)])
S_inf_comm1_age_2 <- run_mod_test[,(2)]
S_inf_comm1_age_3 <- run_mod_test[,(3)]
S_inf_comm1_age_4 <- run_mod_test[,(4)]
S_inf_comm1_age_5 <- run_mod_test[,(5)]
S_inf_comm1_age_6 <- run_mod_test[,(6)]
S_inf_comm1_age_7 <- run_mod_test[,(7)]
S_inf_comm1_age_8 <- run_mod_test[,(8)]
S_inf_comm1_age_9 <- run_mod_test[,(9)]
S_inf_comm1_age_10 <- run_mod_test[,(10)]
S_inf_comm1_age_11 <- run_mod_test[,(11)]
S_inf_comm1_age_12 <- run_mod_test[,(12)]
S_inf_comm1_age_13 <- run_mod_test[,(13)]
S_inf_comm1_age_14 <- run_mod_test[,(14)]
S_inf_comm1_age_15 <- run_mod_test[,(15)]
S_inf_comm1_age_16 <- run_mod_test[,(16)]
S_inf_comm1_age_17 <- run_mod_test[,(17)]
S_inf_comm1_age_18 <- run_mod_test[,(18)]
S_inf_comm1_age_19 <- run_mod_test[,(19)]
S_inf_comm1_age_20 <- run_mod_test[,(20)]
S_inf_comm1_age_21 <- run_mod_test[,(21)]
S_inf_comm1_age_22 <- run_mod_test[,(N_agegroups)]

S_inf_comm2_age1 <- (run_mod_test[,(N_agegroups+1)])
S_inf_comm2_age_2 <- run_mod_test[,(N_agegroups+2)]
S_inf_comm2_age_3 <- run_mod_test[,(N_agegroups+3)]
S_inf_comm2_age_4 <- run_mod_test[,(N_agegroups+4)]
S_inf_comm2_age_5 <- run_mod_test[,(N_agegroups+5)]
S_inf_comm2_age_6 <- run_mod_test[,(N_agegroups+6)]
S_inf_comm2_age_7 <- run_mod_test[,(N_agegroups+7)]
S_inf_comm2_age_8 <- run_mod_test[,(N_agegroups+8)]
S_inf_comm2_age_9 <- run_mod_test[,(N_agegroups+9)]
S_inf_comm2_age_10 <- run_mod_test[,(N_agegroups+10)]
S_inf_comm2_age_11 <- run_mod_test[,(N_agegroups+11)]
S_inf_comm2_age_12 <- run_mod_test[,(N_agegroups+12)]
S_inf_comm2_age_13 <- run_mod_test[,(N_agegroups+13)]
S_inf_comm2_age_14 <- run_mod_test[,(N_agegroups+14)]
S_inf_comm2_age_15 <- run_mod_test[,(N_agegroups+15)]
S_inf_comm2_age_16 <- run_mod_test[,(N_agegroups+16)]
S_inf_comm2_age_17 <- run_mod_test[,(N_agegroups+17)]
S_inf_comm2_age_18 <- run_mod_test[,(N_agegroups+18)]
S_inf_comm2_age_19 <- run_mod_test[,(N_agegroups+19)]
S_inf_comm2_age_20 <- run_mod_test[,(N_agegroups+20)]
S_inf_comm2_age_21 <- run_mod_test[,(N_agegroups+21)]
S_inf_comm2_age_22 <- run_mod_test[,(N_agegroups*2)]


S_inf_comm3_age1 <- (run_mod_test[,(N_agegroups*2+1)])
S_inf_comm3_age_2 <- run_mod_test[,(N_agegroups*2+2)]
S_inf_comm3_age_3 <- run_mod_test[,(N_agegroups*2+3)]
S_inf_comm3_age_4 <- run_mod_test[,(N_agegroups*2+4)]
S_inf_comm3_age_5 <- run_mod_test[,(N_agegroups*2+5)]
S_inf_comm3_age_6 <- run_mod_test[,(N_agegroups*2+6)]
S_inf_comm3_age_7 <- run_mod_test[,(N_agegroups*2+7)]
S_inf_comm3_age_8 <- run_mod_test[,(N_agegroups*2+8)]
S_inf_comm3_age_9 <- run_mod_test[,(N_agegroups*2+9)]
S_inf_comm3_age_10 <- run_mod_test[,(N_agegroups*2+10)]
S_inf_comm3_age_11 <- run_mod_test[,(N_agegroups*2+11)]
S_inf_comm3_age_12 <- run_mod_test[,(N_agegroups*2+12)]
S_inf_comm3_age_13 <- run_mod_test[,(N_agegroups*2+13)]
S_inf_comm3_age_14 <- run_mod_test[,(N_agegroups*2+14)]
S_inf_comm3_age_15 <- run_mod_test[,(N_agegroups*2+15)]
S_inf_comm3_age_16 <- run_mod_test[,(N_agegroups*2+16)]
S_inf_comm3_age_17 <- run_mod_test[,(N_agegroups*2+17)]
S_inf_comm3_age_18 <- run_mod_test[,(N_agegroups*2+18)]
S_inf_comm3_age_19 <- run_mod_test[,(N_agegroups*2+19)]
S_inf_comm3_age_20 <- run_mod_test[,(N_agegroups*2+20)]
S_inf_comm3_age_21 <- run_mod_test[,(N_agegroups*2+21)]
S_inf_comm3_age_22 <- run_mod_test[,(N_agegroups*3)]

par(mfrow=c(4,3))

plot(S_inf_comm1_age1, col= "black", type="l", ylim=c(0,5000), xlim=c(0,15000), main= "Susceptibles", xlab= "Time", ylab="N_S")

points(S_inf_comm1_age_2, col="#FF0000FF", lty=2, type="l")
points(S_inf_comm1_age_3, col= "#FF4600FF", lty=2, type="l")
points(S_inf_comm1_age_4, col="#FF8B00FF", lty=2, type="l")
points(S_inf_comm1_age_5, col="#FFD100FF", lty=2, type="l")
points(S_inf_comm1_age_6, col="#E8FF00FF", lty=2, type="l")
points(S_inf_comm1_age_7, col="#A2FF00FF", lty=2, type="l")
points(S_inf_comm1_age_8, col="#5DFF00FF", lty=2, type="l")
points(S_inf_comm1_age_9, col="#17FF00FF", lty=2, type="l")
points(S_inf_comm1_age_10, col="#00FF2EFF", lty=2, type="l")
points(S_inf_comm1_age_11, col="#00FF74FF", lty=2, type="l")
points(S_inf_comm1_age_12, col="#00FFB9FF", lty=2, type="l")
points(S_inf_comm1_age_13, col="#00FFFFFF", lty=2, type="l")
points(S_inf_comm1_age_14, col= "#00B9FFFF", lty=2, type="l")
points(S_inf_comm1_age_15, col="#0074FFFF", lty=2, type="l")
points(S_inf_comm1_age_16, col="#002EFFFF", lty=2, type="l")
points(S_inf_comm1_age_17, col="#1700FFFF", lty=2, type="l")
points(S_inf_comm1_age_18, col= "#5D00FFFF", lty=2, type="l")
points(S_inf_comm1_age_19, col="#A200FFFF" , lty=2, type="l")
points(S_inf_comm1_age_20, col="#E800FFFF", lty=2, type="l")
points(S_inf_comm1_age_21, col="#FF00D1FF", lty=2, type="l")
points(S_inf_comm1_age_22, col="#FF008BFF", lty=2, type="l")

plot(S_inf_comm2_age1, col= "black", type="l", ylim=c(0,5000), xlim=c(0,15000), main= "Susceptibles", xlab= "Time", ylab="N_S")

points(S_inf_comm2_age_2, col="#FF0000FF", lty=2, type="l")
points(S_inf_comm2_age_3, col= "#FF4600FF", lty=2, type="l")
points(S_inf_comm2_age_4, col="#FF8B00FF", lty=2, type="l")
points(S_inf_comm2_age_5, col="#FFD100FF", lty=2, type="l")
points(S_inf_comm2_age_6, col="#E8FF00FF", lty=2, type="l")
points(S_inf_comm2_age_7, col="#A2FF00FF", lty=2, type="l")
points(S_inf_comm2_age_8, col="#5DFF00FF", lty=2, type="l")
points(S_inf_comm2_age_9, col="#17FF00FF", lty=2, type="l")
points(S_inf_comm2_age_10, col="#00FF2EFF", lty=2, type="l")
points(S_inf_comm2_age_11, col="#00FF74FF", lty=2, type="l")
points(S_inf_comm2_age_12, col="#00FFB9FF", lty=2, type="l")
points(S_inf_comm2_age_13, col="#00FFFFFF", lty=2, type="l")
points(S_inf_comm2_age_14, col= "#00B9FFFF", lty=2, type="l")
points(S_inf_comm2_age_15, col="#0074FFFF", lty=2, type="l")
points(S_inf_comm2_age_16, col="#002EFFFF", lty=2, type="l")
points(S_inf_comm2_age_17, col="#1700FFFF", lty=2, type="l")
points(S_inf_comm2_age_18, col= "#5D00FFFF", lty=2, type="l")
points(S_inf_comm2_age_19, col="#A200FFFF" , lty=2, type="l")
points(S_inf_comm2_age_20, col="#E800FFFF", lty=2, type="l")
points(S_inf_comm2_age_21, col="#FF00D1FF", lty=2, type="l")
points(S_inf_comm2_age_22, col="#FF008BFF", lty=2, type="l")

plot(S_inf_comm3_age1, col= "black", type="l", ylim=c(0,5000), xlim=c(0,15000), main= "Susceptibles", xlab= "Time", ylab="N_S")

points(S_inf_comm3_age_2, col="#FF0000FF", lty=2, type="l")
points(S_inf_comm3_age_3, col= "#FF4600FF", lty=2, type="l")
points(S_inf_comm3_age_4, col="#FF8B00FF", lty=2, type="l")
points(S_inf_comm3_age_5, col="#FFD100FF", lty=2, type="l")
points(S_inf_comm3_age_6, col="#E8FF00FF", lty=2, type="l")
points(S_inf_comm3_age_7, col="#A2FF00FF", lty=2, type="l")
points(S_inf_comm3_age_8, col="#5DFF00FF", lty=2, type="l")
points(S_inf_comm3_age_9, col="#17FF00FF", lty=2, type="l")
points(S_inf_comm3_age_10, col="#00FF2EFF", lty=2, type="l")
points(S_inf_comm3_age_11, col="#00FF74FF", lty=2, type="l")
points(S_inf_comm3_age_12, col="#00FFB9FF", lty=2, type="l")
points(S_inf_comm3_age_13, col="#00FFFFFF", lty=2, type="l")
points(S_inf_comm3_age_14, col= "#00B9FFFF", lty=2, type="l")
points(S_inf_comm3_age_15, col="#0074FFFF", lty=2, type="l")
points(S_inf_comm3_age_16, col="#002EFFFF", lty=2, type="l")
points(S_inf_comm3_age_17, col="#1700FFFF", lty=2, type="l")
points(S_inf_comm3_age_18, col= "#5D00FFFF", lty=2, type="l")
points(S_inf_comm3_age_19, col="#A200FFFF" , lty=2, type="l")
points(S_inf_comm3_age_20, col="#E800FFFF", lty=2, type="l")
points(S_inf_comm3_age_21, col="#FF00D1FF", lty=2, type="l")
points(S_inf_comm3_age_22, col="#FF008BFF", lty=2, type="l")


#############################################################
## Plot VT
##


VT_inf_comm1_age1 <- run_mod_test[,(N_agegroups*3+1)]
VT_inf_comm1_age1 <- (run_mod_test[,(N_agegroups*3+1)])
VT_inf_comm1_age_2 <- run_mod_test[,(N_agegroups*3+2)]
VT_inf_comm1_age_3 <- run_mod_test[,(N_agegroups*3+3)]
VT_inf_comm1_age_4 <- run_mod_test[,(N_agegroups*3+4)]
VT_inf_comm1_age_5 <- run_mod_test[,(N_agegroups*3+5)]
VT_inf_comm1_age_6 <- run_mod_test[,(N_agegroups*3+6)]
VT_inf_comm1_age_7 <- run_mod_test[,(N_agegroups*3+7)]
VT_inf_comm1_age_8 <- run_mod_test[,(N_agegroups*3+8)]
VT_inf_comm1_age_9 <- run_mod_test[,(N_agegroups*3+9)]
VT_inf_comm1_age_10 <- run_mod_test[,(N_agegroups*3+10)]
VT_inf_comm1_age_11 <- run_mod_test[,(N_agegroups*3+11)]
VT_inf_comm1_age_12 <- run_mod_test[,(N_agegroups*3+12)]
VT_inf_comm1_age_13 <- run_mod_test[,(N_agegroups*3+13)]
VT_inf_comm1_age_14 <- run_mod_test[,(N_agegroups*3+14)]
VT_inf_comm1_age_15 <- run_mod_test[,(N_agegroups*3+15)]
VT_inf_comm1_age_16 <- run_mod_test[,(N_agegroups*3+16)]
VT_inf_comm1_age_17 <- run_mod_test[,(N_agegroups*3+17)]
VT_inf_comm1_age_18 <- run_mod_test[,(N_agegroups*3+18)]
VT_inf_comm1_age_19 <- run_mod_test[,(N_agegroups*3+19)]
VT_inf_comm1_age_20 <- run_mod_test[,(N_agegroups*3+20)]
VT_inf_comm1_age_21 <- run_mod_test[,(N_agegroups*3+21)]
VT_inf_comm1_age_22 <- run_mod_test[,(N_agegroups*4)]




VT_inf_comm2_age1 <- run_mod_test[,(N_agegroups*4+1)]
VT_inf_comm2_age1 <- (run_mod_test[,(N_agegroups*4+1)])
VT_inf_comm2_age_2 <- run_mod_test[,(N_agegroups*4+2)]
VT_inf_comm2_age_3 <- run_mod_test[,(N_agegroups*4+3)]
VT_inf_comm2_age_4 <- run_mod_test[,(N_agegroups*4+4)]
VT_inf_comm2_age_5 <- run_mod_test[,(N_agegroups*4+5)]
VT_inf_comm2_age_6 <- run_mod_test[,(N_agegroups*4+6)]
VT_inf_comm2_age_7 <- run_mod_test[,(N_agegroups*4+7)]
VT_inf_comm2_age_8 <- run_mod_test[,(N_agegroups*4+8)]
VT_inf_comm2_age_9 <- run_mod_test[,(N_agegroups*4+9)]
VT_inf_comm2_age_10 <- run_mod_test[,(N_agegroups*4+10)]
VT_inf_comm2_age_11 <- run_mod_test[,(N_agegroups*4+11)]
VT_inf_comm2_age_12 <- run_mod_test[,(N_agegroups*4+12)]
VT_inf_comm2_age_13 <- run_mod_test[,(N_agegroups*4+13)]
VT_inf_comm2_age_14 <- run_mod_test[,(N_agegroups*4+14)]
VT_inf_comm2_age_15 <- run_mod_test[,(N_agegroups*4+15)]
VT_inf_comm2_age_16 <- run_mod_test[,(N_agegroups*4+16)]
VT_inf_comm2_age_17 <- run_mod_test[,(N_agegroups*4+17)]
VT_inf_comm2_age_18 <- run_mod_test[,(N_agegroups*4+18)]
VT_inf_comm2_age_19 <- run_mod_test[,(N_agegroups*4+19)]
VT_inf_comm2_age_20 <- run_mod_test[,(N_agegroups*4+20)]
VT_inf_comm2_age_21 <- run_mod_test[,(N_agegroups*4+21)]
VT_inf_comm2_age_22 <- run_mod_test[,(N_agegroups*5)]


VT_inf_comm3_age1 <- run_mod_test[,(N_agegroups*5+1)]
VT_inf_comm3_age1 <- (run_mod_test[,(N_agegroups*5+1)])
VT_inf_comm3_age_2 <- run_mod_test[,(N_agegroups*5+2)]
VT_inf_comm3_age_3 <- run_mod_test[,(N_agegroups*5+3)]
VT_inf_comm3_age_4 <- run_mod_test[,(N_agegroups*5+4)]
VT_inf_comm3_age_5 <- run_mod_test[,(N_agegroups*5+5)]
VT_inf_comm3_age_6 <- run_mod_test[,(N_agegroups*5+6)]
VT_inf_comm3_age_7 <- run_mod_test[,(N_agegroups*5+7)]
VT_inf_comm3_age_8 <- run_mod_test[,(N_agegroups*5+8)]
VT_inf_comm3_age_9 <- run_mod_test[,(N_agegroups*5+9)]
VT_inf_comm3_age_10 <- run_mod_test[,(N_agegroups*5+10)]
VT_inf_comm3_age_11 <- run_mod_test[,(N_agegroups*5+11)]
VT_inf_comm3_age_12 <- run_mod_test[,(N_agegroups*5+12)]
VT_inf_comm3_age_13 <- run_mod_test[,(N_agegroups*5+13)]
VT_inf_comm3_age_14 <- run_mod_test[,(N_agegroups*5+14)]
VT_inf_comm3_age_15 <- run_mod_test[,(N_agegroups*5+15)]
VT_inf_comm3_age_16 <- run_mod_test[,(N_agegroups*5+16)]
VT_inf_comm3_age_17 <- run_mod_test[,(N_agegroups*5+17)]
VT_inf_comm3_age_18 <- run_mod_test[,(N_agegroups*5+18)]
VT_inf_comm3_age_19 <- run_mod_test[,(N_agegroups*5+19)]
VT_inf_comm3_age_20 <- run_mod_test[,(N_agegroups*5+20)]
VT_inf_comm3_age_21 <- run_mod_test[,(N_agegroups*5+21)]
VT_inf_comm3_age_22 <- run_mod_test[,(N_agegroups*6)]


plot(VT_inf_comm1_age1, col= "black", type="l", ylim=c(0,2500), xlim=c(0,15000))

points(VT_inf_comm1_age_2, col="#FF0000FF", lty=2, type="l")
points(VT_inf_comm1_age_3, col= "#FF4600FF", lty=2, type="l")
points(VT_inf_comm1_age_4, col="#FF8B00FF", lty=2, type="l")
points(VT_inf_comm1_age_5, col="#FFD100FF", lty=2, type="l")
points(VT_inf_comm1_age_6, col="#E8FF00FF", lty=2, type="l")
points(VT_inf_comm1_age_7, col="#A2FF00FF", lty=2, type="l")
points(VT_inf_comm1_age_8, col="#5DFF00FF", lty=2, type="l")
points(VT_inf_comm1_age_9, col="#17FF00FF", lty=2, type="l")
points(VT_inf_comm1_age_10, col="#00FF2EFF", lty=2, type="l")
points(VT_inf_comm1_age_11, col="#00FF74FF", lty=2, type="l")
points(VT_inf_comm1_age_12, col="#00FFB9FF", lty=2, type="l")
points(VT_inf_comm1_age_13, col="#00FFFFFF", lty=2, type="l")
points(VT_inf_comm1_age_14, col= "#00B9FFFF", lty=2, type="l")
points(VT_inf_comm1_age_15, col="#0074FFFF", lty=2, type="l")
points(VT_inf_comm1_age_16, col="#002EFFFF", lty=2, type="l")
points(VT_inf_comm1_age_17, col="#1700FFFF", lty=2, type="l")
points(VT_inf_comm1_age_18, col= "#5D00FFFF", lty=2, type="l")
points(VT_inf_comm1_age_19, col="#A200FFFF" , lty=2, type="l")
points(VT_inf_comm1_age_20, col="#E800FFFF", lty=2, type="l")
points(VT_inf_comm1_age_21, col="#FF00D1FF", lty=2, type="l")
points(VT_inf_comm1_age_22, col="#FF008BFF", lty=2, type="l")


plot(VT_inf_comm2_age1, col= "black", type="l", ylim=c(0,2500), xlim=c(0,15000))

points(VT_inf_comm2_age_2, col="#FF0000FF", lty=2, type="l")
points(VT_inf_comm2_age_3, col= "#FF4600FF", lty=2, type="l")
points(VT_inf_comm2_age_4, col="#FF8B00FF", lty=2, type="l")
points(VT_inf_comm2_age_5, col="#FFD100FF", lty=2, type="l")
points(VT_inf_comm2_age_6, col="#E8FF00FF", lty=2, type="l")
points(VT_inf_comm2_age_7, col="#A2FF00FF", lty=2, type="l")
points(VT_inf_comm2_age_8, col="#5DFF00FF", lty=2, type="l")
points(VT_inf_comm2_age_9, col="#17FF00FF", lty=2, type="l")
points(VT_inf_comm2_age_10, col="#00FF2EFF", lty=2, type="l")
points(VT_inf_comm2_age_11, col="#00FF74FF", lty=2, type="l")
points(VT_inf_comm2_age_12, col="#00FFB9FF", lty=2, type="l")
points(VT_inf_comm2_age_13, col="#00FFFFFF", lty=2, type="l")
points(VT_inf_comm2_age_14, col= "#00B9FFFF", lty=2, type="l")
points(VT_inf_comm2_age_15, col="#0074FFFF", lty=2, type="l")
points(VT_inf_comm2_age_16, col="#002EFFFF", lty=2, type="l")
points(VT_inf_comm2_age_17, col="#1700FFFF", lty=2, type="l")
points(VT_inf_comm2_age_18, col= "#5D00FFFF", lty=2, type="l")
points(VT_inf_comm2_age_19, col="#A200FFFF" , lty=2, type="l")
points(VT_inf_comm2_age_20, col="#E800FFFF", lty=2, type="l")
points(VT_inf_comm2_age_21, col="#FF00D1FF", lty=2, type="l")
points(VT_inf_comm2_age_22, col="#FF008BFF", lty=2, type="l")


plot(VT_inf_comm3_age1, col= "black", type="l", ylim=c(0,2500), xlim=c(0,15000))

points(VT_inf_comm3_age_2, col="#FF0000FF", lty=2, type="l")
points(VT_inf_comm3_age_3, col= "#FF4600FF", lty=2, type="l")
points(VT_inf_comm3_age_4, col="#FF8B00FF", lty=2, type="l")
points(VT_inf_comm3_age_5, col="#FFD100FF", lty=2, type="l")
points(VT_inf_comm3_age_6, col="#E8FF00FF", lty=2, type="l")
points(VT_inf_comm3_age_7, col="#A2FF00FF", lty=2, type="l")
points(VT_inf_comm3_age_8, col="#5DFF00FF", lty=2, type="l")
points(VT_inf_comm3_age_9, col="#17FF00FF", lty=2, type="l")
points(VT_inf_comm3_age_10, col="#00FF2EFF", lty=2, type="l")
points(VT_inf_comm3_age_11, col="#00FF74FF", lty=2, type="l")
points(VT_inf_comm3_age_12, col="#00FFB9FF", lty=2, type="l")
points(VT_inf_comm3_age_13, col="#00FFFFFF", lty=2, type="l")
points(VT_inf_comm3_age_14, col= "#00B9FFFF", lty=2, type="l")
points(VT_inf_comm3_age_15, col="#0074FFFF", lty=2, type="l")
points(VT_inf_comm3_age_16, col="#002EFFFF", lty=2, type="l")
points(VT_inf_comm3_age_17, col="#1700FFFF", lty=2, type="l")
points(VT_inf_comm3_age_18, col= "#5D00FFFF", lty=2, type="l")
points(VT_inf_comm3_age_19, col="#A200FFFF" , lty=2, type="l")
points(VT_inf_comm3_age_20, col="#E800FFFF", lty=2, type="l")
points(VT_inf_comm3_age_21, col="#FF00D1FF", lty=2, type="l")
points(VT_inf_comm3_age_22, col="#FF008BFF", lty=2, type="l")


##################################
## For NVTs
##
NVT_inf_comm1_age1 <- run_mod_test[,(N_agegroups*6+1)]
NVT_inf_comm1_age1 <- (run_mod_test[,(N_agegroups*6+1)])
NVT_inf_comm1_age_2 <- run_mod_test[,(N_agegroups*6+2)]
NVT_inf_comm1_age_3 <- run_mod_test[,(N_agegroups*6+3)]
NVT_inf_comm1_age_4 <- run_mod_test[,(N_agegroups*6+4)]
NVT_inf_comm1_age_5 <- run_mod_test[,(N_agegroups*6+5)]
NVT_inf_comm1_age_6 <- run_mod_test[,(N_agegroups*6+6)]
NVT_inf_comm1_age_7 <- run_mod_test[,(N_agegroups*6+7)]
NVT_inf_comm1_age_8 <- run_mod_test[,(N_agegroups*6+8)]
NVT_inf_comm1_age_9 <- run_mod_test[,(N_agegroups*6+9)]
NVT_inf_comm1_age_10 <- run_mod_test[,(N_agegroups*6+10)]
NVT_inf_comm1_age_11 <- run_mod_test[,(N_agegroups*6+11)]
NVT_inf_comm1_age_12 <- run_mod_test[,(N_agegroups*6+12)]
NVT_inf_comm1_age_13 <- run_mod_test[,(N_agegroups*6+13)]
NVT_inf_comm1_age_14 <- run_mod_test[,(N_agegroups*6+14)]
NVT_inf_comm1_age_15 <- run_mod_test[,(N_agegroups*6+15)]
NVT_inf_comm1_age_16 <- run_mod_test[,(N_agegroups*6+16)]
NVT_inf_comm1_age_17 <- run_mod_test[,(N_agegroups*6+17)]
NVT_inf_comm1_age_18 <- run_mod_test[,(N_agegroups*6+18)]
NVT_inf_comm1_age_19 <- run_mod_test[,(N_agegroups*6+19)]
NVT_inf_comm1_age_20 <- run_mod_test[,(N_agegroups*6+20)]
NVT_inf_comm1_age_21 <- run_mod_test[,(N_agegroups*6+21)]
NVT_inf_comm1_age_22 <- run_mod_test[,(N_agegroups*7)]

NVT_inf_comm2_age1 <- run_mod_test[,(N_agegroups*7+1)]
NVT_inf_comm2_age1 <- (run_mod_test[,(N_agegroups*7+1)])
NVT_inf_comm2_age_2 <- run_mod_test[,(N_agegroups*7+2)]
NVT_inf_comm2_age_3 <- run_mod_test[,(N_agegroups*7+3)]
NVT_inf_comm2_age_4 <- run_mod_test[,(N_agegroups*7+4)]
NVT_inf_comm2_age_5 <- run_mod_test[,(N_agegroups*7+5)]
NVT_inf_comm2_age_6 <- run_mod_test[,(N_agegroups*7+6)]
NVT_inf_comm2_age_7 <- run_mod_test[,(N_agegroups*7+7)]
NVT_inf_comm2_age_8 <- run_mod_test[,(N_agegroups*7+8)]
NVT_inf_comm2_age_9 <- run_mod_test[,(N_agegroups*7+9)]
NVT_inf_comm2_age_10 <- run_mod_test[,(N_agegroups*7+10)]
NVT_inf_comm2_age_11 <- run_mod_test[,(N_agegroups*7+11)]
NVT_inf_comm2_age_12 <- run_mod_test[,(N_agegroups*7+12)]
NVT_inf_comm2_age_13 <- run_mod_test[,(N_agegroups*7+13)]
NVT_inf_comm2_age_14 <- run_mod_test[,(N_agegroups*7+14)]
NVT_inf_comm2_age_15 <- run_mod_test[,(N_agegroups*7+15)]
NVT_inf_comm2_age_16 <- run_mod_test[,(N_agegroups*7+16)]
NVT_inf_comm2_age_17 <- run_mod_test[,(N_agegroups*7+17)]
NVT_inf_comm2_age_18 <- run_mod_test[,(N_agegroups*7+18)]
NVT_inf_comm2_age_19 <- run_mod_test[,(N_agegroups*7+19)]
NVT_inf_comm2_age_20 <- run_mod_test[,(N_agegroups*7+20)]
NVT_inf_comm2_age_21 <- run_mod_test[,(N_agegroups*7+21)]
NVT_inf_comm2_age_22 <- run_mod_test[,(N_agegroups*8)]

NVT_inf_comm3_age1 <- run_mod_test[,(N_agegroups*8+1)]
NVT_inf_comm3_age1 <- (run_mod_test[,(N_agegroups*8+1)])
NVT_inf_comm3_age_2 <- run_mod_test[,(N_agegroups*8+2)]
NVT_inf_comm3_age_3 <- run_mod_test[,(N_agegroups*8+3)]
NVT_inf_comm3_age_4 <- run_mod_test[,(N_agegroups*8+4)]
NVT_inf_comm3_age_5 <- run_mod_test[,(N_agegroups*8+5)]
NVT_inf_comm3_age_6 <- run_mod_test[,(N_agegroups*8+6)]
NVT_inf_comm3_age_7 <- run_mod_test[,(N_agegroups*8+7)]
NVT_inf_comm3_age_8 <- run_mod_test[,(N_agegroups*8+8)]
NVT_inf_comm3_age_9 <- run_mod_test[,(N_agegroups*8+9)]
NVT_inf_comm3_age_10 <- run_mod_test[,(N_agegroups*8+10)]
NVT_inf_comm3_age_11 <- run_mod_test[,(N_agegroups*8+11)]
NVT_inf_comm3_age_12 <- run_mod_test[,(N_agegroups*8+12)]
NVT_inf_comm3_age_13 <- run_mod_test[,(N_agegroups*8+13)]
NVT_inf_comm3_age_14 <- run_mod_test[,(N_agegroups*8+14)]
NVT_inf_comm3_age_15 <- run_mod_test[,(N_agegroups*8+15)]
NVT_inf_comm3_age_16 <- run_mod_test[,(N_agegroups*8+16)]
NVT_inf_comm3_age_17 <- run_mod_test[,(N_agegroups*8+17)]
NVT_inf_comm3_age_18 <- run_mod_test[,(N_agegroups*8+18)]
NVT_inf_comm3_age_19 <- run_mod_test[,(N_agegroups*8+19)]
NVT_inf_comm3_age_20 <- run_mod_test[,(N_agegroups*8+20)]
NVT_inf_comm3_age_21 <- run_mod_test[,(N_agegroups*8+21)]
NVT_inf_comm3_age_22 <- run_mod_test[,(N_agegroups*9)]

plot(NVT_inf_comm1_age1, col= "black", type="l", ylim=c(0,2500), xlim=c(0,15000))

points(NVT_inf_comm1_age_2, col="#FF0000FF", lty=2, type="l")
points(NVT_inf_comm1_age_3, col= "#FF4600FF", lty=2, type="l")
points(NVT_inf_comm1_age_4, col="#FF8B00FF", lty=2, type="l")
points(NVT_inf_comm1_age_5, col="#FFD100FF", lty=2, type="l")
points(NVT_inf_comm1_age_6, col="#E8FF00FF", lty=2, type="l")
points(NVT_inf_comm1_age_7, col="#A2FF00FF", lty=2, type="l")

points(NVT_inf_comm1_age_8, col="#5DFF00FF", lty=2, type="l")
points(NVT_inf_comm1_age_9, col="#17FF00FF", lty=2, type="l")
points(NVT_inf_comm1_age_10, col="#00FF2EFF", lty=2, type="l")
points(NVT_inf_comm1_age_11, col="#00FF74FF", lty=2, type="l")

points(NVT_inf_comm1_age_12, col="#00FFB9FF", lty=2, type="l")
points(NVT_inf_comm1_age_13, col="#00FFFFFF", lty=2, type="l")
points(NVT_inf_comm1_age_14, col= "#00B9FFFF", lty=2, type="l")
points(NVT_inf_comm1_age_15, col="#0074FFFF", lty=2, type="l")
points(NVT_inf_comm1_age_16, col="#002EFFFF", lty=2, type="l")
points(NVT_inf_comm1_age_17, col="#1700FFFF", lty=2, type="l")
points(NVT_inf_comm1_age_18, col= "#5D00FFFF", lty=2, type="l")
points(NVT_inf_comm1_age_19, col="#A200FFFF" , lty=2, type="l")
points(NVT_inf_comm1_age_20, col="#E800FFFF", lty=2, type="l")
points(NVT_inf_comm1_age_21, col="#FF00D1FF", lty=2, type="l")
points(NVT_inf_comm1_age_22, col="#FF008BFF", lty=2, type="l")


plot(NVT_inf_comm2_age1, col= "black", type="l", ylim=c(0,2500), xlim=c(0,15000))

points(NVT_inf_comm2_age_2, col="#FF0000FF", lty=2, type="l")
points(NVT_inf_comm2_age_3, col= "#FF4600FF", lty=2, type="l")
points(NVT_inf_comm2_age_4, col="#FF8B00FF", lty=2, type="l")
points(NVT_inf_comm2_age_5, col="#FFD100FF", lty=2, type="l")
points(NVT_inf_comm2_age_6, col="#E8FF00FF", lty=2, type="l")
points(NVT_inf_comm2_age_7, col="#A2FF00FF", lty=2, type="l")
points(NVT_inf_comm2_age_8, col="#5DFF00FF", lty=2, type="l")
points(NVT_inf_comm2_age_9, col="#17FF00FF", lty=2, type="l")
points(NVT_inf_comm2_age_10, col="#00FF2EFF", lty=2, type="l")
points(NVT_inf_comm2_age_11, col="#00FF74FF", lty=2, type="l")
points(NVT_inf_comm2_age_12, col="#00FFB9FF", lty=2, type="l")
points(NVT_inf_comm2_age_13, col="#00FFFFFF", lty=2, type="l")
points(NVT_inf_comm2_age_14, col= "#00B9FFFF", lty=2, type="l")
points(NVT_inf_comm2_age_15, col="#0074FFFF", lty=2, type="l")
points(NVT_inf_comm2_age_16, col="#002EFFFF", lty=2, type="l")
points(NVT_inf_comm2_age_17, col="#1700FFFF", lty=2, type="l")
points(NVT_inf_comm2_age_18, col= "#5D00FFFF", lty=2, type="l")
points(NVT_inf_comm2_age_19, col="#A200FFFF" , lty=2, type="l")
points(NVT_inf_comm2_age_20, col="#E800FFFF", lty=2, type="l")
points(NVT_inf_comm2_age_21, col="#FF00D1FF", lty=2, type="l")
points(NVT_inf_comm2_age_22, col="#FF008BFF", lty=2, type="l")


plot(NVT_inf_comm3_age1, col= "black", type="l", ylim=c(0,2500), xlim=c(0,15000))

points(NVT_inf_comm3_age_2, col="#FF0000FF", lty=2, type="l")
points(NVT_inf_comm3_age_3, col= "#FF4600FF", lty=2, type="l")
points(NVT_inf_comm3_age_4, col="#FF8B00FF", lty=2, type="l")
points(NVT_inf_comm3_age_5, col="#FFD100FF", lty=2, type="l")
points(NVT_inf_comm3_age_6, col="#E8FF00FF", lty=2, type="l")
points(NVT_inf_comm3_age_7, col="#A2FF00FF", lty=2, type="l")
points(NVT_inf_comm3_age_8, col="#5DFF00FF", lty=2, type="l")
points(NVT_inf_comm3_age_9, col="#17FF00FF", lty=2, type="l")
points(NVT_inf_comm3_age_10, col="#00FF2EFF", lty=2, type="l")
points(NVT_inf_comm3_age_11, col="#00FF74FF", lty=2, type="l")
points(NVT_inf_comm3_age_12, col="#00FFB9FF", lty=2, type="l")
points(NVT_inf_comm3_age_13, col="#00FFFFFF", lty=2, type="l")
points(NVT_inf_comm3_age_14, col= "#00B9FFFF", lty=2, type="l")
points(NVT_inf_comm3_age_15, col="#0074FFFF", lty=2, type="l")
points(NVT_inf_comm3_age_16, col="#002EFFFF", lty=2, type="l")
points(NVT_inf_comm3_age_17, col="#1700FFFF", lty=2, type="l")
points(NVT_inf_comm3_age_18, col= "#5D00FFFF", lty=2, type="l")
points(NVT_inf_comm3_age_19, col="#A200FFFF" , lty=2, type="l")
points(NVT_inf_comm3_age_20, col="#E800FFFF", lty=2, type="l")
points(NVT_inf_comm3_age_21, col="#FF00D1FF", lty=2, type="l")
points(NVT_inf_comm3_age_22, col="#FF008BFF", lty=2, type="l")

########################################
## plot for both
##

B_inf_comm1_age1 <- run_mod_test[,(N_agegroups*9+1)]
B_inf_comm1_age1 <- (run_mod_test[,(N_agegroups*9+1)])
B_inf_comm1_age_2 <- run_mod_test[,(N_agegroups*9+2)]
B_inf_comm1_age_3 <- run_mod_test[,(N_agegroups*9+3)]
B_inf_comm1_age_4 <- run_mod_test[,(N_agegroups*9+4)]
B_inf_comm1_age_5 <- run_mod_test[,(N_agegroups*9+5)]
B_inf_comm1_age_6 <- run_mod_test[,(N_agegroups*9+6)]
B_inf_comm1_age_7 <- run_mod_test[,(N_agegroups*9+7)]
B_inf_comm1_age_8 <- run_mod_test[,(N_agegroups*9+8)]
B_inf_comm1_age_9 <- run_mod_test[,(N_agegroups*9+9)]
B_inf_comm1_age_10 <- run_mod_test[,(N_agegroups*9+10)]
B_inf_comm1_age_11 <- run_mod_test[,(N_agegroups*9+11)]
B_inf_comm1_age_12 <- run_mod_test[,(N_agegroups*9+12)]
B_inf_comm1_age_13 <- run_mod_test[,(N_agegroups*9+13)]
B_inf_comm1_age_14 <- run_mod_test[,(N_agegroups*9+14)]
B_inf_comm1_age_15 <- run_mod_test[,(N_agegroups*9+15)]
B_inf_comm1_age_16 <- run_mod_test[,(N_agegroups*9+16)]
B_inf_comm1_age_17 <- run_mod_test[,(N_agegroups*9+17)]
B_inf_comm1_age_18 <- run_mod_test[,(N_agegroups*9+18)]
B_inf_comm1_age_19 <- run_mod_test[,(N_agegroups*9+19)]
B_inf_comm1_age_20 <- run_mod_test[,(N_agegroups*9+20)]
B_inf_comm1_age_21 <- run_mod_test[,(N_agegroups*9+21)]
B_inf_comm1_age_22 <- run_mod_test[,(N_agegroups*10)]

B_inf_comm2_age1 <- run_mod_test[,(N_agegroups*10+1)]
B_inf_comm2_age1 <- (run_mod_test[,(N_agegroups*10+1)])
B_inf_comm2_age_2 <- run_mod_test[,(N_agegroups*10+2)]
B_inf_comm2_age_3 <- run_mod_test[,(N_agegroups*10+3)]
B_inf_comm2_age_4 <- run_mod_test[,(N_agegroups*10+4)]
B_inf_comm2_age_5 <- run_mod_test[,(N_agegroups*10+5)]
B_inf_comm2_age_6 <- run_mod_test[,(N_agegroups*10+6)]
B_inf_comm2_age_7 <- run_mod_test[,(N_agegroups*10+7)]
B_inf_comm2_age_8 <- run_mod_test[,(N_agegroups*10+8)]
B_inf_comm2_age_9 <- run_mod_test[,(N_agegroups*10+9)]
B_inf_comm2_age_10 <- run_mod_test[,(N_agegroups*10+10)]
B_inf_comm2_age_11 <- run_mod_test[,(N_agegroups*10+11)]
B_inf_comm2_age_12 <- run_mod_test[,(N_agegroups*10+12)]
B_inf_comm2_age_13 <- run_mod_test[,(N_agegroups*10+13)]
B_inf_comm2_age_14 <- run_mod_test[,(N_agegroups*10+14)]
B_inf_comm2_age_15 <- run_mod_test[,(N_agegroups*10+15)]
B_inf_comm2_age_16 <- run_mod_test[,(N_agegroups*10+16)]
B_inf_comm2_age_17 <- run_mod_test[,(N_agegroups*10+17)]
B_inf_comm2_age_18 <- run_mod_test[,(N_agegroups*10+18)]
B_inf_comm2_age_19 <- run_mod_test[,(N_agegroups*10+19)]
B_inf_comm2_age_20 <- run_mod_test[,(N_agegroups*10+20)]
B_inf_comm2_age_21 <- run_mod_test[,(N_agegroups*10+21)]
B_inf_comm2_age_22 <- run_mod_test[,(N_agegroups*11)]


B_inf_comm3_age1 <- run_mod_test[,(N_agegroups*11+1)]
B_inf_comm3_age1 <- (run_mod_test[,(N_agegroups*11+1)])
B_inf_comm3_age_2 <- run_mod_test[,(N_agegroups*11+2)]
B_inf_comm3_age_3 <- run_mod_test[,(N_agegroups*11+3)]
B_inf_comm3_age_4 <- run_mod_test[,(N_agegroups*11+4)]
B_inf_comm3_age_5 <- run_mod_test[,(N_agegroups*11+5)]
B_inf_comm3_age_6 <- run_mod_test[,(N_agegroups*11+6)]
B_inf_comm3_age_7 <- run_mod_test[,(N_agegroups*11+7)]
B_inf_comm3_age_8 <- run_mod_test[,(N_agegroups*11+8)]
B_inf_comm3_age_9 <- run_mod_test[,(N_agegroups*11+9)]
B_inf_comm3_age_10 <- run_mod_test[,(N_agegroups*11+10)]
B_inf_comm3_age_11 <- run_mod_test[,(N_agegroups*11+11)]
B_inf_comm3_age_12 <- run_mod_test[,(N_agegroups*11+12)]
B_inf_comm3_age_13 <- run_mod_test[,(N_agegroups*11+13)]
B_inf_comm3_age_14 <- run_mod_test[,(N_agegroups*11+14)]
B_inf_comm3_age_15 <- run_mod_test[,(N_agegroups*11+15)]
B_inf_comm3_age_16 <- run_mod_test[,(N_agegroups*11+16)]
B_inf_comm3_age_17 <- run_mod_test[,(N_agegroups*11+17)]
B_inf_comm3_age_18 <- run_mod_test[,(N_agegroups*11+18)]
B_inf_comm3_age_19 <- run_mod_test[,(N_agegroups*11+19)]
B_inf_comm3_age_20 <- run_mod_test[,(N_agegroups*11+20)]
B_inf_comm3_age_21 <- run_mod_test[,(N_agegroups*11+21)]
B_inf_comm3_age_22 <- run_mod_test[,(N_agegroups*12)]


plot(B_inf_comm1_age1, col= "black", type="l", ylim=c(0,2500), xlim=c(0,15000))

points(B_inf_comm1_age_2, col="#FF0000FF", lty=2, type="l")
points(B_inf_comm1_age_3, col= "#FF4600FF", lty=2, type="l")
points(B_inf_comm1_age_4, col="#FF8B00FF", lty=2, type="l")
points(B_inf_comm1_age_5, col="#FFD100FF", lty=2, type="l")
points(B_inf_comm1_age_6, col="#E8FF00FF", lty=2, type="l")
points(B_inf_comm1_age_7, col="#A2FF00FF", lty=2, type="l")
points(B_inf_comm1_age_8, col="#5DFF00FF", lty=2, type="l")
points(B_inf_comm1_age_9, col="#17FF00FF", lty=2, type="l")
points(B_inf_comm1_age_10, col="#00FF2EFF", lty=2, type="l")
points(B_inf_comm1_age_11, col="#00FF74FF", lty=2, type="l")
points(B_inf_comm1_age_12, col="#00FFB9FF", lty=2, type="l")
points(B_inf_comm1_age_13, col="#00FFFFFF", lty=2, type="l")
points(B_inf_comm1_age_14, col= "#00B9FFFF", lty=2, type="l")
points(B_inf_comm1_age_15, col="#0074FFFF", lty=2, type="l")
points(B_inf_comm1_age_16, col="#002EFFFF", lty=2, type="l")
points(B_inf_comm1_age_17, col="#1700FFFF", lty=2, type="l")
points(B_inf_comm1_age_18, col= "#5D00FFFF", lty=2, type="l")
points(B_inf_comm1_age_19, col="#A200FFFF" , lty=2, type="l")
points(B_inf_comm1_age_20, col="#E800FFFF", lty=2, type="l")
points(B_inf_comm1_age_21, col="#FF00D1FF", lty=2, type="l")
points(B_inf_comm1_age_22, col="#FF008BFF", lty=2, type="l")


plot(B_inf_comm2_age1, col= "black", type="l", ylim=c(0,2500), xlim=c(0,15000))

points(B_inf_comm2_age_2, col="#FF0000FF", lty=2, type="l")
points(B_inf_comm2_age_3, col= "#FF4600FF", lty=2, type="l")
points(B_inf_comm2_age_4, col="#FF8B00FF", lty=2, type="l")
points(B_inf_comm2_age_5, col="#FFD100FF", lty=2, type="l")
points(B_inf_comm2_age_6, col="#E8FF00FF", lty=2, type="l")
points(B_inf_comm2_age_7, col="#A2FF00FF", lty=2, type="l")
points(B_inf_comm2_age_8, col="#5DFF00FF", lty=2, type="l")
points(B_inf_comm2_age_9, col="#17FF00FF", lty=2, type="l")
points(B_inf_comm2_age_10, col="#00FF2EFF", lty=2, type="l")
points(B_inf_comm2_age_11, col="#00FF74FF", lty=2, type="l")
points(B_inf_comm2_age_12, col="#00FFB9FF", lty=2, type="l")
points(B_inf_comm2_age_13, col="#00FFFFFF", lty=2, type="l")
points(B_inf_comm2_age_14, col= "#00B9FFFF", lty=2, type="l")
points(B_inf_comm2_age_15, col="#0074FFFF", lty=2, type="l")
points(B_inf_comm2_age_16, col="#002EFFFF", lty=2, type="l")
points(B_inf_comm2_age_17, col="#1700FFFF", lty=2, type="l")
points(B_inf_comm2_age_18, col= "#5D00FFFF", lty=2, type="l")
points(B_inf_comm2_age_19, col="#A200FFFF" , lty=2, type="l")
points(B_inf_comm2_age_20, col="#E800FFFF", lty=2, type="l")
points(B_inf_comm2_age_21, col="#FF00D1FF", lty=2, type="l")
points(B_inf_comm2_age_22, col="#FF008BFF", lty=2, type="l")






plot(B_inf_comm3_age1, col= "black", type="l", ylim=c(0,2500), xlim=c(0,15000))

points(B_inf_comm3_age_2, col="#FF0000FF", lty=2, type="l")
points(B_inf_comm3_age_3, col= "#FF4600FF", lty=2, type="l")
points(B_inf_comm3_age_4, col="#FF8B00FF", lty=2, type="l")
points(B_inf_comm3_age_5, col="#FFD100FF", lty=2, type="l")
points(B_inf_comm3_age_6, col="#E8FF00FF", lty=2, type="l")
points(B_inf_comm3_age_7, col="#A2FF00FF", lty=2, type="l")
points(B_inf_comm3_age_8, col="#5DFF00FF", lty=2, type="l")
points(B_inf_comm3_age_9, col="#17FF00FF", lty=2, type="l")
points(B_inf_comm3_age_10, col="#00FF2EFF", lty=2, type="l")
points(B_inf_comm3_age_11, col="#00FF74FF", lty=2, type="l")
points(B_inf_comm3_age_12, col="#00FFB9FF", lty=2, type="l")
points(B_inf_comm3_age_13, col="#00FFFFFF", lty=2, type="l")
points(B_inf_comm3_age_14, col= "#00B9FFFF", lty=2, type="l")
points(B_inf_comm3_age_15, col="#0074FFFF", lty=2, type="l")
points(B_inf_comm3_age_16, col="#002EFFFF", lty=2, type="l")
points(B_inf_comm3_age_17, col="#1700FFFF", lty=2, type="l")
points(B_inf_comm3_age_18, col= "#5D00FFFF", lty=2, type="l")
points(B_inf_comm3_age_19, col="#A200FFFF" , lty=2, type="l")
points(B_inf_comm3_age_20, col="#E800FFFF", lty=2, type="l")
points(B_inf_comm3_age_21, col="#FF00D1FF", lty=2, type="l")
points(B_inf_comm3_age_22, col="#FF008BFF", lty=2, type="l")


