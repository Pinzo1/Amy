prev.function <- function(par_ode){
  
  ###############################
  ## declare init variables using endem equi data
  
  #par_ode <- c(0.5883626, 0.2182841, 0.04778841, 0.5098258, 0.2072802, 0.04003783)
  
  #dS  = matrix  (0, nrow = N_comms, ncol = N_agegroups) ## S
  #dVT  = matrix (0, nrow = N_comms, ncol = N_agegroups) ## I
  #dNVT = matrix (0, nrow = N_comms, ncol = N_agegroups) ## IA
  #dB  = matrix  (0, nrow = N_comms, ncol = N_agegroups) ## A
  #dV  = matrix  (0, nrow = N_comms, ncol = N_agegroups) ## A
  
  #S 	= matrix( x[(0*N_comms*N_agegroups+1) : (1*N_comms*N_agegroups)], nrow = N_comms, ncol = N_agegroups, byrow = TRUE)
  #VT 	= matrix( x[(1*N_comms*N_agegroups+1) : (2*N_comms*N_agegroups)], nrow = N_comms, ncol = N_agegroups, byrow = TRUE)
  #NVT = matrix( x[(2*N_comms*N_agegroups+1) : (3*N_comms*N_agegroups)], nrow = N_comms, ncol = N_agegroups, byrow = TRUE)
  #B 	= matrix( x[(3*N_comms*N_agegroups+1) : (4*N_comms*N_agegroups)], nrow = N_comms, ncol = N_agegroups, byrow = TRUE)
  #V 	= matrix( x[(4*N_comms*N_agegroups+1) : (5*N_comms*N_agegroups)], nrow = N_comms, ncol = N_agegroups, byrow = TRUE)
  
  
  
  betaVT <- as.vector(c(par_ode[1],par_ode[1],par_ode[1],par_ode[1],par_ode[1],par_ode[1],par_ode[1],par_ode[1],par_ode[1],par_ode[1],par_ode[1],par_ode[1],par_ode[2],par_ode[2],par_ode[3],par_ode[3],par_ode[3],par_ode[3],par_ode[3],par_ode[3],par_ode[3],par_ode[3]))
  
  betaNVT <- as.vector(c(par_ode[4],par_ode[4],par_ode[4],par_ode[4],par_ode[4],par_ode[4],par_ode[4],par_ode[4],par_ode[4],par_ode[4],par_ode[4],par_ode[4],par_ode[5],par_ode[5],par_ode[6],par_ode[6],par_ode[6],par_ode[6],par_ode[6],par_ode[6],par_ode[6],par_ode[6]))
  
  
  param.list=list(betaVT = as.matrix(sweep(mixing_matrix,1,betaVT,"*")),
                  betaNVT = as.matrix(sweep(mixing_matrix,1,betaNVT,"*")),
                  clearVT=clearVT ,
                  clearNVT=clearNVT ,
                  ageout=ageout,
                  mixing_matrix=as.matrix(mixing_matrix))
  

  
  #res = runsteady(y=state_prop,
                  #fun=SIS_ode_cpp,
                  #parms=param.list, 
                  #times=c(0,1e5))
  
  #steady_state_S =res$y[1:(N_comms*N_agegroups)] %>% matrix(nrow=N_comms, ncol=N_agegroups, byrow=TRUE)#;  #colnames(steady_state_S)=c("age1","age2","age3") ## we need a matrix for each state variable 
  
  #steady_state_VT =res$y[((N_comms*N_agegroups)+1):((N_comms*N_agegroups)*2)] %>% matrix(nrow=N_comms, ncol=N_agegroups, byrow=TRUE)#;  #colnames(steady_state_VT)=c("age1","age2","age3") ## we need a matrix for each state variable 
  
  #steady_state_NVT =res$y[((N_comms*N_agegroups)*2+1):((N_comms*N_agegroups)*3)] %>% matrix(nrow=N_comms, ncol=N_agegroups, byrow=TRUE)#;  colnames(steady_state_NVT)=c("age1","age2","age3") 
  
  #steady_state_B =res$y[((N_comms*N_agegroups)*3+1):((N_comms*N_agegroups)*4)] %>% matrix(nrow=N_comms, ncol=N_agegroups, byrow=TRUE)#;  colnames(steady_state_B)=c("age1","age2","age3") 
  
  
  ## extract the prevalence from the steady state 
  
  time_days  = seq(from=0, to=2500, by=0.2)
  
  #y_start_endem = c( as.vector(t(steady_state_S)), as.vector(t(steady_state_VT)), as.vector(t(steady_state_NVT)), as.vector(t(steady_state_B)))
  
  ###############################
  ## run simualtion
  

  state_prop=c(rep(c(400), N_agegroups*N_comms),rep(1,N_agegroups*N_comms),rep(1,N_agegroups*N_comms),rep(0,N_agegroups*N_comms)) ## put a certian number of people in each state 
 
   #state_prop=c((c(2200,2100,2000,1900,1800,1700,1600,1500,1400,1300,1200,1100,1000,900,800,700,600,500,400,300,200,100,4200,4100,4000,3900,3800,3700,3600,3500,3400,3300,3200,3100,3000,2900,2800,2700,2600,2500,2400,2300,2200,2100,5200,5100,5000,4900,4800,4700,4600,4500,4400,4300,4200,4100,4000,3900,3800,3700,3600,3500,3400,3300,3200,3100)),
   #            rep(c(22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1), N_comms),
   #           rep(c(22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1), N_comms),rep(0,N_agegroups*N_comms)) ## put a certian number of people in each state 
  
  

  run_mod_test = as.data.frame(lsoda(y=state_prop, times=time_days, func=SIS_ode_cpp, parms=param.list))

  
  run_mod_test =  run_mod_test[,-1]
  
 

  test_var_2 <- vector()
  
  for(check in 1:((N_comms*N_agegroups)*4)){
    
    test_var_2[check] <- sd(run_mod_test[12296:12496,check])
    
    
  }

  
  
  if(all( test_var_2 < 1.5) == "TRUE"){
  
    equi = (as.numeric(as.vector(run_mod_test[12496,1:352])))
           
    
 
  }
    
  
  return(equi)
  
}
