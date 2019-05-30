#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List SIS_ode_cpp(NumericVector times, vec state, List parms) {
  
  double N_comms = 4;
  double N_agegroups = 22;
  /*double times = 10;*/
  double comp = 0.1; 
  mat betaVT = as<mat>(parms["betaVT"]);
  mat betaNVT = as<mat>(parms["betaNVT"]);
  mat clearVT = as<mat>(parms["clearVT"]);
  mat clearNVT = as<mat>(parms["clearNVT"]);
  mat ageout = as<mat>(parms["ageout"]);
  mat mixing_matrix = as<mat>(parms["mixing_matrix"]);
  

  /* we firstly need to define empty matrices for each of the state variable 
   * start by filling everything with zeros*/

  mat Sus(N_comms, N_agegroups, fill::zeros);
  mat VT (N_comms, N_agegroups, fill::zeros);
  mat NVT(N_comms, N_agegroups, fill::zeros);
  mat B  (N_comms, N_agegroups, fill::zeros);
  mat N  (N_comms, N_agegroups, fill::zeros);
  mat VT_prev_by_comm_byage  (N_comms, N_agegroups, fill::zeros);
  mat NVT_prev_by_comm_byage (N_comms, N_agegroups, fill::zeros);
  
  
  /* objects for the mirror of all state variable */
   
  mat Sustmp(N_comms, N_agegroups, fill::zeros); 
  mat VTtmp(N_comms, N_agegroups, fill::zeros);
  mat NVTtmp(N_comms, N_agegroups, fill::zeros);
  mat Btmp(N_comms, N_agegroups, fill::zeros);
 
 mat Susnew(N_comms, N_agegroups, fill::zeros); 
 mat VTnew(N_comms, N_agegroups, fill::zeros);
 mat NVTnew(N_comms, N_agegroups, fill::zeros);
 mat Bnew(N_comms, N_agegroups, fill::zeros);
 
  /* A vector to store all of the simulation output for 1 time step */

  vec All_vec(352);
  
  
  /* having defined the empty matrices above we can go about filling  them
   * with the values that are in the vector vec state which are the initial
   * conditions that we pass from R. The argument state.subvec allows you to
   * subset elements from a longer vector. However this dosen't work for the
   * matrix so we need to loop through and insert all the individual elements 
   */
    
    for(int j=0; j<N_comms; j++){
    for(int i=0; i<N_agegroups; i++){
        Sus(j,i)  = state(i);
        VT (j,i)  = state((N_agegroups*4)+i); 
        NVT (j,i) = state((N_agegroups*8)+i);
        B (j,i)   = state((N_agegroups*12)+i);
    }
  }
   
  

  mat FOIVT(N_comms, N_agegroups, fill::zeros); 
  mat FOINVT(N_comms, N_agegroups, fill::zeros); 
  
  
  
  /* Need to define mirror objectives so we modify the state variables with all the values that they were read in as 
   * at a single time point */
   
  
  Sustmp = Sus;
  VTtmp = VT;
  NVTtmp = NVT;
  Btmp = B;
  
  N = Sustmp + VTtmp + NVTtmp + Btmp;
  
  double test_sum = accu(VTtmp + NVTtmp + Btmp  );
  
  if(test_sum > 2){
  
  
  VT_prev_by_comm_byage = ((VT + B)/N) ;
    
  NVT_prev_by_comm_byage = ((NVT + B)/N) ;
  
  mat VT_calc_c1 (22,1,fill::ones);
    
  mat NVT_calc_c1 (22,1,fill::ones);
  
  mat VT_calc_c2 (22,1,fill::ones);
  
  mat NVT_calc_c2 (22,1,fill::ones);
  
  mat VT_calc_c3 (22,1,fill::ones);
  
  mat NVT_calc_c3 (22,1,fill::ones);
  
  mat VT_calc_c4 (22,1,fill::ones);
  
  mat NVT_calc_c4 (22,1,fill::ones);
  
  mat result_c1(22,1, fill::zeros);
    
  mat result2_c1 (22,1, fill::zeros);
  
  mat result_c2(22,1, fill::zeros);
  
  mat result2_c2 (22,1, fill::zeros);
  
  mat result_c3(22,1, fill::zeros);
  
  mat result2_c3 (22,1, fill::zeros);
  
  mat result_c4(22,1, fill::zeros);
  
  mat result2_c4 (22,1, fill::zeros);
  
  // Calculate the FOI for each age group and commune

  VT_calc_c1 = VT_prev_by_comm_byage.submat(0,0,0,(N_agegroups-1));
  
  VT_calc_c1 = trans(VT_calc_c1);
  
  NVT_calc_c1 = NVT_prev_by_comm_byage.submat(0,0,0,(N_agegroups-1));   
  
  NVT_calc_c1 = trans(NVT_calc_c1);
  
  result_c1 =    betaVT * VT_calc_c1   ;

  result2_c1 =    betaNVT * NVT_calc_c1   ;
  
  FOIVT.submat(0,0,0,(N_agegroups-1)) = trans(result_c1);
  
  FOINVT.submat(0,0,0,(N_agegroups-1)) = trans(result2_c1);
  
  // Calculate the FOI for each age group and commune 2
  
  VT_calc_c2 = VT_prev_by_comm_byage.submat(1,0,1,(N_agegroups-1));
  
  VT_calc_c2 = trans(VT_calc_c2);
  
  NVT_calc_c2 = NVT_prev_by_comm_byage.submat(1,0,1,(N_agegroups-1));   
  
  NVT_calc_c2 = trans(NVT_calc_c2);
  
  result_c2 =    betaVT * VT_calc_c2   ;
  
  result2_c2 =    betaNVT * NVT_calc_c2   ;
  
  FOIVT.submat(1,0,1,(N_agegroups-1)) = trans(result_c2);
  
  FOINVT.submat(1,0,1,(N_agegroups-1)) = trans(result2_c2);
  
  // Calculate the FOI for each age group and commune 3
  
  VT_calc_c3 = VT_prev_by_comm_byage.submat(2,0,2,(N_agegroups-1));
  
  VT_calc_c3 = trans(VT_calc_c3);
  
  NVT_calc_c3 = NVT_prev_by_comm_byage.submat(2,0,2,(N_agegroups-1));   
  
  NVT_calc_c3 = trans(NVT_calc_c3);
  
  result_c3 =    betaVT * VT_calc_c3   ;
  
  result2_c3 =    betaNVT * NVT_calc_c3   ;
  
  FOIVT.submat(2,0,2,(N_agegroups-1)) = trans(result_c3);
  
  FOINVT.submat(2,0,2,(N_agegroups-1)) = trans(result2_c3);
  
  // Calculate the FOI for each age group and commune 4
  
  VT_calc_c4 = VT_prev_by_comm_byage.submat(3,0,3,(N_agegroups-1));
  
  VT_calc_c4 = trans(VT_calc_c4);
  
  NVT_calc_c4 = NVT_prev_by_comm_byage.submat(3,0,3,(N_agegroups-1));   
  
  NVT_calc_c4 = trans(NVT_calc_c4);
  
  result_c4 =    betaVT * VT_calc_c4   ;
  
  result2_c4 =    betaNVT * NVT_calc_c4   ;
  
  FOIVT.submat(3,0,3,(N_agegroups-1)) = trans(result_c4);
  
  FOINVT.submat(3,0,3,(N_agegroups-1)) = trans(result2_c4);
    
  // determine the number of people leaving the compartment at each time point
  
  mat Nexit (1,4, fill::zeros);
    
  Nexit= Sustmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1)) + VTtmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1)) + NVTtmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1)) + Btmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1)); 
    
 
  Susnew.submat(0,0,(N_comms-1),0)	                             =  - (FOIVT.submat(0,0,(N_comms-1),0)                             + FOINVT.submat(0,0,(N_comms-1),0))%Sustmp.submat(0,0,(N_comms-1),0)                                                         + clearVT.submat(0,0,(N_comms-1),0)%VTtmp.submat(0,0,(N_comms-1),0)                                                          + clearNVT.submat(0,0,(N_comms-1),0)%NVTtmp.submat(0,0,(N_comms-1),0)                                                                                                                                                                                                  - ageout.submat(0,0,(N_comms-1),0)%Sustmp.submat(0,0,(N_comms-1),0)                                                           +  ageout.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))%Nexit ; 
    
  Susnew.submat(0,1,(N_comms-1),(N_agegroups-2))                =  - (FOIVT.submat(0,1,(N_comms-1),(N_agegroups-2))               + FOINVT.submat(0,1,(N_comms-1),(N_agegroups-2)))%Sustmp.submat(0,1,(N_comms-1),(N_agegroups-2))                             + clearVT.submat(0,1,(N_comms-1),(N_agegroups-2))%VTtmp.submat(0,1,(N_comms-1),(N_agegroups-2))                              + clearNVT.submat(0,1,(N_comms-1),(N_agegroups-2))%NVTtmp.submat(0,1,(N_comms-1),(N_agegroups-2))                                                                                                                                                                        - ageout.submat(0,1,(N_comms-1),(N_agegroups-2))%Sustmp.submat(0,1,(N_comms-1),(N_agegroups-2))                             + ageout.submat(0,0,(N_comms-1),(N_agegroups-3))%Sustmp.submat(0,0,(N_comms-1),(N_agegroups-3))  ;
  
  Susnew.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))	 =  - (FOIVT.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1)) + FOINVT.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1)))%Sustmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1)) + clearVT.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))%VTtmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))  + clearNVT.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))%NVTtmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))                                                                                                                                           - ageout.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))%Sustmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))    + ageout.submat(0,(N_agegroups-2),(N_comms-1),(N_agegroups-2))%Sustmp.submat(0,(N_agegroups-2),(N_comms-1),(N_agegroups-2))   ;     
  
  
  VTnew.submat(0,0,(N_comms-1),0) 		  	                       =  +	 FOIVT.submat(0,0,(N_comms-1),0)%Sustmp.submat(0,0,(N_comms-1),0) 	                                                                                                                      -  clearVT.submat(0,0,(N_comms-1),0)%VTtmp.submat(0,0,(N_comms-1),0)                                                            -	comp*FOINVT.submat(0,0,(N_comms-1),0)%VTtmp.submat(0,0,(N_comms-1),0)                                                         +  clearNVT.submat(0,0,(N_comms-1),0)%Btmp.submat(0,0,(N_comms-1),0)                                                                  - ageout.submat(0,0,(N_comms-1),0)%VTtmp.submat(0,0,(N_comms-1),0)	 ;        		  
    
  VTnew.submat(0,1,(N_comms-1),(N_agegroups-2))                 =  +	 FOIVT.submat(0,1,(N_comms-1),(N_agegroups-2))%Sustmp.submat(0,1,(N_comms-1),(N_agegroups-2))                                                                                             -  clearVT.submat(0,1,(N_comms-1),(N_agegroups-2))%VTtmp.submat(0,1,(N_comms-1),(N_agegroups-2))                                -  comp*FOINVT.submat(0,1,(N_comms-1),(N_agegroups-2))%VTtmp.submat(0,1,(N_comms-1),(N_agegroups-2))                              + clearNVT.submat(0,1,(N_comms-1),(N_agegroups-2))%Btmp.submat(0,1,(N_comms-1),(N_agegroups-2))                                   - ageout.submat(0,1,(N_comms-1),(N_agegroups-2))%VTtmp.submat(0,1,(N_comms-1),(N_agegroups-2))	                             + ageout.submat(0,0,(N_comms-1),(N_agegroups-3))%VTtmp.submat(0,0,(N_comms-1),(N_agegroups-3)); 
      
  VTnew.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))	 =  +  FOIVT.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))%Sustmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))                                                                 - clearVT.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))%VTtmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))     - comp*FOINVT.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))%VTtmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))   + clearNVT.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))%Btmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))           - ageout.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))%VTtmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))	 + ageout.submat(0,(N_agegroups-2),(N_comms-1),(N_agegroups-2))%VTtmp.submat(0,(N_agegroups-2),(N_comms-1),(N_agegroups-2)) ;
  
  
  NVTnew.submat(0,0,(N_comms-1),0)  		  	                       =  +	FOINVT.submat(0,0,(N_comms-1),0)%Sustmp.submat(0,0,(N_comms-1),0)	                                                                                                                      - clearNVT.submat(0,0,(N_comms-1),0)%NVTtmp.submat(0,0,(N_comms-1),0)                                                           -	comp*FOIVT.submat(0,0,(N_comms-1),0)%NVTtmp.submat(0,0,(N_comms-1),0)                                                          + clearVT.submat(0,0,(N_comms-1),0)%Btmp.submat(0,0,(N_comms-1),0)                                                                  - ageout.submat(0,0,(N_comms-1),0)%NVTtmp.submat(0,0,(N_comms-1),0);             	 	  
    
  NVTnew.submat(0,1,(N_comms-1),(N_agegroups-2))                  =  + FOINVT.submat(0,1,(N_comms-1),(N_agegroups-2))%Sustmp.submat(0,1,(N_comms-1),(N_agegroups-2))                                                                                           - clearNVT.submat(0,1,(N_comms-1),(N_agegroups-2))%NVTtmp.submat(0,1,(N_comms-1),(N_agegroups-2))                               - comp*FOIVT.submat(0,1,(N_comms-1),(N_agegroups-2))%NVTtmp.submat(0,1,(N_comms-1),(N_agegroups-2))                               + clearVT.submat(0,1,(N_comms-1),(N_agegroups-2))%Btmp.submat(0,1,(N_comms-1),(N_agegroups-2))                                          - ageout.submat(0,1,(N_comms-1),(N_agegroups-2))%NVTtmp.submat(0,1,(N_comms-1),(N_agegroups-2))                             + ageout.submat(0,0,(N_comms-1),(N_agegroups-3))%NVTtmp.submat(0,0,(N_comms-1),(N_agegroups-3))     ;           
      
  NVTnew.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))	   =  + FOINVT.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))%Sustmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))                                                               - clearNVT.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))%NVTtmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))   - comp*FOIVT.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))%NVTtmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))   + clearVT.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))%Btmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))             - ageout.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))%NVTtmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))   + ageout.submat(0,(N_agegroups-2),(N_comms-1),(N_agegroups-2))%NVTtmp.submat(0,(N_agegroups-2),(N_comms-1),(N_agegroups-2))	   ;               	
  
 
  Bnew.submat(0,0,(N_comms-1),0)    		  	                       =  + comp*FOINVT.submat(0,0,(N_comms-1),0)%VTtmp.submat(0,0,(N_comms-1),0)                                                                                                                   +	comp*FOIVT.submat(0,0,(N_comms-1),0)%NVTtmp.submat(0,0,(N_comms-1),0)                                                      - (clearNVT.submat(0,0,(N_comms-1),0)+clearVT.submat(0,0,(N_comms-1),0))%Btmp.submat(0,0,(N_comms-1),0)                                                                                                                                                                   - ageout.submat(0,0,(N_comms-1),0)%Btmp.submat(0,0,(N_comms-1),0);	    		  
    
  Bnew.submat(0,1,(N_comms-1),(N_agegroups-2))                    =  + comp*FOINVT.submat(0,1,(N_comms-1),(N_agegroups-2))%VTtmp.submat(0,1,(N_comms-1),(N_agegroups-2))                                                                                       + comp*FOIVT.submat(0,1,(N_comms-1),(N_agegroups-2))%NVTtmp.submat(0,1,(N_comms-1),(N_agegroups-2))                          - clearNVT.submat(0,1,(N_comms-1),(N_agegroups-2))%B.submat(0,1,(N_comms-1),(N_agegroups-2)) - clearVT.submat(0,1,(N_comms-1),(N_agegroups-2))%Btmp.submat(0,1,(N_comms-1),(N_agegroups-2))                                                                                - ageout.submat(0,1,(N_comms-1),(N_agegroups-2))%Btmp.submat(0,1,(N_comms-1),(N_agegroups-2)) 	                               + ageout.submat(0,0,(N_comms-1),(N_agegroups-3))%Btmp.submat(0,0,(N_comms-1),(N_agegroups-3));
                                                                         
  Bnew.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))	     =  + comp*FOINVT.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))%VTtmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))                                                           + comp*FOIVT.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))%NVTtmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1)) - (clearNVT.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1)) +   clearVT.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1)))%Btmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))                     	                                                  - ageout.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))%Btmp.submat(0,(N_agegroups-1),(N_comms-1),(N_agegroups-1))	   + ageout.submat(0,(N_agegroups-2),(N_comms-1),(N_agegroups-2))%Btmp.submat(0,(N_agegroups-2),(N_comms-1),(N_agegroups-2));
                                                                           

  
  All_vec.subvec(0,21) = trans((Susnew.submat(0,0,0,(N_agegroups-1))));
  All_vec.subvec(22,43) = trans((Susnew.submat(1,0,1,(N_agegroups-1))));
  All_vec.subvec(44,65) = trans((Susnew.submat(2,0,2,(N_agegroups-1))));
  All_vec.subvec(66,87) = trans((Susnew.submat(3,0,3,(N_agegroups-1))));
  
  All_vec.subvec(88,109) = trans((VTnew.submat(0,0,0,(N_agegroups-1))));
  All_vec.subvec(110,131) = trans((VTnew.submat(1,0,1,(N_agegroups-1))));
  All_vec.subvec(132,153) = trans((VTnew.submat(2,0,2,(N_agegroups-1))));
  All_vec.subvec(154,175) = trans((VTnew.submat(3,0,3,(N_agegroups-1))));
  
  All_vec.subvec(176,197) = trans((NVTnew.submat(0,0,0,(N_agegroups-1))));
  All_vec.subvec(198,219) = trans((NVTnew.submat(1,0,1,(N_agegroups-1))));
  All_vec.subvec(220,241) = trans((NVTnew.submat(2,0,2,(N_agegroups-1))));
  All_vec.subvec(242,263) = trans((NVTnew.submat(3,0,3,(N_agegroups-1))));
  
  All_vec.subvec(264,285) = trans((Bnew.submat(0,0,0,(N_agegroups-1))));
  All_vec.subvec(286,307) = trans((Bnew.submat(1,0,1,(N_agegroups-1))));
  All_vec.subvec(308,329) = trans((Bnew.submat(2,0,2,(N_agegroups-1))));
  All_vec.subvec(330,351) = trans((Bnew.submat(3,0,3,(N_agegroups-1))));
  
  } else {
    
    
    /* mat Sustmp(N_comms, N_agegroups, fill::zeros); 
    mat VTtmp(N_comms, N_agegroups, fill::zeros);
    mat NVTtmp(N_comms, N_agegroups, fill::zeros);
    mat Btmp(N_comms, N_agegroups, fill::zeros);*/
    
    All_vec.subvec(0,21) = vec(22, fill::zeros);
    All_vec.subvec(22,43) = vec(22, fill::zeros);
    All_vec.subvec(44,65) = vec(22, fill::zeros);
    All_vec.subvec(66,87) = vec(22, fill::zeros);
    
    All_vec.subvec(88,109) = vec(22, fill::zeros);
    All_vec.subvec(110,131) = vec(22, fill::zeros);
    All_vec.subvec(132,153) = vec(22, fill::zeros);
    All_vec.subvec(154,175) = vec(22, fill::zeros);
    
    All_vec.subvec(176,197) = vec(22, fill::zeros);
    All_vec.subvec(198,219) = vec(22, fill::zeros);
    All_vec.subvec(220,241) = vec(22, fill::zeros);
    All_vec.subvec(242,263) = vec(22, fill::zeros);
    
    All_vec.subvec(264,285) = vec(22, fill::zeros);
    All_vec.subvec(286,307) = vec(22, fill::zeros);
    All_vec.subvec(308,329) = vec(22, fill::zeros);
    All_vec.subvec(330,351) = vec(22, fill::zeros);
    
  }
  
  return List::create((All_vec));
   

  
 // return List::create(Sus, VT, NVT, B,FOIVT,FOINVT);

  }


