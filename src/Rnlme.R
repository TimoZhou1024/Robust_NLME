#' @param nlmeObject
#' @param long.data
#' @param idVar

# independent.raneff: "byModel"; "byOne"; FALSE

Rnlme <- function(nlmeObjects, long.data, idVar, 
                  sd.method="None", dispersion.SD=FALSE, independent.raneff=FALSE,
                  sdghsize=4, itertol=1e-3, Ptol=1e-2, iterMax=15, Verbose=FALSE){

  #set.seed(123)
  ##################################### settings for nlme model 
  JReturn <- get_Jloglike(nlmeObjects)
  
  Jloglike <- JReturn$Jloglike # log likelihood conditional on random effects
  
  Jfixed  <-  JReturn$Jfixed   # fixed parameters
  Jraneff <- JReturn$Jraneff # random effects 
  Jdisp   <-  JReturn$Jdisp # dispersion parameters
   
  
  p <- length(Jfixed)  #  dimension of fixed parameters 
  q <- length(Jraneff) # dimension of random effects
  qSIGMA <- JReturn$SIGMA.dim # dimension of SIGMA
  
  SIGMA.block <- JReturn$SIGMA.block
  
  ##################################### initial values
  fixedest0 <- JReturn$str.fixed
  dispest0 <- JReturn$str.disp
  
  invSIGMA0 <- SIGMA0 <- diag(1,qSIGMA,qSIGMA)
  Lval0 <- NULL # re-parameters for COV matrix
  
  #################################### bounds
  lower.fixed <- JReturn$lower.fixed
  upper.fixed <- JReturn$upper.fixed
  lower.disp <- JReturn$lower.disp
  upper.disp <- JReturn$upper.disp
  
  
  #################################### settings for dataset
  group <- long.data[ , idVar]  # grouping variable, e.g patient ID
  uniqueID <- unique(group)   
  n <- length(uniqueID)  # sample size
  ni <- table(group)   # number of repeated measurements for each subject 
  N <- nrow(long.data)
  
  
  #################################### condition for iteration
  likDiff <- Diff <- Diff0 <- 1
  convergence <- 1
  M <- 1

  while(!(likDiff <= itertol | (Diff <= Ptol & Diff0 <= Ptol) |  M > iterMax)) {

    #################################### estimation
    Diff0 <- Diff
    cat("############## Iteration:", M, "###############","\n")
    
    # estimate random effects
    cat("Start estimating random effects ...\n")
    ran.time <- system.time({
      ran.output <- est_raneff(RespLog=Jloglike, long.data, idVar, Jraneff, 
                               fixedest0, dispest0, invSIGMA0,
                               uniqueID, n,ni,q, N, q_split,df.sigma, df.randisp,
                               Verbose=Verbose, scale=TRUE)
    })
    Bi <- ran.output$Bi
    B <- ran.output$B
    cat("done. elapsed=", round(ran.time[["elapsed"]], 2), "s\n")
    
    # estimate fixed parameters
    cat("Start estimating fixed parameters ... \n")
    fixed.time <- system.time({
      fixed.output <- est_fixed(RespLog=Jloglike, long.data,Jfixed,
                                fixedest0, dispest0, invSIGMA0,
                                Bi, B, 
                                lower=lower.fixed, upper=upper.fixed,
                                Verbose=Verbose)
    })
    
    fixedest <- fixed.output$beta
    cat("done. elapsed=", round(fixed.time[["elapsed"]], 2), "s\n")
    
    # estimate dispersion parameters
    cat("Start estimating dispersion parameters ... \n")
    disp.time <- system.time({
      disp.output <- est_disp_ml(RespLog=Jloglike, long.data, Jdisp, Jfixed, Jraneff,
                                    fixedest, dispest0, invSIGMA0, Lval0,
                                    Bi, B,
                                    lower=lower.disp, upper=upper.disp,
                                    independent=independent.raneff,block=SIGMA.block, Verbose=Verbose)
    })
    
    dispest <- disp.output$disp
    invSIGMA <- disp.output$invSIGMA
    SIGMA <- disp.output$SIGMA
    Lval <- disp.output$Lval
    Lmat <- disp.output$Lmat
    
    cat("done. elapsed=", round(disp.time[["elapsed"]], 2), "s\n")
    
    ####################################################    
    ################## update results ##################
    ####################################################
    # calculate approximated log marginal likelihood value
    loglike_value <- get_loglike_value(RespLog=Jloglike, long.data,
                                       fixedest, dispest, invSIGMA, Bi, B)
    
    
    if(M == 1){
      likDiff <- 1
    } else{
      likDiff <- abs(loglike_value-loglike_value0)/abs(loglike_value0)
    }
    
    # calculate relative changes in mean parameters
    
    Diff <- mean(c(abs((fixedest - fixedest0)/(fixedest0 + 1e-6))))
    
    
    ############## print
    
    cat("fixed.par:", round(fixedest, 2), "\n")
    cat("FixedParDiff = ", Diff, '\n')
    cat("likDiff = ", likDiff, '\n')
    if(!is.null(unlist(dispest))){
      cat("dispersion.par:", round(unlist(dispest), 2), "\n")
    }
    cat("loglike:", loglike_value, "\n")
    cat("SIGMA:", as.matrix(SIGMA), "\n")
    cat("##########################################","\n")
    
    fixedest0 <- fixedest
    dispest0 <- dispest
    invSIGMA0 <- invSIGMA
    SIGMA0 <- SIGMA
    Lval0 <- Lval
    loglike_value0 <- loglike_value
    
    M <- M+1  
  }
  
  ## messages about convergence success or failure
  if((likDiff > itertol  & Diff>Ptol ) ){
    message("Iteration limit reached without covergence.")
    convergence <- 1
  }
  if(likDiff <= itertol & likDiff >= 0){
    message("Successful convergence. Iteration stops because likDiff <= itertol.")
    convergence <- 0
  }

  if(Diff0 <= Ptol & Diff <= Ptol){
    message("Successful convergence. Iteration stops because FixedParDiff <= Ptol in consective two iterations.")
    convergence <- 0
  }

  

  
  # estimate sd's of parameter estimates  

  if(sd.method=="HL") {
    cat("Start estimating SD for fixed parameters ...\n ...\n")
    
    sd.time <- system.time({
      sd_output <- get_sd(RespLog=Jloglike, long.data,  idVar,
                               fixedest0, dispest0, invSIGMA0, SIGMA0,
                               Bi, B,
                               Jfixed,Jraneff)
    })
    fixedSD <- sd_output
    cat("done. elapsed=", round(sd.time[["elapsed"]], 2), "s\n")
  
  } else if(sd.method=="aGH"){
    cat("Start estimating SD for fixed parameters ...\n ...\n")
    sd_output <- get_sd_aGH(RespLog=Jloglike, long.data, idVar, 
                                        fixedest0, dispest0, invSIGMA0,Bi, B,
                                        Jfixed, Jraneff,  
                                        ghsize=ghsize, Silent=T, epsilon=10^{-6}, 
                                        parallel=TRUE)
    fixedSD <- sd_output
    cat("done.\n")
   
  } else if(sd.method=="Both"){
    cat("Start estimating SD for fixed parameters ...\n ...\n")
    sd_HL <- get_sd(RespLog=Jloglike, long.data,  idVar,
                        fixedest0, dispest0, invSIGMA0, SIGMA0,
                        Bi, B,
                        Jfixed,Jraneff)
    sd_aGH <-  get_sd_aGH(RespLog=Jloglike, long.data, idVar, 
                          fixedest0, dispest0, invSIGMA0,Bi, B,
                          Jfixed, Jraneff,  
                          ghsize=ghsize, Silent=T, epsilon=10^{-6}, 
                          parallel=TRUE)
    fixedSD=list(HL=sd_HL, aGH=sd_aGH)
    cat("done.\n")
  } else if(sd.method=="None") {
    fixedSD <- NULL
  }
  

  
  if(dispersion.SD==TRUE){
    cat("Start estimating SD for dispersion parameters ...\n ...\n")
    sd.disp.time <- system.time({
      sd_disp <- get_sd_dipsersion(RespLog=Jloglike, long.data, idVar,
                                   fixedest0, dispest0, invSIGMA0,SIGMA0, Lval0, Lmat,
                                   Bi, B,
                                   Jfixed, Jraneff, Jdisp)
    })
    cat("done. elapsed=", round(sd.disp.time[["elapsed"]], 2), "s\n")
  } else sd_disp <- NULL
  #### AIC
  AIC <- 2*length(c(fixedest0, dispest0, Lval0)) -2*loglike_value0
  BIC <- length(c(fixedest0, dispest0, Lval0))*log(nrow(B))-2*loglike_value0
  
  list(fixedest = fixedest0,
       fixedSD  = fixedSD,
       dispersion = dispest0,
       dispSD=sd_disp,
       Bi = Bi, 
       B = B,
       SIGMA = solve(invSIGMA0), 
       convergence = convergence==0,
       loglike_value = loglike_value0,
       AIC=AIC,
       BIC=BIC,
       long.data = long.data,
       #surv.data = surv.data,
       RespLog = Jloglike,Jfixed=Jfixed, Jraneff = Jraneff,
       idVar = idVar, uniqueID = uniqueID
  )
}
