# separate mu.loglike, so we can incorporate missing data 
est_individual_raneff1 <- function(RespLog, data, raneff, 
                                  fixedest, dispest, invSIGMA, 
                                  Verbose=TRUE) {
  
  Yrandisp <- !is.null(RespLog$randisp.loglike)
  Ysigma <- !is.null(RespLog$sigma.loglike)
  
  # ff() returns the value of h-likelihood
  ff <- function(xx){  
    fy <- numeric(1)
    
    # assign values to parameters
    B.val <- make_name(raneff, xx) 
    par.val <- c(B.val, fixedest, dispest) 
    par.val$invSIGMA <- invSIGMA
    
    # values from longitudinal data
    mu.val <- unlist(map(RespLog$mu.loglike, function(t){with(par.val, with(data, eval(parse(text=t))))}))
    
    # values from residual dispersion model
    if(Ysigma) {
      sigma.val <- with(par.val, eval(parse(text=RespLog$sigma.loglike)))
    } else sigma.val <- 0
    
    # values from random eff dispersion model
    if(Yrandisp) {
      randisp.val <- with(par.val, eval(parse(text=RespLog$randisp.loglike)))
    } else randisp.val<-0 
    
    # values from the distribution of random effect
    ran.val <- c(with(par.val,  eval(parse(text=RespLog$ran.loglike))))
    
    fy <- sum(mu.val, na.rm=TRUE)+sigma.val+randisp.val+ ran.val
    # print(fy)
    return(-fy)
  }
  
  k <-  length(raneff)
  
  gr.mu <- map(RespLog$mu.loglike, function(t){get_Deriv_cached(t, raneff)})
  if(Ysigma) gr.sigma <- get_Deriv_cached(RespLog$sigma.loglike, raneff)
  if(Yrandisp) {gr.randisp <- get_Deriv_cached(RespLog$randisp.loglike, raneff)}
  gr.ran <- get_Deriv_cached(RespLog$ran.loglike, raneff)
  
  gr <- function(xx){
    fy <- numeric(k)
    
    # assign values to parameters
    B.val <- make_name(raneff, xx)              
    par.val <- c(B.val, fixedest, dispest) 
    par.val$invSIGMA <- invSIGMA
    
    gr.mu.val <- gr.sigma.val <- gr.randisp.val <- gr.ran.val <- rep(NA, k)
    
    ##
    
    gr.mu.val <- map(gr.mu, function(t){
      val <- c()
      for(i in 1:nrow(data)){
      val <- rbind(val, with(par.val, with(data[i,], eval(parse(text=t)))))
      }
      val <- val[complete.cases(val),]
      as.vector(apply(val, 2, sum, na.rm=TRUE))
    }
    )
    
    gr.mu.val <- apply(do.call(rbind, gr.mu.val), 2, sum)
    
    if(Ysigma) {
      gr.sigma.val <- as.vector(with(par.val, eval(parse(text=gr.sigma))))
    } else gr.sigma.val <- rep(0,k)
    
    if(Yrandisp) {
      gr.randisp.val <- as.vector(with(par.val, eval(parse(text=gr.randisp))))
    } else {
      gr.randisp.val <- rep(0,k)
    }
    gr.ran.val <- as.vector(with(par.val, eval(parse(text=gr.ran))))
    
    fy <-  gr.mu.val + gr.sigma.val + gr.randisp.val + gr.ran.val
    
    
    return(-fy)
  }
  
  # start iteration
  message <- -1
  M <- 0
  bval <- rep(0,k)
  
  while(message!=0 & M<50){
    error_mess="try-error"
    
    bval <- rnorm(k,0,0.01)
    
    result <- try(optim(par=bval,  fn=ff, gr=gr, method="L-BFGS-B", 
                        control=list(maxit = 2000, trace=0)), silent=TRUE)
    
    error_mess <- attr(result, "class")
    
    M <- M+1
    if(length(error_mess)==1) message=-1 else message <- result$convergence
    
    if(Verbose==TRUE) cat(paste0("\n M=",M,", Convergence=", message==0,".\n"))
  }
  
  if(message==0){
    bi <- result$par
  } else {
    stop("Iteration stops because random effects can not be successfully estimated.")  
  }
  
  return(bi)
}
