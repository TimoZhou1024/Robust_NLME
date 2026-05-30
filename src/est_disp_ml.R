
est_disp_ml <- function(RespLog, long.data, Jdisp,Jfixed, Jraneff,
                        fixedest, dispest0, invSIGMA0, Lval0,
                        Bi, B,
                        lower, upper,
                        independent, block, Verbose=TRUE){
  
  n <- nrow(Bi)
  q1 <- length(Jdisp)
  q2 <- ncol(invSIGMA0)
  if(q2>1 & independent==FALSE) {
    Lmat  <- make_strMat(q2) 
  } 
  if(q2>1 & independent=="byModel" & max(block)>1){
    Lmat  <- make_strMat(q2, block)
  }
  if(q2==1 | independent=="byOne"|max(block)==1) {
    M=diag(rep(1,q2))
    M.expr <- paste0("matrix(c(", paste0(c(M), collapse=","), "),nrow=", q2, ",ncol=", q2, ", byrow=TRUE)") 
    Lmat <- list(M=M.expr, Mpar=NULL)
  }
  qL <- length(Lmat$Mpar)
  ran.loglike <- str_replace_all(RespLog$ran.loglike, "invSIGMA", paste0("solve(",Lmat$M, ")"))

  Yrandisp <- !is.null(RespLog$randisp.loglike)
  Ysigma <- !is.null(RespLog$sigma.loglike)
  
  Jdisp_new <- c(Jdisp, Lmat$Mpar)
  
  # ff() returns the negative value of h-likelihood to be optimized.
  ff <- function(xx){  # xx are input values of dispersion parameters
    fy <- numeric(1)
    
    # assign values to the parameters
    par.val <- make_name(Jdisp_new , xx)  
    par.val <- c(fixedest, par.val)
    
    
    # evaluate the h-likelihood
    
    mu.val <- unlist(map(RespLog$mu.loglike, function(t){
      with(long.data, with(par.val, with(B, eval(parse(text=t)))))
    }))
    
    if(Ysigma){
      sigma.val <- with(par.val, with(Bi, eval(parse(text=RespLog$sigma.loglike))))
    } else sigma.val <- 0
    
    if(Yrandisp) {
      randisp.val <- with(par.val, with(Bi, eval(parse(text=RespLog$randisp.loglike))))
    } else {
      randisp.val <- 0
    }
    ran.val <- vector("list", n)
    for(i in 1:n){
      ran.val[[i]] <-  with(par.val, with(Bi[i,], eval(parse(text=ran.loglike))))
    }
    ran.val <- unlist(ran.val) 
    
    # evaluate the adjusted profile h-likelihood
    fy <- sum(mu.val, na.rm=TRUE)+sum(sigma.val)+sum(randisp.val)+sum(ran.val)
    
    return(-fy)
  }
  
  k <- length(Jdisp_new)

  gr.mu <- map(RespLog$mu.loglike, function(t){
     get_Deriv_cached(t, Jdisp)
  })
  gr.ran <- get_Deriv_cached(ran.loglike, Jdisp_new)
  
  ############ gradient function
  
  gr <- function(xx){
    fy <- numeric(k)
    
    # assign values to the parameters
    par.val <- make_name(Jdisp_new , xx)  
    par.val <- c(fixedest, par.val)
    
    
    # derivative of sd parameters 
    gr.mu.val <- map(gr.mu, function(t){
      val <- c()
      for(i in 1:nrow(B)){
      val <- rbind(val, with(par.val, with(long.data[i,], with(B[i,], eval(parse(text=t))))))
      }
      val <- val[complete.cases(val),]
      c(apply(val, 2, sum, na.rm=TRUE), rep(0, qL))
    })

    
    gr.mu.val <- apply(do.call(rbind, gr.mu.val), 2, sum)
    
    if(!is.null(Lmat$Mpar)){
      gr.ran.val <- c()
      for(i in 1:n){
        temp.val <-  with(par.val, with(Bi[i,], eval(parse(text=gr.ran))))
        gr.ran.val <- rbind(gr.ran.val, temp.val)
      }
      gr.ran.val <- as.vector(apply(gr.ran.val, 2, sum))
      
      fy <- gr.mu.val+gr.ran.val
    } else {
      fy <- gr.mu.val
    }
    return(-fy)
  }
  
  str_val00 <- dispest0
  convge <- -1
  M <- 0
  
  if(Verbose==TRUE) check=1  else check=0
  
  # start iteration
  while(convge != 0 & M<20){
    if(is.null(Lval0)){
      Lval00 <- runif( qL, 0, pi)
    } else{
      Lval00 <- Lval0+rnorm(qL,0,0.01)
    }
    str_val0 <-  unlist(c(str_val00+rnorm(q1,0,0.1), Lval00))
    
    
    result <- try(optim(par=str_val0, fn=ff, gr=gr,
                        method="L-BFGS-B",
                        lower=c(lower, rep(0.01,qL)), 
                        upper=c(upper,rep(pi*0.95,qL)),
                        control = list(trace=check,maxit=1000)), 
                  silent=T)
    
    error_mess <- attr(result, "class")
    
    if(is.null(error_mess)){ 
      convge <- result$convergence 
      str_val00 <- result$par[1:q1]
      Lval0 <- result$par[-c(1:q1)]
    } else {
      str_val00 <- dispest0
      convge = -1
      
    }
    M <- M+1
    if(Verbose==TRUE){ 
      cat(paste0("\n M=",M,", Convergence=", convge==0,".\n"))
      print(result)}
  }
  
  if(convge==0){
    if(q1>0){
      output <- result$par[1:q1]
      names(output) <- Jdisp
      Lval <- make_name(Lmat$Mpar, result$par[-(1:q1)])
      SIGMAest <- with(Lval, eval(parse(text=Lmat$M)))
      invSIGMAest <- solve(SIGMAest)
    } else{ 
      output=NULL
      Lval <- make_name(Lmat$Mpar, result$par)
      SIGMAest <- with(Lval, eval(parse(text=Lmat$M)))
      invSIGMAest <- solve(SIGMAest)
    }
    
  } else {
    stop("Iteration stops because dispersion parameters can not be successfully estimated.")  
  }
  
  return(list(disp=output, invSIGMA=invSIGMAest, Lval=Lval, SIGMA=SIGMAest, Lmat=Lmat))
}
