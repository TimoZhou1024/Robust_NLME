# for exp(a0) ~ df/chisq
get_sd_bootstrap2<- function(Rnlme.fit, simdat,at.rep ,k.runs=50, big1=0.1, big2=0.15, 
                             independent.raneff = "byModel",df){
  group <- simdat$patid  # grouping variable, e.g patient ID
  uniqueID <- unique(group)   
  n <- length(uniqueID)  # sample size
  ni <- table(group)   # number of repeated measurements for each subject 
  N <- nrow(simdat)
  t <- simdat$day
  
  # estimates from cd4.fit
  gamma <- Rnlme.fit$fixedest[c(6:8)]
  xi <- Rnlme.fit$dispersion["xi"]
  sigma_b <- Rnlme.fit$dispersion["sigb1"]
  
  # estimates from Rnlme
  d <- Rnlme.fit$dispersion[c(1:2)]
  Mat <- Rnlme.fit$SIGMA
  beta <- Rnlme.fit$fixedest[c(1:3)]
  alpha0 <- Rnlme.fit$fixedest["alpha0"]
  alpha1 <-  Rnlme.fit$fixedest["alpha1"]
  alpha <- c(alpha0, alpha1)


  
  true.fixed <- c(beta, alpha, gamma)
  n.par <- length(true.fixed)
  ################ models
  nf1 <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
  nf2 <- function(p1,p2,p3,t) p1+p2*t+p3*t^2
  
  sigma1.bt <- list(
    model=NULL,
    str.disp=xi,
    lower.disp=NULL,
    upper.disp=NULL,
    parName="xi"
  )
  
  
  lmeObject.bt <- list(
    nf = "nf2" ,
    model= cd4 ~ nf(p1,p2,p3,day),
    var=c("day"),
    fixed = p1+p2+p3 ~1,
    random = p1 ~1,
    family='normal', 
    ran.dist='normal',
    fixName="gamma",
    ranName="b",
    dispName="sigb",
    sigma=sigma1.bt,    # residual dispersion model (include residual random eff)
    ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
    str.fixed=gamma,  # starting value for fixed effect
    str.disp=c(sigma_b),  # starting value for fixed dispersion of random eff
    lower.fixed=NULL, # lower bounds for fixed eff
    upper.fixed=rep(100,3), # upper bounds for fixed eff
    lower.disp=c(0), # lower bounds for fixed dispersion of random eff
    upper.disp=c(Inf) # upper bounds for  fixed dispersion of random eff
  )
  
  # residual dispersion model:  
  sigma2.bt <- list(
    model=~1+cd4.true+(1|patid),
    link='log',
    ran.dist="inverse-Chi",
    str.fixed=c(alpha0, alpha1),
    lower.fixed=NULL,
    upper.fixed=NULL,
    fixName="alpha",
    ranName="a",
    df=df,
    trueVal.model=list(var="cd4.true", model=lmeObject.bt)
  )
  
  nlmeObject.bt <- list(
    nf = "nf1",
    model= lgcopy ~ nf(p1,p2,p3,day),
    var=c("day"),
    fixed = p1+p2+p3 ~1,
    random = p1+p3 ~1,
    family='normal', 
    ran.dist='normal',
    fixName="beta",
    ranName="u",
    dispName="d",
    sigma=sigma2.bt,    # residual dispersion model (include residual random eff)
    ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
    str.fixed=beta,  # starting value for fixed effect
    str.disp=d,  # starting value for fixed dispersion of random eff
    lower.fixed=NULL, # lower bounds for fixed eff
    upper.fixed=rep(100,3), # upper bounds for fixed eff
    lower.disp=c(0,0), # lower bounds for fixed dispersion of random eff
    upper.disp=c(Inf,Inf) # upper bounds for  fixed dispersion of random eff
  )
  
  
  nlmeObjects.bt <- list(nlmeObject.bt, lmeObject.bt)
  
  ## bootstrap
  est <- disp.est <-  c()
  List.Rnlme <- NULL
  max.attempts <- as.integer(Sys.getenv("SIM_BOOTSTRAP_MAX_ATTEMPTS", "100"))
  
  for(k in 1:k.runs){
    
    cat("This is run", at.rep, "--- BT run", k, "\n")
    
    model.fit  <-  0
    class(model.fit) <- "try-error"
    convg <- FALSE
    attempt <- 0
    
    
    while(convg==FALSE | class(model.fit)=="try-error"){
      attempt <- attempt + 1
      if(attempt > max.attempts){
        stop("Exceeded SIM_BOOTSTRAP_MAX_ATTEMPTS while trying to obtain a converged bootstrap fit.")
      }
      
      ## generate random effects
      a0 <- log(df/rchisq(n, df))
      
      D <- diag(c(d, sigma_b)) %*% Mat %*% diag(c(d, sigma_b))
      ran <- rmvnorm(n, sigma=D)
      
      u <- ran[,c(1:2)]
      b1 <- ran[,3]
      
      u <- cbind(u[,1], 0, u[,2])
      simdat.bt <- c()
      
      for(i in 1:n){
        indexi <- simdat$patid==uniqueID[i] 
        nii <- ni[i]
        ti <- t[indexi]
        
        ## simulate CD4
        b1i <- b1[i]
        cd_errori <- rnorm(nii, sd=xi)
        cd_truei <- nf2(gamma[1]+b1i, gamma[2], gamma[3], ti)
        cd_obsi <- cd_truei+cd_errori
        
        ## get time-varying variance
        a0i <- a0[i]
        sdi <- sqrt(exp(alpha0+alpha1*cd_truei+a0i))
        errori <- rnorm(nii, sd=sdi)
        
        ## simulate lgcopy
        ui <- u[i,]
        tolEffi <- beta+ui
        lgcopyi <- nf1(tolEffi[1], tolEffi[2],tolEffi[3],ti)+errori
        
        
        dati <- data.frame(patid=uniqueID[i], day=ti, lgcopy=lgcopyi, cd4=cd_obsi)
        
        simdat.bt <- rbind(simdat.bt, dati)
      }
      
      simdat.bt <- simdat.bt %>% arrange(patid, day)
      
      fit.time <- system.time({
        model.fit <- try(Rnlme(nlmeObject=nlmeObjects.bt, long.data=simdat.bt, idVar="patid", 
                               independent.raneff = independent.raneff ))
      })
      cat("BT Rnlme elapsed=", round(fit.time[["elapsed"]], 2), "s\n")
      if(class(model.fit)!="try-error") convg <- model.fit$convergence
    }
    
    est <- rbind(est, model.fit$fixedest)
    disp.est <- rbind(disp.est, model.fit$dispersion)
    
  }
  
  
  
  drop.index1 <- apply(est,1,FUN=function(t){max(abs((t-true.fixed)/true.fixed))>big1})
  
  
  drop.index2 <- apply(est,1,FUN=function(t){max(abs((t-true.fixed)/true.fixed))>big2})
  
  if((k.runs-sum(drop.index1))==0) {
    se.bt1=rep(NA, n.par)
  } else {
    se.bt1=apply(matrix(est[!drop.index1,], ncol=n.par), 2, sd)
  }
  
  
  if((k.runs-sum(drop.index2))==0) {
    se.bt2=rep(NA, length(true.fixed))
  } else {
    se.bt2=apply(matrix(est[!drop.index2,], ncol=n.par), 2, sd)
  }
  
  return(list(se.bt=apply(est, 2, sd), se.bt1=se.bt1, 
              se.bt2=se.bt2, 
              runs.bt1=k.runs-sum(drop.index1), 
              runs.bt2=k.runs-sum(drop.index2)))
  
  
}
