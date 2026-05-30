##########################################
# Setting IV:

## number of subject: n=100
## repeated measurements: ni=30
## a_i ~ N(0, sigma_a)
## number of simulation repetition: rep=100
## number of bootstrap runs: k.runs=100

# the model for sigma contains cd4 only
# for the Two-step (TS): the model for sigma contains cd4_pred
# for Joint models (JM): the model for sigma contains cd4* (unobserved true value)
# This simulation compare LB, TS and JM
# only p1 and p3 contains random effects in NLME
########################### load packages
library(nlme)
library(dplyr)
library(Deriv)
library(stringr)
library(LaplacesDemon)
library(mvtnorm)
library(purrr)
library(MASS)
library(Matrix)
library(xtable)
library(R.utils)

rm(list=ls())
########################## source all functions
file.sources <- list.files(path=here::here("src"), pattern="*.R$", full.names=TRUE)
file.sources <- file.sources[basename(file.sources)!="00_dependencies.R"]
sapply(file.sources, source)
######################### Simulation setting
set.seed(as.integer(Sys.getenv("SIM_SEED", "1")))
rep <- as.integer(Sys.getenv("SIM_REP", "100"))
k.runs <- as.integer(Sys.getenv("SIM_K_RUNS", "100")) # number of bootstrap runs
skip.bootstrap <- tolower(Sys.getenv("SIM_SKIP_BOOTSTRAP", "false")) %in% c("1", "true", "yes", "y")
max.attempts <- as.integer(Sys.getenv("SIM_MAX_ATTEMPTS", "100"))
output.dir <- Sys.getenv("SIM_OUTPUT_DIR", here::here())
dir.create(output.dir, recursive=TRUE, showWarnings=FALSE)
n <- 100
ni <- 15
N <- n*ni
ti <- seq(0, 1, length.out=ni)

patid <- rep(1:n, each=ni)
day <- rep(ti, n)
uniqueID <- seq(1:n)   

beta <- c(2.5, 3.0, 7.5)
gamma <- c(5.2, 1.6, -1.2)
d <- c(0.4, 1.0)
Mat <- matrix(c(1,0.5,  0.5, 1), ncol=2)

alpha0 <- -8
alpha1 <-  1.5
alpha <- c(alpha0, alpha1)
sigma_a <- 0.7
sigma_b <- 0.4
xi <- 0.2
true.fixed <- c(beta, alpha, gamma)

nf1 <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
nf2 <- function(p1,p2,p3,t) p1+p2*t+p3*t^2

script <- "
nf1 <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
"
writeLines(script, here::here("src","nf1.R"))

script <- "
nf2 <- function(p1,p2,p3,t) p1+p2*t+p3*t^2
"
writeLines(script, here::here("src","nf2.R"))
########################## model specification

####Two step

# residual dispersion model:  
sigmaObject_TS <- list(
  model=~1+cd4.pred+(1|patid),
  link='log',
  ran.dist="normal",
  str.fixed=c(alpha0, alpha1),
  lower.fixed=NULL,
  upper.fixed=NULL,
  fixName="alpha",
  ranName="a",
  dispName="siga",
  str.disp=sigma_a
)


nlmeObject_TS <- list(
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
  sigma=sigmaObject_TS,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=beta,  # starting value for fixed effect
  str.disp=d,  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

nlmeObjects_TS=list(nlmeObject_TS)

#### JM

sigma1 <- list(
  model=NULL,
  str.disp=xi,
  lower.disp=NULL,
  upper.disp=NULL,
  parName="xi"
)

lmeObject_JM <- list(
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
  sigma=sigma1,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=gamma,  # starting value for fixed effect
  str.disp=c(sigma_b),  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf) # upper bounds for  fixed dispersion of random eff
)

# residual dispersion model:  
sigma2 <- list(
  model=~1+cd4.true+(1|patid),
  link='log',
  ran.dist="normal",
  str.fixed=c(alpha0, alpha1),
  lower.fixed=NULL,
  upper.fixed=NULL,
  fixName="alpha",
  ranName="a",
  dispName="siga",
  str.disp=sigma_a,
  trueVal.model=list(var="cd4.true", model=lmeObject_JM)
)

nlmeObject_JM <- list(
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
  sigma=sigma2,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=beta,  # starting value for fixed effect
  str.disp=d,  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

nlmeObjects_JM <- list(nlmeObject_JM, lmeObject_JM)

################################# Simulation runs

est.NLME <- sd.NLME <- est.TS <- sd.TS <-  est.JM <- sd.JM <- c()
alpha.NLME <- c()
sd.bt.JM <- sd.bt1.JM <- sd.bt2.JM <- c()
runs.bt1 <- runs.bt2 <- c()

for(k in 1:rep){
  cat("This is run", k, "\n")
  
  nlme.fit <- cd4.fit <- TS <- JM <- 0
  class(nlme.fit) <- class(cd4.fit) <- class(TS) <- class(JM)<- "try-error"
  TSconvg <- JMconvg <-  FALSE
  attempt <- 0
  
  while(class(nlme.fit)[1]=="try-error"|class(cd4.fit)=="try-error" | class(JM)=="try-error" | class(TS)=="try-error"
        | TSconvg==FALSE| JMconvg==FALSE){
    attempt <- attempt + 1
    if(attempt > max.attempts){
      stop("Exceeded SIM_MAX_ATTEMPTS while trying to obtain converged NLME, TS, and JM fits.")
    }
    
    ##########################  simulate data set
    cat("--Simulating data\n")
    ## generate random effects
    a0 <- rnorm(n, sd=sigma_a)
    
    
    D <- diag(d) %*% Mat %*% diag(d)
    u <- rmvnorm(n, sigma=D)
    u <- cbind(u[,1], 0, u[,2])
    
    b1 <- rnorm(n, sd=sigma_b)
    
    simdat <- c()
    ## data set
    for(i in 1:n){
      
      ## simulate CD4
      b1i <- b1[i]
      cd_errori <- rnorm(ni, sd=xi)
      cdi_true <- nf2(gamma[1]+b1i, gamma[2], gamma[3], ti)
      cdi_obs <- cdi_true+cd_errori
      
      
      ## get time-varying variance
      a0i <- a0[i]
      sdi <- sqrt(exp(alpha0+alpha1*cdi_true+a0i))
      errori <- rnorm(ni, sd=sdi)
      
      ## simulate lgcopy
      ui <- u[i,]
      tolEffi <- beta+ui
      lgcopyi <- nf1(tolEffi[1], tolEffi[2],tolEffi[3],ti)+errori
      
      
      dati <- data.frame(patid=uniqueID[i], day=ti, lgcopy=lgcopyi, cd4=cdi_obs)
      
      simdat <- rbind(simdat, dati)
    }
    
    simdat <- simdat %>% dplyr::arrange(patid, day)
    cat("--Done \n")
    ########################## Modeling simdat
    #### nlme()
    cat("--Fit NLME() \n")
    dat_g <- groupedData(lgcopy~day|patid, data=simdat)
    nlme.fit <- withTimeout({try(nlme(lgcopy~nf1(p1,p2,p3, day),fixed = p1+p2+p3 ~1,random  = p1+p3 ~1,
                     data=dat_g,start=beta))},timeout=1.5,onTimeout="warning")
    
    cat("--Done \n")
    if(class(nlme.fit)[1]!="try-error"){
      
      #### Two step
      
      cd4.fit <- try(lme(cd4~day+I(day^2), data=simdat, random=~1|patid))
      if(class(cd4.fit)!="try-error"){
        simdat$cd4.pred <- fitted(cd4.fit)
        
        cat("--Runing Two-step model\n")
        TS <- try(Rnlme(nlmeObjects=nlmeObjects_TS, long.data=simdat, 
                        idVar="patid", sd.method="HL", dispersion.SD = TRUE,
                        independent.raneff=FALSE))
        cat("--done\n")
        if(class(TS)!="try-error") TSconvg <- TS$convergence
        
        if(class(TS)!="try-error" & TSconvg){
          #### JM
          cat("--Runing Joint model\n")
          JM <- try(Rnlme(nlmeObjects=nlmeObjects_JM, long.data=simdat, 
                          idVar="patid", sd.method="HL", dispersion.SD = TRUE,
                          independent.raneff="byModel"))
          cat("--done\n")
          if(class(JM)!="try-error") JMconvg <- JM$convergence
        }
      }
    }
  }
  ############### Bootstrapping SE #####################
  if(skip.bootstrap){
    cat("--Skipping Bootstrapping SE for Joint model\n\n")
    JM.SD.BT <- list(se.bt=rep(NA_real_, length(true.fixed)),
                     se.bt1=rep(NA_real_, length(true.fixed)),
                     se.bt2=rep(NA_real_, length(true.fixed)),
                     runs.bt1=NA_integer_,
                     runs.bt2=NA_integer_)
  } else {
    cat("--Runing Bootstrapping SE for Joint model\n\n")
    JM.SD.BT <- get_sd_bootstrap1(Rnlme.fit=JM, simdat, at.rep=k ,k.runs=k.runs, independent.raneff = "byModel")
    cat("--done\n")
  }
  
  runs.bt1 <- c(runs.bt1, JM.SD.BT$runs.bt1)
  runs.bt2 <- c(runs.bt2, JM.SD.BT$runs.bt2)
  ############### store output
  #### nlme
  est.NLME <- rbind(est.NLME, fixef(nlme.fit))
  alpha.NLME <- c(alpha.NLME, 2*log(nlme.fit$sigma))
  sd.NLME <- rbind(sd.NLME, summary(nlme.fit)$tTable[,"Std.Error"])
  
  #### TS
  est.TS <- rbind(est.TS, c(TS$fixedest, fixef(cd4.fit)))
  sd.TS <- rbind(sd.TS, c(TS$fixedSD, summary(cd4.fit)$tTable[,"Std.Error"]))
  
  
  
  #### JM
  est.JM <- rbind(est.JM, JM$fixedest)
  sd.JM <- rbind(sd.JM, JM$fixedSD)
  sd.bt.JM <- rbind(sd.bt.JM, JM.SD.BT$se.bt)
  sd.bt1.JM <- rbind(sd.bt1.JM, JM.SD.BT$se.bt1)
  sd.bt2.JM <- rbind(sd.bt2.JM, JM.SD.BT$se.bt2)
  
}
  ##################### organize output 
  
  colnames(est.TS) <- colnames(sd.TS) <-   colnames(sd.JM) <- colnames(est.JM)
  colnames(sd.bt.JM) <- colnames(sd.bt1.JM) <-colnames(sd.bt2.JM) <- colnames(est.JM)

  
  NLME.out <- list(True=beta,Est=est.NLME, SD=sd.NLME)
  TS.out <- list(True=true.fixed,Est=est.TS, SD=sd.TS)
  JM.out <- list(True=true.fixed,Est=est.JM, SD=sd.JM,
                 SD.BT=sd.bt.JM, SD.BT1=sd.bt1.JM, SD.BT2=sd.bt2.JM)
  

big <- as.numeric(Sys.getenv("SIM_RBIAS_REMOVE_CUTOFF", "15"))
rm.big <- function(out, big){
  out.bias <- t(apply(out$Est, 1, FUN=function(t){abs(t-out$True)/abs(out$True)*100}))
  out.rm <-  apply(out.bias,1,max)>big
  cat("\n",  deparse(substitute(out)) ,"output remove", sum(out.rm), "row with rBias >",big ,"%\n")
  out1 <- out
  out1$Est <- out$Est[!out.rm,,drop=FALSE]
  out1$SD <- out$SD[!out.rm,,drop=FALSE]
  if(!is.null(out$SD.BT)){
    out1$SD.BT <- out$SD.BT[!out.rm,,drop=FALSE]
    out1$SD.BT1 <- out$SD.BT1[!out.rm,,drop=FALSE]
    out1$SD.BT2 <- out$SD.BT2[!out.rm,,drop=FALSE]
  }
  
  return(list(out_df=out1, out.rm=out.rm))
}


get_summary<- function(out){
  Est <- apply(out$Est, 2, mean, na.rm=TRUE)
  
  Bais_mat <- t(apply(out$Est, 1, FUN=function(t){t-out$True}))
  
  rBias <- abs(Est-out$True)/abs(out$True)*100
  
  rMSE <-  apply(Bais_mat, 2, FUN=function(t){mean(t^2, na.rm=TRUE)})/abs(out$True)*100
  
  SE.em <- apply(out$Est, 2, sd, na.rm=TRUE)
  
  SE <- apply(out$SD,2, FUN = function(t){sqrt(mean(t^2, na.rm=TRUE))})
  
  
  lower <- out$Est-1.96*out$SD
  upper <- out$Est+1.96*out$SD
  
  Coverage <- c()
  for(i in 1:length(out$True)){
    cov <- (lower[,i] <= out$True[i]) & (out$True[i]<= upper[,i])
    Coverage <- c(Coverage, mean(cov, na.rm=TRUE))
  }
  
  if(!is.null(out$SD.BT)){
    SE.BT <- apply(out$SD.BT,2, FUN = function(t){sqrt(mean(t^2, na.rm=TRUE))})
    SE.BT1 <- apply(out$SD.BT1,2, FUN = function(t){sqrt(mean(t^2, na.rm=TRUE))})
    SE.BT2 <- apply(out$SD.BT2,2, FUN = function(t){sqrt(mean(t^2, na.rm=TRUE))})
    
    lower.bt <- out$Est-1.96*out$SD.BT
    upper.bt <- out$Est+1.96*out$SD.BT
    
    lower.bt1 <- out$Est-1.96*out$SD.BT1
    upper.bt1 <- out$Est+1.96*out$SD.BT1
    
    lower.bt2 <- out$Est-1.96*out$SD.BT2
    upper.bt2 <- out$Est+1.96*out$SD.BT2
    
    Coverage.bt <- Coverage.bt1 <-  Coverage.bt2 <- c()
    for(i in 1:length(out$True)){
      cov <- (lower.bt[,i] <= out$True[i]) & (out$True[i]<= upper.bt[,i])
      Coverage.bt <- c(Coverage.bt, mean(cov, na.rm=TRUE))
      
      cov <- (lower.bt1[,i] <= out$True[i]) & (out$True[i]<= upper.bt1[,i])
      Coverage.bt1 <- c(Coverage.bt1, mean(cov, na.rm=TRUE))
      
      cov <- (lower.bt2[,i] <= out$True[i]) & (out$True[i]<= upper.bt2[,i])
      Coverage.bt2 <- c(Coverage.bt2, mean(cov, na.rm=TRUE))
    }
    
  }
  
  res <- cbind(Est, rBias, rMSE, SE.em, SE, Coverage)
  
  if(!is.null(out$SD.BT)) res <- cbind(Est, rBias, rMSE, SE.em, SE, Coverage, 
                                       SE.BT, Coverage.bt,
                                       SE.BT1, Coverage.bt1,
                                       SE.BT2, Coverage.bt2)
  
  res
}

cat("\n Average runs for BT1 is", mean(runs.bt1), ".\n")
cat("\n Average runs for BT2 is", mean(runs.bt2), ".\n")
(nl <- get_summary(NLME.out))
mean(alpha.NLME)
(ts <- get_summary(TS.out))
(jm <- get_summary(JM.out))

NLME.out1 <- rm.big(NLME.out, big)$out_df
(nl1 <- get_summary(NLME.out1))
mean(alpha.NLME[!rm.big(NLME.out, big)$out.rm])


TS.out1 <- rm.big(TS.out, big)$out_df
(ts1 <- get_summary(TS.out1))

JM.out1 <- rm.big(JM.out, big)$out_df
(jm1 <- get_summary(JM.out1))

cat("\n xtable for original output \n ")
xtable(cbind(rbind(nl,0,0,0,0,0), ts, jm), type = "latex",digits = 3)


cat("\n xtable for output with large rBias removed \n ")
xtable(cbind(rbind(nl1,0,0,0,0,0),ts1, jm1), type = "latex",digits = 3)

saveRDS(list(NLME.out=NLME.out,TS.out=TS.out,JM.out=JM.out, alpha.NLME=alpha.NLME), 
        file.path(output.dir, "s4.rds"))

save.image(file.path(output.dir, "s4.RData"))



