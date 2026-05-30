##########################################
# Supplementary setting
## for revision:

## 1. Use LB estimates as starting value (except alpha's)
## 2. number of simulation repetition: rep=200
## 3. add a one-step model: a model with non-time dependence of the variance  
## 4. add prediction (Bias, MSE, coverage) **please ignore**
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
rep <- as.integer(Sys.getenv("SIM_REP", "200"))
k.runs <- as.integer(Sys.getenv("SIM_K_RUNS", "100")) # number of bootstrap runs
skip.bootstrap <- tolower(Sys.getenv("SIM_SKIP_BOOTSTRAP", "false")) %in% c("1", "true", "yes", "y")
max.attempts <- as.integer(Sys.getenv("SIM_MAX_ATTEMPTS", "100"))
output.dir <- Sys.getenv("SIM_OUTPUT_DIR", here::here())
dir.create(output.dir, recursive=TRUE, showWarnings=FALSE)
n <- 100
ni_train <- 15 # number of training points
test_points <- c(16, 18, 20)
ni_test <- length(test_points) # number of test values
ni <- ni_train+ni_test
 
N <- n*ni
ti <- c(seq(1, ni_train), test_points)/ni_train


patid <- rep(1:n, each=ni)
day <- rep(ti, n)
uniqueID <- seq(1:n)   

beta <- c(2.5, 3.0, 7.5)
gamma <- c(5.2, 1.6, -1.2)
d <- c(0.4, 1.0)
Mat <- matrix(c(1,0.5,  0.5, 1), ncol=2)

alpha0 <- -8
alpha1 <-  2
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


################################# Simulation runs

est.NLME <- sd.NLME <- est.OS <- sd.OS <- est.TS <- sd.TS <-  est.JM <- sd.JM <- c()
pred.cover.LB<- pred.cover.OS <- pred.cover.TS <- pred.cover.JM <- pred.cover.JM.bt <-  c()
pred.bias.LB <- pred.bias.OS <- pred.bias.TS <- pred.bias.JM <- pred.bias.JM.bt <-  c()
pred.mse.LB <- pred.mse.OS <-  pred.mse.TS <- pred.mse.JM <- pred.mse.JM.bt <- c()
alpha.NLME <- c()
sd.bt.JM <- sd.bt1.JM <- sd.bt2.JM <- c()
runs.bt1 <- runs.bt2 <- c()

for(k in 1:rep){
  cat("This is run", k, "\n")
  
  nlme.fit <- cd4.fit <- OS <- TS <- JM <- 0
  class(nlme.fit) <- class(cd4.fit) <- class(OS) <-  class(TS) <- class(JM)<- "try-error"
  OSconvg <- TSconvg <- JMconvg <-  FALSE
  attempt <- 0
  
  while(class(nlme.fit)[1]=="try-error"|
        class(cd4.fit)=="try-error" | 
        class(OS)=="try-error" | class(JM)=="try-error" | class(TS)=="try-error" |
        OSconvg==FALSE | TSconvg==FALSE| JMconvg==FALSE){
    attempt <- attempt + 1
    if(attempt > max.attempts){
      stop("Exceeded SIM_MAX_ATTEMPTS while trying to obtain converged NLME, OS, TS, and JM fits.")
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
      
      
      dati <- data.frame(patid=uniqueID[i], day=ti, lgcopy=lgcopyi, cd4=cdi_obs )
      
      simdat <- rbind(simdat, dati)
    }
    
    simdat <- simdat %>% dplyr::arrange(patid, day)
    
    simdat_train <- simdat %>% filter(day <= 1)
    
    simdat_test <- NULL
    for(m in c(1:ni_test)){
      simdat_test[[m]] <- simdat %>% filter(day==test_points[m]/ni_train)
    }
    
    
    
    cat("--Done \n")
    ########################## Modeling simdat
    #### nlme()
    cat("--Fit NLME() \n")
    dat_g <- groupedData(lgcopy~day|patid, data=simdat_train)
    nlme.fit <- withTimeout({try(nlme(lgcopy~nf1(p1,p2,p3, day),fixed = p1+p2+p3 ~1,random  = p1+p3 ~1, 
                                      control=nlmeControl(msMaxIter = 500),
                                      data=dat_g,start=beta))},timeout=1.5,onTimeout="warning")
    
    cat("--Done \n")
    if(class(nlme.fit)[1]!="try-error"){
      
      #### nlme prediction for test data
      
      for(m in c(1:ni_test)){
      ## coverage
      simdat_test[[m]]$lgcopy.pred.LB <- predict(nlme.fit, simdat_test[[m]], level=1)
      upper <- simdat_test[[m]]$lgcopy.pred.LB+1.96*nlme.fit$sigma
      lower <- simdat_test[[m]]$lgcopy.pred.LB-1.96*nlme.fit$sigma
      simdat_test[[m]]$cover.LB <- lower<=simdat_test[[m]]$lgcopy & simdat_test[[m]]$lgcopy<=upper
      ## bias
      simdat_test[[m]] <- simdat_test[[m]] %>% mutate(bias.LB=(lgcopy.pred.LB-lgcopy)/abs(lgcopy)*100)
      # MSE
      simdat_test[[m]] <- simdat_test[[m]] %>% mutate(mse.LB=abs(lgcopy.pred.LB-lgcopy)/abs(lgcopy)*100)
      }
      
      
      #### One step 
      
      ########################## model specification
      
      # residual dispersion model:  
      sigmaObject_OS <- list(
        model=~1+(1|patid),
        link='log',
        ran.dist="normal",
        str.fixed=2*log(nlme.fit$sigma),
        lower.fixed=NULL,
        upper.fixed=NULL,
        fixName="alpha",
        ranName="a",
        dispName="siga",
        str.disp=sigma_a
      )
      
      
      nlmeObject_OS <- list(
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
        sigma=sigmaObject_OS,    # residual dispersion model (include residual random eff)
        ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
        str.fixed=fixef(nlme.fit),  # starting value for fixed effect
        str.disp=as.numeric(VarCorr(nlme.fit)[,"StdDev"][c("p1", "p3")]),  # starting value for fixed dispersion of random eff
        #str.disp=d,
        lower.fixed=NULL, # lower bounds for fixed eff
        upper.fixed=rep(100,3), # upper bounds for fixed eff
        lower.disp=c(0,0), # lower bounds for fixed dispersion of random eff
        upper.disp=c(Inf,Inf) # upper bounds for  fixed dispersion of random eff
      )
      
      nlmeObjects_OS=list(nlmeObject_OS)
      
      cat("--Runing One-step model\n")
      
      OS <- try(Rnlme(nlmeObjects=nlmeObjects_OS, long.data=simdat_train, 
                      idVar="patid", sd.method="HL", dispersion.SD = TRUE,
                      independent.raneff=FALSE))
      cat("--done\n")
      
      if(class(OS)!="try-error") OSconvg <- OS$convergence
      
      if(class(OS)!="try-error" & OSconvg){
        
        #### Two step prediction for test data
        for(m in c(1:ni_test)){
          simdat_test[[m]]$lgcopy.pred.OS <- nf1(p1=OS$fixedest["beta1"]+OS$Bi["u1"]*OS$dispersion["d1"], 
                                                 p2=OS$fixedest["beta2"], 
                                                 p3=OS$fixedest["beta3"]+OS$Bi["u3"]*OS$dispersion["d3"],
                                                 t=simdat_test[[m]]$day)[[1]]
          
          OS.sigma <- sqrt(exp(OS$fixedest["alpha0"]+OS$Bi["a0"]*OS$dispersion["siga0"]))[[1]]
          
          upper <- simdat_test[[m]]$lgcopy.pred.OS+1.96*OS.sigma
          lower <- simdat_test[[m]]$lgcopy.pred.OS-1.96*OS.sigma
          
          simdat_test[[m]]$cover.OS <- lower<=simdat_test[[m]]$lgcopy & simdat_test[[m]]$lgcopy<=upper
          
          ## bias
          simdat_test[[m]] <- simdat_test[[m]] %>% mutate(bias.OS=(lgcopy.pred.OS-lgcopy)/abs(lgcopy)*100)
          # MSE
          simdat_test[[m]] <- simdat_test[[m]] %>% mutate(mse.OS=abs(lgcopy.pred.OS-lgcopy)/abs(lgcopy)*100)
        }
        
      
      #### Two step
      
      cd4.fit <- try(lme(cd4~day+I(day^2), data=simdat_train, random=~1|patid))
      
      if(class(cd4.fit)!="try-error"){
        
        simdat_train$cd4.pred <- fitted(cd4.fit)
        
        cat("--Runing Two-step model\n")
        
        
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
          str.fixed=fixef(nlme.fit),  # starting value for fixed effect
          str.disp= as.numeric(VarCorr(nlme.fit)[,"StdDev"][c("p1", "p3")]),  # starting value for fixed dispersion of random eff
          #str.disp=d,
          lower.fixed=NULL, # lower bounds for fixed eff
          upper.fixed=rep(100,3), # upper bounds for fixed eff
          lower.disp=c(0,0), # lower bounds for fixed dispersion of random eff
          upper.disp=c(Inf,Inf) # upper bounds for  fixed dispersion of random eff
        )
        
        nlmeObjects_TS=list(nlmeObject_TS)
        
        TS <- try(Rnlme(nlmeObjects=nlmeObjects_TS, long.data=simdat_train, 
                        idVar="patid", sd.method="HL", dispersion.SD = TRUE,
                        independent.raneff=FALSE))
        cat("--done\n")
        if(class(TS)!="try-error") TSconvg <- TS$convergence
        
        if(class(TS)!="try-error" & TSconvg){
          
          #### Two step prediction for test data
          for(m in c(1:ni_test)){
          simdat_test[[m]]$lgcopy.pred.TS <- nf1(p1=TS$fixedest["beta1"]+TS$Bi["u1"]*TS$dispersion["d1"], 
                                                 p2=TS$fixedest["beta2"], 
                                                 p3=TS$fixedest["beta3"]+TS$Bi["u3"]*TS$dispersion["d3"],
                                                 t=simdat_test[[m]]$day)[[1]]
          
          simdat_test[[m]]$cd4.pred.TS <- predict(cd4.fit, simdat_test[[m]], level=1)
          TS.sigma <- sqrt(exp(TS$fixedest["alpha0"]+TS$fixedest["alpha1"]*simdat_test[[m]]$cd4.pred.TS+
                                 TS$Bi["a0"]*TS$dispersion["siga0"]))[[1]]
          
          
          upper <- simdat_test[[m]]$lgcopy.pred.TS+1.96*TS.sigma
          lower <- simdat_test[[m]]$lgcopy.pred.TS-1.96*TS.sigma
          
          simdat_test[[m]]$cover.TS <- lower<=simdat_test[[m]]$lgcopy & simdat_test[[m]]$lgcopy<=upper
          
          ## bias
          simdat_test[[m]] <- simdat_test[[m]] %>% mutate(bias.TS=(lgcopy.pred.TS-lgcopy)/abs(lgcopy)*100)
          # MSE
          simdat_test[[m]] <- simdat_test[[m]] %>% mutate(mse.TS=abs(lgcopy.pred.TS-lgcopy)/abs(lgcopy)*100)
          }
          
          #### JM
          cat("--Runing Joint model\n")
          
          ########################## model specification
          
          #### JM
          
          sigma1 <- list(
            model=NULL,
            str.disp=as.numeric(VarCorr(cd4.fit)[,"StdDev"][2]),
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
            str.fixed=fixef(cd4.fit),  # starting value for fixed effect
            str.disp=as.numeric(VarCorr(cd4.fit)[,"StdDev"][1]),  # starting value for fixed dispersion of random eff
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
            str.fixed=fixef(nlme.fit),  # starting value for fixed effect
            str.disp=as.numeric(VarCorr(nlme.fit)[,"StdDev"][c("p1", "p3")]),  # starting value for fixed dispersion of random eff
            #str.disp=d,  # starting value for fixed dispersion of random eff
            lower.fixed=NULL, # lower bounds for fixed eff
            upper.fixed=rep(100,3), # upper bounds for fixed eff
            lower.disp=c(0,0), # lower bounds for fixed dispersion of random eff
            upper.disp=c(Inf,Inf) # upper bounds for  fixed dispersion of random eff
          )
          
          nlmeObjects_JM <- list(nlmeObject_JM, lmeObject_JM)
          
          JM <- try(Rnlme(nlmeObjects=nlmeObjects_JM, long.data=simdat_train, 
                          idVar="patid", sd.method="HL", dispersion.SD = TRUE,
                          independent.raneff="byModel"))
          cat("--done\n")
          if(class(JM)!="try-error") JMconvg <- JM$convergence
          
          if(class(TS)!="try-error" & TSconvg){
            
            #### Joint model prediction for test data
            for(m in c(1:ni_test)){
            simdat_test[[m]]$lgcopy.pred.JM <- nf1(p1=JM$fixedest["beta1"]+JM$Bi["u1"]*JM$dispersion["d1"], 
                                                   p2=JM$fixedest["beta2"], 
                                                   p3=JM$fixedest["beta3"]+JM$Bi["u3"]*JM$dispersion["d3"],
                                              t=simdat_test[[m]]$day)[[1]]
          
            simdat_test[[m]]$cd4.pred.JM <- nf2(p1=JM$fixedest["gamma1"]+JM$Bi["b1"]*JM$dispersion["sigb1"], 
                                                p2=JM$fixedest["gamma2"],
                                                p3=JM$fixedest["gamma3"], 
                                                t=simdat_test[[m]]$day)[[1]]
            
            JM.sigma <- sqrt(exp(JM$fixedest["alpha0"]+
                                   JM$fixedest["alpha1"]*simdat_test[[m]]$cd4.pred.JM+
                                   JM$Bi["a0"]*JM$dispersion["siga0"]))[[1]]
            
            upper <- simdat_test[[m]]$lgcopy.pred.JM+1.96*JM.sigma
            lower <- simdat_test[[m]]$lgcopy.pred.JM-1.96*JM.sigma
            
            simdat_test[[m]]$cover.JM <- lower<=simdat_test[[m]]$lgcopy & simdat_test[[m]]$lgcopy<=upper
            
            ## bias
            simdat_test[[m]] <- simdat_test[[m]] %>% mutate(bias.JM=(lgcopy.pred.JM-lgcopy)/abs(lgcopy)*100)
            # MSE
            simdat_test[[m]] <- simdat_test[[m]] %>% mutate(mse.JM=abs(lgcopy.pred.JM-lgcopy)/abs(lgcopy)*100)
            }
          }
        }
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
    simdat_test <- map(simdat_test, function(t){
      t$cover.JM.bt <- NA_real_
      t$bias.JM.bt <- NA_real_
      t$mse.JM.bt <- NA_real_
      t
    })
  } else {
    cat("--Runing Bootstrapping SE for Joint model\n\n")
    JM.SD.BT <- get_sd_bootstrap_with_pred(Rnlme.fit=JM, simdat_train, simdat_test, a0_dist=sigma2$ran.dist, a0_df=sigma2$df,
                                  at.rep=k ,k.runs=k.runs, independent.raneff = "byModel")
    cat("--done\n")
    simdat_test <- JM.SD.BT$simdat_test
  }
  
  runs.bt1 <- c(runs.bt1, JM.SD.BT$runs.bt1)
  runs.bt2 <- c(runs.bt2, JM.SD.BT$runs.bt2)
  
  ############### store output
  #### nlme
  est.NLME <- rbind(est.NLME, fixef(nlme.fit))
  alpha.NLME <- c(alpha.NLME, 2*log(nlme.fit$sigma))
  sd.NLME <- rbind(sd.NLME, summary(nlme.fit)$tTable[,"Std.Error"])
  pred.cover.LB <- rbind(pred.cover.LB,  map_dbl(simdat_test, function(t){mean(t$cover.LB)}))
  pred.bias.LB <- rbind(pred.bias.LB,  map_dbl(simdat_test, function(t){mean(t$bias.LB)}))
  pred.mse.LB <- rbind(pred.mse.LB,  map_dbl(simdat_test, function(t){mean(t$mse.LB)}))
  
  #### OS
  est.OS <- rbind(est.OS, OS$fixedest)
  sd.OS <- rbind(sd.OS, OS$fixedSD)
  pred.cover.OS <- rbind(pred.cover.OS,  map_dbl(simdat_test, function(t){mean(t$cover.OS)}))
  pred.bias.OS <- rbind(pred.bias.OS,  map_dbl(simdat_test, function(t){mean(t$bias.OS)}))
  pred.mse.OS <- rbind(pred.mse.OS,  map_dbl(simdat_test, function(t){mean(t$mse.OS)}))
  
  #### TS
  est.TS <- rbind(est.TS, c(TS$fixedest, fixef(cd4.fit)))
  sd.TS <- rbind(sd.TS, c(TS$fixedSD, summary(cd4.fit)$tTable[,"Std.Error"]))
  pred.cover.TS <- rbind(pred.cover.TS,  map_dbl(simdat_test, function(t){mean(t$cover.TS)}))
  pred.bias.TS <- rbind(pred.bias.TS,  map_dbl(simdat_test, function(t){mean(t$bias.TS)}))
  pred.mse.TS <- rbind(pred.mse.TS,  map_dbl(simdat_test, function(t){mean(t$mse.TS)}))
  
  #### JM
  est.JM <- rbind(est.JM, JM$fixedest)
  sd.JM <- rbind(sd.JM, JM$fixedSD)
  sd.bt.JM <- rbind(sd.bt.JM, JM.SD.BT$se.bt)
  sd.bt1.JM <- rbind(sd.bt1.JM, JM.SD.BT$se.bt1)
  sd.bt2.JM <- rbind(sd.bt2.JM, JM.SD.BT$se.bt2)
  
  pred.cover.JM <- rbind(pred.cover.JM,  map_dbl(simdat_test, function(t){mean(t$cover.JM)}))
  pred.bias.JM <- rbind(pred.bias.JM,  map_dbl(simdat_test, function(t){mean(t$bias.JM)}))
  pred.mse.JM <- rbind(pred.mse.JM,  map_dbl(simdat_test, function(t){mean(t$mse.JM)}))
  
  pred.cover.JM.bt <- rbind(pred.cover.JM.bt,  map_dbl(simdat_test, function(t){mean(t$cover.JM.bt)}))
  pred.bias.JM.bt <- rbind(pred.bias.JM.bt,  map_dbl(simdat_test, function(t){mean(t$bias.JM.bt)}))
  pred.mse.JM.bt <- rbind(pred.mse.JM.bt,  map_dbl(simdat_test, function(t){mean(t$mse.JM.bt)}))
}
##################### organize output 

colnames(est.OS) <- colnames(sd.OS) <- colnames(est.JM)[c(1:4)]
colnames(est.TS) <- colnames(sd.TS) <-   colnames(sd.JM) <- colnames(est.JM)
colnames(sd.bt.JM) <- colnames(sd.bt1.JM) <-colnames(sd.bt2.JM) <- colnames(est.JM)
colnames(pred.cover.LB) <- colnames(pred.cover.OS) <- colnames(pred.cover.TS) <- colnames(pred.cover.JM) <- colnames(pred.cover.JM.bt) <- test_points
colnames(pred.bias.LB) <- colnames(pred.bias.OS) <- colnames(pred.bias.TS) <- colnames(pred.bias.JM) <- colnames(pred.bias.JM.bt)  <- test_points
colnames(pred.mse.LB) <- colnames(pred.mse.OS) <- colnames(pred.mse.TS) <- colnames(pred.mse.JM) <- colnames(pred.mse.JM.bt) <- test_points

NLME.out <- list(True=beta,Est=est.NLME, SD=sd.NLME, 
                 Pred.Bias=pred.bias.LB,
                 Pred.MSE=pred.mse.LB,
                 Pred.Coverage=pred.cover.LB
                 )
OS.out <- list(True=c(beta, alpha0),Est=est.OS, SD=sd.OS, 
               Pred.Bias=pred.bias.OS,
               Pred.MSE=pred.mse.OS,
               Pred.Coverage=pred.cover.OS
)
TS.out <- list(True=true.fixed,Est=est.TS, SD=sd.TS, 
               Pred.Bias=pred.bias.TS,
               Pred.MSE=pred.mse.TS,
               Pred.Coverage=pred.cover.TS
               )
JM.out <- list(True=true.fixed,Est=est.JM, SD=sd.JM, SD.BT=sd.bt.JM, 
               Pred.Bias=pred.bias.JM,
               Pred.MSE=pred.mse.JM,
               Pred.Coverage=pred.cover.JM,
               Pred.Bias.BT=pred.bias.JM.bt,
               Pred.MSE.BT=pred.mse.JM.bt,
               Pred.Coverage.BT=pred.cover.JM.bt,
               SD.BT1=sd.bt1.JM, SD.BT2=sd.bt2.JM)



saveRDS(list(NLME.out=NLME.out,OS.out=OS.out, TS.out=TS.out,JM.out=JM.out, alpha.NLME=alpha.NLME), 
        file.path(output.dir, "s5.rds"))

#save.image(here::here("s5.RData"))



