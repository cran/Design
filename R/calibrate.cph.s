#Resampling optimism of reliability of a Cox survival model
#For predicting survival at a fixed time u, getting grouped K-M estimates
#with avg. of m subjects in a group, or using cutpoints cuts if present,
#e.g. cuts=c(0,.2,.4,.6,.8,1).
#B: # reps  method=see predab.resample
#bw=T to incorporate backward stepdown (using fastbw) with params rule,type,sls
#pr=T to print results of each rep
#what="observed" to get optimism in observed (Kaplan-Meier) survival for
#groups by predicted survival
#what="observed-predicted" to get optimism in KM - Cox - more suitable if
#distributions of predicted survival vary greatly withing quantile groups
#defined from original sample's predicted survival

calibrate.cph <- function(fit,method="boot",u,m=150,cuts,B=40,
                          bw=FALSE,rule="aic",
                          type="residual",sls=.05,aics=0,
                          pr=FALSE,what="observed-predicted",tol=1e-12, ...)	{

  call <- match.call()
                                        #.Options$digits <- 3  14Sep00
  oldopt <- options(digits=3)
  on.exit(options(oldopt))
  unit <- fit$units
  if(unit=="") unit <- "Day"
  ssum <- fit$surv.summary  ##14Nov00
  if(!length(ssum)) stop('did not use surv=T for cph( )')
  cat("Using Cox survival estimates at ", dimnames(ssum)[[1]][2],
      " ", unit,"s\n",sep="")
  surv.by.strata <- ssum[2,,1] #2nd time= at u, all strata
  xb <- fit$linear.predictors
  if(length(stra <- Surv.strata(xb))) 
    surv.by.strata <- surv.by.strata[stra]
  survival <- surv.by.strata^exp(xb)
  if(missing(cuts)) {
    g <- max(1,floor(length(xb)/m))
    cuts <- quantile(c(0,1,survival), seq(0,1,length=g+1),na.rm=TRUE)
  }


  distance <- function(x,y,fit,iter,u,fit.orig,what="observed",
                       orig.cuts, ...) {
    ##Assumes y is matrix with 1st col=time, 2nd=event indicator, 3rd=strata
    ## This is now invalid for Surv objects.  strata info stored in attribute.

    if(sum(y[,2])<5)return(NA)
    surv.by.strata <- fit$surv.summary[2,,1]
    ##2 means to use estimate at first time past t=0 (i.e., at u)
    if(length(stra <- Surv.strata(y)))
      surv.by.strata <- surv.by.strata[stra] #Get for each stratum in data
    cox <- surv.by.strata^exp(x - fit$center)
    ##Assumes x really= x * beta
    pred.obs <- groupkm(cox,Surv(y[,1],y[,2]),u=u,cuts=orig.cuts)
    if(what=="observed") dist <- pred.obs[,"KM"]	else
    dist <- pred.obs[,"KM"] - pred.obs[,"x"]
    if(iter==0) {
      print(pred.obs)
      ##		assign("pred.obs", pred.obs, where=1)   18Apr01
      storeTemp(pred.obs, "pred.obs")  #Store externally for plotting
    }
    
    dist
  }

  coxfit <- function(x,y,u,iter=0, ...) {
    etime <- y[,1]
    e <- y[,2]
    stra <- Surv.strata(y)
      
    ##	attr(x,"strata") <- strata
    if(sum(e)<5)return(list(fail=TRUE))
    ##	f <- coxph(x,etime,e,surv="summary",time.inc=u,pr=F)
    x <- x	#Get around lazy evaluation creating complex expression
    f <- if(length(x)) {
#      cph(Surv(etime,e) ~ x + strat(stra), surv=TRUE, time.inc=u)
      if(length(stra)) {
        cph(Surv(etime,e) ~ x + strat(stra), surv=TRUE, time.inc=u)
      } else {
        cph(Surv(etime,e) ~ x, surv=TRUE, time.inc=u)
      }
    }else{
      cph(Surv(etime,e) ~ strat(stra), surv=TRUE, time.inc=u)
    }
    ## length(x)==0 case 25apr03
    ##Get predicted survival at times 0, u, 2u, 3u, ...
    attr(f$terms, "Design") <- NULL
    ##	f$non.slopes <- f$assume <- f$assume.code <- f$assign <- f$name <- NULL
    ##Don't fool fastbw called from predab.resample
    f
  }

  b <- min(10,B)
  overall.reps <- max(1,round(B/b)) #Bug in S prevents>10 loops in predab.resample
  cat("\nAveraging ", overall.reps," repetitions of B=",b,"\n\n")
  rel <- 0
  opt <- 0
  nrel <- 0
  B <- 0

  for(i in 1:overall.reps) {

    reliability <-
      predab.resample(fit, method=method,
                      fit=coxfit,measure=distance,
                      pr=pr, B=b, bw=bw, rule=rule, type=type,  
                      u=u, m=m, what=what, sls=sls, aics=aics, strata=TRUE,
                      orig.cuts=cuts, tol=tol, ...)
    n <- reliability[,"n"]
    rel <- rel + n * reliability[,"index.corrected"]
    opt <- opt + n * reliability[,"optimism"]
    nrel <- nrel + n
    B <- B + max(n)	
    ##	cat("Memory used after ",i," overall reps:",memory.size(),"\n")
  }

  mean.corrected <- rel/nrel
  mean.opt <- opt/nrel
  rel <- cbind(mean.optimism=mean.opt,mean.corrected=mean.corrected,n=nrel)
  cat("\nMean over ",overall.reps," overall replications\n\n")
  print(rel)

  pred <- pred.obs[,"x"]
  KM <- pred.obs[,"KM"]
  obs.corrected <- KM - mean.opt

  e <- fit$y[,2]

  structure(cbind(reliability[,c("index.orig","training","test"),drop=FALSE],
                  rel,mean.predicted=pred,KM=KM,
                  KM.corrected=obs.corrected,std.err=pred.obs[,"std.err",drop=FALSE]),
            class="calibrate", u=u, units=unit, n=length(e), d=sum(e),
            p=length(fit$coef), m=m, B=B, what=what, call=call)
}

