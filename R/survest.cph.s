##Use x= if input is a design matrix, newdata= if a data frame or data matrix
##or vector.  Can specify (centered) linear predictor values instead (linear.predictors).
##Strata is attached to linear.predictors or x as "strata" attribute.
##data matrix assumes that categorical variables are coded with integer codes
##7Jun99: took away checks on consistency between type and type used
##in fit
##18Sep99: added what='parallel' for val.surv


survest.cph <- function(fit, newdata, linear.predictors, x, times, fun,
                        loglog=FALSE, conf.int=.95,
                        type=NULL, vartype=NULL,
                        conf.type=c("log-log","log","plain","none"),
                        se.fit=TRUE, what=c("survival","parallel"),
                        individual=FALSE, ...) {

  at <- fit$Design
  if(!length(at)) at <- getOldDesign(fit)

  f <- sum(at$assume.code!=8)		#non-strata factors
  nf <- length(at$name)-f
  num.strata <- if(nf==0)1 else length(fit$strata)
  strata.levels <- fit$strata  ## 7may02

#  N <- if(!missing(newdata)) nrow(newdata) else
#  if(!missing(linear.predictors))length(linear.predictors) else
#  if(!missing(x)) nrow(x) else
#  if(ll <- length(fit$linear.predictors)) ll else
#  if(length(fit$x)) nrow(fit$x)

  conf.type <- match.arg(conf.type)
  what <- match.arg(what)
  if(what=='parallel') {conf.int <- FALSE; conf.type <- 'none'}

  inputData <- !(missing(newdata) && missing(linear.predictors) &&
                 missing(x))
    
  if(!se.fit) conf.int <- 0
  ##stype <- attr(fit$surv,"type")
  ##if(length(stype)==0)stype <- "tsiatis"

  if(individual && (length(fit$x)==0 || length(fit$y)==0 || 
                    attr(fit$y,'type')!='counting')) 
    stop('must specify x=T, y=T, and start and stop time to cph when individual=T')

  if(missing(fun)) fun <-
    if(loglog) function(x) logb(-logb(ifelse(x==0|x==1,NA,x)))
    else function(x) x

  naa <- fit$na.action

  ##First see if use method that depends on x and y being stored in fit

  if(!missing(linear.predictors) && length(fit$surv)==0)
    stop('when using linear.predictors= must have specified surv=T to cph')

  if(length(fit$y) && (f==0 || length(fit$x)) &&
     ((conf.int>0 && f>0) | length(fit$surv)==0) &
     (!missing(newdata) | (missing(linear.predictors) && missing(x))))  {

    if(conf.int>0 && conf.type!="log-log" && length(fit$surv))
      warning(paste("Using conf.type=",conf.type,
                    ".\n Survival estimates stored with fit used conf.type=log-log",sep=""))
    if(!missing(linear.predictors) | !missing(x))
      stop(paste("may not specify linear.predictors or x when survival estimation",
                 "is not using underlying survival estimates stored with surv=T"))

    sf <- function(..., type=NULL, vartype=NULL, cphnull=FALSE) {
      g <- list(...)
      if(length(type))    g$type <- type
      if(length(vartype)) g$vartype <- vartype
      do.call(if(cphnull)'survfit.cph.null' else 'survfit.cph', g)
    }
    if(f==0) {
      g <- sf(fit, se.fit=se.fit, conf.int=conf.int, conf.type=conf.type,
              type=type, vartype=vartype, cphnull=TRUE)
      if(missing(newdata))sreq <- 1 else sreq <- 
        attr(predict(fit, newdata, type="lp", expand.na=FALSE),"strata")
    }
    else {
      if(missing(newdata)) {
        g <- sf(fit, se.fit=se.fit, conf.int=conf.int, conf.type=conf.type,
                type=type, vartype=vartype)
        sreq <- 1
      }
      else {
        if(nrow(newdata)>1 && !individual && missing(times))
          stop("must specify times= if predicting for >1 observation")
        g <- sf(fit, newdata=newdata, se.fit=se.fit, conf.int=conf.int,
                conf.type=conf.type, individual=individual,
                type=type, vartype=vartype)
        sreq <- g$requested.strata
      }
      naa <- g$na.action
    }
    sreq <- oldUnclass(sreq)

    if(missing(times)) {
      ##delete extra S(t) curves added by survfit for all strata
      ##No newdata -> requested underlying survival for all strata
      if(missing(newdata)) return(g)
      else {
        if(nf==0) j <- TRUE
        else {
          stemp <- rep(1:num.strata, g$strata)
          j <- stemp==sreq
        }

        tim <- c(0,g$time[j]); nr <- c(g$n.risk[j][1],g$n.risk[j])
        ne <- c(0,g$n.event[j]); surv <- c(1, g$surv[j])
        se <- c(NA, g$std.err[j])
        upper <- c(NA, g$upper[j]); lower <- c(NA,g$lower[j])

        yy <- fit$y; ny <- ncol(yy)
        str <- oldUnclass(Surv.strata(yy))
        if(length(str)) yy <- yy[str==sreq,ny-1] else yy <- yy[,ny-1]
        maxt <- max(yy)
        if(maxt>tim[length(tim)]) {
          tim <- c(tim,maxt); nr <- c(nr, sum(yy>=maxt-1e-6))
          ne <- c(ne, 0); surv <- c(surv, surv[length(surv)])
          se <- c(se, NA); upper <- c(upper, NA); lower <- c(lower, NA)
        }
        se <- -se/logb(ifelse(surv==0|surv==1,NA,surv))

        surv  <- fun(surv);  surv[is.infinite(surv)] <- NA
        lower <- fun(lower); lower[is.infinite(lower)] <- NA
        upper <- fun(upper); upper[is.infinite(upper)] <- NA
        retlist <- list(time=tim,n.risk=nr,
                        n.event=ne,
                        surv=surv, 
                        std.err=se,
                        upper=upper, lower=lower, 
                        conf.type=g$conf.type,
                        conf.int=g$conf.int, call=g$call)
        if(nf>0) retlist$strata <- sreq  #was g$strata[sreq]
        return(retlist)
      }
    }
    else {
      g <- summary(g, print.it=FALSE, times=times)
      if(!individual && nf>0) { #delete extra cells added by survfit for strat
        if(length(g$time) != length(times)*num.strata)
          stop('summary.survfit could not compute estimates for all strata at all times requested.\nYou probably requested times where data are limited.')
          
        d <- dim(g$surv); if(length(d)==0) d <- c(length(g$surv),1)
        strata.col <- matrix(rep(sreq,d[1]),ncol=d[2],byrow=TRUE)
        gs <- factor(g$strata, strata.levels)
        strata.row <- matrix(rep(oldUnclass(gs),d[2]),ncol=d[2]) # was g$strata
        m <- strata.col==strata.row
        g$surv    <- matrix(g$surv[m],   ncol=d[2])[,,drop=TRUE]
        g$lower   <- matrix(g$lower[m],  ncol=d[2])[,,drop=TRUE]
        g$upper   <- matrix(g$upper[m],  ncol=d[2])[,,drop=TRUE]
        g$std.err <- matrix(g$std.err[m],ncol=d[2])[,,drop=TRUE]
      }
    }
    tim  <- g$time
    nr   <- g$n.risk
    ne   <- g$n.event
    surv <- g$surv
    se   <- g$std.err
    low  <- g$lower
    up   <- g$upper
    tim <- unique(tim)
    if(is.matrix(surv)) {
      surv <- t(surv); se <- t(se); low <- t(low); up <- t(up)
      dn <- list(row.names(newdata),format(tim))
      dimnames(surv) <- dn; dimnames(se) <- dn; dimnames(low) <- dn;
      dimnames(up) <- dn
    }

    se   <- -se/logb(ifelse(surv==0|surv==1,NA,surv))
    if(!.R.) {
      storage.mode(surv) <- "single"; storage.mode(se) <- "single"
      storage.mode(low)  <- "single"; storage.mode(up) <- "single"
    }
    surv <- fun(surv); low <- fun(low); up <- fun(up)
    surv[is.infinite(surv)] <- NA
    low[is.infinite(low)] <- NA
    up[is.infinite(up)] <- NA
    retlist <- list(time=tim, surv=naresid(naa,surv), std.err=naresid(naa,se)
                    ,lower=naresid(naa,low), upper=naresid(naa,up))
    if(nf>0) retlist$strata <- naresid(naa,sreq)
    return(retlist)
  }
  
#  strata.levels <- fit$strata  7may02
  asnum.strata <- function(str, strata.levels) {
    if(!length(str)) return(NULL)  # 4Aug01; thanks Mike Kattan
    if(is.numeric(str) && any(str<1 | str>length(strata.levels)))
      stop('illegal stratum number')

    if(is.category(str) || is.numeric(str)) return(as.integer(str))

    i <- match(str, strata.levels, nomatch=0)
    if(any(i==0)) stop(paste('illegal strata:',
             paste(str[i==0],collapse=' ')))
    i
  }

  ##Instead use the baseline survival computed at fit time with cph(...,surv=T)

  nt <- if(missing(times))0 else length(times)
  if(conf.int>0 && f>0)
    warning(paste("S.E. and confidence intervals are approximate except",
                  "at predictor means.\nUse cph(...,x=T,y=T) (and don't use linear.predictors=) for better estimates."))
#  if(!missing(type) && method!=stype)
#    warning(paste("method=",stype,
#                  " used since\nusing survival estimates stored with fit",sep=""))
  if(conf.type!="log-log" && conf.type!="none")
	stop("only conf.type=log-log can be used since using survival estimates stored with fit")

  if(missing(linear.predictors)) {
    if(missing(x) && missing(newdata)) {
      linear.predictors <- fit$linear.predictors	#assume was centered
      rnam <- names(linear.predictors)
      if(length(linear.predictors)==0) {
        if(length(fit$x)==0)
          stop("newdata, x, linear.predictors not given but x nor linear.predictors stored in fit")
        linear.predictors <- matxv(fit$x, fit$coef)-fit$center
        strata <- fit$strata
        rnam <- dimnames(fit$x)[[1]]
      }
      else strata <- attr(linear.predictors,"strata")
    }
    else {
      if(missing(x)) 
        {x <- predict(fit,newdata,type="x",expand.na=FALSE)
         naa <- attr(x,"na.action")}
      strata <- attr(x,"strata")
      if(f>0) linear.predictors <- matxv(x,fit$coef) - fit$center 
      else linear.predictors <- 0
      rnam <- dimnames(x)[[1]]
    }
  }
  else {
    strata <- asnum.strata(attr(linear.predictors, "strata"),strata.levels)
    rnam <- names(linear.predictors)
  }
  if(length(strata)==0 && nf>0) 
    stop("strata not stored in x or linear.predictors")
  attr(strata, "class") <- NULL

  if(length(fit$surv)==0 && length(fit$x)==0 && length(fit$y)==0)
    stop("you did not specify surv=T or x=T, y=T in cph")

  if(conf.int>0) zcrit <- qnorm((conf.int+1)/2)
  if(length(strata)==0) {
	n <- length(linear.predictors)
	strata <- rep(1,n)
	ns <- 1
  } else {
    ns <- max(strata, na.rm=TRUE)   #na.rm added 8Jun94
    n <- length(strata)
  }

  if(what=='parallel') {
    if(length(times)>1 && length(times) != n)
      stop('length of times must equal 1 or number of subjects being predicted')
    if(!length(fit$surv))stop('must specify surv=T to cph')
    if(diff(range(strata))==0) {
      estsurv <- approx(fit$time, fit$surv, xout=times,
                        method="constant", f=0)$y
      return(estsurv ^ exp(linear.predictors))
    }
    est.surv <- if(.R.)double(n) else single(n)
    for(zs in unique(strata)) {
      this <- strata==zs
      estsurv <- approx(fit$time[[zs]], fit$surv[[zs]],
                        xout=if(length(times)==1)times else times[this],
                        method='constant', f=0)$y
      est.surv[this] <-
        estsurv ^ exp(if(length(linear.predictors)==1)
                      linear.predictors else
                      linear.predictors[this])
    }
    return(est.surv)
  }
  
  if(n>1 && nt==0)
    stop("must specify times if getting predictions for >1 obs.")
  if(nt==0)	{  #Get est for 1 obs
	if(!is.list(fit$time)) {
      times <- fit$time
      surv <- fit$surv^exp(linear.predictors)
      std.err <- fit$std.err
    }	else {
      times <- fit$time[[strata]]
      surv <- fit$surv[[strata]]^exp(linear.predictors)
      std.err <- fit$std.err[[strata]]
    }
	if(conf.int>0) {
      lower <- surv^exp(zcrit*std.err)
      upper <- surv^exp(-zcrit*std.err)
      lower[1] <- 1
      upper[1] <- 1
      attr(lower,"type") <- NULL
      attr(upper,"type") <- NULL
    }
    surv <- fun(surv); surv[is.infinite(surv)] <- NA
    if(conf.int>0) {
      lower <- fun(lower); lower[is.infinite(lower)] <- NA
      upper <- fun(upper); upper[is.infinite(upper)] <- NA
    }
    if(!.R.) {
      storage.mode(times) <- "single"
      storage.mode(surv) <- "single"
    }
	if(conf.int>0 && !.R.) {      
      storage.mode(std.err) <- "single"
      storage.mode(lower) <- "single"
      storage.mode(upper) <- "single"
    }
	if(!.R.) storage.mode(linear.predictors) <- "single"
	if(nf==0) strata <- NULL
	retlist <- list(time=times,surv=surv,linear.predictors=linear.predictors)
	if(conf.int>0) retlist <- 
      c(retlist,list(lower=lower,upper=upper,std.err=std.err))
	if(nf>0) {
      retlist$strata <- strata
      retlist$requested.strata <- oldUnclass(strata)
    }
	return(retlist)
  }

  ##Selected times for >=1 obs
  ##if(conf.int>0)
  ## warning("for >1 observation being predicted, does not store confidence intervals")
  ##First get survival at times "times" for each stratum

  surv <- matrix(if(.R.)double(1) else single(1),nrow=ns,ncol=nt)
  serr <- matrix(if(.R.)double(1) else single(1),nrow=ns,ncol=nt)
  for(i in 1:ns) {
	if(!is.list(fit$time)) {
      tim <- fit$time
      se  <- fit$std.err
      srv <- fit$surv
    }	else {
      tim <- fit$time[[i]]
      se  <- fit$std.err[[i]]
      srv <- fit$surv[[i]]
    }
	m <- length(tim)
	j <- 0
	for(u in times)	{
      j <- j + 1
      ##		tm <- max((1:length(tim))[max(tim[tim<=u])==tim])
      tm <- max((1:length(tim))[tim<=u+1e-6])
      s <- srv[tm]; Se <- se[tm]
      ##		if(s != min(srv[tim<=u])) warning()
      if(u > tim[m] && srv[m]>0) {s <- NA; Se <- NA}
      surv[i,j] <- s; serr[i,j] <- Se
    }
  }
  srv <- surv[strata,]^exp(linear.predictors)		#was surv[strata,1:nt]
  ft <- format(times)
  if(is.matrix(srv)) {
    dn <- list(rnam, ft)
    dimnames(srv) <- dn	} else names(srv) <- if(n==1) ft else rnam
  if(conf.int>0) {
    serr <- serr[strata,]
    lower <- srv^exp(zcrit*serr); upper <- srv^exp(-zcrit*serr)
    if(is.matrix(lower)) {
      dimnames(serr) <- dn; dimnames(lower) <- dn; dimnames(upper) <- dn
    }
    else {
      names(serr) <- names(lower) <- names(upper) <- if(n==1) ft else rnam
    }
    lower <- fun(lower); lower[is.infinite(lower)] <- NA
    upper <- fun(upper); upper[is.infinite(upper)] <- NA
    if(!.R.) {
      storage.mode(serr)  <- "single"
      storage.mode(lower) <- "single"
      storage.mode(upper) <- "single"
    }
  }

  srv <- fun(srv); srv[is.infinite(srv)] <- NA
  if(!.R.) storage.mode(srv) <- "single"
  nar <- if(inputData) function(naa,w) w else function(...) naresid(...)
  ## 15mar04 nar            
  if(conf.int==0)
    return(list(time=times, surv=nar(naa,srv))) #was return(srv)
  retlist <- list(time=times, surv=nar(naa,srv), lower=nar(naa,lower),
                  upper=nar(naa,upper), std.err=nar(naa,serr))
  if(nf>0) retlist$requested.strata <- nar(naa,oldUnclass(strata))
  retlist
}
