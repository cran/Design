cph <- function(formula=formula(data),
                data=if(.R.)parent.frame() else sys.parent(),
                weights, subset, na.action=na.delete, 
                method=c("efron","breslow","exact",
                  "model.frame", "model.matrix"),
                singular.ok=FALSE, robust=FALSE,
                model=FALSE, x=FALSE, y=FALSE, se.fit=FALSE,
                eps=1e-4, init, iter.max=10, tol=1e-9,
                surv=FALSE, time.inc, type, vartype, conf.type, ...) {

  if(.R.) require('survival')
  getN <- if(.R.) function(obj) {
    if(existsFunction(obj)) get(obj) else
       getFromNamespace(obj,'survival')
 } else function(obj) get(obj)

  method <- match.arg(method)
  call <- match.call()
  m <- match.call(expand=FALSE)
  m$na.action <- na.action    ## 30apr02
  m$method <- m$model <- m$x <- m$y <- m$... <- m$se.fit <-
    m$type <- m$vartype <-
    m$surv <- m$time.inc <- m$eps <- m$init <- m$iter.max <- m$tol <-
      m$weights <- m$singular.ok <- m$robust <- NULL
  m$na.action <- na.action
  if(.R.) m$drop.unused.levels <- TRUE  ## 31jul02
  m[[1]] <- as.name("model.frame")

  if (!inherits(formula,"formula")) {
	## I allow a formula with no right hand side
	## The dummy function stops an annoying warning message "Looking for
	##  'formula' of mode function, ignored one of mode ..."
	if (inherits(formula,"Surv")) {
      xx <- function(x) formula(x)
      formula <- xx(paste(deparse(substitute(formula)), 1, sep="~"))
    }
	else stop("Invalid formula")
  }
  m$formula <- formula

  nstrata <- 0; Strata <- NULL

  if(!missing(data) || (length(z <- attr(terms(formula),"term.labels"))>0 &&
                        any(z!="."))) { #X's present
    if(.R.) {
      dul <- .Options$drop.unused.levels
      if(!length(dul) || dul) {
        on.exit(options(drop.unused.levels=dul))
        options(drop.unused.levels=FALSE)
      }
    }
    X <- Design(eval(m, if(.R.) parent.frame() else sys.parent()))
    atrx <- attributes(X)
    atr  <- atrx$Design
    nact <- atrx$na.action
    if(method=="model.frame") return(X)

#    atr <- attr(attr(X,'terms'),'Design')
    Terms <- if(missing(data)) terms(formula, 
                                     specials=c("strat","cluster")) else
    terms(formula, 
          specials=c("strat","cluster"), data=data)

    asm   <- atr$assume.code
    name  <- atr$name

    cluster <- attr(Terms, "specials")$cluster
    stra <- attr(Terms, "specials")$strat
    Terms.ns <- Terms
    if(length(cluster)) {
      if(missing(robust)) robust <- TRUE
      tempc <- untangle.specials(Terms.ns, "cluster", 1:10)
      ord <- attr(Terms, 'order')[tempc$terms]
      if(any(ord>1)) stop("Cluster can not be used in an interaction")
      cluster <- strata(X[,tempc$vars], shortlabel=TRUE)  #allow multiples
      Terms.ns <- Terms.ns[-tempc$terms]
    }
    if(length(stra)) {
      temp <- untangle.specials(Terms.ns, "strat", 1)
      Terms.ns <- Terms.ns[-temp$terms]	#uses [.terms function
      ##  Set all factors=2 (-> interaction effect not appearing in main effect
      ##  that was deleted strata effect
      if(!.R.) {
        tfac <- attr(Terms,'factors')
        ## For some reason attr(...) <- pmin(attr(...)) changed a detail
        ## in factors attribute in R but doesn't seem to be needed in
        ## R or SV4 anyway 30apr02
        if(length(tfac) && any(tfac > 1))
          attr(Terms,'factors') <- pmin(tfac, 1)
        tfac <- attr(Terms.ns,'factors')
        if(length(tfac) && any(tfac > 1))
          attr(Terms.ns,'factors') <-  pmin(tfac, 1)
      }

      ##added 21Apr94, if( ) 28Apr02
      Strata <- list()
      for(i in (1:length(asm))[asm==8])	{
        nstrata <- nstrata+1
        xi <- X[[i+1]]
        levels(xi) <- paste(name[i],"=",levels(xi),sep="")
        Strata[[nstrata]] <- xi			}
      names(Strata) <- paste("S",1:nstrata,sep="")
      Strata <- interaction(as.data.frame(Strata),drop=TRUE)
    }
    offs <- offset<- attr(Terms, "offset")  ## offs 23nov02 moved up 6dec02
    ## 23nov02
    if(length(nact$nmiss)) {
      jia <- grep('*%ia%*',names(nact$nmiss))  ## 8feb03
      if(length(jia)) nact$nmiss <- nact$nmiss[-jia]
      s <- if(length(offs)) names(nact$nmiss) !=  atrx$names[offs] else TRUE
      names(nact$nmiss)[s] <- 
        c(as.character(formula[2]), atr$name[atr$assume.code!=9])
    }
    ## [s] 23nov02
    xpres <- length(asm) && any(asm!=8)
    Y <- model.extract(X, 'response')
    if(!inherits(Y,"Surv"))stop("response variable should be a Surv object")
    weights <- model.extract(X, 'weights')
    tt <- length(offset)
    offset <- if(tt == 0) rep(0, nrow(Y))
    else if(tt == 1)
      X[[offset]]
    else {
      ff <- X[[offset[1]]]
      for(i in 2:tt)   # for case with multiple offset terms
        ff <- ff + X[[offset[i]]]
      ff
    }
    
    if(model) m <- X

    ##No mf if only strata factors

    if(!xpres)
      {
        X <- if(.R.)matrix(nrow=0,ncol=0) else NULL
        assign <- NULL
      }
    else
      {
        X <- model.matrix(Terms.ns, X)[,-1,drop=FALSE]
        assign <- attr(X, "assign")
        assign[[1]] <- NULL  # remove intercept position, renumber
      }
    nullmod <- FALSE
  }

  else	# model with no right-hand side
    {
      X <- NULL
      Terms <- terms(formula)
      yy <- attr(terms(formula),"variables")[1]
      Y <- eval(yy,data=data)     #sys.parent(2))
      if(!inherits(Y,"Surv"))stop("response variable should be a Surv object")
      Y <- Y[!is.na(Y)]
      assign <- NULL
      xpres <- FALSE
      nullmod <- TRUE
      nact <- NULL
    }
  ny <- ncol(Y)
  time.units <- attr(Y, "units")
  maxtime <- max(Y[,ny-1])

  rnam <- dimnames(Y)[[1]]
  if(xpres) dimnames(X) <- list(rnam, atr$colnames)

  if(method=="model.matrix") return(X)
  if(!length(time.units)) time.units <- "Day"
  if(missing(time.inc)) {
    time.inc <- switch(time.units,Day=30,Month=1,Year=1,maxtime/10)
    if(time.inc>=maxtime | maxtime/time.inc>25) time.inc <- max(pretty(c(0,maxtime)))/10
  }
  if(nullmod) f <- NULL
  else
    {
      ytype <- attr(Y, "type")
      if( method=="breslow" || method =="efron") {
        if (ytype== 'right')  fitter <- getN("coxph.fit")
        else if (ytype=='counting') fitter <- getN("agreg.fit")
        else stop(paste("Cox model doesn't support \"", ytype,
                        "\" survival data", sep=''))
      }
      else if (method=='exact') fitter <- getN("agexact.fit")
      else stop(paste ("Unknown method", method))
      if (missing(init)) init <- NULL

      ## S-Plus 6 has control parameter.  S-Plus 5 has toler.chol.
      ## Previous to 5 has neither.  R has names(fitter)=NULL
      nf <- names(fitter)
      if(any(nf=='toler.chol'))
        f <- fitter(X, Y, strata=Strata, offset=offset, iter.max=iter.max,
                    eps=eps, weights=weights, init=init,
                    method=method, rownames=rnam, toler.chol=tol) else
      if(.R. || any(nf=='control'))
        f <- fitter(X, Y, strata=Strata, offset=offset,
                    weights=weights, init=init,
                    method=method, rownames=rnam,
                    control=getN('coxph.control')(eps=eps, toler.chol=tol,
                      toler.inf=1, iter.max=iter.max)) else
      f <- fitter(X, Y, strata=Strata, offset=offset, iter.max=iter.max,
                  eps=eps, weights=weights, init=init,
                  method=method, rownames=rnam)
    }
  if (is.character(f)) {
    cat("Failure in cph:\n",f,"\n")
    if(.SV4.)return(structure(list(fail=TRUE,fitFunction='cph'),
                              class='Design')) else 
    return(structure(list(fail=TRUE),class="cph"))
  }   else {
    if(length(f$coef) && any(is.na(f$coef))) {  ## length 1may02
      vars <- names(f$coef)[is.na(f$coef)]
      msg <- paste("X matrix deemed to be singular; variable",
                   paste(vars, collapse=" "))
      if(singular.ok) warning(msg)
	  else {
		cat(msg,"\n")
        if(.SV4.)return(structure(list(fail=TRUE,fitFunction='cph'),
                                  class='Design')) else 
        return(structure(list(fail=TRUE),class="cph")) ##13Nov00
	  }
    }
  }
#  attr(Terms, 'Design') <- atr
  f$terms <- Terms

  if(robust) {
    f$naive.var <- f$var
    ## Terry gets a little tricky here, calling resid before adding
    ## na.action method to avoid re-inserting NAs.  Also makes sure
    ## X and Y are there
    if(missing(cluster)) cluster <- FALSE
    fit2 <- c(f, list(x=X, y=Y, method=method))
    if(length(stra)) fit2$strata <- Strata
    temp <- residuals.coxph(fit2, type='dfbeta', collapse=cluster)
    f$var <- t(temp) %*% temp
  }
  
  if(length(weights) && any(weights!=1)) f$weights <- weights

  nvar <-  length(f$coefficients)

  temp <- factor(Y[,ny], levels=0:1, labels=c("No Event","Event"))
  ## was category 30apr02
  n.table <- if(.R.) {
    if(!length(Strata)) table(temp,dnn='Status') else
    table(Strata, temp, dnn=c('Stratum','Status'))
  } else {
    if(!length(Strata)) table(temp) else table(Strata, temp)
  }
  f$n <- n.table
  nnn <- nrow(Y)
  nevent <- sum(Y[,ny])
  if(xpres)	{
    logtest <- -2 * (f$loglik[1] - f$loglik[2])
    R2.max <- 1 - exp(2*f$loglik[1]/nnn)
    R2 <- (1 - exp(-logtest/nnn))/R2.max
    P <- 1-pchisq(logtest,nvar)
    stats <- c(nnn, nevent, logtest, nvar, P, f$score, 
               1-pchisq(f$score,nvar), R2)
    names(stats) <- c("Obs", "Events", "Model L.R.", "d.f.", "P", 
                      "Score", "Score P","R2")
  }
  else { stats <- c(nnn, nevent); names(stats) <- c("Obs","Events") }
  f$method <- NULL
  if(xpres) dimnames(f$var) <- list(atr$colnames, atr$colnames)  #2Apr95
  f <- c(f, list(call=call, Design=atr, 
                 assign=DesignAssign(atr, 0, atrx$terms), # was =assign 1may02
                 na.action=nact,  #terms=Terms
                 fail = FALSE, non.slopes = 0, stats = stats, method=method,
                 maxtime = maxtime, time.inc = time.inc,
                 units = time.units, fitFunction=c('cph','coxph')))

  if(xpres) {
    f$center <- sum(f$means*f$coef)
    f$scale.pred <- c("log Relative Hazard","Hazard Ratio")
    attr(f$linear.predictors,"strata") <- Strata
    names(f$linear.predictors) <- rnam	#23Feb93
    if(se.fit) {
      XX <- X - rep(f$means,rep.int(nnn,nvar))   # see scale() function
      ##  XX <- sweep(X, 2, f$means)	# center   (slower)
      se.fit <- drop(((XX %*% f$var) * XX) %*% rep(1,ncol(XX)))^.5
      if(!.R.) storage.mode(se.fit) <- "single"
      names(se.fit) <- rnam
      f$se.fit <- se.fit  	
    }
  }
  if(model) f$model <- m
  if(nstrata > 0) {
    attr(X, "strata") <- attr(Y, "strata") <- Strata
    f$strata <- levels(Strata)
  }
  if(x) f$x <- X
  if(y) f$y <- Y

  if(is.character(surv) || surv) {
    if(!length(Strata)) Strata <- rep(1, nnn)
    else Strata <- oldUnclass(Strata)
    nstr <- max(Strata, na.rm=TRUE)
    srv <- NULL
    tim <- NULL
    s.e. <- NULL
    timepts <- seq(0, maxtime, by=time.inc)
    s.sum <- array(if(.R.)double(1) else single(1),c(length(timepts),nstr,3),
                   list(format(timepts),paste("Stratum",1:nstr),
                        c("Survival","n.risk","std.err")))
    
    if(xpres) {
      g <- f; if(!x) g$x <- X; if(!y) g$y <- Y
      fname <- 'survfit.cph'
    } else {
      g <- f; if(!y) g$y <- Y
      fname <- 'survfit.cph.null'
      }
    g <- list(g)
    if(!missing(type))      g$type      <- type
    if(!missing(vartype))   g$vartype   <- vartype
    if(!missing(conf.type)) g$conf.type <- conf.type
    g <- do.call(fname, g)

    if(nstr==1) stemp <- rep(1, length(g$time))
    else stemp <- rep(1:nstr,g$strata)
    i <- 0
    for(k in 1:nstr) {
      j <- stemp==k; i <- i+1
      yy <- Y[Strata==i,ny-1]
      maxt <- max(yy)
      ##n.risk from surv.fit does not have usual meaning if not Kaplan-Meier
      
      tt <- c(0,g$time[j])
      su <- c(1,g$surv[j])
      se <- c(NA,-g$std.err[j]/logb(g$surv[j]))
      if(!.R.) {
        storage.mode(tt) <- 'single'
        storage.mode(su) <- 'single'
        storage.mode(se) <- 'single'
      }
      if(maxt>tt[length(tt)]) {
        tt <- c(tt, maxt); su <- c(su, su[length(su)]); se <- c(se, NA)
      }
      kk <- 0
      for(tp in timepts) {
        kk <- kk + 1
        ##	  t.choice <- max((1:length(tt))[max(tt[tt<=tp])==tt])
        t.choice <- max((1:length(tt))[tt<=tp+1e-6])
        if(tp > max(tt)+1e-6 & su[length(su)]>0) {
          Su <- NA
          Se <- NA			}	else {
            Su <- su[t.choice]
            Se <- se[t.choice]		 }
        n.risk <- sum(yy>=tp)
        s.sum[kk,i,1:3] <- c(Su, n.risk, Se)		}
      if(!is.character(surv)) {
        if(nstr==1)	{
          tim <- tt
          srv <- su
          s.e. <- se		} else {
            tim <- c(tim, list(tt))
            srv <- c(srv, list(su))
            s.e. <- c(s.e., list(se))		}	}
    }
    if(is.character(surv)) f$surv.summary <- s.sum
    else {
      attr(srv, "type") <- if(missing(type)) method else type #7Jun99
      if(nstr>1) { names(srv) <- names(tim) <- names(s.e.) <- f$strata }
      f <- c(f, list(time=tim, surv=srv,
                     std.err=s.e., surv.summary=s.sum))		
    }
  }
  oldClass(f) <- if(.SV4.) 'Design' else c("cph", "Design", "coxph") ##13Nov00
  f
}

## Define a version of coxph.fit that works in S-Plus post version
## 4.5 as well as in earlier versions, as toler.chol argument was
## added in survival5 for S-Plus 2000 and Unix S-Plus 5.0
## coxphFit also handles the case where toler.chol and 4 other
## arguments were forgotten.  In S-Plus 6 control= is used.

coxphFit <- if(.R. || .SV4.)
  function(..., strata=NULL, rownames=NULL, offset=NULL,
           init=NULL, toler.chol=1e-9, eps=.0001, iter.max=10) {
    if(!existsFunction('coxph.fit'))
      coxph.fit <- getFromNamespace('coxph.fit','survival')
    if(!existsFunction('coxph.control'))
      coxph.control <- getFromNamespace('coxph.control', 'survival')
    res <- coxph.fit(..., strata=strata, rownames=rownames,
                     offset=offset, init=init,
                     control=coxph.control(toler.chol=toler.chol, toler.inf=1,
                       eps=eps, iter.max=iter.max))
    if(is.character(res)) return(list(fail=TRUE))
    if(iter.max > 1 && res$iter >= iter.max) return(list(fail=TRUE))
    res$fail <- FALSE
    res
} else  function(..., strata=NULL, rownames=NULL, offset=NULL,
                     init=NULL, toler.chol=1e-9, eps=.0001, iter.max=10) {

    nf <- names(coxph.fit)  ##13Nov00
    res <-
      if(any(nf=='control'))
        coxph.fit(..., strata=strata, rownames=rownames,
                  offset=offset, init=init,
                  control=coxph.control(toler.chol=toler.chol, toler.inf=1,
                    eps=eps, iter.max=iter.max))
      else if(all(c('toler.chol','eps','iter.max') %in% nf))
        coxph.fit(..., strata=strata, rownames=rownames,
                  offset=offset, init=init, toler.chol=toler.chol,
                  eps=eps, iter.max=iter.max)
      else if(all(c('iter.max','eps') %in% nf))
        coxph.fit(..., strata=strata, rownames=rownames,
                  offset=offset, init=init, eps=eps, iter.max=iter.max)
      else coxph.fit(..., strata=strata, rownames=rownames,
                     offset=offset, init=init)
    if(is.character(res)) return(list(fail=TRUE))
    if(length(res$iter) && iter.max > 1 && res$iter >= iter.max)
      return(list(fail=TRUE))
    res$fail <- FALSE
    res
  }

Survival.cph <- function(object, ...) {
if(!length(object$time) || !length(object$surv))
  stop("did not specify surv=T with cph")
f <- function(times, lp=0, stratum=1, type=c("step","polygon"),
              time, surv) {
  type <- match.arg(type)
  if(length(stratum)>1) stop("does not handle vector stratum")
  if(length(times)==0) {
    if(length(lp)>1) stop("lp must be of length 1 if times=NULL")
    return(surv[[stratum]]^exp(lp))
  }
  s <- matrix(NA, nrow=length(lp), ncol=length(times),
              dimnames=list(names(lp), format(times)))
  if(is.list(time)) {time <- time[[stratum]]; surv <- surv[[stratum]]}
  if(type=="polygon") {
    if(length(lp)>1 && length(times)>1)
      stop('may not have length(lp)>1 & length(times>1) when type="polygon"')
    su <- approx(time, surv, times)$y
    return(su ^ exp(lp))
  }
  for(i in 1:length(times)) {
    tm <- max((1:length(time))[time <= times[i]+1e-6])
    su <- surv[tm]
    if(times[i] > max(time)+1e-6) su <- NA
    s[,i] <- su^exp(lp)
  }
  drop(s)
}
formals(f) <- list(times=NULL, lp=0, stratum=1,
                   type=c("step","polygon"),
                   time=object$time, surv=object$surv)
f
}

Quantile.cph <- function(object, ...) {
if(!length(object$time) || !length(object$surv))
  stop("did not specify surv=T with cph")
f <- function(q=.5, lp=0, stratum=1, type=c("step","polygon"), time, surv) {
  type <- match.arg(type)
  if(length(stratum)>1) stop("does not handle vector stratum")
  if(is.list(time)) {time <- time[[stratum]]; surv <- surv[[stratum]]}
  Q <- matrix(NA, nrow=length(lp), ncol=length(q),
              dimnames=list(names(lp), format(q)))
  for(j in 1:length(lp)) {
    s <- surv^exp(lp[j])
    if(type=="polygon") Q[j,] <- approx(s, time, q)$y
    else for(i in 1:length(q))
      if(any(s <= q[i])) Q[j,i] <- min(time[s<=q[i]])  #is NA if none
    ## if(any()) 20may02
  }
  drop(Q)
}
formals(f) <- list(q=.5, lp=0, stratum=1,
                   type=c('step','polygon'),
                   time=object$time, surv=object$surv)
f
}


Mean.cph <- function(object, method=c("exact","approximate"),
type=c("step","polygon"), n=75, tmax, ...) {
method <- match.arg(method)
type   <- match.arg(type)

if(!length(object$time) || !length(object$surv))
  stop("did not specify surv=T with cph")

if(method=="exact") {
  f <- function(lp=0, stratum=1, type=c("step","polygon"),
                tmax=NULL, time, surv) {
    type <- match.arg(type)
    if(length(stratum)>1) stop("does not handle vector stratum")
    if(is.list(time)) {time <- time[[stratum]]; surv <- surv[[stratum]]}
    Q <- lp
    if(!length(tmax)) {
      if(min(surv)>1e-3)
        warning(paste("Computing mean when survival curve only defined down to",
                      format(min(surv)),"\n Mean is only a lower limit"))
      k <- rep(TRUE,length(time))
    } else {
      if(tmax>max(time)) stop(paste("tmax=",format(tmax),
                                    "> max follow-up time=",format(max(time))))
      k <- (1:length(time))[time<=tmax]
    }
    for(j in 1:length(lp)) {
      s <- surv^exp(lp[j])
      Q[j] <- if(type=="step") sum(c(diff(time[k]),0) * s[k]) else 
      trap.rule(time[k], s[k])
    }
    Q
  }
  formals(f) <- list(lp=0, stratum=1,
                     type=if(!missing(type))type else c("step","polygon"),
                     tmax=tmax,
                     time=object$time, surv=object$surv)
##    if(!missing(tmax)) f$tmax <- tmax ??
} else {
  lp     <- object$linear.predictors
  lp.seq <- if(length(lp)) lp.seq <- seq(min(lp), max(lp), length=n) else 0
  
  time   <- object$time
  surv   <- object$surv
  nstrat <- if(is.list(time)) length(time) else 1
  areas  <- list()

  for(is in 1:nstrat) {
    tim <- if(nstrat==1) time else time[[is]]
    srv <- if(nstrat==1) surv else surv[[is]]
    if(!length(tmax)) {
      if(min(srv)>1e-3)
        warning(paste("Computing mean when survival curve only defined down to",
                      format(min(srv)),"\n Mean is only a lower limit"))
      k <- rep(TRUE,length(tim))
    } else {
      if(tmax>max(tim)) stop(paste("tmax=",format(tmax),
                                   "> max follow-up time=",format(max(tim))))
      k <- (1:length(tim))[tim<=tmax]
    }
    ymean <- lp.seq
    for(j in 1:length(lp.seq)) {
      s <- srv^exp(lp.seq[j])
      ymean[j] <- if(type=="step") sum(c(diff(tim[k]),0) * s[k]) else 
      trap.rule(tim[k], s[k])
    }
    if(!.R.) storage.mode(ymean) <- "single"
    areas[[is]] <- ymean
  }
  if(nstrat>1) names(areas) <- names(time)

  f <- function(lp=0, stratum=1, lp.seq, areas) {

    if(length(stratum)>1) stop("does not handle vector stratum")
    area <- areas[[stratum]]
    if(length(lp.seq)==1 && all(lp==lp.seq)) ymean <- rep(area,length(lp)) else
    ymean <- approx(lp.seq, area, xout=lp)$y
    if(any(is.na(ymean)))
      warning("means requested for linear predictor values outside range of linear\npredictor values in original fit")
    names(ymean) <- names(lp)
    ymean
  }
if(!.R.) storage.mode(lp.seq) <- "single"
formals(f) <- list(lp=0, stratum=1, lp.seq=lp.seq, areas=areas)
}
f
}

## cox.zph demands that the fit object inherit 'coxph'  14Nov00
## The following slightly changed cox.zph also explicitly invokes
## residuals.cph
cox.zph <- function(fit, transform = "km", global = TRUE)
{
	call <- match.call()
    clas <- c(oldClass(fit), fit$fitFunction)  ##FEH
	if(!any(c('coxph','cph') %in% clas))       ##FEH
		stop("Argument must be the result of coxph or cph")
      if('coxph.null' %in% clas)               ##FEH + next 5
		stop("The are no score residuals for a Null model")
	sresid <- if('cph' %in% clas) residuals.cph(fit, "schoenfeld") else
      if('coxph.penal' %in% clas)
      residuals.coxph.penal(fit,'schoenfeld') else
      residuals.coxph(fit, 'schoenfeld') ## coxph.penal 18mar04
	varnames <- names(fit$coef)
	nvar <- length(varnames)
	ndead <- length(sresid)/nvar
	if(nvar == 1)
		times <- as.numeric(names(sresid))
	else times <- as.numeric(dimnames(sresid)[[1]])
	if(missing(transform) && attr(fit$y, "type") != "right")
		transform <- "identity"
	if(is.character(transform)) {
		tname <- transform
		ttimes <- switch(transform,
			identity = times,
			rank = rank(times),
			log = log(times),
			km = {
              if(.R. && !existsFunction('survfit.km'))
                survfit.km <-
                  getFromNamespace('survfit.km','survival')
              temp <- survfit.km(factor(rep(1, nrow(fit$
					y))), fit$y, se.fit = FALSE)
				# A nuisance to do left cont KM
				t1 <- temp$surv[temp$n.event > 0]
				t2 <- temp$n.event[temp$n.event > 0]
				km <- rep(c(1, t1), c(t2, 0))
				if(!length(attr(sresid, "strata")))
					1 - km
				else (1 - km[sort.list(sort.list(times))])
			}
			,
			stop("Unrecognized transform"))
	}
	else {
		tname <- deparse(substitute(transform))
		ttimes <- transform(times)
	}
	xx <- ttimes - mean(ttimes)
	r2 <- sresid %*% fit$var * ndead
	test <- xx %*% r2
	# time weighted col sums
	corel <- c(cor(xx, r2))
	z <- c(test^2/(diag(fit$var) * ndead * sum(xx^2)))
	Z.ph <- cbind(corel, z, 1 - pchisq(z, 1))
	if(global && nvar > 1) {
		test <- c(xx %*% sresid)
		z <- (c(test %*% fit$var %*% test) * ndead)/sum(xx^2)
		Z.ph <- rbind(Z.ph, c(NA, z, 1 - pchisq(z, ncol(sresid))))
		dimnames(Z.ph) <- list(c(varnames, "GLOBAL"), c("rho", "chisq",
			"p"))
	}
	else dimnames(Z.ph) <- list(varnames, c("rho", "chisq", "p"))
	dimnames(r2) <- list(times, names(fit$coef))
	temp <- list(table = Z.ph, x = ttimes, y = r2 + outer(rep(1, ndead),
		fit$coef), var = fit$var, call = call, transform = tname)
	oldClass(temp) <- "cox.zph"
	temp
}

predict.cph <- 
  function(object, newdata=NULL,
           type=c("lp","x","data.frame","terms","adjto","adjto.data.frame",
             "model.frame"),
           se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
           incl.non.slopes=NULL, non.slopes=NULL, kint=1,
           na.action=na.keep, expand.na=TRUE, center.terms=TRUE, ...)
  predictDesign(object, newdata, type, se.fit, conf.int, conf.type,
                incl.non.slopes, non.slopes, kint,
                na.action, expand.na, center.terms, ...)

