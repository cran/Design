# psm for SV4/R is a modification of Therneau's survreg from survival5
# (survReg in S-Plus 6)   FEH 17Apr02

# SCCS @(#)survreg.s	5.8 07/10/00
#  The newest version of survreg, that accepts penalties and strata
#

# .newSurvival. <- .R. || existsFunction('survReg')

psm <- if(.newSurvival.) 
  function(formula=formula(data),
           data=if(.R.)parent.frame() else sys.parent(),
           weights, subset, na.action=na.delete, dist='weibull', 
           init=NULL,  scale=0,
           control=if(!.R.)survReg.control() else survreg.control(),
           parms=NULL, model=FALSE, x=FALSE, y=TRUE, time.inc, ...) {

    if(.R.) {
      require('survival')
      if(!existsFunction('survreg.fit'))
        survreg.fit <- getFromNamespace('survreg.fit','survival')
    }
    call <- match.call()
    m <- match.call(expand=FALSE)
    if(dist=='extreme')
warning('Unlike earlier versions of survreg, dist="extreme" does not fit\na Weibull distribution as it uses an identity link.  To fit the Weibull\ndistribution use the default for dist or specify dist="weibull".')
    m$na.action <- na.action  ## FEH
    temp <- c("", "formula", "data", "weights", "subset", "na.action")
    m <- m[ match(temp, names(m), nomatch=0)]
    if(.R.) m$drop.unused.levels <- TRUE  ## 31jul02
    m[[1]] <- as.name("model.frame")
    special <- c("strata", "cluster")
    Terms <- if(missing(data)) terms(formula, special)
             else              terms(formula, special, data=data)
    m$formula <- Terms
    ## Start FEH
    offs <- offset <- attr(Terms, "offset")  ## offs 23nov02 moved 6dec02
    if(!.R.) survreg.distributions <- survReg.distributions
    if(.R.) {
      dul <- .Options$drop.unused.levels
      if(!length(dul) || dul) {
        on.exit(options(drop.unused.levels=dul))
        options(drop.unused.levels=FALSE)
      }
    }
    m <- Design(eval(m, if(.R.)parent.frame() else sys.parent()))
    atrx <- attributes(m)
    nact <- atrx$na.action
    Terms <- atrx$terms
    atr   <- atrx$Design
    if(length(nact$nmiss)) {
      jia <- grep('*%ia%*',names(nact$nmiss))  ## 8feb03
      if(length(jia)) nact$nmiss <- nact$nmiss[-jia]
      s <- if(length(offs)) names(nact$nmiss) !=  atrx$names[offs] else TRUE
      ## 23nov02
      names(nact$nmiss)[s] <- 
        c(as.character(formula[2]), atr$name[atr$assume.code!=9])
    }
    ## End FEH    [s] 23nov02
    
    weights <- model.extract(m, 'weights')
    Y <- model.extract(m, "response")

    ## Start FEH
    atY <- attributes(Y)
    ncy <- ncol(Y)
    maxtime <- max(Y[,-ncy])
    nnn <- c(nrow(Y),sum(Y[,ncy]))
    time.units <- attr(Y,'units')
    if(!length(time.units)) time.units <- "Day"
    if(missing(time.inc))	{
      time.inc <- switch(time.units,Day=30,Month=1,Year=1,maxtime/10)
      if(time.inc>=maxtime | maxtime/time.inc>25)
        time.inc <- max(pretty(c(0,maxtime)))/10
    }
    ## End FEH

    if (!inherits(Y, "Surv")) stop("Response must be a survival object")

    strats <- attr(Terms, "specials")$strata
    cluster<- attr(Terms, "specials")$cluster
    dropx <- NULL
    if (length(cluster)) {
        if (missing(robust)) robust <- TRUE
        tempc <- untangle.specials(Terms, 'cluster', 1:10)
        ord <- attr(Terms, 'order')[tempc$terms]
        if (any(ord>1)) stop ("Cluster can not be used in an interaction")
        cluster <- strata(m[,tempc$vars], shortlabel=TRUE)  #allow multiples
        dropx <- tempc$terms
        }
    if (length(strats)) {
        temp <- untangle.specials(Terms, 'strata', 1)
        dropx <- c(dropx, temp$terms)
        if (length(temp$vars)==1) strata.keep <- m[[temp$vars]]
        else strata.keep <- strata(m[,temp$vars], shortlabel=TRUE)
        strata <- as.numeric(strata.keep)
	nstrata <- max(strata)
        }
    else {
	nstrata <- 1
	strata <- 0
	}

    if (length(dropx)) newTerms<-Terms[-dropx]
    else               newTerms<-Terms
    X <- model.matrix(newTerms,m)

    ## Start FEH
    rnam <- dimnames(Y)[[1]]
    dimnames(X) <- list(rnam, c("(Intercept)",atr$colnames))
    ## End FEH except for 23nov02 and later changes
    
    n <- nrow(X)
    nvar <- ncol(X)

    if (length(offset)) offset <- as.numeric(m[[offset]])  ## length 23nov02
    else                offset <- rep(0, n)

    if (is.character(dist)) {
	dlist <- survreg.distributions[[dist]]
	if (is.null(dlist)) stop(paste(dist, ": distribution not found"))
	}
    else if (is.list(dist)) dlist <- dist
    else stop("Invalid distribution object")
    if (is.null(dlist$dist)) {
	if (is.character(dlist$name) && is.function(dlist$init) &&
	    is.function(dlist$deviance)) {}
	else stop("Invalid distribution object")
	}
    else {
	if (!is.character(dlist$name) || !is.character(dlist$dist) ||
	    !is.function(dlist$trans) || !is.function(dlist$dtrans))
		stop("Invalid distribution object")
	}	

    type <- attr(Y, "type")
    if (type== 'counting') stop ("Invalid survival type")
    
    logcorrect <- 0   #correction to the loglik due to transformations
    if (!is.null(dlist$trans)) {
	tranfun <- dlist$trans
	exactsurv <- Y[,ncol(Y)] ==1
	if (any(exactsurv)) logcorrect <- sum(logb(dlist$dtrans(Y[exactsurv,1])))

	if (type=='interval') {
	    if (any(Y[,3]==3))
		    Y <- cbind(tranfun(Y[,1:2]), Y[,3])
	    else Y <- cbind(tranfun(Y[,1]), Y[,3])
	    }
	else if (type=='left')
	     Y <- cbind(tranfun(Y[,1]), 2-Y[,2])
	else     Y <- cbind(tranfun(Y[,1]), Y[,2])
	if (!all(is.finite(Y))) 
	    stop("Invalid survival times for this distribution")
	}
    else {
	if (type=='left') Y[,2] <- 2- Y[,2]
	else if (type=='interval' && all(Y[,3]<3)) Y <- Y[,c(1,3)]
	}

    if (is.null(dlist$itrans)) itrans <- function(x) x
    else itrans <- dlist$itrans

    if (!is.null(dlist$scale)) {
	if (!missing(scale)) warning(paste(dlist$name, 
			   "has a fixed scale, user specified value ignored"))
	scale <- dlist$scale
	}
    if (!is.null(dlist$dist)) dlist <- survreg.distributions[[dlist$dist]]

    if (missing(control)) control <- if(!.R.)
      survReg.control(...) else survreg.control(...)

    if (scale < 0) stop("Invalid scale value")
    if (scale >0 && nstrata >1) 
	    stop("Cannot have multiple strata with a fixed scale")

    # Check for penalized terms
    pterms <- sapply(m, inherits, 'coxph.penalty')
    if (any(pterms)) {
	pattr <- lapply(m[pterms], attributes)
	# 
	# the 'order' attribute has the same components as 'term.labels'
	#   pterms always has 1 more (response), sometimes 2 (offset)
	# drop the extra parts from pterms
	temp <- c(attr(Terms, 'response'), attr(Terms, 'offset'))
	if (length(dropx)) temp <- c(temp, dropx+1)
	pterms <- pterms[-temp]
	temp <- match((names(pterms))[pterms], attr(Terms, 'term.labels'))
	ord <- attr(Terms, 'order')[temp]
	if (any(ord>1)) stop ('Penalty terms cannot be in an interaction')
	##pcols <- (attr(X, 'assign')[-1])[pterms]
        assign<-attrassign(X,newTerms)
        pcols<-assign[-1][pterms]
  
        fit <- survpenal.fit(X, Y, weights, offset, init=init,
				controlvals = control,
			        dist= dlist, scale=scale,
			        strata=strata, nstrat=nstrata,
				pcols, pattr,assign, parms=parms)
	}
    else fit <- if(!.R.) survReg.fit(X, Y, weights, offset, 
			    init=init, controlvals=control,
			    dist= dlist, scale=scale, nstrat=nstrata, 
			    strata, parms=parms) else
      survreg.fit(X, Y, weights, offset, 
			    init=init, controlvals=control,
			    dist= dlist, scale=scale, nstrat=nstrata, 
			    strata, parms=parms)

    # Next line: FEH added fitFunction='psm'
    if (is.character(fit))
      fit <- list(fail=fit, fitFunction='psm')  #error message
    else {
	if (scale==0) {
	    nvar <- length(fit$coef) - nstrata
	    fit$scale <- exp(fit$coef[-(1:nvar)])
	    if (nstrata==1) names(fit$scale) <- NULL
	    else names(fit$scale) <- levels(strata.keep)
	    fit$coefficients  <- fit$coefficients[1:nvar]
	    fit$idf  <- 1 + nstrata
	    }
	else {
	    fit$scale <- scale
	    fit$idf  <- 1
	    }
	fit$loglik <- fit$loglik + logcorrect
	}

    ## na.action <- attr(m, "na.action")
    ##if (length(na.action)) fit$na.action <- na.action
    if(length(nact)) fit$na.action <- nact  ## FEH
    fit$df.residual <- n - sum(fit$df)
#   fit$fitted.values <- itrans(fit$linear.predictors)
    fit$terms <- Terms
    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$means <- apply(X,2, mean)
    fit$call <- call
    fit$dist <- dist
    fit$df.resid <- n-sum(fit$df) ##used for anova.survreg
    if (model) fit$model <- m
    if (x)     fit$x <- X
##    if (y)     fit$y <- Y   FEH
    if (length(parms)) fit$parms <- parms
    ## Start FEH
    ##    if (any(pterms)) class(fit)<- c('survreg.penal', 'survreg')
    ##    else	     class(fit) <- 'survreg'
    fit$assign <- DesignAssign(atr, 1, Terms)
    fit$formula <- formula
    if(y) {
      oldClass(Y) <- 'Surv'
      attr(Y,'type') <- atY$type
      fit$y <- Y
    }
    if(.newSurvival.) scale.pred <-
      if(dist %in% c('weibull','exponential','lognormal','loglogistic'))
      c('log(T)','Survival Time Ratio') else 'T' else
    scale.pred <- if(substring(dist,1,3)=='log')
      c("log(T)","Survival Time Ratio") else "T"

    logtest <- 2*diff(fit$loglik)
    Nn <- if(length(weights)) sum(weights) else nnn[1]  ## 5jun02
    R2.max <- 1 - exp(2*fit$loglik[1]/Nn)
    R2 <- (1 - exp(-logtest/Nn))/R2.max
    df <- length(fit$coef)-1
    P <- if(df==0) NA else 1-pchisq(logtest,df)
    stats <- c(nnn, logtest, df, P, R2)
    names(stats) <- c("Obs", "Events", "Model L.R.", "d.f.", "P",
                      "R2")
    if(length(weights)) stats <- c(stats, 'Sum of Weights'=sum(weights))
    fit <- c(fit, list(stats=stats, maxtime=maxtime, units=time.units,
                       time.inc=time.inc, scale.pred=scale.pred,
                       non.slopes=1, Design=atr, fail=FALSE,
                       fitFunction=c("psm", "survreg", "glm", "lm")))
    if (any(pterms)) oldClass(fit) <- if(.SV4.)'Design' else
                        c('psm','Design','survreg.penal','survreg')
    else oldClass(fit) <- if(.SV4.)'Design' else c('psm','Design','survreg')
    ## End FEH
                       
    fit
    
  } else function(formula=formula(data), data,
	subset, na.action=na.delete, method="fit",
	link="log",
	dist=c("extreme", "logistic", "gaussian", "exponential","rayleigh","t"),
	init=NULL, fixed=list(), control,
	model=FALSE, x=FALSE, y=FALSE, time.inc, ...) {

    call <- match.call()
    dist <- match.arg(dist)
    
    m <- match.call(expand=FALSE)
    m$dist <- m$link <- m$model <- m$x <- m$y <- m$... <-  NULL
    m$init <- m$fixed <- m$control <- m$time.inc <- NULL
    m$na.action <- na.action
    m[[1]] <- as.name("model.frame")
    X <- Design(eval(m, sys.parent()))    # 24Apr01
    atrx <- attributes(X)                 # 24Apr01 + next 4
    nact <- atrx$na.action
    if(method=="model.frame") return(X)
    Terms <- atrx$terms
    offs <- offset <- attr(Terms, "offset")  ## offs 23nov02 moved 6dec02
    atr   <- atrx$Design
    s <- if(length(offs)) names(nact$nmiss) !=  atrx$names[offs] else TRUE
    ## 23nov02
    if(length(nact$nmiss))
      names(nact$nmiss)[s] <- 
        c(as.character(formula[2]), atr$name[atr$assume.code!=9])
    ## [s] 23nov02
    lnames <- if(.R.) c("logit","probit","cloglog","identity","log","sqrt",
      "1/mu^2","inverse") else dimnames(glm.links)[[2]]
    link <- pmatch(link, lnames, 0)
    if(link==0) stop("invalid link function")
    link <- lnames[link]
    Y <- model.extract(X, "response")
    atY <- attributes(Y)    # 1 Apr 95
    ncy <- ncol(Y)
    maxtime <- max(Y[,-ncy])
    nnn <- c(nrow(Y),sum(Y[,ncy]))
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    if(model) m <- X

    if(length(offset)) offset <- as.numeric(X[[offset]])
    else offset <- rep(0, nrow(X))

    X <- model.matrix(Terms, X)

    time.units <- attr(Y, "units")
    if(!length(time.units)) time.units <- "Day"
    if(missing(time.inc))	{
      time.inc <- switch(time.units,Day=30,Month=1,Year=1,maxtime/10)
      if(time.inc>=maxtime | maxtime/time.inc>25) time.inc <- max(pretty(c(0,maxtime)))/10
      				}
    rnam <- dimnames(Y)[[1]]
    dimnames(X) <- list(rnam, c("(Intercept)",atr$colnames))
    if(method=="model.matrix") return(X)
    n <- nrow(X)
    nvar <- ncol(X)

    type <- attr(Y, "type")
    linkfun <- if(.R.) make.link(link)$linkfun else
                       glm.links["link", link][[1]]
    if (type== 'counting') stop ("Invalid survival type")
    else if (type=='interval') {
	if (any(Y[,3]==3))
	     Y <- cbind(linkfun(Y[,1:2]), Y[,3])
	else Y <- cbind(linkfun(Y[,1]), Y[,3])
	}
    else if (type=='left')
	     Y <- cbind(linkfun(Y[,1]), 2-Y[,2])
    else     Y <- cbind(linkfun(Y[,1]), Y[,2])

    controlvals <- survreg.control()
    if(!missing(control)) controlvals[names(control)] <- control

    if(dist=="exponential") { fixed$scale <- 1; dist <- "extreme" }
    else if (dist=="rayleigh") { fixed$scale <- .5; dist <- "extreme" }

    sd <- survreg.distributions[[dist]]
    if (length(fixed)>0) {
	ifix <- match(names(fixed), names(sd$parms), nomatch=0)
	if (any(ifix==0))
	    stop (paste("Parameter(s)", paste(names(fixed)[ifix==0]),
			"in the fixed list not valid for this dist"))
	}
    if (is.list(init) && length(init)>0) {
	ifix <- match(names(init), c('eta',names(sd$parms)), nomatch=0)
	if (any(ifix==0))
	    stop (paste("Parameter(s)", paste(names(init)[ifix==0]),
			"in the init list not valid for this dist"))
	}

    sfit <- if(.R.) survreg.fit(X, Y, offset=offset, init=init,
                                controlvals=controlvals,
                                dist=dist, parms=fixed) else
    survreg.fit(X, Y, offset, init=init, controlvals=controlvals,
                dist= dist, fixed=fixed)

    if (is.character(sfit))  	{
	cat("Failure in psm:\n",sfit,"\n")
	fit <- list(fail=TRUE, fitFunction='psm')  #error message; 14Nov00
	oldClass(fit) <- if(.SV4.)'Design' else "psm" ##14Nov00
	return(fit)
				}
    else {
	# There may be more clever ways to do this, but ....
	#  In order to make it look like IRLS, and get all the components
	#  that I need for glm inheritance, do one step of weighted least
	#  squares.
	eta <- c(X %*% sfit$coef[1:nvar]) + offset
	wt<- -sfit$deriv[,3]
	fit <- lm.wfit(X, eta + sfit$deriv[,2]/wt, wt, "qr", ...)

    ifun <- if(.R.) make.link(link)$linkinv else
                    glm.links["inverse",link][[1]]
    fit$fitted.values <- ifun(fit$fitted.values)
	fit$family <- c(name=dist, link=link, "")
	fit$linear.predictors <- eta
	fit$iter <- sfit$iter
    fit$parms <- sfit$parms
    fit$df.residual <- fit$df.residual-sum(!sfit$fixed)

	# If singular.ok=T, there may be NA coefs.  The var matrix should
	#   be an inversion of the "non NA" portion.

	var <- 0*sfit$imat
	good <- c(!is.na(fit$coef), rep(TRUE, ncol(var)-nvar))
	var[good,good] <- solve(qr(sfit$imat[good,good], tol=1e-12))
	fit$var <- var
	fit$fixed <- sfit$fixed
	dev <- sd$deviance(Y, fit$parms, sfit$deriv[,1])
	fit$dresiduals <- sign(fit$residuals)*sqrt(dev)
	fit$deviance <- sum(dev)
	fit$null.deviance <- fit$deviance +2*(sfit$loglik[2]- sfit$ndev[2])
	fit$loglik <- c(sfit$ndev[2], sfit$loglik[2])
	}

    if (length(nact)) fit$na.action <- nact

    i <- 1:nvar
    var <- var[i,i,drop=FALSE]	#omit scale row and col.
    fit$se.fit <-  drop(sqrt(((X %*% var) * X) %*% rep(1,nvar)))

    logtest <- fit$null.deviance - fit$deviance
    R2.max <- 1 - exp(-fit$null.deviance/nnn[1])
    R2 <- (1 - exp(-logtest/nnn[1]))/R2.max
    df <- length(fit$coef)-1
    P <- if(df==0) NA else 1-pchisq(logtest,df)
    stats <- c(nnn, logtest, df, P, R2)
    names(stats) <- c("Obs", "Events", "Model L.R.", "d.f.", "P", "R2")
    scale.pred <- if(link=="log") c("log(T)","Survival Time Ratio") else "T"
    fit <- c(fit, list(maxtime=maxtime, units=time.units,
                       time.inc=time.inc,scale.pred=scale.pred,
                       non.slopes=1,
                       fitFunction=c("psm", "survreg", "glm", "lm"))) ##13Nov00
    fit$Design <- atr
    fit$stats  <- stats
    oldClass(fit) <- if(.SV4.)'Design' else
     c("psm", "Design", "survreg", "glm", "lm")  ##14Nov00
    fit$terms <- Terms
    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$call <- call
    fit$fail <- FALSE
    if (model) fit$model <- m
    if (x)     fit$x <- X
    if (y) {
      oldClass(Y) <- 'Surv'    # 1 Apr 95
      attr(Y,'type') <- atY$type
      fit$y <- Y
    }
    fit
    }

Hazard   <- function(object, ...) UseMethod("Hazard")
Survival <- function(object, ...) UseMethod("Survival")

Hazard.psm <- if(.newSurvival.) function(object, ...) {
 dist <- object$dist
 g <- survreg.auxinfo[[dist]]$hazard
 formals(g) <- list(times=NA, lp=NULL, parms=logb(object$scale))
 g
} else function(object) {
  fam <- object$family
  dist <- fam["name"]
  transform <- fam[2]
  g <- survreg.auxinfo[[dist]]$hazard
  formals(g) <- list(times=NULL, lp=NULL,
                     parms=object$parms, transform=transform)
  g
}

Survival.psm <- if(.newSurvival.) function(object, ...) {
 dist <- object$dist
 g <- survreg.auxinfo[[dist]]$survival
 formals(g) <- list(times=NULL, lp=NULL, parms=logb(object$scale))
 g
} else function(object) {
  fam <- object$family
  dist <- fam["name"]
  transform <- fam[2]
  g <- survreg.auxinfo[[dist]]$survival
  formals(g) <- list(times=NULL, lp=NULL,
                     parms=object$parms, transform=transform)
  g
}

Quantile.psm <- if(.newSurvival.) function(object, ...) {
  dist <- object$dist
  g <- survreg.auxinfo[[dist]]$Quantile
  formals(g) <- list(q=.5, lp=NULL, parms=logb(object$scale))
  g
} else function(object, ...) {
  fam <- object$family
  dist <- fam["name"]
  transform <- fam[2]
  g <- survreg.auxinfo[[dist]]$quantile
  formals(g) <- list(q=.5, lp=NULL,
                     parms=object$parms, transform=transform)
  g
}

Mean.psm <- if(.newSurvival.) function(object, ...) {
 dist <- object$dist
 g <- survreg.auxinfo[[dist]]$mean
 formals(g) <- list(lp=NULL, parms=logb(object$scale))
 g
} else function(object, ...) {
  fam <- object$family
  dist <- fam["name"]
  transform <- fam[2]
  g <- survreg.auxinfo[[dist]]$mean
  formals(g) <- list(lp=NULL, parms=object$parms,
                     transform=transform)
  g
}

predict.psm <- 
  function(object, newdata,
           type=c("lp","x","data.frame","terms","adjto","adjto.data.frame",
             "model.frame"),
           se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
           incl.non.slopes, non.slopes, kint=1,
           na.action=na.keep, expand.na=TRUE, center.terms=TRUE, ...)
  predictDesign(object, newdata, type, se.fit, conf.int, conf.type,
                incl.non.slopes, non.slopes, kint,
                na.action, expand.na, center.terms, ...)



residuals.psm <- function(object, type = "censored.normalized", ...)
{
    type <- match.arg(type)
    if(type!='censored.normalized') {
      if(type=='score' && (.newSurvival.))
        stop('score residuals not implemented')
      ## TODO
      return(if(!.R.)residuals.survReg(object, type=type) else
        residuals.survreg(object, type=type))
    }

    y <- object$y
    aty <- attributes(y)
    if(length(y)==0) stop('did not use y=T with fit')
    ncy <- ncol(y)
    if(.newSurvival.) {  ## 17Apr02
      scale <- object$scale
      dist  <- object$dist
    } else {
      scale <- exp(object$parms)
      dist  <- object$family[1]
    }
    r <- (y[,-ncy,drop=FALSE]-object$linear.predictors)/scale
    r <- cbind(r, y[,ncy])
    ## Moved the following line here from bottom
    if(length(object$na.action)) r <- naresid(object$na.action, r)
    attr(r,'dist') <- dist
    attr(r,'type') <- aty$type
    attr(r,'units') <- ' '
    attr(r,'time.label') <- 'Normalized Residual'
    attr(r,'event.label') <- aty$event.label
    oldClass(r) <- c('residuals.psm.censored.normalized','Surv')
    g <- survreg.auxinfo[[dist]]$survival
    formals(g) <- if(.newSurvival.) list(times=NULL, lp=0, parms=0)
    else list(times=NULL, lp=0, parms=0, transform='identify')  ## 17Apr02
    attr(r,'theoretical') <- g
    r
    }

lines.residuals.psm.censored.normalized <- 
  function(x, n=100, lty=1, xlim=range(r[,-ncol(r)],na.rm=TRUE),
           lwd=3, ...) {
    r <- x
	x <- seq(xlim[1], xlim[2], length=n)
    tx <- x
    if(.newSurvival.) {
      dist <- attr(r, 'dist')
      if(dist %in% c('weibull','loglogistic','lognormal')) tx <- exp(x)
      ## $survival functions log x
    }
	lines(x, attr(r,'theoretical')(tx), lwd=lwd, lty=lty)
	invisible()
  }

survplot.residuals.psm.censored.normalized <- 
  function(fit, x, g=4, col, main, ...) {
    r <- fit
    
  if(missing(x)) {
	survplot(survfit(r), conf='none', xlab='Residual', 
			 col=if(missing(col))par('col') else col, ...)
    ## was survfit(r, data=list(r)) 28apr02
	if(!missing(main)) title(main)
  } else {
	if(is.character(x)) x <- as.factor(x)
	if(!is.category(x) && length(unique(x))>5) x <- cut2(x, g=g)
	s <- is.na(r[,1]) | is.na(x)
	if(any(s)) {r <- r[!s,]; x <- x[!s,drop=TRUE]}
	survplot(survfit(r ~ x, data=data.frame(x,r)),  xlab='Residual',
			 conf='none',
			 col=if(missing(col))1:length(levels(x)) else par('col'), ...)
	if(missing(main)) main <- if(length(lab <- attr(x,'label'))) lab 
	  else if(.R.) '' else deparse(substitute(x))
	if(main != '') title(main)
  }
  lines(r, lty=1, lwd=3)
  invisible()
}



