if(.newSurvival.) {
validate.psm <- function(fit,method="boot",B=40,
		   bw=FALSE,rule="aic",type="residual",sls=.05,aics=0,pr=FALSE,
		   dxy=FALSE,tol=1e-12, rel.tolerance=1e-5, maxiter=15, ...) {

xb <- fit$linear.predictors
ny <- dim(fit$y)
nevents <- sum(fit$y[,ny[2]])

#Note: fit$y already has been transformed by the link function by psm

dist  <- fit$dist
scale <- fit$scale
parms <- fit$parms

distance <- function(x,y,fit,iter,evalfit=FALSE,fit.orig,dxy=FALSE,dist,parms,
	tol=1e-12, maxiter=15, rel.tolerance=1e-5, ...){
	#Assumes y is matrix with 1st col=time, 2nd=event indicator
if(evalfit)		{	#Fit was for training sample
  lr <- 2*diff(fit$loglik)
  ll0 <- -2*fit$loglik[1]
  R2.max <- 1 - exp(-ll0/length(x))
  R2 <- (1 - exp(-lr/length(x)))/R2.max
  intercept <- 0
  slope <- 1
  D <- (lr-1)/ll0
  U <- -2/ll0		}
else			{
  f <- survreg.fit2(cbind(1,x),y,iter=iter,dist=dist,parms=parms,tol=tol,
	maxiter=maxiter,rel.tolerance=rel.tolerance)
  if(f$fail) stop("survreg.fit2 failed in distance")
  lr <- 2*diff(f$loglik)
  ll0 <- -2*f$loglik[1]
  R2.max <- 1 - exp(-ll0/length(x))
  R2 <- (1 - exp(-lr/length(x)))/R2.max
  intercept <- f$coefficients[1]
  slope <- f$coefficients[2]
  D <- (lr-1)/ll0
  init <- c(0, 1, if(length(f$scale)) log(f$scale) else NULL)
  f.frozen <- survreg.fit2(cbind(1,x),y,dist=dist,parms=parms,tol=tol,
        maxiter=0, init=init)
  if(f.frozen$fail) stop('survreg.fit2 failed for frozen coefficient re-fit')
  ll0 <- -2*f.frozen$loglik[1]
  frozen.lr <- 2*diff(f.frozen$loglik)
  U <- (frozen.lr - lr)/ll0
}

Q <- D - U
z <- c(R2, intercept, slope, D, U, Q)
nam <- c("R2", "Intercept", "Slope", "D", "U", "Q")
if(dxy)		{
  Dxy <- rcorr.cens(x,y)["Dxy"]; z <- c(Dxy, z); nam <- c("Dxy", nam)   }
names(z) <- nam
z
}

predab.resample(fit, method=method,
		fit=survreg.fit2, measure=distance,
		pr=pr, B=B, bw=bw, rule=rule, type=type,  
		dxy=dxy,
		dist=dist, inverse=inverse, parms=parms,
		sls=sls, aics=aics, strata=FALSE, tol=tol, maxiter=maxiter, 
		rel.tolerance=rel.tolerance, ...)
}


survreg.fit2 <- function(x,y,iter=0,dist,parms=NULL,tol,maxiter=15, 
	init=NULL, rel.tolerance=1e-5, fixed=NULL, ...) {
	e <- y[,2]
	if(sum(e)<5)return(list(fail=TRUE))
	x <- x	#Get around lazy evaluation creating complex expression
    if(!.R.) survreg.distributions <- survReg.distributions
    dlist <- survreg.distributions[[dist]]
    logcorrect <- 0
    if (length(dlist$trans)) {
      exactsurv <- y[,ncol(y)] ==1
      if(any(exactsurv)) {
        ytrans <- if(length(dlist$itrans)) dlist$itrans(y[exactsurv,1]) else
                 y[exactsurv,1]
        logcorrect <- sum(logb(dlist$dtrans(ytrans)))
      }
    }
    if (length(dlist$dist)) dlist <- survreg.distributions[[dlist$dist]]

	f <- if(!.R.)
      survReg.fit(as.matrix(x),y,dist=dlist,parms=parms,
                  controlvals=survReg.control(maxiter=maxiter,
                    rel.tolerance=rel.tolerance,
                    failure=2),
                  offset=rep(0,length(e)),init=init) else {
                    if(.R. && !existsFunction('survreg.fit'))
                      survreg.fit <- getFromNamespace('survreg.fit','survival')
    survreg.fit(as.matrix(x),y,dist=dlist,parms=parms,
                controlvals=survreg.control(maxiter=maxiter,
                  rel.tolerance=rel.tolerance,
                  failure=2),
                offset=rep(0,length(e)),init=init)
                  }
	if(is.character(f)) { warning(f); return(list(fail=TRUE)) }
	f$fail <- FALSE
    
    # TODO: fetch scale properly if fixed
    nstrata <- length(f$icoef) - 1
    if (nstrata > 0) {
      nvar <- length(f$coef) - nstrata
      f$scale <- exp(f$coef[-(1:nvar)])
      names(f$scale) <- NULL  # get rid of log( ) in names
      f$coefficients  <- f$coefficients[1:nvar]
    }
	else f$scale <- scale
    f$loglik <- f$loglik + logcorrect
    
#	f$var <- solvet(f$imat, tol=tol)
#	sd <- survreg.distributions[[dist]]
#	f$deviance <- sum(sd$deviance(y,f$parms, f$deriv[,1]))
#	f$null.deviance <- f$deviance + 2*(f$loglik[2] - f$ndev[2])
#	f$loglik <- c(f$ndev[2], f$loglik[2])
    f
  }

NULL
} else {

  validate.psm <- function(fit,method="boot",dxy=FALSE,B=40,
		   bw=FALSE,rule="aic",type="residual",sls=.05,aics=0,pr=FALSE,
		   tol=1e-12, rel.tolerance=1e-5, maxiter=15, ...) {

## f$coef -> f$coefficients 14Aug01
#call <- match.call()
xb <- fit$linear.predictors
ny <- dim(fit$y)
nevents <- sum(fit$y[,ny[2]])

#Note: fit$y already has been transformed by the link function by psm

link <- fit$family["link"]
inverse <- glm.links["inverse",link][[1]]
family <- fit$family
fixed <- fit$fixed
fixed <- if(length(fixed)==1 && is.logical(fixed) && !fixed) list() else
   list(scale=TRUE)
dist <- fit$family["name"]


distance <- function(x,y,fit,iter,evalfit=FALSE,fit.orig,dxy=FALSE,dist,fixed,
	family, tol=1e-12, maxiter=15, rel.tolerance=1e-5, ...){
	#Assumes y is matrix with 1st col=time, 2nd=event indicator
if(evalfit)		{	#Fit was for training sample
  ll0 <- fit$null.deviance
  lr <- ll0-fit$deviance
  R2.max <- 1 - exp(-ll0/length(x))
  R2 <- (1 - exp(-lr/length(x)))/R2.max
  intercept <- 0
  slope <- 1
  D <- (lr-1)/ll0
  U <- -2/ll0		}
else			{
  f <- survreg.fit2(cbind(1,x),y,iter=iter,dist=dist,fixed=fixed,tol=tol,
	family=family,maxiter=maxiter,rel.tolerance=rel.tolerance)
  if(f$fail) stop("survreg.fit2 failed in distance")
  ll0 <- f$null.deviance
  lr <- ll0-f$deviance
  R2.max <- 1 - exp(-ll0/length(x))
  R2 <- (1 - exp(-lr/length(x)))/R2.max
  intercept <- f$coefficients[1]
  slope <- f$coefficients[2]
  D <- (lr-1)/ll0
  init <- if(length(fixed)==0 || (is.logical(fixed) && !fixed)) 
	c(0,1,f$parms) else c(0, 1)

  f.frozen <- survreg.fit2(cbind(1,x),y,dist=dist,fixed=fixed,tol=tol,
        maxiter=0, init=init,family=family)
  if(f.frozen$fail) stop('survreg.fit2 failed for frozen coefficient re-fit')
  ll0 <- f.frozen$null.deviance
  U <- (f.frozen$deviance-f$deviance)/ll0
			}

Q <- D - U
z <- c(R2, intercept, slope, D, U, Q)
nam <- c("R2", "Intercept", "Slope", "D", "U", "Q")
if(dxy)		{
  Dxy <- rcorr.cens(x,y)["Dxy"]; z <- c(Dxy, z); nam <- c("Dxy", nam)   }
names(z) <- nam
z								}

predab.resample(fit, method=method,
		fit=survreg.fit2, measure=distance,
		pr=pr, B=B, bw=bw, rule=rule, type=type,  
		dxy=dxy,
		dist=dist, inverse=inverse, fixed=fixed, family=family,
		sls=sls, aics=aics, strata=FALSE, tol=tol, maxiter=maxiter, 
		rel.tolerance=rel.tolerance, ...)
								}


survreg.fit2 <- function(x,y,iter=0,dist,fixed=NULL,family,tol,maxiter=15, 
	init=NULL, rel.tolerance=1e-5, parms=NULL, ...) {
	e <- y[,2]
	if(sum(e)<5)return(list(fail=TRUE))
	x <- x	#Get around lazy evaluation creating complex expression
    if(.R. && !existsFunction('survreg.fit'))
      survreg.fit <- getFromNamespace('survreg.fit','survival')
    
	f <- survreg.fit(as.matrix(x),y,dist=dist,
		fixed=fixed,
		controlvals=survreg.control(maxiter=maxiter,
		  rel.tolerance=rel.tolerance,
		  failure=2),
		offset=rep(0,length(e)),init=init)
	if(is.character(f)) { warning(f); return(list(fail=TRUE)) }
    if(!length(f$coefficients)) {
      f$coefficients <- f$coef # 14Aug01
      f$coef <- NULL
    }
	if(length(fixed)==0 || !fixed$scale)
      f$coefficients <- f$coefficients[-length(f$coefficients)]
	f$family <- family
	f$fail <- FALSE
	f$var <- solvet(f$imat, tol=tol)
	sd <- survreg.distributions[[dist]]
	f$deviance <- sum(sd$deviance(y,f$parms, f$deriv[,1]))
	f$null.deviance <- f$deviance + 2*(f$loglik[2] - f$ndev[2])
	f$loglik <- c(f$ndev[2], f$loglik[2])
	f
  }
NULL
}

