#SCCS @(#)residuals.coxph.s	4.5 7/14/92
residuals.cph <-
  function(object, type=c("martingale", "deviance", "score", "schoenfeld",
	"dfbeta","dfbetas", "scaledsch"), collapse=FALSE, weighted=FALSE, ...)
    {
    type <- match.arg(type)
    if(.R.) require('survival')
    otype <- type
    if(type=="dfbeta" | type=="dfbetas") type <- "score"
    if(type=="scaledsch") type <- "schoenfeld"
    n <- length(object$residuals)
    rr <- object$residual
    y <- object$y
    x <- object$x
	weights <- object$weights
    strat <- attr(x,"strata")
    method <- object$method
    if (method=='exact' && (type=='score' || type=='schoenfeld'))
	stop(paste(type, 'residuals are not available for the exact method'))
    if(type!="martingale" & is.null(x))
      stop("you must specify x=T in the fit")
    if(type!="deviance" & type!="martingale" & is.null(y))
      stop("you must specify y=T in the fit")

    if(type!="martingale")				{

	ny <- ncol(y)
	status <- y[,ny,drop=TRUE]

	if (type != 'deviance') {
	    nvar <- ncol(x)
	    if (is.null(strat)) {
		ord <- order(y[,ny-1], -status)
		newstrat <- rep(0,n)
		}
	    else {
		ord <- order(strat, y[,ny-1], -status)
		newstrat <- c(diff(as.numeric(strat[ord]))!=0 ,1)
		}
	    newstrat[n] <- 1

	    # sort the data
	    x <- x[ord,]
	    y <- y[ord,]
	    score <- exp(object$linear.predictor)[ord]
		weights <- if(length(weights)) weights[ord] else rep(1,n)

	    if (ny ==3) subs <- paste("agres", 1:2, sep='')
	    else        subs <- paste("coxres",1:2, sep='')
	    }
	}

    #
    # Now I have gotten the data that I need-- do the work
    #
    if (type=='schoenfeld') {
	if (ny==2)  {
      mintime <- min(y[, 1])
      if(mintime < 0) y <- cbind(2 * mintime - 1, y)
      else y <- cbind(-1, y)
    }
	temp <- if(.R.)
      .C("coxscho",
         n=as.integer(n),
         as.integer(nvar),
         as.double(y),
         resid= x,
         score * weights,
         as.integer(newstrat),
         as.integer(method=='efron'),
         double(3*nvar), PACKAGE="survival") else
      .C(if(.SV4.)'S_coxscho' else "coxscho",  ##14Nov00
         n=as.integer(n),
         as.integer(nvar),
         as.double(y),
         resid= x,
         score * weights,
         as.integer(newstrat),
         as.integer(method=='efron'),
         double(3*nvar))

	deaths <- y[,3]==1

	if (nvar==1) rr <- temp$resid[deaths]
	else rr <- matrix(temp$resid[deaths,], ncol=nvar) #pick rows, and kill attr

	if (length(object$strata)) 
	  attr(rr, "strata")  <- table((strat[ord])[deaths])
	time <- c(y[deaths,2])  # 'c' kills all of the attributes
	if (is.matrix(rr)) dimnames(rr)<- list(time, names(object$coef))
	else               names(rr) <- time

	if (otype=='scaledsch') {
	  ndead <- sum(deaths)
	  coef <- ifelse(is.na(object$coef), 0, object$coef)
	    if (nvar==1) rr <- rr*object$var * ndead + coef
	    else         rr <- rr %*% object$var * ndead + 
		  outer(rep(1,nrow(rr)), coef)
	    }
	return(rr)
	}

    if (type=='score') {
	if (ny==2) {
	    resid <- if(.R.)
          .C("coxscore",
             as.integer(n),
             as.integer(nvar),
             as.double(y),
             x=as.double(x),
             as.integer(newstrat),
             as.double(score),
             as.double(weights),
             as.integer(method=='efron'),
             resid= double(n*nvar),
             double(2*nvar), PACKAGE="survival")$resid else
        .C(if(.SV4.)'S_coxscore' else "coxscore",
           as.integer(n),
           as.integer(nvar),
           as.double(y),
           x=as.double(x),
           as.integer(newstrat),
           as.double(score),
           as.double(weights),
           as.integer(method=='efron'),
           resid= double(n*nvar),
           double(2*nvar))$resid
	    }
	else {
	    resid <- if(.R.)
          .C("agscore",
             as.integer(n),
             as.integer(nvar),
             as.double(y),
             as.double(x),
             as.integer(newstrat),
             as.double(score),
             as.double(weights),
             as.integer(method=='efron'),
             resid=double(n*nvar),
             double(nvar*6), PACKAGE="survival")$resid
        .C(if(.SV4.)'S_agscore' else "agscore",
           as.integer(n),
           as.integer(nvar),
           as.double(y),
           as.double(x),
           as.integer(newstrat),
           as.double(score),
           as.double(weights),
           as.integer(method=='efron'),
           resid=double(n*nvar),
           double(nvar*6))$resid
      }
	if (nvar >1) {
	    rr <- matrix(0, n, nvar)
	    rr[ord,] <- matrix(resid, ncol=nvar)
	    dimnames(rr) <- list(names(object$resid), names(object$coef))
	    }
	else rr[ord] <- resid
	}

    #Expand out the missing values in the result
    if (!is.null(object$na.action)) {
	rr <- naresid(object$na.action, rr)
	if (is.matrix(rr)) n <- nrow(rr)
	else               n <- length(rr)
	if (type=='deviance') status <- naresid(object$na.action, status)
	}

    # Collapse if desired
    if (!missing(collapse)) {
	if (length(collapse) !=n) stop("Wrong length for 'collapse'")
	rr <- rowsum(rr, collapse)
	}

    # Deviance residuals are computed after collapsing occurs
    if (type=='deviance')
	rr <- sign(rr) *sqrt(-2* (rr+
			      ifelse(status==0, 0, status*logb(status-rr))))

    if      (otype=='dfbeta') rr %*% object$var
    else if (otype=='dfbetas') (rr %*% object$var) %*% 
		diag(sqrt(1/diag(object$var)))
    else  rr
    }
