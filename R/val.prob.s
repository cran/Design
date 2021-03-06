#Compute various measures of reliability and discrimination for a set
#of predicted probabilities p or predicted logits logit.
#If pl=T, the following apply:
#  Plots reliability curve, for which xlab is optional label.
#  If smooth=T and pl=T, plots lowess(p,y,iter=0)
#  lim is x-axis and y-axis range, default=c(0,1)
#  If m or g is specified, also computes and plots proportions of y=1
#  by quantile groups of p (or 1/(1+exp(-logit))).  If m is given,
#  groups are constructed to have m observations each on the average.
#  Otherwise, if g is given, g quantile groups will be constructed.
#  If instead cuts is given, proportions will be computed based on the
#  cut points in the vector cuts, e.g. cuts<-seq(0,1,by=.2). 
#  If legendloc is given, a legend will be plotted there
#  Otherwise, it is placed at (.6, .38)
#  Use legendloc=locator(1) to use the mouse for legend positioning.
#  Use legendloc="none" to suppress legend.
#  If statloc is given, some statistics will be plotted there
#  Use statloc=locator(1) to use the mouse.  This is done after the legend.
#  legendloc and statloc can be lists as returned by locator() or they
#  can be vectors, e.g. c(x,y).
#
#Frank Harrell 1 Jun 91
#
val.prob <- function(p,y,logit,group,weights=rep(1,length(y)),normwt=FALSE,
					 pl=TRUE,smooth=TRUE,logistic.cal=TRUE,
					 xlab="Predicted Probability", ylab="Actual Probability",
					 lim=c(0,1),m,g,cuts,emax.lim=c(0,1),
					 legendloc=lim[1]+c(.55*diff(lim),.27*diff(lim)),
					 statloc=c(0,.9),riskdist="calibrated",cex=.75, mkh=.02,
					 connect.group=FALSE, connect.smooth=TRUE, 
					 g.group=4, evaluate=100, nmin=0) {

if(missing(p)) p <- 1/(1+exp(-logit))	else logit <- logb(p/(1-p))
if(length(p)!=length(y))stop("lengths of p or logit and y do not agree")
names(p) <- names(y) <- names(logit) <- NULL

if(!missing(group)) {
  if(length(group)==1 && is.logical(group) && group)
	group <- rep('',length(y))
  if(!is.factor(group)) group <- 
				 if(is.logical(group) || is.character(group)) 
				 as.factor(group) else cut2(group, g=g.group)
  names(group) <- NULL
  nma <- !(is.na(p + y + weights) | is.na(group))
  ng  <- length(levels(group))
} else {
  nma <- !is.na(p + y + weights)
  ng <- 0
}

logit <- logit[nma]
y <- y[nma]
p <- p[nma]
if(ng > 0) {
  group <- group[nma]
  weights <- weights[nma]
  return(val.probg(p, y, group, evaluate, weights, normwt, nmin) )
}

if(length(unique(p))==1) {   #22Sep94
  P <- mean(y)
  Intc <- logb(P/(1-P))
  n <- length(y)
  D <- -1/n
  L01 <- -2 * sum(y * logit - logb(1 + exp(logit)), na.rm=TRUE)
  L.cal <- -2 * sum(y * Intc - logb(1 + exp(Intc)), na.rm=TRUE)
  U.chisq <- L01 - L.cal
  U.p <- 1 - pchisq(U.chisq, 1)
  U <- (U.chisq - 1)/n
  Q <- D - U
  stats <- c(0, .5, 0, D, 0, 1, U, U.chisq, U.p, Q, mean((y-p[1])^2), Intc, 0,
             rep(abs(p[1]-P),2)) 
  names(stats) <- c("Dxy","C (ROC)", 
	"R2","D","D:Chi-sq","D:p","U","U:Chi-sq","U:p","Q",
	"Brier","Intercept","Slope","Emax","Eavg")
  return(stats)
}

i <- !is.infinite(logit)
nm <- sum(!i)
if(nm>0) warning(paste(nm,
   "observations deleted from logistic calibration due to probs. of 0 or 1"))
f <- lrm.fit(logit[i], y[i])
stats <- f$stats
n <- stats["Obs"]
predprob <- seq(emax.lim[1],emax.lim[2],by=.0005)
lt <- f$coef[1]+f$coef[2]*logb(predprob/(1-predprob))
calp <- 1/(1+exp(-lt))
emax <- max(abs(predprob-calp))

if(pl)	{
	plot(.5,.5,xlim=lim, ylim=lim, type="n", xlab=xlab, ylab=ylab)
	abline(0,1,lty=2)
	lt <- 2; leg <- "Ideal"; marks <- -1
	if(logistic.cal)	{
	   lt <- c(lt, 1); leg <- c(leg, "Logistic calibration")
	   marks <- c(marks,-1)	}
	if(smooth)	{
		Sm <- lowess(p,y,iter=0)
		if(connect.smooth)	{ 
 		  lines(Sm,lty=3)
		  lt <- c(lt, 3)
		  marks <- c(marks, -1)	}	else	{
		  points(Sm)
		  lt <- c(lt, 0)
		  marks <- c(marks, 1)	}
		leg <- c(leg, "Nonparametric")
		cal.smooth <- approx(Sm, xout=p)$y
		eavg <- mean(abs(p-cal.smooth))
			}
	if(!missing(m) | !missing(g) | !missing(cuts))	{
		if(!missing(m)) q <- cut2(p, m=m, levels.mean=TRUE, digits=7)
		else if(!missing(g)) q <- cut2(p, g=g, levels.mean=TRUE, digits=7)
		else if(!missing(cuts)) q <- cut2(p, cuts=cuts, levels.mean=TRUE,
			digits=7)
#		means <- tapply(p, q, function(x)mean(x,na.rm=TRUE))
		means <- if(.R.)as.double(levels(q)) else as.single(levels(q))
		prop <- tapply(y, q, function(x)mean(x,na.rm=TRUE))
		points(means, prop, pch=2)
		if(connect.group) {lines(means, prop); lt <- c(lt, 1)}
		else lt <- c(lt, 0)
		leg <- c(leg, "Grouped observations")
		marks <- c(marks, 2)
							}
	}
lr <- stats["Model L.R."]
p.lr <- stats["P"]
D <- (lr-1)/n
L01 <- -2 * sum(y * logit - logb(1 + exp(logit)), na.rm=TRUE)
U.chisq <- L01 - f$deviance[2]
p.U <- 1 - pchisq(U.chisq, 2)
U <- (U.chisq - 2)/n
Q <- D - U
Dxy <- stats["Dxy"]
C <- stats["C"]
R2 <- stats["R2"]
B <- sum((p-y)^2)/n
stats <- c(Dxy, C, R2, D, lr, p.lr, U, U.chisq, p.U, Q, B, f$coef, emax)
names(stats) <- c("Dxy","C (ROC)", 
	"R2","D","D:Chi-sq","D:p","U","U:Chi-sq","U:p","Q",
	"Brier","Intercept","Slope","Emax")
if(smooth) stats <- c(stats, c(Eavg=eavg))
if(pl)			{
	logit <- seq(-7, 7, length=200)
	prob <- 1/(1+exp(-logit))
	pred.prob <- f$coef[1] + f$coef[2] * logit
	pred.prob <- 1/(1+exp(-pred.prob))
	if(logistic.cal)lines(prob, pred.prob, lty=1)
#	pc <- rep(" ", length(lt))
#	pc[lt==0] <- "."
	lp <- legendloc
	if(!is.logical(lp))	{
		if(!is.list(lp)) lp <- list(x=lp[1],y=lp[2])
        if(.R.) legend(lp, leg, lty=lt, pch=marks, cex=cex, bty="n") else
		legend(lp, leg, lty=lt, marks=marks, mkh=mkh, cex=cex,
               bty="n")
				}
	if(!is.logical(statloc))	{
		dostats <- c(1,2,3,4,7,10,11,12,13,14)
		leg <- format(names(stats)[dostats]) #constant length
		leg <- paste(leg, ":", format(stats[dostats]),sep="")
		if(!is.list(statloc)) statloc <- list(x=statloc[1],y=statloc[2])
		text(statloc,paste(format(names(stats[dostats])),collapse="\n"),
		   adj=0,cex=cex)
		text(statloc$x+.225*diff(lim),statloc$y,
                     paste(format(round(stats[dostats],3)),
	                   collapse="\n"),adj=1,cex=cex)
	#	legend(statloc, leg, lty=rep(0, length(dostats)))
      }
	if(is.character(riskdist))		{
      if(riskdist=="calibrated")   {
        x <- f$coef[1]+f$coef[2]*logb(p/(1-p))
        x <- 1/(1+exp(-x))
        x[p==0] <- 0; x[p==1] <- 1}
      else x <- p
      bins <- seq(lim[1],lim[2],length=101)
      x <- x[x>=lim[1] & x<=lim[2]]
      f <- table(if(version$major < 5)cut(x,bins) else oldCut(x,bins))
      j <- f>0
      bins <- (bins[-101])[j] ; f <- f[j]; f <- lim[1]+.15*diff(lim)*f/max(f)
      segments(bins,0,bins,f)
    }
  }	
stats
}


val.probg <- function(p, y, group, evaluate=100, weights, normwt, nmin) {
  if(normwt) weights <- length(y)*weights/sum(weights)
  ng <- length(lg <- levels(group))
  if(ng==1) {ng <- 0; lg <- character(0)}
  stats <- matrix(NA, nrow=ng+1, ncol=12,
				  dimnames=list(nn <- c(lg,'Overall'), 
					c('n','Pavg','Obs','ChiSq','ChiSq2','Eavg',
					  'Eavg/P90','Med OR','C','B','B ChiSq','B cal')))
  curves <- vector('list',ng+1)
  names(curves) <- nn
  q.limits <- c(.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99)
  limits <- matrix(NA, nrow=ng+1, ncol=length(q.limits),
				   dimnames=list(nn, as.character(q.limits)))
  for(i in 1:(ng+1)) {
	s <- if(i==(ng+1)) 1:length(p) else group==lg[i]
	P <- p[s]
	Y <- y[s]
	wt <- weights[s]
	lims <- wtd.quantile(P, wt, q.limits, na.rm=FALSE, normwt=FALSE)
	limits[i,] <- lims
	n <- sum(wt)
	n1 <- sum(wt[Y == 1])
	c.index <- (mean(wtd.rank(P,wt,na.rm=FALSE,normwt=FALSE)[Y == 1]) - 
				(n1 + 1)/2)/(n - n1)
	# c.index <- somers2(P,Y,wt,normwt=FALSE,na.rm=FALSE)['C']
	sm <- wtd.loess.noiter(P, Y, wt, na.rm=FALSE, type='all')  
	##all -> return all points
	curve <- if(length(sm$x) > evaluate)
	  approx(sm, xout=seq(min(P),max(P),length=evaluate)) else {
		o <- order(sm$x)
		nd <- !duplicated(sm$x[o])
		list(x=sm$x[o][nd], y=sm$y[o][nd])
	  }
	if(nmin > 0) {
	  cuts <- wtd.quantile(P, wt, c(nmin, n-nmin)/n, normwt=FALSE, na.rm=FALSE)
	  keep <- curve$x >= cuts[1] & curve$x <= cuts[2]
	  curve <- list(x=curve$x[keep], y=curve$y[keep])
	}
	curves[[i]] <- curve
	cal.smooth <- sm$y
	eavg <- sum(wt*abs(P-cal.smooth))/n
	b    <- sum(wt*((P-Y)^2))/n
	E0b  <- sum(wt*P*(1-P))/n
	Vb   <- sum(wt*((1-2*P)^2)*P*(1-P))/n/n
	bchisq <- (b - E0b)^2 / Vb
	b.cal  <- sum(wt*((cal.smooth-Y)^2))/n

	pred  <- sum(wt*P)/n
	obs   <- sum(wt*Y)/n
	L <- ifelse(P==0 | P==1, NA, logb(P/(1-P)))
	w <- !is.na(L)
	del <- matrix(c(sum((wt*(Y-P))[w]),sum((wt*L*(Y-P))[w])),ncol=2)
	v <- rbind(c(sum((wt*P*(1-P))[w]), sum((wt*L*P*(1-P))[w])),
			   c(NA, sum((wt*L*L*P*(1-P))[w])))
	v[2,1] <- v[1,2]
	chisq  <- (sum(wt*(P-Y))^2) / sum(wt*P*(1-P))
	chisq2 <- del %*% solve(v) %*% t(del)
	p90    <- diff(lims[c(3,9)])
	Lcal   <- ifelse(cal.smooth <= 0 | cal.smooth >= 1, NA,
				   logb(cal.smooth/(1-cal.smooth)))
	or <- exp(wtd.quantile(abs(L - Lcal), wt, .5, na.rm=TRUE, normwt=FALSE))
	stats[i,] <- c(n,pred,obs,chisq,chisq2,eavg,eavg/p90,or,c.index,
				   b,bchisq,b.cal)
  }
  structure(list(stats=stats, cal.curves=curves, quantiles=limits), 
			class='val.prob')
}

print.val.prob <- function(x, ...) {
  print(round(x$stats,3))
  cat('\nQuantiles of Predicted Probabilities\n\n')
  print(round(x$quantiles,3))
  invisible()
}

plot.val.prob <- function(x, 
						  xlab="Predicted Probability", 
						  ylab="Actual Probability",
						  lim=c(0,1), statloc=lim, stats=1:12, cex=.5, 
						  lwd.overall=4, quantiles=c(0.05,0.95),
						  flag=function(stats) ifelse(
						   stats[,'ChiSq2'] > qchisq(.99,2) |
						   stats[,'B ChiSq'] > qchisq(.99,1),'*',' '), ...) {

  stats <- x$stats[,stats,drop=FALSE]
  lwd <- rep(par('lwd'), nrow(stats))
  lwd[dimnames(stats)[[1]]=='Overall'] <- lwd.overall
  curves <- x$cal.curves

  labcurve(curves, pl=TRUE, xlim=lim, ylim=lim, 
		   xlab=xlab, ylab=ylab, cex=cex, lwd=lwd, ...)
  abline(a=0,b=1,lty=2)
  if(is.logical(statloc) && !statloc) return(invisible())

  if(length(quantiles)) {
	limits <- x$quantiles
	quant <- round(as.numeric(dimnames(limits)[[2]]),3)
	w <- quant %in% round(quantiles,3)
	if(any(w)) for(j in 1:nrow(limits)) {
	  qu <- limits[j,w]
	  scat1d(qu, y=approx(curves[[j]], xout=qu)$y)
	}
  }

  xx <- statloc[1]; y <- statloc[2]
  for(i in 0:ncol(stats)) {
	column.text <- if(i==0) c('Group',
						paste(flag(stats),dimnames(stats)[[1]],sep='')) else
	  c(dimnames(stats)[[2]][i], 
		format(round(stats[,i],if(i %in% c(4:5,11))1 else 3)))
    cat(column.text,'\n')
	text(xx, y, paste(column.text,collapse='\n'), adj=0, cex=cex)
	xx <- xx + (1+.8*max(nchar(column.text)))*cex*par('cxy')[1]
  }
invisible()
}

val.surv <- function(fit, newdata, S, est.surv, censor) {
  if(missing(S)) {
    S <- fit$y
    if(!length(S)) stop('when S is omitted you must use y=T in the fit')
    itrans <- if(.R.)survreg.distributions[[fit$dist]]$itrans else
              if(.SV4.)survReg.distributions[[fit$dist]]$itrans else
              glm.links["inverse", fit$family[2]]$inverse
    S[,1] <- itrans(S[,1])
  }
  if(!any(attr(S,'class')=='Surv')) stop('S must be a Surv object')
  if(ncol(S)!=2) stop('S must be a right-censored Surv object')
  if(missing(est.surv))
  est.surv <- if(missing(newdata))
    survest(fit, times=S[,1], what='parallel') else
    survest(fit, newdata, times=S[,1], what='parallel')

  n <- nrow(S)
  nac <- if(!missing(fit)) fit$na.action
  if(!missing(censor) && length(censor)>1 && !missing(fit)) {
    if(length(censor) > n && length(nac)) {
      ## Missing observations were deleted during fit
      j <- !is.na(naresid(nac, censor))
      censor <- censor[j]
    }
  if(length(censor) != n)
    stop("length of censor does not match # rows used in fit")
  }

  est.surv.censor <- NULL; lp <- NULL
  if(!missing(censor)) {
    if(missing(fit))
      stop('fit must be specified when censor is specified')
    est.surv.censor <- if(missing(newdata))
      survest(fit, times=censor, what='parallel') else
      survest(fit, newdata, times=censor, what='parallel')
    if(mc <- sum(is.na(est.surv.censor)))
      warning(paste(mc,'observations had missing survival estimates at censoring time'))
    lp <- if(missing(newdata)) predict(fit, type='lp') else
      predict(fit, newdata, type='lp')
  }

  if(length(est.surv) != n)
    stop('length of est.surv must equal number of rows in S')
  if(!.R.) storage.mode(S) <- 'single'
  as <- if(.R.)function(x)x else as.single
  structure(list(S=S, est.surv=as(est.surv),
                 censor.est.surv=if(length(est.surv.censor))
                   as(est.surv.censor),
                 lp=if(length(lp))as(lp),
                 na.action=nac), class='val.surv')
}

plot.val.surv <- function(x, group, g.group=4,
                          what=c('difference','ratio'),
                          type=c('l','b','p'),
                          xlab, ylab, xlim, ylim, datadensity=TRUE,
                          ...) {
  S <- x$S
  est.surv <- x$est.surv
  censor.est.surv <- x$censor.est.surv
  what <- match.arg(what)
  type <- match.arg(type)

  n <- length(est.surv)
  nac <- x$na.action
  
  if(!missing(group)) {
    if(length(group) > n && length(nac)) {
      ## Missing observations were deleted during fit
      j <- !is.na(naresid(nac, est.surv))
      group <- group[j]
    }
    if(length(group) != n)
      stop("length of group does not match # rows used in fit")
    if(!is.factor(group)) group <- 
      if(is.logical(group) || is.character(group)) 
        as.factor(group) else cut2(group, g=g.group)
  }
 
  if(length(censor.est.surv)) {
    if(missing(group)) group <- rep(1, length(censor.est.surv))
    i <- S[,2]==1
    group <- group[i]
    if(sum(i)<2) stop('fewer than 2 uncensored observations')
    y <- switch(what,
                difference=1-est.surv - .5*(1-censor.est.surv),
                ratio=(1 - est.surv) / (.5 * (1 - censor.est.surv)))
    meanF <- tapply(1 - est.surv[i], group, mean, na.rm=TRUE)
    meanE <- tapply(.5*(1-censor.est.surv[i]), group, mean, na.rm=TRUE)
    res <- matrix(cbind(meanF,meanE), ncol=2,
                  dimnames=list(levels(group),
                    c('Mean F(T|T<C,X)','Expected')))
    lp <- x$lp
    lp <- lp[i]; y <- y[i]
    if(missing(xlab)) xlab <- 'Linear Predictor'
    if(missing(ylab)) ylab <-
      switch(what,
             difference='F(T|X,T<=C) = .5F(C|X)',
             ratio='F(T|X,T<=C)/.5F(C|X), C=Censoring Time')
    if(missing(xlim)) xlim <- range(lp, na.rm=TRUE)
    if(missing(ylim)) ylim <- range(y,  na.rm=TRUE)
    if(type=='l') plsmo(lp, y, group=group, xlab=xlab, ylab=ylab,
         xlim=xlim, ylim=ylim,
         datadensity=datadensity, trim=0, ...) else {
      plot(lp, y, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
      if(type!='p') plsmo(lp, y, group=group, add=TRUE,
           xlab=xlab, ylab=ylab, datadensity=datadensity, trim=0, ...)
    }
    abline(h=if(what=='difference')0 else 1, lty=2)
    return(res)
  }

  if(missing(xlab)) xlab <- 'Predicted Pr[T <= observed T]'
  if(missing(ylab)) ylab <- 'Fraction <= x'
  if(missing(xlim)) xlim <- 0:1
  if(missing(ylim)) ylim <- 0:1
  
  if(missing(group)) {
    nma <- !is.na(est.surv + S[,2])
    est.surv <- est.surv[nma]
    S <- S[nma,,drop=FALSE]
    f <- survfitKM(factor(rep(1,length(est.surv))),
                    Surv(1-est.surv,S[,2]),
                    se.fit = FALSE, conf.type = "none")
	tt <- c(0, f$time)
	ss <- c(1, f$surv)
    plot(tt, 1-ss, xlab=xlab, ylab=ylab, type='s',
         xlim=xlim, ylim=ylim, ...)
    abline(a=0, b=1, lty=2)
    return(invisible())
  }

 
  nma <- !(is.na(est.surv + S[,1] + S[,2]) | is.na(group))
  S <- S[nma,,drop=FALSE]
  est.surv <- est.surv[nma]

  ng <- length(lg <- levels(group))
  curves <- vector('list',ng+1)
  names(curves) <- c(lg, 'Overall')
  for(i in 1:(ng+1)) {
	s <- if(i==(ng+1)) rep(TRUE,length(est.surv)) else group==lg[i]
    f <- survfitKM(factor(rep(1,sum(s))),
                    Surv(1-est.surv[s],S[s,2]),
                    se.fit = FALSE, conf.type = "none")
    curves[[i]] <- list(x=c(0,f$time), y=1-c(1,f$surv))
}  
  labcurve(curves, pl=TRUE, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
  abline(a=0,b=1,lty=2)
  invisible()
}



