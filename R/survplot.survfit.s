survplot.survfit <- function(fit, xlim, 
							 ylim, xlab, ylab, time.inc,
							 conf=c("bars","bands","none"), add=FALSE, 
							 label.curves=TRUE, abbrev.label=FALSE,
							 lty,lwd=par('lwd'),col=1,
							 loglog=FALSE,fun,n.risk=FALSE,logt=FALSE,
							 dots=FALSE,dotsize=.003,
							 grid=FALSE,
							 srt.n.risk=0,sep.n.risk=.056,adj.n.risk=1,
							 y.n.risk,cex.n.risk=.6, pr=FALSE, ...) {

  conf <- match.arg(conf)
  conf.int <- fit$conf.int
  if(!length(conf.int) | conf=="none") conf.int <- 0

  ##??Mgp <- mgp.axis.labels(type='x and y') #7Feb98 - see Hmisc Misc.s, gs.slide

  units <- fit$units
  if(!length(units)) units <- "Day"
  maxtime <- fit$maxtime
  if(!length(maxtime)) maxtime <- max(fit$time)
  mintime <- min(fit$time,0)
  pret <- pretty(c(mintime,maxtime))  #1Apr95
  maxtime <- max(pret)
  mintime <- min(pret)
  if(missing(time.inc)) {
	time.inc <- switch(units,Day=30,Month=1,Year=1,(maxtime-mintime)/10)
	if(time.inc>maxtime) time.inc <- (maxtime-mintime)/10
  }
  if(n.risk && !length(fit$n.risk)) {
	n.risk <- FALSE
	warning("fit does not have number at risk\nIs probably from a parametric model\nn.risk set to F")
  }
  trans <- loglog | !missing(fun)
  if(missing(ylab))	 {
	if(loglog) ylab <- "log(-log Survival Probability)"
	  else if(trans) ylab <- ""
		else ylab <- "Survival Probability"
  }
  if(loglog) fun <- function(w) logb(-logb(ifelse(w==0|w==1,NA,w)))
	else if(!trans) fun <- function(w) w

  if(missing(xlab)) 		 {
	if(logt) xlab <- paste("log Survival Time in ",units,"s",sep="")
	  else xlab <- if(units==' ') '' else paste(units,"s",sep="")	}

  if(missing(xlim)) 
	xlim <- if(logt)logb(c(maxtime/100,maxtime)) else c(mintime,maxtime)

  if(trans)		 {
	fit$surv <- fun(fit$surv)
	fit$surv[is.infinite(fit$surv)] <- NA
										#  handle e.g. logit(1) - Inf would mess up ylim in plot()
	if(conf.int>0)	 {
      fit$lower <- fun(fit$lower)
      fit$upper <- fun(fit$upper)
      fit$lower[is.infinite(fit$lower)] <- NA
      fit$upper[is.infinite(fit$upper)] <- NA
      if(missing(ylim)) ylim <- range(c(fit$lower,fit$upper),na.rm=TRUE)
	}
	  else if(missing(ylim)) ylim <- range(fit$surv,na.rm=TRUE)
  }
	else if(missing(ylim)) ylim <- c(0,1)

  if(grid) {dots <- FALSE; if(is.logical(grid)) grid <- .05}
  if(logt|trans) { dots <- FALSE; grid <- FALSE }

  slev <- names(fit$strata)
  sleva <- if(abbrev.label) abbreviate(slev) else slev
  ns <- length(slev)
  slevp <- ns>0

  labelc <- is.list(label.curves) || label.curves
  if(!slevp) labelc <- FALSE

  ns <- max(ns,1)


  y <- 1:ns
  if(ns==1) stemp <- rep(1, length(fit$time))
	else stemp <- rep(1:ns, fit$strata)

  if(n.risk | (conf.int>0 & conf=="bars"))  {
	stime <- seq(mintime,maxtime,time.inc)
	v <- summary(fit, times=stime, print.it=FALSE)
	vs <- if(ns==1) rep(1, length(v$time)) else oldUnclass(v$strata)
  }
  xd <- xlim[2]-xlim[1]
  yd <- ylim[2]-ylim[1]

  if(n.risk && !add)	 {
	mar <- par()$mar
	if(mar[4]<4) {mar[4] <- mar[4]+2; par(mar=mar)}
  }
  ## One curve for each value of y, excl style used for C.L.
  lty <- if(missing(lty)) seq(ns+1)[-2] else rep(lty, length=ns)
  lwd <- rep(lwd, length=ns)
  col <- rep(col, length=ns)

  if(labelc) curves <- vector('list',ns)

  ##if(n.risk && .R. && !missing(y.n.risk)) {  ## 3nov02  24apr03
    if(.R.) {
      oxpd <- par('xpd')
      par(xpd=NA)
      on.exit(par(xpd=oxpd))
    }
  for(i in 1:ns) {
	st <- stemp==i
	time <- fit$time[st]
	surv <- fit$surv[st]
	if(logt) time <- logb(time)
	s <- !is.na(time) & (time>=xlim[1])  #& (time<=xlim[2])
	if(i==1 & !add) {
	  plot(time,surv,xlab=xlab,xlim=xlim,
		   ylab=ylab,ylim=ylim,type="n",axes=FALSE)
	  mgp.axis(1,at=if(logt)pretty(xlim) else
               seq(xlim[1],max(pretty(xlim)),time.inc),
               labels=TRUE)  #7Feb98, 2Jun99

	  mgp.axis(2, at=pretty(ylim))  #2Jun99
	  if(dots|grid) {
		xlm <- pretty(xlim)
		xlm <- c(xlm[1],xlm[length(xlm)])
		xp <- seq(xlm[1],xlm[2],by=time.inc)
		yd <- ylim[2]-ylim[1]
		if(yd<=.1)yi <- .01
		  else if(yd<=.2) yi <- .025
			else if(yd<=.4) yi <- .05
			  else yi <- .1
		yp <- seq(ylim[2],ylim[1]+if(n.risk && missing(y.n.risk))yi else 0,
				  by=-yi)
		if(dots)
		  for(tt in xp)symbols(rep(tt,length(yp)),yp,
							   circle=rep(dotsize,length(yp)),
							   inches=dotsize,add=TRUE)
		  else abline(h=yp, v=xp, col=grid)
	  }
	}
	tim <- time[s]; srv <- surv[s]
    if(conf.int > 0 && conf=="bands") {  ## 5Apr02
      blower <- fit$lower[st][s]
      bupper <- fit$upper[st][s]
    }
	##don't let step function go beyond x-axis -
	##this cuts it off but allows step to proceed axis end
	if(max(tim)>xlim[2])	 {
	  srvl <- srv[tim<=xlim[2]+1e-6]
	  ##	  s.last <- min(srv[tim<=xlim[2]+1e-6]) #not work w/fun
	  s.last <- srvl[length(srvl)]
	  k <- tim<xlim[2]
	  tim <- c(tim[k], xlim[2]); srv <- c(srv[k], s.last)
      if(conf.int > 0 && conf=="bands") {  ## 5Apr02
          low.last <- blower[time <= xlim[2] + 1e-6]
          low.last <- low.last[length(low.last)]
          up.last  <- bupper[time <= xlim[2] + 1e-6]
          up.last  <- up.last[length(up.last)]
          blower <- c(blower[k],low.last)
          bupper <- c(bupper[k],up.last)
      }
	}
	if(logt) {
	  lines(tim,srv,type="s",lty=lty[i],col=col[i],lwd=lwd[i])
	  if(labelc) curves[[i]] <- list(tim,srv)
	}
	  else {
		xxx <- c(mintime,tim);  yyy <- c(fun(1),srv)
		lines(xxx,yyy,type="s",lty=lty[i],col=col[i],lwd=lwd[i])
		if(labelc) curves[[i]] <- list(xxx,yyy)
	  }
	if(pr){zest <- rbind(time[s],surv[s])
		   dimnames(zest) <- list(c("Time","Survival"),
								  rep("",sum(s)))
		   if(slevp)cat("\nEstimates for ", slev[i],"\n\n")
		   print(zest, digits=3)  }
	if(conf.int>0) {
	  if(conf=="bands") {
		if(logt) {
		  lines(tim,blower,type="s",lty=2,col=col[i],lwd=1)  ## 5Apr02
		  lines(tim,bupper,type="s",lty=2,col=col[i],lwd=1)
		}
		  else {
			lines(c(mintime,tim),c(fun(1),blower),	#temp 5Apr02
				  type="s",lty=2,col=col[i],lwd=1)
			lines(c(0,tim),c(fun(1),bupper),	#temp  5Apr02
				  type="s",lty=2,col=col[i],lwd=1)
		  }
	  } else {
		j <- vs==i
		tt <- v$time[j]  #may not get predictions at all t
		ss <- v$surv[j]; lower <- v$lower[j]; 
		upper <- v$upper[j]
		if(logt) tt <- logb(ifelse(tt==0,NA,tt))
		tt <- tt + xd*(i-1)*.01
		errbar(tt, ss, upper, lower, add=TRUE, lty=lty[i],
			   col=col[i])
	  }
	}

	if(n.risk) {
	  j <- vs==i; tt <- v$time[j]; nrisk <- v$n.risk[j]
	  tt[1] <- xlim[1]  #was xd*.015, .030, .035
	  if(missing(y.n.risk)) y.n.risk <- ylim[1]   ## 11Oct96
	  yy <- y.n.risk+yd*(ns-i)*sep.n.risk
										#was .029, .038, .049
	  nri <- nrisk; nri[tt>xlim[2]] <- NA # added 2Sep94	
	  text(tt[1],yy,nri[1],cex=cex.n.risk,
		   adj=adj.n.risk,srt=srt.n.risk)
	  text(tt[-1],yy,nri[-1],cex=cex.n.risk,adj=1)
	  if(slevp)text(xlim[2]+xd*.025,
					yy,adj=0,sleva[i],cex=cex.n.risk)
	}
  }


  if(labelc) labcurve(curves, sleva, type='s', lty=lty, lwd=lwd,
					  opts=label.curves, col=col)
  
  invisible(slev)
}


