survplot.Design <-
  function(fit, ..., xlim, 
           ylim=if(loglog)c(-5,1.5) else 
           if(what=="survival" & missing(fun))c(0,1),
           xlab, ylab, time.inc,
           what=c("survival","hazard"),
           type=c("tsiatis","kaplan-meier"),
           conf.type=c("log-log","log","plain","none"),
           conf.int=FALSE, conf=c("bars","bands"), add=FALSE, 
           label.curves=TRUE, abbrev.label=FALSE,
           lty,lwd=par('lwd'),col=1,
           adj.subtitle,loglog=FALSE,fun,n.risk=FALSE,logt=FALSE,
           dots=FALSE,dotsize=.003,grid=FALSE,
           srt.n.risk=0,sep.n.risk=.056,adj.n.risk=1,
           y.n.risk,cex.n.risk=.6, pr=FALSE) {

  what <- match.arg(what)

  if(.R.) ylim <- ylim  ## before R changes missing(fun)

type <- match.arg(type)
conf.type <- match.arg(conf.type)
conf <- match.arg(conf)

psmfit <- inherits(fit,'psm') || (length(fit$fitFunction) &&
                                  any(fit$fitFunction == 'psm'))
##14Nov00 22May01
if(what=="hazard" && !psmfit)
   stop('what="hazard" may only be used for fits from psm')
if(what=="hazard" & conf.int>0) {
   warning('conf.int may only be used with what="survival"')
   conf.int <- FALSE
}

if(loglog) {fun <- function(x) logb(-logb(ifelse(x==0|x==1,NA,x))); use.fun <- TRUE}
else if(!missing(fun)) {
  use.fun <- TRUE
  if(loglog) stop("cannot specify loglog=T with fun")
}
else { fun <- function(x) x; use.fun <- FALSE }

if(what=="hazard" & loglog) stop('may not specify loglog=T with what="hazard"')

if(use.fun | logt | what=="hazard") { dots <- FALSE; grid <- FALSE }

cox <- inherits(fit,"cph") || (length(fit$fitFunction) &&
                               any(fit$fitFunction == 'cph'))
                               ##14Nov00 22May01
if(cox)	{
  if(n.risk | conf.int>0) surv.sum <- fit$surv.summary
  exactci <- !(is.null(fit$x)|is.null(fit$y))
  ltype <- "s"	#step functions for cph
}

else 	{
#   if(n.risk | loglog)
#    stop("the n.risk and loglog options only apply to fits from cph")
   if(n.risk) stop("the n.risk option applies only to fits from cph")
   exactci <- TRUE
   ltype <- "l"
}

##  if(n.risk && .R. && !missing(y.n.risk)) {  ## 3nov02  24apr03
  if(.R.) {
    oxpd <- par('xpd')
    par(xpd=NA)
    on.exit(par(xpd=oxpd))
  }
 
  
labelc <- is.list(label.curves) || label.curves

#atr <- attr(fit$terms, "Design")   17Apr01
atr <- fit$Design
if(!length(atr)) atr <- getOldDesign(fit)

Limval <- Getlim(atr, allow.null=TRUE, need.all=FALSE)
values <- Limval$values
assume <- atr$assume.code
if(is.null(assume))stop("fit does not have design information")
non.ia <- assume!=9	#limit list to main effects factors
f <- sum(non.ia)
name <- atr$name	#interactions are placed at end by design
label <- atr$label
parms <- atr$parms
units <- fit$units
if(missing(ylab)) {
  if(loglog) ylab <- "log(-log Survival Probability)"
  else if(use.fun) ylab <- ""
  else if(what=="hazard") ylab <- "Hazard Function"
  else ylab <- "Survival Probability"
}
if(missing(xlab)) {
  if(logt) xlab <- paste("log Survival Time in ",units,"s",sep="")
  else xlab <- paste(units,"s",sep="")
  }

maxtime <- fit$maxtime
maxtime <- max(pretty(c(0,maxtime)))
if(missing(time.inc)) time.inc <- fit$time.inc

if(missing(xlim)) 
  xlim <- if(logt)logb(c(maxtime/100,maxtime)) else c(0,maxtime)


if(grid) {dots <- FALSE; if(is.logical(grid)) grid <- .05}

factors <- list(...)
nf <- length(factors)
if(nf<1)stop("must specify 1 factor to plot")
which <- charmatch(names(factors), name, 0)
if(any(which==0))stop(paste("factor(s) not in design:",
	paste(names(factors)[which==0],collapse=" ")))
if(any(assume[which]==9 | assume[which]==10))
   stop("may not plot interaction or matrix effects")

#Number of non-slopes
nrp <- num.intercepts(fit)

if(is.logical(conf.int)) {
  if(conf.int) conf.int <- .95	else conf.int <- 0
}
zcrit <- qnorm((1+conf.int)/2)

xadj <- Limval$limits[2,assume!=9,drop=FALSE]
#for(i in 1:length(xadj)) if(is.factor(xadj[[i]]))
#  xadj[[i]] <- as.character(xadj[[i]])  #preserves class data.frame
#commented out 14Feb95
#Use default adjusted to, replace some later

nv <- 1
y <- factors[[1]]
iy <- which[1]
iyt <- assume[iy]
y <- value.chk(atr, iy, y, 0, Limval)
nc <- length(y)

if(missing(adj.subtitle)) {
  if(add)adj.subtitle <- FALSE else adj.subtitle <- f-nv <= 4
}
						
jf <- nv
if(nf>nv) for(i in which[(nv+1):nf]) {
  jf <- jf+1
  z <- factors[[jf]]
  if(!is.na(z)) z <- value.chk(atr, i, z, 0, Limval)
  if(is.na(z) | length(z)!=1)
	   stop("must specify single value for adjustment factors")
  if(!is.na(z)) xadj[[name[i]]] <- z
}
#Put a legal value in for factor that's moving - will be ignored later.
#Just don't want it excluded from model frame
niy <- name[iy]
if(assume[iy]==5 | assume[iy]==8) xadj[[niy]] <- parms[[niy]][1] else
if(assume[iy]<8 && !is.null(V <- values[[niy]]) && is.character(V))
	xadj[[niy]] <- V[1] else
xadj[[niy]] <- 1
isna <- sapply(xadj, is.na)
if(any(isna)) stop(
   paste("settings not defined here or with datadist for",
	 paste(name[isna],collapse=" ")))

beta <- fit$coef
#if(conf.int>0) cov <- fit$var

curve.labels <- NULL
xd <- xlim[2]-xlim[1]
if(n.risk & !add) {
  mar <- par()$mar
  if(mar[4]<4){mar[4] <- mar[4]+2; par(mar=mar)}
}

# One curve for each value of y, excl style used for C.L.
if(missing(lty)) lty <- seq(nc+1)[-2] else
lty <- rep(lty, length=nc)
col <- rep(col, length=nc)
lwd <- rep(lwd, length=nc)

i <- 0
abbrevy <- if(abbrev.label) abbreviate(y) else y
abbrevy <- if(is.factor(abbrevy)) as.character(abbrevy) else format(abbrevy)

if(labelc) curves <- vector('list',nc)

for(ay in y) {
  i <- i+1
  adj <- xadj
  adj[[name[iy]]] <- ay
  adj <- data.frame(adj)

  index.y <- i
  ci <- conf.int
  w <- survest(fit,newdata=adj,fun=fun,what=what,
		conf.int=ci,type=type,conf.type=conf.type)
  time <- w$time
  if(logt) time <- logb(time)
  s <- !is.na(time) & (time>=xlim[1])  #& (time<=xlim[2])
  surv <- w$surv
  if(is.null(ylim)) ylim <- range(surv, na.rm=TRUE)
  stratum <- w$strata
  if(is.null(stratum)) stratum <- 1
  if(!is.na(stratum)) {
    #can be NA if illegal strata combinations requested
    if(is.factor(ay)) cl <- as.character(ay)
    else cl <- format(ay)
    curve.labels <- c(curve.labels, abbrevy[i])
    if(i==1 & !add) {				
      plot(time,surv,xlab=xlab,xlim=xlim,
           ylab=ylab,ylim=ylim,type="n",axes=FALSE)	
      mgp.axis(1,at=if(logt)pretty(xlim) else
               seq(xlim[1],max(pretty(xlim)),time.inc),labels=TRUE)  #2Jun99

      mgp.axis(2, at=pretty(ylim))  #7Feb98, 2Jun99
      if(!logt & (dots|grid)) {
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
        if(dots) for(tt in xp)symbols(rep(tt,length(yp)),yp,
									  circle=rep(dotsize,length(yp)),
									  inches=dotsize,add=TRUE)
		  else abline(h=yp, v=xp, col=grid)
      }
    }
    tim <- time[s]; srv <- surv[s]
    if(conf.int>0 && conf=='bands') {  ## 5Apr02
      blower <- w$lower[s]
      bupper <- w$upper[s]
    }
    if(max(tim)>xlim[2]) {
      if(ltype=="s") {
        ##Get estimate at last permissible point to plot
        ## s.last <- min(srv[tim<=xlim[2]+1e-6])  #not work with function
        s.last <- srv[tim <= xlim[2] + 1e-6]
        s.last <- s.last[length(s.last)]
        k <- tim < xlim[2]
        tim <- c(tim[k], xlim[2]); srv <- c(srv[k], s.last)
        if(conf.int>0 && conf=='bands') {  ## 5Apr02
          low.last <- blower[time <= xlim[2] + 1e-6]
          low.last <- low.last[length(low.last)]
          up.last  <- bupper[time <= xlim[2] + 1e-6]
          up.last  <- up.last[length(up.last)]
          blower <- c(blower[k],low.last)
          bupper <- c(bupper[k],up.last)
        }
      }
      else tim[tim>xlim[2]] <- NA
    }
		   
    ##don't let step function go beyond x-axis -
    ##this cuts it off but allows step to proceed axis end

    lines(tim,srv,type=ltype,lty=lty[i],col=col[i],lwd=lwd[i])
	if(labelc) curves[[i]] <- list(tim, srv)

    if(pr) {
      zest <- rbind(tim,srv)
      dimnames(zest) <- list(c("Time","Survival"), rep("",length(srv)))
      cat("\nEstimates for ", cl,"\n\n")
      print(zest, digits=3)
    }
    if(conf.int>0) {
      if(conf=="bands") {
        lines(tim,blower,type=ltype,lty=2,col=col[i],lwd=1)
        lines(tim,bupper,type=ltype,lty=2,col=col[i],lwd=1)
        ## was w$lower[s], w$upper[s] 5Apr02
      } else {
	if(exactci) { # not from cph(surv=T)
      tt <- seq(0, maxtime, time.inc)
      v <- survest(fit, newdata=adj,times=tt, what=what, fun=fun,
                   conf.int=ci, type=type,conf.type=conf.type)
      tt <- v$time   #may not get predictions at all t
      ss <- v$surv; lower <- v$lower; upper <- v$upper
      if(is.null(ylim)) ylim <- range(ss, na.rm=TRUE)
      if(logt) tt <- logb(ifelse(tt==0,NA,tt))
    }
	else {
      tt <- as.single(dimnames(surv.sum)[[1]])
      if(logt) tt <- logb(tt)
      ss <- surv.sum[,stratum,1]^exp(w$linear.predictors)
      se <- surv.sum[,stratum,3]
      ss <- fun(ss)
      lower <- fun(ss^exp(zcrit*se))
      upper <- fun(ss^exp(-zcrit*se))
	  ss[is.infinite(ss)] <- NA; lower[is.infinite(lower)] <- NA
      upper[is.infinite(upper)] <- NA
    }
    tt <- tt + xd*(i-1)*.01
    errbar(tt, ss, upper, lower, add=TRUE, lty=lty[i], col=col[i])
  }
    }
    if(n.risk) {
      if(!is.null(Y <- fit$y)) {
        tt <- seq(max(0,xlim[1]),min(maxtime,xlim[2]),by=time.inc)
        ny <- ncol(Y)
        if(is.null(str <- attr(Y,"strata"))) Y <- Y[,ny-1]
        else Y <- Y[oldUnclass(str)==oldUnclass(stratum),ny-1]
        nrisk <- rev(cumsum(table(
                                  if(version$major < 5)
                     cut(-Y,sort(unique(-c(tt,range(Y)+c(-1,1)))))
      else oldCut(-Y,sort(unique(-c(tt,range(Y)+c(-1,1)))))
                                  )[-length(tt)-1]))	
      }
      else {
        if(is.null(surv.sum))
          stop("you must use surv=T or y=T in fit to use n.risk=T")
        tt <- as.single(dimnames(surv.sum)[[1]])
        l <- (tt >= xlim[1]) & (tt <= xlim[2])
        tt <- tt[l]; nrisk <- surv.sum[l,stratum,2]
      }
      tt[1] <- xlim[1]  #was xd*.015, .030, .035
      yd <- ylim[2]-ylim[1]
	  if(missing(y.n.risk)) y.n.risk <- ylim[1]
      yy <- y.n.risk+yd*(nc-i)*sep.n.risk #was .029, .038, .049
	  # Generalized 11Oct96 to y.n.risk from ylim[1]
      nri <- nrisk; nri[tt>xlim[2]] <- NA   # added 2Sep94	
      text(tt[1],yy,nri[1],cex=cex.n.risk,
		 adj=adj.n.risk,srt=srt.n.risk)
      text(tt[-1],yy,nri[-1],cex=cex.n.risk,adj=1)
      text(xlim[2]+xd*.025,yy,adj=0,curve.labels[i],cex=cex.n.risk)
    }
  }
}

if(labelc) labcurve(curves, curve.labels, type=ltype, lty=lty, col=col,
					lwd=lwd, opts=label.curves)

adjust <- ""
if(f>nv) for(i in 1:f) {
  if(!any(i==which[1:nv])) {
    if(is.numeric(xadj[,i])) fv <- format(xadj[,i])
    #format won't work with factors
    else fv <- as.character(xadj[,i])
    adjust <- paste(adjust, name[i], "=", fv, " ",sep="")}
  }

if(adjust!="" & adj.subtitle) title(sub=paste("Adjusted to:",adjust),
									adj=0,cex=.6)

invisible(list(adjust=adjust, curve.labels=curve.labels))
}

