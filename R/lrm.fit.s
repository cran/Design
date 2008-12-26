#If est is specified and it is not 1:ncol(x), user will have to negate
# $var[est,est] before running matinv on non.slopes+(1:nx)[-est] in 
# obtaining score
# statistics.
# Use est=NULL to compute score stat components for all vars
# fit$non.slopes = # intercepts
# trace to print loglik at each iteration
# Set tol=0 to turn off singularity checking
# tol is only used during iteractions, not for final inversion (since
# solve does not pass the tolerance argument).
#Mod 1-2-91 : change !is.matrix to is.vector
#Sent to statlib : ??/??/??
#Mod 5-24-91: if maxiter=1, does not compute p-values
#Mod 6-11-91: model with no covariables return proper fail, added is.null to
#             missing(x), is.vector to !is.matrix(x)
#Mod 10-8-91: Changed missing data routines to tstna, naset (see na.fortran.f),
#	      added specialsok=T to .Fortran calls
#Mod 10-11-91:Added class attribute "logist" to fit object, improved call
#Mod 10-30-91:Changed to lrm.fit for use with lrm function,
#             removed subset, print.iter->trace,maxiter->maxit,dlike->eps,
#             eps->tol, f$coef->f$coefficients
#             Remove attributes(x) from fit object
#Mod 3-5-92  :Use solvet instead of solve, to pass tol argument
#    6-9-92  :Change to Nagelkerke R2
#    9-27.92 :Remove dyn load commands (using .First.lib now)
#    5-23-94 :Check for 0 length as well as NULL
#   11-28-94 :added Brier score,
#	      return linear predictor, get rid of "nused", improve NA logic
#    1-17-95 :added penalty, penalty.matrix
#    9-30-95 :changed penalty matrix to be self-contained
#    5-06-96 :return information matrix
#    6-06-02 :added back weights, normwt like SAS PROC LOGIST
#    1-17-03 :made all versions use weights, double precision for x,y

lrm.fit <- function(x,y,offset,initial,est,
                    maxit=12,eps=.025,tol=1E-7,trace=FALSE,
                    penalty.matrix=NULL,weights=NULL,normwt=FALSE) {	
cal <- match.call()
opts <- double(12)
opts[1:4] <- c(tol, eps, maxit, trace)
len.penmat <- length(penalty.matrix)

n <- length(y)

wtpres <- TRUE
if(!length(weights)) {
  wtpres <- FALSE
  normwt <- FALSE
  weights <- rep(1, n)
}
if(length(weights) != n) stop('length of wt must equal length of y')
if(normwt) weights <- weights*n/sum(weights)
storage.mode(weights) <- 'double'
opts[12] <- normwt

initial.there <- !missing(initial)
if(missing(x) || length(x)==0) {
	nx <- 0
	xname <- NULL
	if(!missing(est))stop("est was given without x")
	est <- NULL
	x <- 0
  }	else {
    if(!is.matrix(x)) x <- as.matrix(x)
	storage.mode(x) <- "double"   # 17jan03
	dx <- dim(x)
	nx <- dx[2]
	if(dx[1]!=n)stop("x and y must have same length")
	if(missing(est)) est <- 1:nx	else if(length(est)) {
		estr <- range(est)
		if(estr[1]<1 | estr[2]>nx)
			stop("est has illegal column number for x")
		if(any(duplicated(est)))stop("est has duplicates")
		storage.mode(est) <- "integer"
      }
	xname <- dimnames(x)[[2]]
	if(length(xname)==0) xname <- paste("x[",1:nx,"]",sep="")
  }

nxin <- length(est)

if(!is.category(y)) y <- as.category(y)
y <- oldUnclass(y)   # in case is.factor
ylevels <- levels(y)

ofpres <- !missing(offset)
opts[5] <- ofpres
if(ofpres)	{
  if(length(offset)!=n)stop("offset and y must have same length")
  storage.mode(offset) <- "double"  ## 17jan03
} else offset <- 0

if(n<3)stop("must have >=3 non-missing observations")
kint <- as.integer(length(ylevels)-1)
ftable <- integer(501*(kint+1))
levels(y) <- ylevels
numy <- table(y)
names(numy) <- ylevels   ## 9Apr02  needed for R
y <- as.integer(y-1)
nvi <- as.integer(nxin+kint)

sumwty <- tapply(weights, y, sum)
sumwt  <- sum(sumwty)
if(!wtpres && any(numy != sumwty)) stop('program logic error 1')
sumw <- if(normwt) numy else as.integer(round(sumwty))

if(missing(initial)) {
  ## ncum <- rev(cumsum(rev(numy)))[2:(kint+1)]
  ## pp <- ncum/n
  ncum <- rev(cumsum(rev(sumwty)))[2:(kint+1)]
  pp <- ncum/sumwt
  initial <-logb(pp/(1-pp))
  if(ofpres) initial <- initial-mean(offset)
}
if(length(initial) < nvi) initial <- c(initial,rep(0,nvi-length(initial)))
storage.mode(initial) <- "double"

loglik <- -2 * sum(sumwty*logb(sumwty/sum(sumwty)))
## loglik <-  -2 * sum(numy * logb(numy/n))

if(nxin>0) {
  if(len.penmat==0) penalty.matrix <- matrix(0,nrow=nx,ncol=nx)
  if(nrow(penalty.matrix)!=nx || ncol(penalty.matrix)!=nx) 
    stop(paste("penalty.matrix does not have",nx,"rows and columns"))
  penmat <- rbind(
	matrix(0,ncol=kint+nx,nrow=kint),
	cbind(matrix(0,ncol=kint,nrow=nx),penalty.matrix))
} else penmat <- matrix(0, ncol=kint, nrow=kint)
storage.mode(penmat) <- 'double'

if(nxin==0 & !ofpres) {
	loglik <- rep(loglik,2)
	z <- list(coef=initial,u=rep(0,kint),opts=c(rep(0,7),.5,0,0,0))
  }

if(ofpres) {
#Fit model with only intercept(s) and offset
	z <- if(.R.)
      .Fortran("lrmfit",coef=initial,as.integer(0),0,x,y,offset,
               u=double(kint),
               double(kint*(kint+1)/2),loglik=double(1),n,as.integer(0),
               sumw,kint,
               v=double(kint*kint),double(kint),double(kint),
               double(kint),pivot=integer(kint),opts=opts,ftable,
               penmat, weights, PACKAGE="Design") else
      .Fortran("lrmfit",coef=initial,as.integer(0),0,x,y,offset,
               u=double(kint),
               double(kint*(kint+1)/2),loglik=double(1),n,as.integer(0),
               sumw,kint,
               v=double(kint*kint),double(kint),double(kint),
               double(kint),pivot=integer(kint),opts=opts,ftable,
               penmat, weights, PACKAGE="Design")
## 17jan03
##      .Fortran("lrmfit",coef=initial,as.integer(0),0,x,y,offset,
##               u=double(kint),
##               double(kint*(kint+1)/2),loglik=double(1),n,as.integer(0),
##               sumw,kint,
##               v=double(kint*kint),double(kint),double(kint),
##               double(kint),pivot=integer(kint),opts=opts,ftable,
##               penmat, PACKAGE="Design")
	loglik <- c(loglik,z$loglik)
	if(z$opts[6] | z$opts[7]<kint) {
      if(.SV4.) return(structure(list(fail=TRUE,fitFunction='lrm'),
                                 class='Design')) else
      return(structure(list(fail=TRUE),class="lrm"))  ##13Nov00
    }
	initial <- z$coef
			}

if(nxin>0)	{
#Fit model with intercept(s), offset, and any fitted covariables
	z <- if(.R.)
      .Fortran("lrmfit",coef=initial,nxin,est,x,y,offset,
               u=double(nvi),
               double(nvi*(nvi+1)/2),loglik=double(1),n,nx,sumw,nvi,
               v=double(nvi*nvi),double(nvi),double(2*nvi),double(nvi),
               pivot=integer(nvi),opts=opts,ftable,penmat,weights,
               PACKAGE="Design") else
      .Fortran("lrmfit",coef=initial,nxin,est,x,y,offset,
               u=double(nvi),
               double(nvi*(nvi+1)/2),loglik=double(1),n,nx,sumw,nvi,
               v=double(nvi*nvi),double(nvi),double(2*nvi),double(nvi),
               pivot=integer(nvi),opts=opts,ftable,penmat,weights,
               PACKAGE="Design")
    ## 17jan03
##      .Fortran("lrmfit",coef=initial,nxin,est,x,y,offset,
##               u=double(nvi),
##               double(nvi*(nvi+1)/2),loglik=double(1),n,nx,sumw,nvi,
##               v=double(nvi*nvi),double(nvi),double(2*nvi),double(nvi),
##               pivot=integer(nvi),opts=opts,ftable,penmat, PACKAGE="Design")   #2*nvi 28Jul95
	irank <- z$opts[7]
	if(irank < nvi) {
      cat("singular information matrix in lrm.fit (rank=",irank,
          ").  Offending variable(s):\n")
      cat(paste(xname[est[z$pivot[nvi:(irank+1)]-kint]],
                collapse=" "),"\n")
      if(.SV4.) return(structure(list(fail=TRUE,fitFunction='lrm'),
                                 class='Design')) else
      return(structure(list(fail=TRUE),class="lrm"))  ##13Nov00
    }
	loglik <- c(loglik, z$loglik)
  }

dvrg <- z$opts[6] > 0

if(nxin!=nx) {
  ##Set up for score statistics - last model is not refitted but derivatives
  ##with respect to all other columns of x are evaluated
  initial <- rep(0,nx)
  if(length(est))initial[est] <- z$coef[(kint+1):nvi]
  initial <- c(z$coef[1:kint], initial)
  nvi <- as.integer(kint+nx)
  opts[3] <- 1	#Max no. iterations 
  z <- if(.R.)
    .Fortran("lrmfit",coef=initial,nx,1:nx,x,y,offset,
             u=double(nvi),double(nvi*(nvi+1)),double(1),n,nx,
             sumw,nvi,v=double(nvi*nvi),double(nvi),double(nvi),
             double(nvi),integer(nvi),opts=opts,ftable,penmat,weights,
             PACKAGE="Design") else
    .Fortran("lrmfit",coef=initial,nx,1:nx,x,y,offset,
           u=double(nvi),double(nvi*(nvi+1)),double(1),n,nx,
           sumw,nvi,v=double(nvi*nvi),double(nvi),double(nvi),
           double(nvi),integer(nvi),opts=opts,ftable,penmat,weights,
             PACKAGE="Design")
  ## 17jan03
##    .Fortran("lrmfit",coef=initial,nx,1:nx,x,y,offset,
##           u=double(nvi),double(nvi*(nvi+1)),double(1),n,nx,
##           sumw,nvi,v=double(nvi*nvi),double(nvi),double(nvi),
##           double(nvi),integer(nvi),opts=opts,ftable,penmat, PACKAGE="Design")
}

##Invert v with respect to fitted variables
if(nxin==0) elements <- 1:kint	else elements <- c(1:kint,kint+est)
if(nx==0 & !ofpres)	{
  v <- NULL; info.matrix <- NULL   ## info.matrix added 21Aug97
  irank <- kint
} else {
  if(nxin==nx) { 
    info.matrix <- matrix(z$v,nrow=nvi,ncol=nvi)
    v <- solvet(info.matrix, tol=tol)
    irank <- nvi
  }  else {
	info.matrix <- matrix(z$v, nrow=nvi, ncol=nvi)
	v <- matinv(info.matrix,elements,negate=TRUE,eps=tol)
	info.matrix <- info.matrix[elements,elements]
	usc <- z$u[-elements]
	resid.chi2 <- usc %*% solve(v[-elements,-elements],tol=tol) %*% usc
	resid.df <- nx-nxin
	irank <- attr(v,"rank")
	attr(v,"rank") <- NULL
  }
}

#if(kint==1) name <- "Intercept" else name <- paste("Alpha",1:kint,sep="")
if(kint==1) name <- "Intercept" else 
	    name <- paste("y>=",ylevels[2:(kint+1)],sep="")
name <- c(name, xname)
names(z$coef) <- name
names(z$u) <- name
if(length(v)) dimnames(v) <- list(name,name)

llnull <- loglik[length(loglik)-1]
model.lr <- llnull-loglik[length(loglik)]
model.df <- irank-kint
if(initial.there) model.p <- NA else {
  if(model.df>0) model.p <- 1-pchisq(model.lr,model.df) else model.p <- 1
}
r2 <- 1-exp(-model.lr/sumwt)
r2.max <- 1-exp(-llnull/sumwt)
r2 <- r2/r2.max
lp <- matxv(x, z$coef)
B <- mean((1/(1+exp(-lp)) - (y>0))^2)
B <- sum(weights*(plogis(lp) - (y>0))^2)/sum(weights)

stats <- c(n,max(abs(z$u[elements])),model.lr,model.df,
           model.p,z$opts[8],z$opts[9],
           z$opts[10], z$opts[11], r2, B)
## was stats <- c(n,max(abs(z$u[elements])),model.lr,model.df,  21Aug97
##	model.p,round(z$opts[8],3),round(z$opts[9],3),
##	round(z$opts[10],3), round(z$opts[11],3), r2, B)
nam <- c("Obs","Max Deriv",
         "Model L.R.","d.f.","P","C","Dxy",
         "Gamma","Tau-a","R2","Brier")

if(nxin!=nx) {
  stats <- c(stats,resid.chi2,resid.df,1-pchisq(resid.chi2,resid.df))
  nam <- c(nam, "Residual Score","d.f.","P")
}
names(stats) <- nam

if(wtpres) stats <- c(stats, 'Sum of Weights'=sumwt)

retlist <- list(call=cal,freq=numy,sumwty=if(wtpres)sumwty else NULL,
                stats=stats,fail=dvrg,coefficients=z$coef,
                var=v,u=z$u,
                deviance=loglik,
                est=est,non.slopes=kint, linear.predictors=lp,
                penalty.matrix=if(nxin>0 && any(penalty.matrix!=0))
                 penalty.matrix else NULL,
                info.matrix=info.matrix,
                weights=if(wtpres)weights else NULL)

oldClass(retlist) <- 'lrm' # was c("lrm","lm") 17Jul01
retlist
}

#---------------------------------------------------------------------

lrm.fit.strat <- function(x,y,strata,offset,initial,
						  maxit=25,eps=.025,tol=1E-7,trace=FALSE,
						  penalty.matrix=NULL,strata.penalty=0,
                          weights=NULL,normwt) {
cal <- match.call()
opts <- double(12)
len.penmat <- length(penalty.matrix)
lev    <- levels(strata)
nstrat <- length(lev)
strata <- oldUnclass(strata)
n <- length(y)

if(!length(weights)) {
  normwt <- FALSE
  weights <- rep(1,n)
}
if(length(weights) != n) stop('weights and y must be same length')
storage.mode(weights) <- 'double'
opts[12] <- normwt
## weights not implemented for stratified models yet

initial.there <- !missing(initial)
if(missing(x) || length(x)==0) {
	nx <- 0
	xname <- NULL
	x <- 0
  }	else {
    if(!is.matrix(x)) x <- as.matrix(x)
	storage.mode(x) <- "double"  ## 17jan03
	dx <- dim(x)
	nx <- dx[2]
	if(dx[1]!=n)stop("x and y must have same length")
	xname <- dimnames(x)[[2]]
	if(length(xname)==0) xname <- paste("x[",1:nx,"]",sep="")
  }

nxin <- nx

if(!is.category(y)) y <- as.category(y)
y <- oldUnclass(y)   # in case is.factor
ylevels <- levels(y)

ofpres <- !missing(offset)
if(ofpres)	{
  if(length(offset)!=n)stop("offset and y must have same length")
  storage.mode(offset) <- "double"  ## 17jan03
} else offset <- 0

if(n<3)stop("must have >=3 non-missing observations")
kint <- as.integer(length(ylevels)-1)
if(kint!=1) stop('only works for binary y')
ftable <- integer(501*(kint+1))
levels(y) <- ylevels
numy <- table(y)
y <- as.integer(y-1)
nvi <- as.integer(nxin+kint+nstrat-1)
if(missing(initial)) {
	ncum <- rev(cumsum(rev(numy)))[2:(kint+1)]
	pp <- ncum/n
	initial <-logb(pp/(1-pp))
	if(ofpres) initial <- initial-mean(offset)	}
if(length(initial)<nvi) initial <- c(initial,rep(0,nvi-length(initial)))
storage.mode(initial) <- "double"
loglik <- -2 * sum(numy * logb(numy/n))

if(nxin>0) {
  if(len.penmat==0) penalty.matrix <- matrix(0,nrow=nx,ncol=nx)
  if(nrow(penalty.matrix)!=nx || ncol(penalty.matrix)!=nx) 
    stop(paste("penalty.matrix does not have",nx,"rows and columns"))
  penmat <- rbind(
	matrix(0,ncol=kint+nx,nrow=kint),
	cbind(matrix(0,ncol=kint,nrow=nx),penalty.matrix))
} else penmat <- matrix(0, ncol=kint, nrow=kint)
storage.mode(penmat) <- 'double'

if(nxin==0 & !ofpres) {
  loglik <- rep(loglik,2)
  z <- list(coef=initial,u=rep(0,kint),opts=c(rep(0,7),.5,0,0,0))
}

if(ofpres)	{
#Fit model with only intercept(s) and offset
  z <- if(.R.)
    .Fortran("lrmfit",coef=initial,as.integer(0),0,x,y,offset,
             u=double(kint),
             double(kint*(kint+1)/2),loglik=double(1),n,as.integer(0),
             numy,kint,
             v=double(kint*kint),double(kint),double(kint),
             double(kint),pivot=integer(kint),opts=opts,ftable,
             penmat,weights, PACKAGE="Design") else
    .Fortran("lrmfit",coef=initial,as.integer(0),0,x,y,offset,
           u=double(kint),
           double(kint*(kint+1)/2),loglik=double(1),n,as.integer(0),
           numy,kint,
           v=double(kint*kint),double(kint),double(kint),
           double(kint),pivot=integer(kint),opts=opts,ftable,
           penmat, weights, PACKAGE="Design")
  ## 17jan03
##    .Fortran("lrmfit",coef=initial,as.integer(0),0,x,y,offset,
##           u=double(kint),
##           double(kint*(kint+1)/2),loglik=double(1),n,as.integer(0),
##           numy,kint,
##           v=double(kint*kint),double(kint),double(kint),
##           double(kint),pivot=integer(kint),opts=opts,ftable,
##           penmat, PACKAGE="Design")
  loglik <- c(loglik,z$loglik)
  if(z$opts[6] | z$opts[7]<kint) return(list(fail=TRUE,class="lrm"))
  initial <- z$coef
}

		
## Newton-Raphson iterations with patterned matrix inversion for speed
theta  <- initial
iter   <- 0
oldobj <- 1e10
x      <- cbind(1,x)
nns    <- nx+1  ## no. of non-strata parameters

while(iter <= maxit) {
  iter <- iter + 1
  beta <- as.matrix(theta[1:nns])
  tau  <- c(0,theta[-(1:nns)])
  logit <- drop(x %*% beta + tau[strata] + offset)
  pred <- 1/(1+exp(-logit))
  obj  <- -2*sum(y*logb(pred) + (1-y)*logb(1-pred)) + 
	 t(beta) %*% penmat %*% beta + 
	 strata.penalty*sum((tau-mean(tau))^2)
  if(trace)cat('-2logL=',format(obj),'')
  d <- y - pred
  u <-  c(drop(matrix(d,nrow=1) %*% x)-drop(penmat %*% beta), 
             tapply(d,strata,sum)[-1] - strata.penalty*(tau[-1]-mean(tau)))
  u <- as.matrix(u)
  if(trace)cat(format(u),'\n')
## Created patterned information matrix  A B / B' C
## Inverse is AA BB / BB' CC
  pq <- pred*(1-pred)
  A <- crossprod(pq * x, x) + penmat
  B <- t(rowsum(pq * x, strata))[,-1,drop=FALSE]  
## above won't work if a stratum not represented
  dd <- tapply(pq, strata, sum)[-1]
  v <- 1/(dd + strata.penalty)
  vm <- as.matrix(v)
  k <- (strata.penalty/nstrat)/(1 - (strata.penalty/nstrat)*sum(v))
# BCi <- t(v * t(B)) + k * (B %*% vm) %*% t(vm)
# C <- diag(dd + strata.penalty) - strata.penalty/nstrat
  BCi <- B*rep(v,rep(nns,nstrat-1)) + k * (B %*% vm) %*% t(vm)
  AA <- solvet(A - BCi %*% t(B), tol=tol)
  BB <- -AA %*% BCi
  u1 <- u[1:nns,,drop=FALSE]
  u2 <- u[(nns+1):nvi,,drop=FALSE]
  theta <- theta + c(AA %*% u1 + BB %*% u2,
					 t(BB) %*% u1 + vm * u2 + k * vm %*% (t(vm) %*% u2) -
					   t(BB) %*% (BCi %*% u2))
#CC <- diag(drop(vm))+k*vm%*%t(vm)-t(BB) %*% BCi
#					 t(BB) %*% u1 + u2/dd)             FAILS
#theta <- theta + c(AA %*% u1 + BB %*% u2,
#					 t(BB) %*% u1 + (diag(1/dd) - t(BB) %*% BCi) %*% u2)
# theta <- theta + c(solve(A) %*% u1, u2/dd)           SLOW
#  theta <- theta + c(diag(1/diag(A)) %*% u1, u2/dd)   FAILS
  if(abs(obj - oldobj) < eps) break
  oldobj <- obj
}
if(iter > maxit) return(list(fail=TRUE, class='lrm'))

  
xname <- c(xname, lev[-1])

if(kint==1) name <- "Intercept" else 
 name <- paste("y>=",ylevels[2:(kint+1)],sep="")
name <- c(name, xname)
theta <- drop(theta)
names(theta) <- name

loglik <- c(loglik, obj)

dimnames(AA) <- list(name[1:nns],name[1:nns])
dimnames(BB) <- dimnames(BCi) <- list(name[1:nns],name[(nns+1):nvi])
names(BCi)   <- NULL


llnull <- loglik[length(loglik)-1]
model.lr <- llnull-loglik[length(loglik)]
model.df <- nvi - kint
if(initial.there) model.p <- NA 	else 				{
if(model.df>0) model.p <- 1-pchisq(model.lr,model.df) else model.p <- 1	}
r2 <- 1-exp(-model.lr/n)
r2.max <- 1-exp(-llnull/n)
r2 <- r2/r2.max
Brier <- mean((pred - (y>0))^2)

stats <- c(n,max(abs(u)),model.lr,model.df,model.p,
	## z$opts[8],z$opts[9],z$opts[10], z$opts[11], 
		   r2, Brier)
nam <- c("Obs","Max Deriv",	"Model L.R.","d.f.","P",
		 ##"C","Dxy","Gamma","Tau-a",
		 "R2","Brier")
names(stats) <- nam

Varcov <- function(fit,which=c('strata.var','var','strata.var.diag')) {
  which <- match.arg(which)
  strata.penalty <- fit$strata.penalty
  v <- 1 / (fit$strata.unpen.diag.info + strata.penalty)
  nstrat <- fit$nstrat
  k <- (strata.penalty/nstrat)/(1 - (strata.penalty/nstrat)*sum(v))
  sname <- fit$strata.levels[-1]
  CC <- diag(v) + k * v %*% t(v) -t(fit$cov.nonstrata.strata) %*% fit$BCi
  switch(which,
		 strata.var = structure(CC, dimnames=list(sname,sname)),
		 strata.var.diag = structure(diag(CC), names=sname),
		 var = structure(rbind(cbind(fit$var,fit$cov.nonstrata.strata),
		   cbind(t(fit$cov.nonstrata.strata),CC)),
		   dimnames=list(nn <- names(fit$coef),nn)))
}

retlist <- list(call=cal,freq=numy,
				stats=stats,fail=FALSE,coefficients=theta[1:nns],
				non.slopes=1,est=1:(nvi-kint),
				var=AA,u=u,
				deviance=loglik,
				linear.predictors=logit,
				penalty.matrix=if(nxin>0 && any(penalty.matrix!=0)) 
				  penalty.matrix else NULL,
				nstrat=nstrat, strata.levels=lev,
				strata.coefficients=theta[(nns+1):nvi],
				strata.penalty=strata.penalty, 
 				strata.unpen.diag.info=dd,
				cov.nonstrata.strata=BB,
				BCi=BCi,
				Varcov=Varcov,
#				info.matrix=rbind(cbind(A,B),cbind(t(B),diag(dd))))
				info.matrix=A)
oldClass(retlist) <- c("lrm","lm")
retlist
}

