lrm.fit.strat <- function(x,y,strata,offset,initial,
						  maxit=25,eps=.025,tol=1E-7,trace=FALSE,
						  penalty.matrix=NULL,strata.penalty=0){	
cal <- match.call()
opts <- double(11)
len.penmat <- length(penalty.matrix)
lev    <- levels(strata)
nstrat <- length(lev)
strata <- oldUnclass(strata)

n <- length(y)
initial.there <- !missing(initial)
if(missing(x) || length(x)==0)	{
	nx <- 0
	xname <- NULL
	x <- 0	}	else				{
    if(!is.matrix(x)) x <- as.matrix(x)
	storage.mode(x) <- if(.R.)"double" else "single"
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
	storage.mode(offset) <- if(.R.)"double" else "single"
		} else offset <- 0

if(n<3)stop("must have >=3 non-missing observations")
kint <- as.integer(length(ylevels)-1)
if(kint!=1) stop('only works for binary y')
ftable <- integer(501*(kint+1))
levels(y) <- ylevels
numy <- table(y)
y <- as.integer(y-1)
nvi <- as.integer(nxin+kint+nstrat-1)
if(missing(initial))	{
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

if(nxin==0 & !ofpres)	{
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
               penmat, PACKAGE="Design") else
      .Fortran("lrmfit",coef=initial,as.integer(0),0,x,y,offset,
               u=double(kint),
               double(kint*(kint+1)/2),loglik=double(1),n,as.integer(0),
               numy,kint,
               v=double(kint*kint),double(kint),double(kint),
               double(kint),pivot=integer(kint),opts=opts,ftable,
               penmat, PACKAGE="Design")
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

##  Ci <- diag(1/dd)
##  BCi <- B %*% Ci
  BCi <- B %*% diag(1/dd)
  v <- 1/(dd + strata.penalty)
  vm <- as.matrix(v)
  k <- (strata.penalty/nstrat)/(1 - (strata.penalty/nstrat)*sum(v))
  BCi <- t(v * t(B)) - k * (B %*% vm) %*% t(vm)
  AA <- solvet(A - BCi %*% t(B), tol=tol)
  BB <- -AA %*% BCi
#  info <- rbind(cbind(A,B),cbind(t(B),diag(dd)))
#  info.inv <- matrix(NA, nrow=nvi,ncol=nvi)
#  info.inv[(nns+1):nvi,(nns+1):nvi] <- Ci - t(BB) %*% BCi
#  info.inv[1:nns,1:nns] <- AA
#  info.inv[1:nns,(nns+1):nvi] <- BB
#  info.inv[(nns+1):nvi,1:nns] <- t(BB)
#  theta <- theta + solvet(info, u)
#  theta <- theta + info.inv %*% u
  u1 <- u[1:nns,,drop=FALSE]
  u2 <- u[(nns+1):nvi,,drop=FALSE]
  theta <- theta + c(AA %*% u1 + BB %*% u2,
					 t(BB) %*% u1 + vm * u2 - k * vm %*% (t(vm) %*% u2) -
					   t(BB) %*% (BCi %*% u2))
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
v <- solvet(A - BCi %*% t(B), tol=tol)
dimnames(v) <- list(name[1:nns],name[1:nns])
llnull <- loglik[length(loglik)-1]
model.lr <- llnull-loglik[length(loglik)]
model.df <- nvi - kint
if(initial.there) model.p <- NA 	else 				{
if(model.df>0) model.p <- 1-pchisq(model.lr,model.df) else model.p <- 1	}
r2 <- 1-exp(-model.lr/n)
r2.max <- 1-exp(-llnull/n)
r2 <- r2/r2.max
##lp <- matxv(x, z$parameters[1:(nx+kint)]) + 
##  c(0,z$parameters[-(1:(nx+kint))])[strata]
Brier <- mean((pred - (y>0))^2)

stats <- c(n,max(abs(u)),model.lr,model.df,model.p,
	## z$opts[8],z$opts[9],z$opts[10], z$opts[11], 
		   r2, Brier)
nam <- c("Obs","Max Deriv",	"Model L.R.","d.f.","P",
		 ##"C","Dxy","Gamma","Tau-a",
		 "R2","Brier")

names(stats) <- nam
retlist <- list(call=cal,freq=numy,
				stats=stats,fail=FALSE,coefficients=theta,
				non.slopes=1,est=1:(nvi-kint),
				var=v,u=u,
				deviance=loglik,
				linear.predictors=logit,
				penalty.matrix=if(nxin>0 && any(penalty.matrix!=0)) 
				  penalty.matrix else NULL,
				strata.penalty=strata.penalty, 
#				info.matrix=rbind(cbind(A,B),cbind(t(B),diag(dd))))
				info.matrix=A)
oldClass(retlist) <- 'lrm' # was c("lrm","lm") 17Jul01
retlist
}

