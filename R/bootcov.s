bootcov <- function(fit, cluster, B=200, fitter, coef.reps=FALSE, 
					loglik=coef.reps, pr=FALSE, maxit=15, group=NULL) {

  coxcph <- inherits(fit,'coxph') || inherits(fit,'cph') ||
  (length(fit$fitFunction) && any(c('cph','coxph') %in%
                                  fit$fitFunction))  ##14Nov00 22May01

  if(length(fit$weights) && coxcph)
    stop('does not handle weights') ##14Nov00

  if(!length(X <- fit$x) | !length(Y <- fit$y))
    stop("you did not specify x=T and y=T in the fit")

  nfit <- fit$fitFunction[1]
  if(!length(nfit)) nfit <- setdiff(oldClass(fit),'Design')[1]  ##14Nov00

  sc.pres <- if(.newSurvival.) match('scale',names(fit),0)>0 else
                               match("parms",names(fit),0)>0

  if(nfit=='psm') {
    fixed <- fit$fixed   #psm only
    fixed <- if(length(fixed)==1 && is.logical(fixed) && !fixed) list() else
             list(scale=TRUE)
    if(.R.) {
      fixed <-  NULL
      dist <-   fit$dist
      parms <-  fit$parms
      Strata <- NULL
    } else if(.newSurvival.) {
      storeTemp(NULL,      'fixed')            
      storeTemp(fit$parms, 'parms')
      storeTemp(fit$dist,  'dist')
      storeTemp(NULL,      'Strata')
    } else {
      storeTemp(NULL,  'parms')
      storeTemp(fixed, 'fixed')
      storeTemp(fit$family['name'], 'dist')
      storeTemp(NULL,  'Strata')
    }
  }

  if(nfit=='cph') {
   if(.R.) Strata <- NULL else storeTemp(NULL, 'Strata')
 }

  if(nfit=='glmD') {  ## was glm 2dec02
    if(.R.) fitFamily <- fit$family else 
    storeTemp(origGlmFamily(fit$family), 'fitFamily')
  }
  
  penalty.matrix <- fit$penalty.matrix


  if(missing(fitter))
    fitter <-
      switch(nfit,
             ols=if(length(penalty.matrix))
             function(x,y,penalty.matrix,...)
             lm.pfit(x,y,penalty.matrix=penalty.matrix,tol=1e-11,
                     regcoef.only=TRUE) else
             function(x,y,...)
             lm.fit.qr.bare(x,y,tolerance=1e-11,intercept=FALSE), 
             lrm=function(x,y,maxit=15,penalty.matrix,...)
             lrm.fit(x,y,maxit=maxit,tol=1E-11,
                     penalty.matrix=penalty.matrix), 
             cph=function(x,y,maxit=15,...)coxphFit(x,y,
                                strata=Strata, iter.max=maxit, 
                                eps=.0001, method="efron", toler.chol=1e-11),
             psm= function(x,y,maxit=15,...) survreg.fit2(x, y,
                                 dist=dist, parms=parms, fixed=fixed,
                                 offset=NULL,
                                 init=NULL, maxiter=maxit),
             bj=function(x,y,maxit=15,eps=.0001,...)
             bj.fit(x,y,
                    control=list(iter.max=maxit,eps=1e-4)),
             glmD=function(x,y,...) glm.fit(x,as.vector(y),
               family=fitFamily)  ## 24nov02 02dec02
             )

  ## psm previously used (20Apr02)
  ##function(x,y,maxit=15,...) survreg.fit(x, y, 
  ## dist=dist, fixed=fixed, offset=NULL, init=NULL,  controlvals=
  ## survreg.control(maxiter=maxit,rel.tol=.0001,failure=2))
  ## survreg.fit2 is defined in validate.psm.s.  It ignores unneeded
  ## parms vs fixed
	
if(!length(fitter))stop("fitter not valid")

if(loglik) {
  oosl <- switch(nfit,   ##14Nov00
				 ols=oos.loglik.ols,
				 lrm=oos.loglik.lrm,
				 cph=oos.loglik.cph,
				 psm=oos.loglik.psm,
                 glmD=oos.loglik.glmD)  ## 6dec02
  if(!length(oosl)) 
	stop('loglik=T but no oos.loglik method for model in Design.Misc')
  Loglik <- if(.R.) double(B+1) else single(B+1)
  Loglik[B+1] <- oosl(fit)
} else Loglik <- NULL


n <- nrow(X)
p <- length(fit$coef)
vname <- names(fit$coef)
if(sc.pres) {p <- p+1; vname <- c(vname, "log scale")}
bar <- rep(0, p)
cov <- matrix(0, nrow=p, ncol=p, dimnames=list(vname,vname))
if(coef.reps) coefs <- matrix(NA, nrow=B, ncol=p, dimnames=list(NULL,vname))

Y <- as.matrix(if(is.category(Y)) oldUnclass(Y) else Y)  ##25Mar98
ny <- ncol(Y)

str.pres <- FALSE
if(inherits(fit,"cph") || (length(fit$fitFunction) &&
                           any(fit$fitFunction=='cph')))		{
  ##14Nov00 22May01
   str.pres <- TRUE
   str <- attr(X, "strata")
   str.pres <- length(str)	}

nac <- fit$na.action

if(length(group)) {
##  if(!missing(cluster))  15nov02
##    stop('group is currently allowed only when cluster is not given')
  if(length(group) > n) {
	## Missing observations were deleted during fit
	if(length(nac)) {
	  j <- !is.na(naresid(nac, Y) %*% rep(1,ny))
	  group <- group[j]
	}
  }
  if(length(group) != n)
	stop('length of group does not match # rows used in fit')
  group.inds <- split(1:n, group)  ## see bootstrap()
  ngroup <- length(group.inds)
} else ngroup <- 0
 
if(missing(cluster))						{
  b <- 0
  for(i in 1:B)				{
     if(pr)cat(i,"")
     if(ngroup) {
        j <- integer(n)
        for(si in 1:ngroup) {
          gi <- group.inds[[si]]
          j[gi] <- sample(gi, length(gi), replace=TRUE)
        }
      } else j <- sample(1:n, n, replace=TRUE)
     if(str.pres) {
       if(.R.) Strata <- str[j] else assign("Strata", str[j], 1)
     }
     f <- fitter(X[j,,drop=FALSE], Y[j,,drop=FALSE], maxit=maxit, 
				 penalty.matrix=penalty.matrix)
     if(length(f$fail) && f$fail) next
     b <- b+1
     cof <- as.vector(f$coef)
#    if(sc.pres) cof <- c(cof, f$parms[1])   #survreg.fit already does this
     if(sc.pres && (.newSurvival.)) cof <- c(cof, log(f$scale))
     if(coef.reps) coefs[b,] <- cof
     bar <- bar + cof
     cof <- as.matrix(cof)
     cov <- cov + cof %*% t(cof)
	 if(loglik) Loglik[b] <- oosl(f, matxv(X,cof), Y)
   }
								}

else								{
   if(length(cluster) > n)		{
      #Missing obs were deleted during fit
      if(length(nac))	{
        j <- !is.na(naresid(nac, Y) %*% rep(1,ny))
        cluster <- cluster[j]
				}	}

   if(length(cluster)!=n) stop("length of cluster does not match # rows used in fit")
   if(any(is.na(cluster))) stop("cluster contains NAs")

   cluster <- as.character(cluster)

   clusters <- unique(cluster)
   nc <- length(clusters)
##   nx <- ncol(X)
   Obsno <- split(1:n, cluster)

## Spread row names along with every column
##   q <- as.vector(t(matrix(rep(cluster, nx), ncol=nx)))
##   X <- split(as.vector(t(X)), q)
##   q <- as.vector(t(matrix(rep(cluster, ny), ncol=ny)))
##   Y <- split(as.vector(t(Y)), q)
##   if(str.pres) str <- split(str, cluster)

   b <- 0
   for(i in 1:B)			{
      if(pr)cat(i,"")
      ## Begin addition Bill Pikounis 15nov02 (work done 1nov99)
      if(ngroup) {
        j <- integer(0)
        for(si in 1:ngroup) {
          gi <- group.inds[[si]]
          cluster.gi <- cluster[gi]
          clusters.gi <- unique(cluster.gi)
          nc.gi <- length(clusters.gi)
          Obsno.gci <- split(gi, cluster.gi)
          j.gci <- sample(clusters.gi, nc.gi, replace = T)
          obs.gci <- unlist(Obsno.gci[j.gci])
          jadd <- c(j, jadd)
        }
        obs <- j
      }
      else {
      ## End addition Bill Pikounis (except for closing brace below)
        j <- sample(clusters, nc, replace=TRUE)
        ##      if(str.pres) assign("Strata", unlist(str[j]), 1)
        ##      x <- matrix(unlist(X[j]), ncol=nx, byrow=T)
        ##      y <- matrix(unlist(Y[j]), ncol=ny, byrow=T)
        obs <- unlist(Obsno[j])
      }   ## 15nov02
      if(str.pres) {
        if(.R.) Strata <- str[obs] else assign("Strata", str[obs], 1)
      }
      ##    f <- fitter(x,y,maxit=maxit)
      f <- fitter(X[obs,,drop=FALSE], Y[obs,,drop=FALSE], 
                  maxit=maxit, penalty.matrix=penalty.matrix)
      ## added ,drop 9Oct01
      if(length(f$fail) && f$fail) next
      b <- b+1
      cof <- as.vector(f$coef)
      ##     if(sc.pres) cof <- c(cof, f$parms[1])
      if(sc.pres && (.newSurvival.)) cof <- c(cof, log(f$scale))
      if(coef.reps) coefs[b,] <- cof
      bar <- bar + cof
      cof <- as.matrix(cof)
      cov <- cov + cof %*% t(cof)
      if(loglik) Loglik[b] <- oosl(f, matxv(X,cof), Y)
	}
 }
  if(b < B) {  # 21Apr02
    warning(paste('fit failure in',B-b,
                  'resamples.  Might try increasing maxit'))
    if(coef.reps) coefs <- coefs[1:b,,drop=FALSE]
    Loglik <- Loglik[1:b]
  }
  bar <- bar/b
  names(bar) <- vname
  fit$boot.coef <- bar
  if(coef.reps) fit$boot.Coef <- coefs
  bar <- as.matrix(bar)
  cov <- (cov - b * bar %*% t(bar))/(b-1)
  fit$orig.var <- fit$var
  fit$var <- cov
  fit$boot.loglik <- Loglik
  ##oldClass(fit) <- c("bootcov",oldClass(fit))    14Nov00
  fit
}

#bootplot <- function(obj, ...) UseMethod('bootplot')  14Nov00
#confplot <- function(obj, ...) UseMethod('confplot')

bootplot <- function(obj, which, X,
                     conf.int=c(.9,.95,.99),
                     what=c('density','qqnorm'),
                     fun=function(x)x,
                     labels., ...) {

what <- match.arg(what)
Coef <- obj$boot.Coef
if(length(Coef)==0) stop('did not specify "coef.reps=T" to bootcov')
if(missing(which)) {
  if(!is.matrix(X)) X <- matrix(X, nrow=1)
  qoi <- X %*% t(Coef)   ##nxp pxB = nxB
  if(missing(labels.)) {
	labels. <- dimnames(X)[[1]]
	if(length(labels.)==0) labels. <- as.character(1:nrow(X))
  }
} else {
  qoi <- t(Coef[,which,drop=FALSE])
  nns <- num.intercepts(obj)
  if(missing(labels.))
	labels. <- paste(ifelse(which > nns, 'Coefficient of ', ''), 
					 dimnames(Coef)[[2]][which], sep='')
}
nq <- nrow(qoi)
qoi <- fun(qoi)
quan <- NULL
if(what=='density') {
  probs <- (1+conf.int)/2
  probs <- c(1-probs, probs)
  quan <- matrix(NA, nrow=nq, ncol=2*length(conf.int),
				 dimnames=list(labels., format(probs)))
  for(j in 1:nq) {
	histdensity(qoi[j,], xlab=labels.[j], ...)
	quan[j,] <- quantile(qoi[j,],probs)
	abline(v=quan[j,], lty=2)
	title(sub=paste('Fraction of effects>',fun(0),' = ',
			format(mean(qoi[j,]>fun(0))),sep=''),adj=0)
  }
} else {
  for(j in 1:nq) {
	qqnorm(qoi[j,], ylab=labels.[j])
	qqline(qoi[j,])
  }
}

invisible(list(qoi=drop(qoi), quantiles=drop(quan)))
}


## histdensity runs hist() and density(), using twice the number of
## class than the default for hist, and 1.5 times the width than the default
## for density

histdensity <- function(y, xlab, nclass, width, mult.width=1, ...) {
  y <- y[is.finite(y)]
  if(missing(xlab)) {
	xlab <- label(y)
	if(xlab=='') xlab <- as.character(sys.call())[-1]
  }
  if(missing(nclass)) nclass <- (logb(length(y),base=2)+1)*2
  hist(y, nclass=nclass, xlab=xlab, probability=TRUE, ...)
  if(missing(width)) {
	nbar <- logb(length(y), base = 2) + 1
	width <- diff(range(y))/nbar*.75*mult.width
  }
  lines(density(y,width=width,n=200))
  invisible()
}


confplot <- function(obj, X, against, 
                     method=c('simultaneous','pointwise'),
                     conf.int=0.95,
                     fun=function(x)x, 
                     add=FALSE, lty.conf=2, ...) {

method <- match.arg(method)
if(length(conf.int)>1) stop('may not specify more than one conf.int value')

boot.Coef <- obj$boot.Coef
if(length(boot.Coef)==0) stop('did not specify "coef.reps=T" to bootcov')

if(!is.matrix(X)) X <- matrix(X, nrow=1)
fitted <- fun(X %*% obj$coefficients)

if(method=='pointwise') {
  pred <- X %*% t(boot.Coef)   ## n x B
  p <- fun(apply(pred, 1, quantile, probs=c((1-conf.int)/2,1-(1-conf.int)/2)))
  lower <- p[1,]
  upper <- p[2,]
} else {
  boot.Coef <- rbind(boot.Coef, obj$coefficients)
  loglik    <- obj$boot.loglik
  if(length(loglik)==0) stop('did not specify "loglik=T" to bootcov')
  crit <- quantile(loglik, conf.int)
  qual <- loglik <= crit
  boot.Coef <- boot.Coef[qual,,drop=FALSE]
  pred   <- X %*% t(boot.Coef)  ## n x B
  upper  <- fun(apply(pred, 1, max))
  lower  <- fun(apply(pred, 1, min))
  pred   <- fun(pred)
}

if(!missing(against)) {
  lab <- label(against)
  if(lab=='') lab <- (as.character(sys.call())[-1])[3]
  if(add) lines(against, fitted, ...) else
  plot(against, fitted, xlab=lab, type='l', ...)
  lines(against, lower, lty=lty.conf)
  lines(against, upper, lty=lty.conf)
}
if(missing(against))list(fitted=fitted, upper=upper, lower=lower) else
invisible(list(fitted=fitted, upper=upper, lower=lower))
}

## 24nov02
if(!.R.) { 
  origGlmFamily <- function(glmfitFamily) {
    ## S-Plus glm.fit stores only first component of
    ## as.family(family) in fit object, but
    ## it hands all of as.family(family) to glm.fit
    familyname <- casefold(glmfitFamily['name'])
    link <- casefold(unpaste(glmfitFamily['link'],':')[[1]])
    eval(parse(text=paste(familyname,'(',link,')',sep='')))
  }
  NULL
}


