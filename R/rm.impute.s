# Start with first visit containing any NAs
		
rm.impute <- function(pformula, y, last, rformula, fitter=ols,
					  which=c('last','mean','auc'),
					  data=sys.parent(1), n.impute=10, g=5, 
					  nk=0, rinteraction, rint.with=c('all','recent'),
					  pr=FALSE, pra=FALSE, npr=1, keep.prop=FALSE, keep.pfits=FALSE) {

rint.with <- match.arg(rint.with)
warning('This is an experimental procedure.  It should only be used for testing, as results are incorrect.')

d <- dim(y)
if(length(d) < 3) {
  y <- if(length(d) < 2) array(y, c(length(y),1,1),
							   list(names(y),NULL,NULL)) else
  array(y, c(d,1), list(dimnames(y)[[1]],dimnames(y)[[2]],NULL))
}
yat <- attributes(y)
d <- yat$dim
nr <- d[3]
nt <- d[2]
n  <- d[1]

if(length(last) != n) stop('length of last != nrow(y)')
which <- match.arg(which)
dimn  <- yat$dimnames
dimn[[4]] <- paste('Imputation',1:n.impute)
Y <- array(NA, c(n, nt, nr, n.impute), dimn)
Rvar <- dimn[[3]]
if(!length(Rvar)) Rvar <- if(nr==1) 'y' else paste('y',1:nr,sep='')

Propensity <- if(keep.prop) array(NA, c(n, nt, n.impute), dimn[c(1,2,4)])

pfits <- if(keep.pfits) vector('list',n.impute*nt)
if(!missing(rformula)) {
  fr    <- vector('list', nr)
  vavg <- bar <- cov <- vector('list',nr)
}

if(!missing(rinteraction) && length(rinteraction)>1)
  rinteraction <- paste('(',paste(rinteraction,collapse='+'),')',sep='')

kk <- 0
for(imp in 1:n.impute) {
  cat(if(pr | imp==1)'\n\nImputation',imp,if(pr)'\n-------------\n')

  form     <- update(pformula, in.period.i ~ ., evaluate=FALSE)  
  ## add/change response variable
  formbase <- form

  for(i in 1:nt) {
	kk <- kk + 1
	in.period.i <- last >= i
	if(imp==1) {
	  for(k in 1:nr) {
		w <- !is.na(y[,i,k]) != in.period.i
		if(any(w)) {
		  cat('Value of last disagrees with missingness of ',
			  Rvar[k],' in period ',i, ' for ', sum(w),
			  ' subjects\n\n', sep='')
		  print(table(ifelse(is.na(y[,i,k]),'Response NA','not NA'),
					  ifelse(in.period.i,'Subject in study','dropped out')))
		}
	  }
	}

	if(!.R.) storeTemp(in.period.i)
    # assign('in.period.i', in.period.i, frame=0) 17Apr01
	frm <- form
	if(i > 1) {
	  trvar <- rvar
	  if(nk > 0) trvar <- paste('rcs(',rvar,',',nk,')',sep='')
	  if(nr > 1) trvar <- paste('(',paste(trvar,collapse='+'),')',sep='')
	  if(!missing(rinteraction)) {
		trvar <- paste(trvar,'*',rinteraction)
		if(i > 2 && rint.with=='all') {
		  trvar <- paste(trvar,'+',rinteraction,'*',if(nr>1 | i>3)'(')
		  for(k in 1:nr) trvar <- paste(trvar, if(k>1) '+',
	   if(nk==0) paste(paste(Rvar[k],'.',1:(i-2),sep=''),collapse='+') else
				 paste(w1 <- paste('rcs(',Rvar[k],'.',1:(i-2),',',nk,')',sep=''),
					   collapse='+'))
		  if(nr>1 | i>3) trvar <- paste(trvar,')')
		  frm <- formbase
		}
	  }
	form <- update(frm, paste('~. +', trvar), evaluate=FALSE)
  }
	if(all(in.period.i)) {
	  if(imp==1) cat('\nTime period',i,': no dropouts\n')
	} else {
		prop.fit      <- lrm(form, data=data)
		if(prop.fit$fail) 
		  stop(paste('propensity model failed to converge for response',i,
					 'imputation',imp))
		if(keep.pfits) pfits[[kk]] <- prop.fit
		if(pr && imp <= npr) {
		  cat('\nTime period',i,'propensity model\n\n')
		  dput(form); cat('\n')
		  print(prop.fit)
		  if(pra) print(anova(prop.fit))
		}
		propensity    <- predict(prop.fit, type='fitted')
		propensity[propensity < 1e-10] <- 0 
		## only needed because bug in cut bombs cut2 below
		if(keep.prop) Propensity[,i,imp] <- propensity
		prop.quantile <- cut2(propensity, g=g)
		if(pr && imp <= npr) {
		  cat('\nFrequencies of propensity quantile groups by dropout status\n\n')
		  print(table(c('In study','Dropped out')[2-in.period.i], 
					  prop.quantile))
		}
		## Fill-in y values for current time corresponding to subjects
		## who dropped out before the current time, by sampling with
		## replacement from a sample with replacement (Rubin approx.
		## Bayes bootstrap) from non-dropout values within the same
		## propensity quantile
		for(pg in levels(prop.quantile)) {
		  s <- prop.quantile==pg
		  needed <- s & !in.period.i
		  avail  <- s & in.period.i
		  if(any(needed)) {
			if(!any(avail))
			  stop(paste(sum(needed),'imputations needed in propensity group',
						 k,'but no responses available; response',i,'imputation',
						 imp))
			indices <- sample(sample((1:n)[avail],sum(avail),rep=TRUE),
							  sum(needed),rep=TRUE)
			for(k in 1:nr) y[needed,i,k] <- y[indices,i,k]
		  }
		}
	  }
	rvar <- character(nr)
	for(k in 1:nr) {
	  rvar[k] <- paste(Rvar[k], '.', i, sep='')
	  #assign(rvar[k], y[,i,k,drop=TRUE], frame=0) 17Apr01
      storeTemp(y[,i,k,drop=TRUE], rvar[k])
	}
  }
  Y[,,,imp] <- y
  if(!missing(rformula)) {
	if(which=='auc') {
	  times <- as.numeric(dimnames(y)[[2]])
	  if(length(times)==0 | any(is.na(times)))
		stop('To use which="auc" y must have column names containing times')
	  mult <- double(nt)
	  for(j in 1:nt) mult[j] <- if(j==1) times[2]-times[1] else
	  if(j==nt) times[nt]-times[nt-1] else times[j+1]-times[j-1]
	}
	for(k in 1:nr) {
	  y. <- switch(which,
				   last = y[,nt,k,drop=TRUE],
				   mean = apply(y[,,k,drop=TRUE], 2, mean, na.rm=TRUE),
				   auc  = y[,,k,drop=TRUE] %*% mult/2)
	  
#	  assign(yy <- paste(Rvar[k],'.',sep=''), y., frame=0)  17Apr01
      yy <- paste(Rvar[k],'.',sep='')
      storeTemp(y., yy)
	  if(imp==1) rformula <- update(rformula, paste(yy,' ~ .'), evaluate=FALSE)
	  frk <- if(is.list(fitter)) fitter[[k]](rformula, data=data) else
	  fitter(rformula, data=data)
	  fr[[k]] <- frk
	  cof <- frk$coef
	  v <- Varcov(frk)
	  if(imp==1) {
		vavg[[k]] <- 0*v
		p <- length(frk$coef)
		bar[[k]] <- rep(0, p)
		vname <- names(frk$coef)
		cov[[k]] <- matrix(0, nrow=p, ncol=p, dimnames=list(vname,vname))
		if(!inherits(frk,'Design'))
		  warning('Not using a Design fitting function; summary(fit) will use\nstandard errors, t, P from last imputation only.  Use Varcov(fit) to get the\ncorrect covariance matrix, sqrt(diag(Varcov(fit))) to get s.e.\n\n')
	  }
	  vavg[[k]] <- vavg[[k]] + v
	  bar[[k]] <- bar[[k]] + cof
	  cof <- as.matrix(cof)
	  cov[[k]] <- cov[[k]] + cof %*% t(cof)
	}
  }
}

if(keep.prop) mode(Propensity) <- 'single'
if(keep.pfits) {
  dim(pfits) <- c(nt, n.impute)
  dimnames(pfits) <- dimnames(Y)[c(2,4)]
}

if(!missing(rformula)) {
  for(k in 1:nr) {
	vavgk <- vavg[[k]] / n.impute
	bark  <- bar[[k]]/n.impute
	bark  <- as.matrix(bark)
	covk  <- (cov[[k]] - n.impute * bark %*% t(bark))/(n.impute-1)
	covk  <- vavgk + (n.impute+1)/n.impute * covk
	r <- diag(covk) / diag(vavgk)
	names(r) <- vname
	cat('\nVariance Inflation Factors Due to Imputation for ',Rvar[k],
		':\n\n')
	print(round(r,2))
	frk <- fr[[k]]
	frk$coefficients <- drop(bark)
	frk$var <- covk
	frk$variance.inflation.impute <- r
	oldClass(frk) <- c('fit.mult.impute',oldClass(fr[[k]]))
	fr[[k]] <- frk
  }
  list(fit=if(nr==1)fr[[1]] else fr, 
	   Y=drop(Y), propensity=Propensity, pfits=pfits)
} else list(Y=drop(Y), propensity=Propensity, pfits=pfits)
}

pbind <- function(...) {
  dotlist <- list(...)
  m1 <- as.matrix(dotlist[[1]])
  d <- dim(m1)
  nam <- names(dotlist)
  if(!length(nam)) nam <- as.character(sys.call())[-1]
  array(unlist(dotlist), c(d,length(dotlist)),
		list(dimnames(m1)[[1]], dimnames(m1)[[2]], nam))
}



