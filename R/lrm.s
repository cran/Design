lrm <- function(formula,data,subset,na.action=na.delete,
				method="lrm.fit",model=FALSE, x=FALSE, y=FALSE, 
				linear.predictors=TRUE, se.fit=FALSE, 
				penalty=0, penalty.matrix, tol=1e-7, strata.penalty=0,
                var.penalty=c('simple','sandwich'),
                weights, normwt=FALSE, ...) {

  call <- match.call()
  var.penalty <- match.arg(var.penalty)
  m <- match.call(expand=FALSE)
  m$method <- m$model <- m$x <- m$y <- m$... <- m$linear.predictors <- 
	m$se.fit <-	m$penalty <- m$penalty.matrix <- m$strata.penalty <- 
      m$tol <- m$var.penalty <- m$normwt <- NULL

  m$na.action <- na.action
  m$formula <- formula
  if(.R.) m$drop.unused.levels <- TRUE
  m[[1]] <- as.name("model.frame")
  nact <- NULL
  if(missing(data)) {
    data <- NULL
  }

  tform <- terms(formula, specials='strat', data=data)
  offs <- attr(tform, "offset")
  nstrata <- 1
  if(!missing(data) || (
						length(atl <- attr(tform,"term.labels")) && 
						any(atl!=".")))	{ ##X's present
    if(.R.) {
      dul <- .Options$drop.unused.levels
      if(!length(dul) || dul) {
        on.exit(options(drop.unused.levels=dul))
        options(drop.unused.levels=FALSE)
      }
    }
    X <- Design(eval(m, sys.parent()))
    atrx <- attributes(X)
    nact <- atrx$na.action
    if(method=="model.frame") return(X)
    Terms <- atrx$terms
    attr(Terms, "formula") <- formula
    atr <- atrx$Design
    if(length(nact$nmiss)) {
      jia <- grep('%ia%',names(nact$nmiss),fixed=TRUE)
      if(length(jia)) nact$nmiss <- nact$nmiss[-jia]
      s <- names(nact$nmiss) %nin% c(atrx$names[offs],'(weights)')
      names(nact$nmiss)[s] <-
        c(as.character(formula[2]), atr$name[atr$assume.code!=9])
    }

    Y <- model.extract(X, response)
    weights <- wt <- model.extract(X, 'weights')
    if(length(weights))
      warning('currently weights are ignored in model validation and bootstrapping lrm fits')
    ##   offset <- attr(X,"offset")
    ofpres <- length(offs) > 0
    if(ofpres) offs <- X[[offs]]  else offs <- 0
    if(!.R.)storage.mode(offs) <- "single"
    if(model) m <- X
    stra <- attr(tform,'specials')$strat
    Strata <- NULL
    Terms.ns <- Terms
    if(length(stra)) {
      temp <- untangle.specials(Terms.ns, 'strat', 1)
      Terms.ns <- Terms.ns[-temp$terms]
      attr(Terms,   "factors") <- pmin(attr(Terms,"factors"),1)
      attr(Terms.ns,"factors") <- pmin(attr(Terms.ns,"factors"),1)
      Strata <- X[[stra]]
      nstrata <- length(levels(Strata))
    }
    X <- model.matrix(Terms.ns, X)
    ##      ass <- attr(X, 'assign')
    ##	  X <- model.matrix(Terms.ns, X)[,-1,drop=FALSE]
    X <- X[,-1,drop=FALSE]
    if(!.R.)storage.mode(X) <- "single"
    dimnames(X)[[2]] <- atr$colnames
    ##	  oldClass(X) <- c("model.matrix", "Design")
    xpres <- length(X)>0

    p <- length(atr$colnames)
    n <- length(Y)

    penpres <- !(missing(penalty) && missing(penalty.matrix))
    if(penpres && missing(var.penalty))
      warning('default for var.penalty has changed to "simple"')

    if(!penpres) penalty.matrix <- matrix(0,ncol=p,nrow=p) else { 
      if(missing(penalty.matrix)) penalty.matrix <- Penalty.matrix(atr, X) else
      if(nrow(penalty.matrix)!=p || ncol(penalty.matrix)!=p) stop(
             paste("penalty.matrix does not have",p,"rows and columns"))
      psetup <- Penalty.setup(atr, penalty)
      penalty <- psetup$penalty
      multiplier <- psetup$multiplier
      if(length(multiplier)==1) penalty.matrix <- multiplier*penalty.matrix else {
        a <- diag(sqrt(multiplier))
        penalty.matrix <- a %*% penalty.matrix %*% a
      }
    }
  }
  else {
    X <- eval(m, sys.parent())
    ofpres <- FALSE
    if(length(offs)) {
      ofpres <- TRUE
      offs <- X[[offs]]
    }
    Y <- model.extract(X, response)
    Y <- Y[!is.na(Y)]
    Terms <- X <- NULL
    xpres <- FALSE
    penpres <- FALSE
    penalty.matrix <- NULL
  }  ##Model: y~. without data= -> no predictors
  
  if(method=="model.matrix") return(X)

  if(nstrata > 1) {  ## 17Dec97
	f <- if(ofpres) lrm.fit.strat(X,Y,Strata,offset=offs,
               penalty.matrix=penalty.matrix,strata.penalty=strata.penalty,
								  tol=tol,
                                  weights=weights,normwt=normwt, ...) else
	lrm.fit.strat(X,Y,Strata,penalty.matrix=penalty.matrix,
				  strata.penalty=strata.penalty,tol=tol,
                  weights=weights, normwt=normwt, ...)
  }		 else {
	if(existsFunction(method)) {
	  fitter <- getFunction(method)
	  if(ofpres) f <- fitter(X,Y,offset=offs,
							 penalty.matrix=penalty.matrix,tol=tol,
                             weights=weights, normwt=normwt, ...)
		else f <- fitter(X,Y,
						 penalty.matrix=penalty.matrix,tol=tol,
                         weights=weights, normwt=normwt, ...)
	}
	  else stop(paste("unimplemented method:", method))
  }
  f$call <- NULL
  if(model) f$model <- m
  if(x) {
	f$x <- X
	f$strata <- Strata
	}
  if(y) f$y <- Y
  nrp <- f$non.slopes
  if(penpres) {
	f$penalty <- penalty
	if(nstrata==1) {
	## Get improved covariance matrix
	v <- f$var
	if(var.penalty=='sandwich') f$var.from.info.matrix <- v
	f.nopenalty <- if(ofpres) fitter(X,Y,offset=offs,initial=f$coef, maxit=1, tol=tol) else
	fitter(X,Y,initial=f$coef, maxit=1, tol=tol)
	##  info.matrix.unpenalized <- solvet(f.nopenalty$var, tol=tol)
	info.matrix.unpenalized <- f.nopenalty$info.matrix
	dag <- diag(info.matrix.unpenalized %*% v)
	f$effective.df.diagonal <- dag
	f$var <- if(var.penalty=='simple')v else
             v %*% info.matrix.unpenalized %*% v
	df <- sum(dag[-(1:nrp)])
	lr <- f.nopenalty$stats["Model L.R."]
	pval <- 1-pchisq(lr, df)
	f$stats[c('d.f.','Model L.R.','P')] <- c(df, lr, pval)  
}
  }
    ass <- if(xpres) DesignAssign(atr, nrp, Terms) else list()

  if(xpres) {
	if(linear.predictors) names(f$linear.predictors) <- names(Y) else
	f$linear.predictors <- NULL

	if(se.fit) 		{
	if(nstrata > 1) stop('se.fit=T not implemented for strat')
	  xint <- matrix(0,nrow=length(Y),ncol=f$non.slopes)
	  xint[,1] <- 1
	  X <- cbind(xint, X)
	  se <- drop((((X %*% f$var) * X) %*% rep(1, ncol(X)))^.5)
	  if(!.R.)storage.mode(se) <- "single"
	  names(se) <- names(Y)
	  f$se.fit <- se	}
  }
  f <- c(f, list(call=call, Design=if(xpres)atr,
                 scale.pred=c("log odds","Odds Ratio"),
				 terms=Terms, assign=ass, na.action=nact, fail=FALSE,
                 nstrata=nstrata, fitFunction=c('lrm','glm')))
  
  oldClass(f) <- if(.SV4.)'Design' else c("lrm","Design","glm")
  f
}
