ols <- function(formula, data, weights, subset, na.action=na.delete, 
                method = "qr", model = FALSE, x = FALSE, y = FALSE,
                se.fit=FALSE, linear.predictors=TRUE,
                penalty=0, penalty.matrix, tol=1e-7, sigma=NULL,
                var.penalty=c('simple','sandwich'), ...){
  call <- match.call()
  var.penalty <- match.arg(var.penalty)
  m <- match.call(expand = FALSE)
  m$method <- m$model <- m$x <- m$y <- m$se.fit <-
    m$linear.predictors <- m$penalty <- 
      m$penalty.matrix <- m$tol <- m$sigma <- m$var.penalty <- m$... <- NULL
  m$na.action <- na.action
  if(.R.) m$drop.unused.levels <- TRUE
  m[[1]] <- as.name("model.frame")
  ##X's present) 
  if(length(attr(terms(formula),"term.labels"))) {
    ## R's model.frame.default gives wrong model frame if [.factor
    ## removes unused factor levels
    if(.R.) {
      dul <- .Options$drop.unused.levels
      if(!length(dul) || dul) {
        on.exit(options(drop.unused.levels=dul))
        options(drop.unused.levels=FALSE)
      }
    }
    X <- Design(eval(m, sys.parent()))
    if(.R.) options(drop.unused.levels=dul)
    atrx <- attributes(X)
    atr <- atrx$Design
	nact <- atrx$na.action
	Terms <- atrx$terms
    assig <- DesignAssign(atr, 1, Terms)
    
    penpres <- FALSE
    if(!missing(penalty) && any(unlist(penalty) != 0)) penpres <- TRUE
    if(!missing(penalty.matrix) && any(penalty.matrix != 0)) penpres <- TRUE
    
    if(penpres && missing(var.penalty))
      warning('default for var.penalty has changed to "simple"')
    
	if(method == "model.frame") return(X)
	scale <- as.character(formula[2])
	attr(Terms, "formula") <- formula
	if(length(nact$nmiss)) {
      jia <- grep('%ia%',names(nact$nmiss), fixed=TRUE)
      if(length(jia)) nact$nmiss <- nact$nmiss[-jia]
      names(nact$nmiss)[names(nact$nmiss)!='(weights)'] <- 
		c(scale,atr$name[atr$assume.code!=9])
    }
	weights <- model.extract(X, weights)
	if(length(weights) && penpres)
	  stop('may not specify penalty with weights')

	Y <- model.extract(X, response)
	n <- length(Y)
	if(model) m <- X
	X <- model.matrix(Terms, X)
	if(length(atr$colnames)) 
      dimnames(X)[[2]] <- c("Intercept",atr$colnames)
	else dimnames(X)[[2]] <- c("Intercept",dimnames(X)[[2]][-1])
	if(method=="model.matrix") return(X)				   }

  ##Model with no covariables:

  else {
    if(length(weights))
      stop('weights not implemented when no covariables are present')
    assig <- NULL
    yy <- attr(terms(formula),"variables")[1]
    Y <- eval(yy,sys.parent(2))
    nmiss <- sum(is.na(Y))
    if(nmiss==0) nmiss <- NULL else names(nmiss) <- as.character(yy)
    Y <- Y[!is.na(Y)]
    yest <- mean(Y)
    coef <- yest
    n <- length(Y)
    if(!length(sigma)) sigma <- sqrt(sum((Y-yest)^2)/(n-1))
    cov <- matrix(sigma*sigma/n, nrow=1, ncol=1,
                  dimnames=list("Intercept","Intercept"))
    fit <- list(coefficients=coef, var=cov,
                non.slopes=1, fail=FALSE, residuals=Y-yest,
                df.residual=n-1, intercept=TRUE)
    if(linear.predictors) {
      fit$linear.predictors <- rep(yest,n); 
      names(fit$linear.predictors) <- names(Y)
    }
    if(model) fit$model <- m
    if(x) fit$x <- matrix(1, ncol=1, nrow=n, 
                          dimnames=list(NULL,"Intercept"))
    if(y) fit$y <- Y
    fit$fitFunction <- c('ols','lm')
    oldClass(fit) <- if(.SV4.)"Design" else c("ols","Design","lm")
    return(fit)
  }

  if(!penpres) {
    fit <- if(length(weights))
      lm.wfit(X, Y, weights, method=method,  ...) else 
    lm.fit(X, Y, method=method, ...)
    if(.R.) cov.unscaled <- chol2inv(fit$qr$qr) else {
      rinv <- solve(fit$R, diag(length(fit$coefficients)))
      cov.unscaled <- rinv %*% t(rinv)
    }
    r <- fit$residuals
    if(length(weights)) { ## see summary.lm
      sse <- sum(weights * r^2)
      yhat <- Y - r
      m <- sum(weights * yhat / sum(weights))
      ssr <- sum(weights * (yhat - m)^2)
      r2 <- ssr / (ssr + sse)
      if(!length(sigma)) sigma <- sqrt(sse/fit$df.residual)
    } else {
      sse <- sum(fit$residuals^2)
      if(!length(sigma)) sigma <- sqrt(sse/fit$df.residual)
      r2 <- 1-sse/sum((Y-mean(Y))^2)
    }
    fit$var <- sigma*sigma*cov.unscaled
    cnam <- dimnames(X)[[2]]
    dimnames(fit$var) <- list(cnam, cnam)
    fit$stats <- c(n=n,'Model L.R.'=-n*logb(1-r2),
                   'd.f.'=length(fit$coef)-1,R2=r2,Sigma=sigma)
  } else {
    p <- length(atr$colnames)
    if(missing(penalty.matrix)) penalty.matrix <- Penalty.matrix(atr,X)
    if(nrow(penalty.matrix)!=p || ncol(penalty.matrix)!=p) 
      stop('penalty matrix does not have',p,'rows and columns')
    psetup <- Penalty.setup(atr, penalty)
    penalty <- psetup$penalty
    multiplier <- psetup$multiplier
    if(length(multiplier)==1) penalty.matrix <- multiplier*penalty.matrix
    else {
      a <- diag(sqrt(multiplier))
      penalty.matrix <- a %*% penalty.matrix %*% a
    }
    fit <- lm.pfit(X, Y,
                   penalty.matrix=penalty.matrix, tol=tol,
                   var.penalty=var.penalty)
    fit$penalty <- penalty
  }
  
  if(model)
    fit$model <- m
  if(linear.predictors) {
    fit$linear.predictors <- Y-fit$residuals
    names(fit$linear.predictors) <- names(Y)
  }
  if(x)
    fit$x <- X
  if(y)
    fit$y <- Y
  if(se.fit) {
    se <- drop((((X %*% fit$var) * X) %*% rep(1, ncol(X)))^0.5)
    if(!.R.) storage.mode(se) <- "single"
    names(se) <- names(Y)
    fit$se.fit <- se
  }
  fit <- c(fit, list(call=call, terms=Terms, Design=atr,
                     non.slopes=1, na.action=nact,
                     scale.pred=scale, fail=FALSE,
                     fitFunction=c('ols','lm')))
  fit$assign <- assig
  oldClass(fit) <- if(.SV4.)'Design' else c("ols","Design","lm")
  fit
}


lm.pfit <- function(X, Y, penalty.matrix, tol=1e-7, regcoef.only=FALSE,
                    var.penalty=c('simple','sandwich')) {

  var.penalty <- match.arg(var.penalty)
  p <- ncol(X)-1
  pm <- rbind(matrix(0, ncol=p+1, nrow=1),
              cbind(matrix(0, ncol=1, nrow=p), penalty.matrix))
  xpx <- t(X) %*% X
  Z <- solvet(xpx+pm, tol=tol)
  coef <- Z %*% t(X) %*% Y
  if(regcoef.only) return(list(coefficients=coef))
  res  <- drop(Y - X %*% coef)
  n <- length(Y)
  sse <- sum(res^2)
  s2 <- drop( (sse + t(coef) %*% pm %*% coef) / n )
  var <- if(var.penalty=='simple') s2 * Z else s2 * Z %*% xpx %*% Z
  cnam <- dimnames(X)[[2]]
  dimnames(var) <- list(cnam, cnam)
  sst <- sum((Y-mean(Y))^2)
  lr <- n*(1+logb(sst/n))-n*logb(s2)-sse/s2
  s2.unpen <- sse/n
  dag <- diag((xpx / s2.unpen) %*% (s2 * Z))
  df <- sum(dag) - 1
  stats <- c(n=n, 'Model L.R.'=lr, 'd.f.'=df, R2=1-sse/sst, Sigma=sqrt(s2))
  
  list(coefficients=drop(coef), var=var, residuals=res, df.residual=n-1,
       penalty.matrix=penalty.matrix, 
       stats=stats, effective.df.diagonal=dag)
}


predict.ols <- 
  function(object, newdata,
           type=c("lp","x","data.frame","terms","adjto","adjto.data.frame",
             "model.frame"),
           se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
           incl.non.slopes, non.slopes, kint=1,
           na.action=na.keep, expand.na=TRUE, center.terms=TRUE, ...)
  predictDesign(object, newdata, type, se.fit, conf.int, conf.type,
                incl.non.slopes, non.slopes, kint,
                na.action, expand.na, center.terms, ...)
