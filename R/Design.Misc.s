#Miscellaneous functions to retrieve characteristics of design
Surv.strata <- function(surv){
  if(.R.) {
    attr(surv, 'strata')
  } else {
    surv[,3]
  }
}

Surv.time <- function(surv) {
  surv[,1]
}

Surv.event <- function(surv) {
  surv[,2]
}


#Function to get the number of intercepts in front of the slope coefficients
#ols - one intercept
#lrm - one or more intercepts (>1 if ordinal model)
#cph - no intercepts,   etc.

num.intercepts <- function(fit)
{
   nrp <- fit$non.slopes
   ## changed is.null(nrp) to below, fit$coefficients to fit$coef 14Aug01
   if(!length(nrp))
   {
	nm1 <- names(fit$coef)[1]  # 14Sep00
	nrp <- 1*(nm1=="Intercept" | nm1=="(Intercept)")
   }
   nrp
}

if(FALSE) DesignAssign <- function(atr, non.slopes=0) {
  ## Given Design attributes atr creates an assign list for R
  ## that is in S-Plus format
  p <- length(atr$name)
  ns <- atr$non.slopes
  a <- vector('list',p+non.slopes)
  names(p) <- c(if(non.slopes>0)'Intercept' else NULL,atr$name)
  if(non.slopes > 0) {
    a[[1]] <- 1:non.slopes
    j <- 1
  } else j <- 0
  for(i in 1:p) {
    j <- j+1
    a[[j]] <- 1:length(atr)
  }
 
}

DesignAssign <- function(atr, non.slopes, Terms) {
  ## Given Design attributes and number of intercepts creates S-Plus
  ## format assign list (needed for R, intercept correction needed for
  ## S-Plus anyway).  If formula is given, names assign using
  ## terms(formul) term.labels, otherwise uses Design predictor names
  ## 23feb03: No, term.labels not useful if "." in formula
  ## formula argument no longer used

  ## ll <- if(missing(formula)) atr$name else attr(terms(formula),'term.labels')
##  ll <- atr$name   ## 22feb03
  ## Changed 24feb03 to pass terms instead of formula, us it
  ll <- if(missing(Terms)) atr$name else attr(Terms,'term.labels')
  if(!length(ll)) return(list())
  nv <- length(ll)
  params <- sapply(atr$nonlinear, length)  ## d.f. per predictor
  asc <- atr$assume.code
  assign <- list() #vector('list', nv+(non.slopes > 0)-sum(asc==8))
  j <- non.slopes + 1
  for(i in 1:length(ll)) {
    if(asc[i]==8) next
    assign[[ll[i]]] <- j:(j+params[i]-1)
    j <- j + params[i]
  }
  assign
}
  

## Functions for Out Of Sample computation of -2 log likelihood
## evaluated at parameter estimates of a given fit

oos.loglik <- function(fit, ...) UseMethod("oos.loglik")

oos.loglik.ols <- function(fit, lp, y, ...) {
  sigma2 <- sum(fit$residuals^2)/length(fit$residuals)
  if(missing(lp)) {
	n <- length(fit$residuals)
	n*logb(2*pi*sigma2)+n
  } else {
	s <- !is.na(lp + y)
	lp <- lp[s]; y <- y[s]
	n <- length(lp)
	sse <- sum((y - lp)^2)
	n*logb(2*pi*sigma2) + sse/sigma2
  }
}

oos.loglik.lrm <- function(fit, lp, y, ...) {
  if(missing(lp)) return(fit$deviance[length(fit$deviance)])
  ns <- fit$non.slopes
  if(ns > 1) stop('ordinal y case not implemented')
  y <- as.integer(as.category(y)) - 1
  s <- !is.na(lp + y)
  lp <- lp[s];  y <- y[s]
  p <- plogis(lp)
  -2*sum(ifelse(y==1, logb(p), logb(1-p)))
}
  
oos.loglik.cph <- function(fit, lp, y, ...) {
  if(missing(lp)) return(-2*fit$loglik[2])
  else stop('not implemented for cph models')
}

oos.loglik.psm <- function(fit, lp, y, ...) {
  if(missing(lp)) return(-2*fit$loglik[2])
  else stop('not implemented for psm models')
}

oos.loglik.glmD <-
  if(.R.) function(fit, lp, y, ...) {
  if(missing(lp)) return(deviance(fit))
  glm.fit.null(x=NULL,y=as.vector(y),offset=lp,family=fit$family)$deviance
} else function(fit, lp, y, ...) {
  if(missing(lp)) return(deviance(fit))
  family <- origGlmFamily(fit$family)
  glm.null(x=NULL, y=as.vector(y),
          offset=lp, family=family, w=rep(1,length(lp)))$deviance
}
  
#Function to retrieve limits and values, from fit (if they are there)
#or from a datadist object.  If need.all=F and input is coming from datadist,
#insert columns with NAs for variables not defined
#at is attr(fit$terms,"Design") (now fit$Design)

Getlim <- function(at, allow.null=FALSE, need.all=TRUE)
{
nam <- at$name[at$assume!="interaction"]
limits <- at$limits
values <- at$values

XDATADIST <- .Options$datadist
X <- lims <- vals <- NULL
if(!is.null(XDATADIST) && exists(XDATADIST)) {
   X <- if(.R.) eval(as.name(XDATADIST)) else
     eval(as.name(XDATADIST),local=FALSE)  #27May99  9Apr02
   lims <- X$limits
   if(is.null(lims)) stop(paste("options(datadist=",XDATADIST,
	") not created with datadist"))
   vals <- X$values
}

if((length(X)+length(limits))==0) {
  if(allow.null) {
    lims <- list()
    for(nn in nam) lims[[nn]] <- rep(NA,7)
    lims <- structure(lims, class="data.frame", 
      row.names=c("Low:effect","Adjust to", "High:effect", "Low:prediction",
		  "High:prediction","Low","High"))
    return(list(limits=lims, values=values))
  }
  stop("no datadist in effect now or during model fit")
}

na <- if(length(limits))
  sapply(limits, function(x) all(is.na(x))) else rep(TRUE, length(nam))
if(length(lims) && any(na)) for(n in nam[na]) { #if() assumes NA stored in fit
						# for missing vars
  z <- limits[[n]]
  u <- if(match(n, names(lims), 0) > 0) lims[[n]] else NULL
  # This requires exact name match, not substring match
  if(is.null(u)) {
    if(need.all) stop(paste("variable",n,
	"does not have limits defined in fit or with datadist"))
    else limits[[n]] <- rep(NA,7)    # Added 28 Jul 94
  }
  else limits[[n]] <- u
}
limits <- structure(limits, class="data.frame", 
   row.names=c("Low:effect","Adjust to", "High:effect", "Low:prediction",
		"High:prediction","Low","High"))

if(length(vals)) values <- c(values, 
	vals[match(names(vals),nam,0)>0 & match(names(vals),names(values),0)==0]
	)   # add in values from datadist corresponding to vars in model
            # not already defined for model

list(limits=limits, values=values)
}

#Function to return limits for an individual variable, given an object
#created by Getlim

Getlimi <- function(name, Limval, need.all=TRUE)
{
   lim <- if(match(name, names(Limval$limits), 0) > 0) 
     Limval$limits[[name]] else NULL
   if(is.null(Limval) || is.null(lim) || all(is.na(lim))) {
      if(need.all) stop(paste("no limits defined by datadist for variable",
			name))
      return(rep(NA,7))
   }
lim
}

#Function to return a list whose ith element contains indexes
#of all predictors related, indirectly or directly, to predictor i
#Predictor i and j are related indirectly if they are related to
#any predictors that interact
#Set type="direct" to only include factors interacting with i
#This function is used by nomogram.Design

related.predictors <- function(at, type=c("all","direct"))
{
  type <- match.arg(type)
  f <- sum(at$assume.code < 9)
  if(any(at$assume.code == 10)) stop("does not work with matrix factors")
  ia <- at$interactions
  x <- rep(NA,f)
  names(x) <- at$name[at$assume.code < 9]
  mode(x) <- "list"
  if(length(ia)==0)
    {
      for(i in 1:f) x[[i]] <- integer(0)
      return(x)
    }
  for(i in 1:f)
    {
      r <- integer(0)
      for(j in 1:ncol(ia))
        {
          w <- ia[,j]
          if(any(w==i)) r <- c(r,w[w>0 & w!=i])
        }
      x[[i]] <- r
    }
  if(type=="direct") return(x)
  
  while(TRUE)
    {
      bigger <- FALSE
      for(j in 1:f)
        {
          xj <- x[[j]]
          y <- unlist(x[xj])
          y <- y[y != j]
          new <- unique(c(y, xj))
          bigger <- bigger | length(new) > length(xj)
          x[[j]] <- new
        }
      if(!bigger) break
    }
  x
}


#Function to list all interaction term numbers that include predictor
#pred as one of the interaction components

interactions.containing <- function(at, pred) {
ia <- at$interactions
if(length(ia)==0) return(NULL)
name <- at$name
parms <- at$parms
ic <- NULL
for(i in (1:length(at$assume.code))[at$assume.code==9]) {
    terms.involved <- parms[[name[i]]][,1]
    if(any(terms.involved==pred)) ic <- c(ic, i)
}
ic
}

#Function to return a vector of logical values corresponding to
#non-intercepts, indicating if the parameter is one of the following types:
# term.order  Meaning
# ----------  -----------------
#     1       all parameters
#     2       all nonlinear or interaction parameters
#     3       all nonlinear parameters (main effects or interactions)
#     4       all interaction parameters
#     5       all nonlinear interaction parameters

param.order <- function(at, term.order) {	#at=Design attributes
if(term.order==1) return(rep(TRUE,length(at$colnames)))
nonlin <- unlist(at$nonlinear[at$name[at$assume!="strata"]]) # omit strat
ia <- NULL
for(i in (1:length(at$name))[at$assume!="strata"])
  ia <- c(ia, rep(at$assume[i]=="interaction",length(at$nonlinear[[i]])))
if(term.order==5) nonlin & ia else if(term.order==4) ia else
if(term.order==3) nonlin else nonlin | ia
}


#	Design.levels
#		Make each variable in an input data frame that is a
#		factor variable in the model be a factor variable with
#		the levels that were used in the model.  This is primarily
#		so that row insertion will work right with <-[.data.frame
#	
#at=Design attributes

Design.levels <- function(df, at) {
ac <- at$assume.code
for(nn in names(df)) {
  j <- match(nn, at$name, 0)
  if(j>0) {
    if((ac[j]==5 | ac[j]==8) & length(lev <- at$parms[[nn]]))
      df[[nn]] <- factor(df[[nn]], lev)
  }
}
df
}


#Function to return a default penalty matrix for penalized MLE,
#according to the design attributes and a design matrix X

Penalty.matrix <- function(at, X) {

d1 <- dimnames(X)[[2]][1]
if(d1=='Intercept' || d1=='(Intercept)') X <- X[,-1,drop=FALSE]

d <- dim(X)
n <- d[1]; p <- d[2]
center <- as.vector(rep(1/n,n) %*% X)   # see scale() function
v <- if(.R.) as.vector(rep(1/(n-1),n) %*%
                       (X - rep(center,rep(n,p)))^2) else
as.vector(rep(1/(n-1),n) %*% (X - rep(center,rep.int(n,p)))^2)

pen <- if(p==1) as.matrix(v) else as.matrix(diag(v))    
# works even if X one column

is <- 1
ac <- at$assume
for(i in (1:length(at$name))[ac!="strata"]) {
  len <- length(at$nonlinear[[i]])
  ie <- is + len - 1
  if(ac[i] == "category") pen[is:ie,is:ie] <- diag(len) - 1/(len+1)
  is <- ie+1
}
pen
}

#Function to take as input a penalty specification of the form
#penalty=constant or penalty=list(simple=,nonlinear=,interaction=,
#nonlinear.interaction=) where higher order terms in the latter notation
#may be omitted, in which case their penalty factors are taken from lower-
#ordered terms.  Returns a new penalty object in full list form along
#with a full vector of penalty factors corresponding to the elements
#in regression coefficient vectors to be estimated

Penalty.setup <- function(at, penalty) {
if(!is.list(penalty)) penalty <- list(simple=penalty, nonlinear=penalty,
  interaction=penalty, nonlinear.interaction=penalty)
tsimple <- penalty$simple
if(!length(tsimple)) tsimple <- 0
tnonlinear <- penalty$nonlinear
if(!length(tnonlinear)) tnonlinear <- tsimple
tinteraction <- penalty$interaction
if(!length(tinteraction)) tinteraction <- tnonlinear
tnonlinear.interaction <- penalty$nonlinear.interaction
if(!length(tnonlinear.interaction)) tnonlinear.interaction <- tinteraction

nonlin <- unlist(at$nonlinear[at$name[at$assume!='strata']])
ia <- NULL
for(i in (1:length(at$name))[at$assume!='strata']) ia <- c(ia,
  rep(at$assume[i]=='interaction',length(at$nonlinear[[i]])))
nonlin.ia <- nonlin & ia
nonlin[nonlin.ia] <- FALSE
ia[nonlin.ia] <- FALSE
simple <- rep(TRUE, length(ia))
simple[nonlin | ia | nonlin.ia] <- FALSE
penfact <- tsimple*simple + tnonlinear*nonlin + tinteraction*ia +
  tnonlinear.interaction*nonlin.ia
list(penalty=list(simple=tsimple, nonlinear=tnonlinear,
  interaction=tinteraction,nonlinear.interaction=tnonlinear.interaction),
  multiplier=penfact)
}

#Function to do likelihood ratio tests from two models that are
# (1) nested and (2) have 'Model L.R.' components of the stats
# component of the fit objects
# For models with scale parameters, it is also assumed that the
# scale estimate for the sub-model was fixed at that from the larger model

lrtest <- function(fit1, fit2) {
  if(length(fit1$fail) && fit1$fail)
    stop('fit1 had failed')
  if(length(fit2$fail) && fit2$fail)
    stop('fit2 had failed')
  
  s1 <- fit1$stats
  s2 <- fit2$stats

  if(!length(s1))
    s1 <- c('Model L.R.'=fit1$null.deviance - fit1$deviance,
            'd.f.'=fit1$rank - (any(names(coef(fit1))=='(Intercept)')))
  if(!length(s2))
    s2 <- c('Model L.R.'=fit2$null.deviance - fit2$deviance,
            'd.f.'=fit2$rank - (any(names(coef(fit2))=='(Intercept)')))
  ## 26nov02
  
  chisq1 <- s1['Model L.R.']
  chisq2 <- s2['Model L.R.']
  if(length(chisq1)==0 || length(chisq2)==2) 
    stop('fits do not have stats component with "Model L.R." or deviance component')
  df1 <- s1['d.f.']
  df2 <- s2['d.f.']
  if(df1==df2) stop('models are not nested')

  lp1 <- length(fit1$parms);  lp2 <- length(fit2$parms)
  if(lp1!=lp2) warning('fits do not have same number of scale parameters') else 
  if(lp1==1 && abs(fit1$parms-fit2$parms)>1e-6)
    warning('fits do not have same values of scale parameters.\nConsider fixing the scale parameter for the reduced model to that from the larger model.')

  chisq <- abs(chisq1-chisq2)
  dof   <- abs(df1-df2)
  p     <- 1-pchisq(chisq,dof)

  r     <- c(chisq,dof,p)
  names(r) <- c('L.R. Chisq','d.f.','P')
  structure(list(stats=r,
                 formula1=formula(fit1),
                 formula2=formula(fit2)),
            class='lrtest')
}

print.lrtest <- function(x, ...) {
  cat('\nModel 1: '); print(x$formula1)
  cat('Model 2: '); print(x$formula2); cat('\n')
  print(x$stats)
  cat('\n')
  invisible()
}


Newlabels <- function(fit, ...) UseMethod('Newlabels')

Newlabels.Design <- function(fit, labels, ...) {

at <- fit$Design
if(!length(at)) at <- getOldDesign(fit)   # 13Apr01
nam <- names(labels)
if(length(nam)==0) {
  if(length(labels)!=length(at$name))
	stop('labels is not a named vector and its length is not equal to the number of variables in the fit')
  nam <- at$name
} 
i <- match(nam, at$name, nomatch=0)

if(any(i==0)) {
  warning(paste('the following variables were not in the fit and are ignored:\n',
			 paste(nam[i==0],collapse=' ')))
  labels <- labels[i>0]
  i <- i[i>0]
}

at$label[i] <- labels

#attr(fit$terms, 'Design') <- at
fit$Design <- at   # 13Apr01
fit
}

Newlevels <- function(fit, ...) UseMethod('Newlevels')

Newlevels.Design <- function(fit, levels, ...) {

at <- fit$Design
if(!length(at)) at <- getOldDesign(fit)  # 13Apr01
nam <- names(levels)
if(length(nam)==0) stop('levels must have names')

i <- match(nam, at$name, nomatch=0)

if(any(i==0)) {
  warning(paste('the following variables were not in the fit and are ignored:\n',
			 paste(nam[i==0],collapse=' ')))
  nam <- nam[i>0]
}

for(n in nam) {
  prm <- at$parms[[n]]
  if(length(prm)!=length(levels[[n]])) stop(paste('new levels for variable',
			 n,'has the wrong length'))

  levs <- levels[[n]]
  if(length(at$values[[n]])) at$values[[n]] <- levs
  if(length(at$limits)) {
	m <- match(at$limits[[n]], at$parms[[n]])
	if(is.category(at$limits[[n]]))
	  attr(at$limits[[n]],'levels') <- levs else
	at$limits[[n]] <- levs[m]
  }
  at$parms[[n]] <- levs
}

fit$Design <- at   # 13Apr01
fit
}


##13Nov00
DesignFit <- function(fit) {
  cl <- oldClass(fit)
  if(cl[1]=='Design') return(fit)
  if(length(fit$Design)==0) getOldDesign(fit)  # 13Apr01
  fit$fitFunction <- cl
  oldClass(fit) <- 'Design'
  fit
}

  
if(.SV4.) {
  predict.Design <- function(object, ...) {
    fitter <- object$fitFunction
    if(!length(fitter)) stop(
      "fit's main class is 'Design' but no fitFunction element is present")
    oldClass(object) <- fitter[1]
    predict(object, ...)
  }

  print.Design <- function(x, ...) {
    fitter <- x$fitFunction
    if(!length(fitter)) stop(
      "fit's main class is 'Design' but no fitFunction element is present")
    oldClass(x) <- fitter[1]
    print(x, ...)
#    f <- paste('print',fitter[1],sep='.')
#    if(!existsFunction(f))print.default(x, ...) else
#    do.call(f, list(x, ...))
  }

  residuals.Design <- function(object, ...) {
    fitter <- object$fitFunction
    if(!length(fitter)) stop(
      "fit's main class is 'Design' but no fitFunction element is present")
    oldClass(object) <- fitter[1]
    residuals(object, ...)
#    f <- paste('residuals',fitter[1],sep='.')
#    oldClass(object) <- 
#    if(!existsFunction(f))residuals.default(object, ...) else
#    do.call(f, list(object, ...))
  }

  validate.Design <- function(fit, ...) {
    fitter <- fit$fitFunction
    if(!length(fitter)) stop(
      "fit's main class is 'Design' but no fitFunction element is present")
    oldClass(fit) <- fitter[1]
    validate(fit, ...)
#    f <- paste('validate',fitter[1],sep='.')
#    if(!existsFunction(f))
#      stop(paste('no validate method exists for fit of type',fitter[1]))
#    do.call(f, list(fit, ...))
  }

  calibrate.Design <- function(fit, ...) {
    fitter <- fit$fitFunction
    if(!length(fitter)) stop(
      "fit's main class is 'Design' but no fitFunction element is present")
    oldClass(fit) <- fitter[1]
    calibrate(fit, ...)
#    f <- paste('calibrate',fitter[1],sep='.')
#    if(!existsFunction(f))calibrate.default(fit, ...) else
#    do.call(f, list(fit, ...))
  }

  Survival.Design <- function(object, ...) {
    fitter <- object$fitFunction
    if(!length(fitter)) stop(
      "fit's main class is 'Design' but no fitFunction element is present")
    oldClass(object) <- fitter[1]
    Survival(object, ...)
  }

  Quantile.Design <- function(object, ...) {
    fitter <- object$fitFunction
    if(!length(fitter)) stop(
      "fit's main class is 'Design' but no fitFunction element is present")
    oldClass(object) <- fitter[1]
    Quantile(object, ...)
  }

  Mean.Design <- function(object, ...) {
    fitter <- object$fitFunction
    if(!length(fitter)) stop(
      "fit's main class is 'Design' but no fitFunction element is present")
    oldClass(object) <- fitter[1]
    Mean(object, ...)
  }

  Hazard.Design <- function(object, ...) {
    fitter <- object$fitFunction
    if(!length(fitter)) stop(
      "fit's main class is 'Design' but no fitFunction element is present")
    oldClass(object) <- fitter[1]
    Hazard(object, ...)
  }

  latex.Design <- function(object, title,
                           file=paste(first.word(deparse(substitute(object))),
                             'tex',sep='.'), ...) {
    fitter <- object$fitFunction
    if(!length(fitter)) stop(
      "fit's main class is 'Design' but no fitFunction element is present")
    oldClass(object) <- fitter[1]
    ## Need to brute-force dispatch because of SV4 problem in latex in Hmisc
    if(existsFunction(p <- paste('latex',fitter[1],sep='.')))
      do.call(p, list(object, file=file, ...)) else latexDesign(object, file=file, ...)
  }

  survest.Design <- function(fit, ...) {
    fitter <- fit$fitFunction
    if(!length(fitter)) stop(
      "fit's main class is 'Design' but no fitFunction element is present")
    f <- paste('survest',fitter[1],sep='.')
    do.call(f, list(fit,...))
  }

  Varcov.Design <- function(object, regcoef.only=FALSE, ...) {
     fitter <- object$fitFunction
    if(!length(fitter)) stop(
      "fit's main class is 'Design' but no fitFunction element is present")
    f <- paste('Varcov',fitter[1],sep='.')
    do.call(f, list(object, regcoef.only, ...))
  }

  oos.loglik.Design <- function(fit, ...) {   ## 6dec02
    fitter <- fit$fitFunction
    if(!length(fitter)) stop(
      "fit's main class is 'Design' but no fitFunction element is present")
    f <- paste('oos.loglik',fitter[1],sep='.')
    do.call(f, list(fit,...))
  }
   
  NULL
}

#3Apr01
getOldDesign <- function(fit) {
  at <- attr(fit$terms,'Design')
  if(is.null(at))
    stop('fit was not created by a Design library fitting function')
  at
}


## make.link used by survreg, survreg.distributions in R
## Design for SV4 uses R version of survreg.distributions
if(!.R.) {
make.link <- function (link) 
{
    if (is.character(link) && length(grep("^power", link) > 0)) 
        return(eval(parse(text = link)))
    else if (!is.character(link) && !is.na(lambda <- as.numeric(link))) {
        linkfun <- function(mu) mu^lambda
        linkinv <- function(eta) pmax(.Machine$double.eps, eta^(1/lambda))
        mu.eta <- function(eta) pmax(.Machine$double.eps, (1/lambda) * 
            eta^(1/lambda - 1))
        valideta <- function(eta) all(eta > 0)
    }
    else switch(link, logit = {
        linkfun <- function(mu) log(mu/(1 - mu))
        linkinv <- function(eta) {
            thresh <- -log(.Machine$double.eps)
            eta <- pmin(thresh, pmax(eta, -thresh))
            exp(eta)/(1 + exp(eta))
        }
        mu.eta <- function(eta) {
            thresh <- -log(.Machine$double.eps)
            res <- rep(.Machine$double.eps, length(eta))
            res[abs(eta) < thresh] <- (exp(eta)/(1 + exp(eta))^2)[abs(eta) < 
                thresh]
            res
        }
        valideta <- function(eta) TRUE
    }, probit = {
        linkfun <- function(mu) qnorm(mu)
        linkinv <- function(eta) {
            thresh <- -qnorm(.Machine$double.eps)
            eta <- pmin(thresh, pmax(eta, -thresh))
            pnorm(eta)
        }
        mu.eta <- function(eta) pmax(dnorm(eta), .Machine$double.eps)
        valideta <- function(eta) TRUE
    }, cloglog = {
        linkfun <- function(mu) log(-log(1 - mu))
        linkinv <- function(eta) pmax(.Machine$double.eps, pmin(1 - 
            .Machine$double.eps, 1 - exp(-exp(eta))))
        mu.eta <- function(eta) {
            eta <- pmin(eta, 700)
            pmax(.Machine$double.eps, exp(eta) * exp(-exp(eta)))
        }
        valideta <- function(eta) TRUE
    }, identity = {
        linkfun <- function(mu) mu
        linkinv <- function(eta) eta
        mu.eta <- function(eta) rep(1, length(eta))
        valideta <- function(eta) TRUE
    }, log = {
        linkfun <- function(mu) log(mu)
        linkinv <- function(eta) pmax(.Machine$double.eps, exp(eta))
        mu.eta <- function(eta) pmax(.Machine$double.eps, exp(eta))
        valideta <- function(eta) TRUE
    }, sqrt = {
        linkfun <- function(mu) mu^0.5
        linkinv <- function(eta) eta^2
        mu.eta <- function(eta) 2 * eta
        valideta <- function(eta) all(eta > 0)
    }, "1/mu^2" = {
        linkfun <- function(mu) 1/mu^2
        linkinv <- function(eta) 1/eta^0.5
        mu.eta <- function(eta) -1/(2 * eta^1.5)
        valideta <- function(eta) all(eta > 0)
    }, inverse = {
        linkfun <- function(mu) 1/mu
        linkinv <- function(eta) 1/eta
        mu.eta <- function(eta) -1/(eta^2)
        valideta <- function(eta) all(eta != 0)
    }, stop(paste(link, "link not recognised")))
    list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
        valideta = valideta)
}
NULL
}

# The following makes formula(Design fit) work - ols was using formula.lm
# if(.R.) formula.Design <- formula.default  16dec03

oldDesignFit2R <- function(f) {
  for(u in Cs(coefficients,residuals,fitted.values,linear.predictors,
              stats,freq,u)) f[[u]] <- Names2names(f[[u]])
  if(inherits(f,'ols') && ('Sigma' %nin% names(f$stats)))
    f$stats <- c(f$stats, c(Sigma=sqrt(sum(f$residuals^2)/f$df.residual)))
  f
}

univarLR <- function(fit) {
## Computes all univariable LR chi-square statistics
  w <- as.character(attr(fit$terms,'variables'))
  if(.R.) w <- w[-1]
  p <- length(w)-1
  stats <- P <- double(p)
  dof <- nobs <- integer(p)
  for(i in 1:p) {
    stat <- update(fit, as.formula(paste(w[1],w[i+1],sep='~')))$stats
    stats[i] <- stat['Model L.R.']
    dof[i]   <- stat['d.f.']
    P[i]     <- stat['P']
    nobs[i]  <- stat['Obs']
  }
data.frame(LR=stats, 'd.f.'=dof, P=P, N=nobs,
           row.names=w[-1], check.names=FALSE)
}
