## Copyright (C) 2001 Frank E Harrell Jr
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation; either version 2, or (at your option) any
## later version.
##
## These functions are distributed in the hope that they will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## The text of the GNU General Public License, version 2, is available
## as http://www.gnu.org/copyleft or by writing to the Free Software
## Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
##
#under.unix   <- !(version$os=='Microsoft Windows' ||
#                  version$os=='Win32' || version$os=='mingw32')
.R.          <- TRUE
.SV4.        <- FALSE
.newSurvival. <- TRUE

.noGenenerics <- TRUE  # faster loading as new methods not used

.First.lib <- function(lib, pkg) {
  cat("Design library by Frank E Harrell Jr\n\n",
      "Type library(help='Design'), ?Overview, or ?Design.Overview')\n",
      "to see overall documentation.\n\n",
      sep='')
  library.dynam("Design", pkg, lib)
  require('Hmisc')
  invisible()
}

#  requirePos('survival',pos=3)
#  p <- .packages()
#  if(match('Design',p) > match('survival',p))
#    warning('By not specifying library(Design,pos=2), functions in the\nsurvival package such as survfit will override those in Design.')
#Miscellaneous functions to retrieve characteristics of design

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
  
#Function to return variance-covariance matrix, optionally deleting
#rows and columns corresponding to parameters such as scale parameters
#in parametric survival models

#Varcov <- function(object, ...) UseMethod("Varcov")  in Hmisc

Varcov.lrm <- function(object, regcoef.only=FALSE, ...)
  Varcov.default(object, regcoef.only, ...)  # for fastbw etc.
Varcov.ols <- function(object, regcoef.only=FALSE, ...)
  Varcov.default(object, regcoef.only, ...)
Varcov.cph <- function(object, regcoef.only=FALSE, ...)
  Varcov.default(object, regcoef.only, ...)
Varcov.psm <- function(object, regcoef.only=FALSE, ...)
  Varcov.default(object, regcoef.only, ...)

#Varcov.default <- function(fit, regcoef.only=F)  Defined in Hmisc
#{
#  vc <- fit$Varcov
#  if(length(vc)) {
#	if(regcoef.only) return(fit$var) else
#	return(vc(fit,which='var'))
#  }
#   cov <- fit$var
#   if(is.null(cov)) stop("fit does not have variance-covariance matrix")
#   if(regcoef.only)
#   {
#	p <- length(fit$coefficients)  # 14Sep00
#	cov <- cov[1:p, 1:p, drop=F]
#   }
#   cov
#}

#Varcov.lm <- function(fit, ...)
#{
#   cof <- fit$coefficients
#   rinv <- solve(fit$R, diag(length(cof)))
#   cov <- rinv %*% t(rinv)
#   cov <- sum(fit$residuals^2)*cov/fit$df.residual
#   nm  <- names(cof)
#   dimnames(cov) <- list(nm, nm)
#   cov
#}

# Varcov.glm was erroneously defined as follows   30Jul99
#   p <- fit$rank
#   if(is.null(p)) p <- sum(!is.na(fit$coefficients))  # 14Sep00
#   R <- fit$R
#   if(p<max(dim(R))) stop("Varcov does not handle singular matrices")
#   rinv <- backsolve(R, diag(p))
#   rinv %*% t(rinv)

#Varcov.glm <- function(fit, ...)
#{
#  s <- summary.glm(fit)
#  s$cov.unscaled * s$dispersion
#}

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

related.predictors <- function(at, type=c("all","direct")) {
type <- match.arg(type)
f <- sum(at$assume.code<9)
if(any(at$assume.code==10))stop("does not work with matrix factors")
ia <- at$interactions
x <- rep(NA,f)
mode(x) <- "list"
if(length(ia)==0) {
  for(i in 1:f) x[[i]] <- integer(0)
  return(x)
  }
for(i in 1:f) {
  r <- integer(0)
  for(j in 1:ncol(ia)) {
    w <- ia[,j]
    if(any(w==i)) r <- c(r,w[w>0 & w!=i])
  }
  x[[i]] <- r
}
if(type=="direct") return(x)

while(TRUE) {
  bigger <- FALSE
  for(j in 1:f) {
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

  Survival.Design <- function(fit, ...) {
    fitter <- fit$fitFunction
    if(!length(fitter)) stop(
      "fit's main class is 'Design' but no fitFunction element is present")
    oldClass(fit) <- fitter[1]
    Survival(fit, ...)
  }

  Quantile.Design <- function(fit, ...) {
    fitter <- fit$fitFunction
    if(!length(fitter)) stop(
      "fit's main class is 'Design' but no fitFunction element is present")
    oldClass(fit) <- fitter[1]
    Quantile(fit, ...)
  }

  Mean.Design <- function(fit, ...) {
    fitter <- fit$fitFunction
    if(!length(fitter)) stop(
      "fit's main class is 'Design' but no fitFunction element is present")
    oldClass(fit) <- fitter[1]
    Mean(fit, ...)
  }

  Hazard.Design <- function(fit, ...) {
    fitter <- fit$fitFunction
    if(!length(fitter)) stop(
      "fit's main class is 'Design' but no fitFunction element is present")
    oldClass(fit) <- fitter[1]
    Hazard(fit, ...)
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
if(.R.) formula.Design <- formula.default

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

  
# Design  FEH 1Aug90, re-written 21Oct91
#
# Augments S formula language to include:
# 
#	name	- name[i] = name of ith original variable in x
#	label	- label[i] = label of ith original variable (=name if none)
#	assume	- assume(original x)
#	assume.code - coded version of assume (1-9, 9=added interaction)
#	parms	- parms(original x)
#		  for interaction effects parms[[i]] is a matrix with dim
#		  3 x (1+# interaction terms).  First element in pair
#		  is 1 if first factor is represented as an expanded
#		  non-linear term, 0 otherwise (this applies to polynomial,
#		  lspline, rcspline, scored).  Second element applies to
#		  second factor in interaction effect.  Third element
#		  applies to third factor (0 if second order interaction)
#		  First column contains factor numbers involved in interaction.
#	limits  - limits(original x)
#	values	- For continuous variables with <=10 unique values, is
#		  vector of values.  NULL otherwise.
#	interactions - 3 x k matrix of factor numbers
#
# Cannot have interactions between two stratification factors. 
#
#
#23Feb95 - added logic to remove cluster() factors
#
Design <- function(mf, allow.offset=TRUE, intercept=1)   			{

Terms <- attr(mf, "terms")
Term.labels <- attr(Terms,'term.labels')
iscluster <- if(length(Term.labels)) substring(Term.labels,1,8)=='cluster('
	else FALSE

#For some reason, model frame sometimes has a blank name if using %ia%

namx <- names(mf)
if(any(namx=="")) {
   namx <- names(mf) <- c(namx[1],Term.labels)
   dimnames(mf)[[2]] <- namx
   dimnames(attr(Terms,"factors"))[[1]] <- namx
 }

wts <- if(any(namx=='(weights)'))(1:length(namx))[namx=='(weights)']
else 0  # 4jun02

if(length(Terms)==0) inner.name <- NULL else {
	#Handle case where a function has two arguments that are names,
	#e.g. rcs(x,knots) -> want x only
	inner.name <- unique(var.inner(Terms))
	# var.inner is stripped down version of terms.inner (see Design.trans)
        #Note: these exclude interaction terms and %ia% terms
  }

response.pres <- attr(Terms, 'response') > 0   # 3Jun99

offs <- attr(Terms, "offset")
if(!length(offs)) offs <- 0
if(offs>0 & !allow.offset)
	stop("offset variable not allowed in formula")


factors <- attr(Terms, "factors")
if(length(factors) && response.pres) factors <- factors[-1,,drop=FALSE]

attr(Terms, "intercept") <- intercept
fname <- flabel <- name <- strt <- asm <- len <- 
	 fname.incl.dup <- ia <- funits <- NULL  # funits 20May99
parm <- nonlinear <- limits <- values <- list()

scol<-1
colnam <- list()

#Corrected 23Jun95 - if user used name 'dist' get would return 'dist'
XDATADIST <- .Options$datadist
if(length(XDATADIST))
{
   if(!exists(XDATADIST)) stop(paste("dataset",XDATADIST,
	"not found for options(datadist=)"))
   datadist <- if(.R.) eval(as.name(XDATADIST)) else
     eval(as.name(XDATADIST), local=FALSE)  #27May99  9Apr02
   Limits <- datadist$limits
   Limnames <- dimnames(Limits)[[2]]
}

nc <- 0

options(Design.attr=NULL, TEMPORARY=FALSE)	#Used internally by asis, rcs, ...
anyfactors <- ncol(mf) > 1*response.pres #3Jun99
i1.noia <- 0
if(anyfactors)for(i in (response.pres+1):ncol(mf)) {
	if(i!=offs && i!=wts && !iscluster[i-response.pres]) {  #3Jun99
	i1 <- i - response.pres   #3Jun99
	xi <- mf[[i]]
	z <- attributes(xi)
	assu <- z$assume.code
	if(!length(assu) || assu!=9) i1.noia <- i1.noia+1
	if(!length(assu))			{#Not processed w/asis,et
	   nam <- inner.name[i1.noia]
	   lab <- attr(xi, "label")
       ord <- is.ordered(xi) && all.is.numeric(levels(xi))   #21Jun99
	   if(!length(lab) || lab=="") lab <- nam
	   if(ord)	{
	      xi <- scored(xi, name=nam, label=lab)
	      attr(mf[,i],"contrasts") <- attr(xi,"contrasts")
			}
	   else if(is.character(xi) | is.category(xi)) {
         ## | is.factor(xi)) 21Feb02
         if(is.ordered(xi) &&
            .Options$contrasts[2]!='contr.treatment')
           warning(paste('Variable',nam,'is an ordered factor.\n',
                         'You should set options(contrasts=c("contr.treatment","contr.treatment"))\nor Design will not work properly.'))
         ## warning 6may03
         xi <- catg(xi, name=nam, label=lab)
       }
	   else if(is.matrix(xi)) xi <- matrx(xi, name=nam, label=lab)
	   else xi <- asis(xi, name=nam, label=lab)
	   z <- c(z,attributes(xi))
						}

	za <- z$assume.code
	zname <- z$name

	fname.incl.dup <- c(fname.incl.dup, zname)
	if(!length(fname) || !any(fname==zname))	{ #unique factor
	  nc <- nc+1
	  fname <- c(fname,zname)
	  flabel <- c(flabel,z$label)
      ##funits <- was here   9Jun99 (see 5 down)
	  asm <- c(asm,za)
	  colnam[[i1]] <- z$colnames
	  if(za!=8) name <- c(name, colnam[[i1]])
	  if(za!=9) 			{
        funits <- c(funits, if(length(z$units))z$units else '')
        if(length(z$parms)) parm[[zname]] <- z$parms
        if(length(XDATADIST))
          {
		limits[[zname]] <- if(any(Limnames==zname)) {
                  j <- match(zname, Limnames, 0) #require EXACT match
                  Limits[,j[j>0]]
                }
			else rep(NA,7)
		j <- match(zname, names(datadist$values), 0)
		if(j>0) {
          values[[zname]] <- datadist$values[[j]]
          l1 <- levels(xi); l2 <- datadist$values[[j]]  #20May99
          if(length(l1) && ((length(l1) != length(l2)) ||
                            any(sort(l1) != sort(l2))))
            warning(paste('Variable',zname,'has levels',paste(l1,collapse=' '),
                          'which do not match levels given to datadist (',
                          paste(l2,collapse=' '),'). datadist values ignored.'))
          values[[zname]] <- l1
        }
      }
       }

	  if(length(nonl <- z$nonlinear)) nonlinear[[zname]] <- nonl

	  if(za==9)			{
		iia <- match(z$ia, fname)
		if(any(is.na(iia)))stop(paste(paste(z$ia,collapse=" "),
	"cannot be used in %ia% since not listed as main effect"))
		ia <- cbind(ia, c(iia,0))
		parms <- rbind(z$parms,0)
		parms[,1] <- c(iia,0)
		if(length(parms)) parm[[zname]] <- parms
      }
    }
	nrows <- if(is.matrix(xi))nrow(xi) else length(xi)
  }
  }

#Save list of which factors where %ia% interactions (before adding automatic ias
which.ia <- (1:length(asm))[asm==9]

#Add automatically created interaction terms
if(anyfactors) {
  if((nrow(factors)-(offs>0)-sum(iscluster))!=length(fname.incl.dup))
	stop("program logic error 1")
for(i in 1:ncol(factors)) {
	f <- factors[,i]  #3Jun99 was -1,i
	j <- (1:length(f))[f>0]
	nia <- length(j)
	if(nia>1) {
	   fn <- fname.incl.dup[j]
	   jf <- match(fn,fname.incl.dup)
	   if(any(is.na(jf))) stop("program logic error 2")
	   nc <- nc + 1
	   asm <- c(asm,9)
	   if(nia==2)ialab <- paste(fn[1],"*",fn[2])
	   else if(nia==3)ialab <- paste(fn[1],"*",fn[2],"*",fn[3])
	   else stop("interaction term not second or third order")
	   fname <- c(fname, ialab)
	   flabel <- c(flabel, ialab)
	   if(sum(asm[jf]==8)>1)
		stop("cannot have interaction between two strata factors")
	   nn <- list()
	   for(k in 1:nia)			{
	      if(asm[jf[k]]==5 | asm[jf[k]]==8)
	         nn[[k]] <- paste(fn[k],"=",parm[[fname[jf[k]]]][-1],sep="")
	      else if(asm[jf[k]]==7)	{
	         nn[[k]] <- c(fn[k], paste(fn[k],"=",
	            parm[[fname[jf[k]]]][c(-1,-2)],sep=""))
					}
	      else nn[[k]] <- colnam[[jf[k]]]
						}
	      if(nia==2) nn[[3]] <- ""
	   parms <- jf
	   if(length(jf)==2) parms <- c(parms, 0)
	   nonlin <- NULL
	   nl1 <- nonlinear[[fname[jf[1]]]]
	   nl2 <- nonlinear[[fname[jf[2]]]]
#Strata factors don't have nonlinear duplicated for # levels - 1
	   if(asm[jf[1]]==8) nl1 <- rep(FALSE, length(parm[[fname[jf[1]]]])-1)
	   if(asm[jf[2]]==8) nl2 <- rep(FALSE, length(parm[[fname[jf[2]]]])-1)
	   if(nia==2) nl3 <- FALSE
	   else if(asm[jf[3]]==8) nl3 <- rep(FALSE, length(parm[[fname[jf[3]]]])-1)
	   else nl3 <- nonlinear[[fname[jf[3]]]]
	   n1 <- nn[[1]]
	   n2 <- nn[[2]]
	   n3 <- nn[[3]]
#model.matrix makes auto-products move first variable fastest, etc.
	   for(j3 in 1:length(n3))				 {
	        for(j2 in 1:length(n2))				{
	          for(j1 in 1:length(n1))		       {
		   parms <- cbind(parms,c(nl1[j1],nl2[j2],nl3[j3]))
		   nonlin <- c(nonlin, nl1[j1] | nl2[j2] | nl3[j3])
		   if(nia==2)name <- c(name, paste(n1[j1],"*",n2[j2]))
		   else name <- c(name, paste(n1[j1],"*",n2[j2],"*",n3[j3]))
							       }}}

#If was 2-way interaction and one of the factors was restricted %ia%,
#adjust indicators
	k <- match(jf, which.ia, 0)
	if(any(k>0)) {
	  if(nia==3) stop("cannot have 2-way interaction with an %ia% interaction")
	  k <- jf[k>0]
	  wparm <- parms[,1]==k; wparm[3] <- TRUE
	  parms[wparm,] <- parm[[fname[k]]][1:2,,drop=FALSE]
	  jf <- parms[,1]
	  nonlin <- apply(parms, 2, any)[-1]
	}

	   if(length(jf)==2) {jf <- c(jf, 0); parms[3,] <- 0}
	   ia <- cbind(ia, jf)
	   if(length(parms)) parm[[ialab]] <- parms
	   if(length(nonlin)) nonlinear[[ialab]] <- nonlin

									}}}

if(!.R.) attr(mf,"names") <- NULL   ## was needed at all?  8Apr02
if(anyfactors) {
if(length(XDATADIST))
   limits <- structure(limits, row.names=c("Low:effect","Adjust to",
	"High:effect",
	"Low:prediction","High:prediction","Low","High"),class="data.frame")
#data.frame converts variables always NA to factor!

if(length(funits) != sum(asm!=9)) warning('program logic warning 1')
else names(funits) <- fname[asm!=9]

atr <- list(name=fname, label=flabel, units=funits, colnames=name,
            assume=c("asis","polynomial","lspline","rcspline","category",
              "","scored","strata","interaction","matrix")[asm],
            assume.code=as.integer(asm), parms=parm, limits=limits,
            values=values,nonlinear=nonlinear,
            interactions=structure(ia,dimnames=NULL))
}

else atr <- list(name=NULL, assume=NULL, assume.code=NULL, parms=NULL)

#attr(Terms,"Design") <- atr  13Apr01 (plus next 4)
#attr(mf,"terms") <- Terms
#attr(mf,"offset") <- offs
attr(mf, 'Design') <- atr
mf
}

# design.trans   FEH 4 Oct 90
# Contains individual functions for creating sub-design matrices from
# vectors, for use with design().
# code  name
# 1	asis	leave variable coded as-is, get default name, label,
#		limits, values
# 2	pol	polynomial expansion
# 3	lsp	linear spline
# 4	rcs	restricted cubic spline
# 5	catg	category
# 7	scored	scored ordinal variable
# 8	strat	stratification factor
#10	matrx	matrix factor - used to keep groups of variables together
#		as one factor
#
#	des.args generic function for retrieving arguments
#	set.atr generic function to set attributes of sub design matrix
#	options sets default options
#	[.Design subsets variables, keeping attributes
#	gparms	retrieve parms for design or fit object.  Not used by any
#		of these routines, but used by analyst to force a new
#		fit to use same parms as previous fit for a given factor.
#	value.chk
#		Check a given list of values for a factor for validity,
#		or if list is NA, return list of possible values
#	var.inner - stripped down terms.inner, returns character strings
#	
# Default label is attr(x,"label") or argument name if label= omitted.
# First argument can be as follows, using asis as an example:
#	asis(x, ...)		name="x", label=attr(x,"label") or "x"
#				if NULL
#	asis(w=abs(q), ...)	name="w", label=attr(x,"label") or "w"
#	asis(age=xx)		name="age", label=label attr or "age"
#	asis(x,label="Age, yr")	name="x", label="Age, yr"
#	asis(age=x,label=	name="age", label="Age in Years"
#		"Age in Years")
#	matrx(dx=cbind(dx1=dx1,dx2=dx2))	name="dx", individual names
#						dx1 and dx2
# For matrx, default label is list of column names.
# An additional argument, name, can be used to instead specify the name of the
# variable.  This is used when the functions are implicitly called from within
# design().
#
# The routines define dimnames for the returned object with column
# names = expanded list of names based on original name.
# assume.code is added to attributes of returned matrix.  Is 1-8
# corresponding to transformation routines asis-strat above, 10 for matrx.
# Adds attribute nonlinear, one element/column of expanded design matrix.
# nonlinear=T if column is a nonlinear expansion of original variable,
# F if linear part or not applicable
# (e.g. dummy variable for category -> F).  For matrx, all are linear.
#
# System options used: nknots for default number of knots in restr. cubic spline
# and poly.degree, default degree of polynomials 
# Second argument to routines is the parameters (parms) of the
# transformation (except for asis), defined as follows:
#
#	poly	order of polynomial, e.g. 2 for quadratic
#	lsp	list of knots
#	rcs	number of knots if parms=1 element (-> compute default
#		knot locations), actual knot locations if >2 elements
#		(2 knots not allowed for restr. cubic spline)
#	catg	list of value labels corresponding to values 1,2,3,...
#	scored	list of unique values of the variable
#	strat	list of value labels corresponding to values 1,2,3
#
# For catg and strat, parms are omitted if the variable is character or
# is already an S category variable.
#
# Argument retrieval: After variable and optional parms, other variables
# may be named or positional, in the following order: label, name.
# For matrx, parms are not allowed.
#
# Function to return list with elements name, parms, label. 
# corresponding to arguments in call to asis, etc.  parms=NULL if 
# parms.allowed=F.  Reason for going to this trouble is that first arg to
# asis, etc. is allowed to be a named argument to set a new name for it.
# With ordinary argument fetching, remaining arguments would have to be
# named.  This logic allows them to be named or positional in order:
# parms (if allowed), label.
#
# If options(Design.attr) is non-null, looks up attributes in elements
# in Design.attr corresponding to the name of the current variable.
# This is used to get predicted values when the original fitting
# function (e.g., rcs) derived parms of the transformation from the data.
#
des.args <- function(x,parms.allowed,call.args) {
  nam <- names(x)
  if(!length(nam)) nam <- rep("",5)
  name <- nam[1]
  if(name=="") {
    form <- formula(call("~",as.name("...y..."),call.args[[2]]))
    name <- var.inner(form)
  }
  pa <- parms.allowed
  argu <- function(x,karg, arg.name, parms.all, nm)	{
	if(!parms.all) karg <- karg-1
	k <- charmatch(arg.name,nm,0)	#k>0 : named arg found
    ## Added karg <= length(x) 9Apr02 for R; R doesn't return NULL
    ## like S+
	if(k>0) x[[k]] else 
	if(nm[karg]!="") NULL else if(karg <= length(x)) x[[karg]] else NULL
  }
  if(parms.allowed) parms <- argu(x,2,"parms",pa,nam) else {
	parms <- NULL
	if(charmatch("parms",nam,0)>0)
      stop(paste("parms not allowed for",as.character(call.args[1])))
  }
 
  nm <- argu(x,5,"name",pa,nam)
  if(!is.null(nm)) name <- nm
  if(!is.null(.Options$Design.attr)) {
	atr <- .Options$Design.attr
	i <- charmatch(name, atr$name, 0)
	if(is.null(i))stop("program logic error for options(factor.number)")
	parmi <- atr$parms[[name]]
	return(list(name=atr$name[i],parms=parmi,label=atr$label[i],
                units=atr$units[i]))		# added units 9Jun99
  }

  label <- argu(x,3,"label",pa,nam)
  atx <- attributes(x[[1]])  # 9Jun99
  if(is.null(label)) label <- atx$label   # 9Jun99 attr(x[[1]],"label")
  if(is.null(label)) label <- name

  list(name=name,parms=parms,label=label,units=atx$units)  #9Jun99
  
}

# Function to list all attributes of new sub-design matrix
set.atr <- function(xd,x,z,colnames,assume,code,parms,nonlinear)	{
#Note: x argument isn't used
  ## Added z$units 9Jun99
	if(is.matrix(xd))
	   list(dim=dim(xd),dimnames=list(NULL,colnames),class="Design",
		name=z$name, label=z$label, assume=assume, assume.code=code,
		parms=parms, 
		nonlinear=nonlinear,colnames=colnames,units=z$units)
	else list(dim=dim(xd),class="Design",
		name=z$name, label=z$label, assume=assume, assume.code=code,
		parms=parms, 
		nonlinear=nonlinear,colnames=colnames,units=z$units)
									}


# asis transformation - no transformation	
asis<-function(...)							{

cal <- sys.call()
xx <- list(...)
z <- des.args(xx, FALSE, cal)
xd <- xx[[1]]
if(is.factor(xd)) attr(xd,"class") <- NULL
if(!(is.numeric(xd) | is.logical(xd))) stop(paste(z$name,"is not numeric"))

storage.mode(xd) <- "single"
attributes(xd) <- set.atr(xd,xd,z,z$name,"asis",1,NULL,FALSE)
xd
								}

# matrx transformation - no transformation, keep original vars as matrix
# column names as parameter names, parms=column medians
matrx <- function(...)							{

cal <- sys.call()
xx <- list(...)
z <- des.args(xx, FALSE, cal)

xd <- xx[[1]]
nc <- ncol(xd)
if(!is.matrix(xd)) stop(paste(z$name,"is not a matrix"))
storage.mode(xd) <- "single"
colname <- dimnames(xd)[[2]]
if(length(colname)==0) colname <- paste(z$name,'[',1:nc,']',sep="") else
if(z$label==z$name) z$label <- paste(colname,collapse=",")
parms <- single(nc)
for(i in 1:nc)parms[i] <- median(xd[,i], na.rm=TRUE)

attributes(xd) <- set.atr(xd,NULL,z,colname,"matrix",10,parms,rep(FALSE,nc))
xd
									}

# Polynomial expansion

pol <- function(...)							{

cal <- sys.call()
xx <- list(...)
z <- des.args(xx,TRUE,cal)
x <- xx[[1]]
if(!is.numeric(x)) stop(paste(z$name,"is not numeric"))
poly.degree <- .Options$poly.degree
if(is.null(poly.degree)) poly.degree <- 2
if(!is.null(z$parms)) poly.degree <- z$parms
if(poly.degree<2)stop("order for polynomial must be 2,3,...")
xd <- matrix(single(1),nrow=length(x),ncol=poly.degree)
nam <- z$name
name <- character(poly.degree)
name[1] <- nam
xd[,1] <- x

for(j in 2:poly.degree)	{
	name[j] <- paste(nam,"^",j,sep="")
	xd[,j] <- x^j	}

attributes(xd) <- set.atr(xd,x,z,name,"polynomial",2,poly.degree,
	c(FALSE,rep(TRUE,poly.degree-1)))
xd
									}


# Linear spline expansion

lsp <- function(...)						{

cal <- sys.call()
xx <- list(...)
z <- des.args(xx,TRUE,cal)
x <- xx[[1]]
if(!is.numeric(x)) stop(paste(z$name,"is not numeric"))
parms <- z$parms
if(is.null(parms) || any(is.na(parms))) 
	stop("must specify knots for linear spline")
suffix <- NULL
nam <- z$name
lp <- length(parms)
xd <- matrix(single(1),nrow=length(x),ncol=lp+1)
name <- character(lp+1)
xd[,1] <- x
name[1] <- nam
for(j in 1:lp)	{
	suffix <- paste(suffix,"'",sep="")
	name[j+1] <- paste(nam,suffix,sep="")
	xd[,j+1] <- pmax(x-parms[j],0)
		}

attributes(xd) <- set.atr(xd,x,z,name,"lspline",3,parms,c(FALSE,rep(TRUE,lp)))
xd
									}

# Restricted cubic spline expansion
rcs <- function(...)							{

cal <- sys.call()
xx <- list(...)
z <- des.args(xx,TRUE,cal)
x <- xx[[1]]
if(!is.numeric(x)) stop(paste(z$name,"is not numeric"))
nknots <- .Options$nknots
if(is.null(nknots)) nknots <- 5
parms <- z$parms
if(is.null(parms)) parms <- nknots
if(length(parms)==1)	{
	nknots <- parms
	knots <- NULL	}
else {	nknots <- length(parms)
	knots <- parms	}
if(is.null(knots))	{
	xd <- rcspline.eval(x,nk=nknots,inclx=TRUE)
	knots <- attr(xd,"knots")
			} else
xd <- rcspline.eval(x,knots=knots,inclx=TRUE)
parms <- knots
nknots <- length(parms)
name <- character(nknots-1)
nam <- z$name
name[1] <- nam
suffix <- NULL
for(j in 1:(nknots-2))	{
	suffix <- paste(suffix,"'",sep="")
	name[j+1] <- paste(nam,suffix,sep="")
			}

attributes(xd) <- set.atr(xd,x,z,name,"rcspline",4,parms,
	c(FALSE,rep(TRUE,nknots-2)))
xd
									}

# Category variable
catg <- function(...)							{
cal <- sys.call()
xx <- list(...)
z <- des.args(xx,TRUE,cal)
nam <- z$name
y <- xx[[1]]
parms <- z$parms
if(is.category(y) & !is.factor(y)) oldClass(y) <- "factor"
if(is.null(parms) & is.category(y)) parms <- levels(y)
if(is.null(parms)) 						{
      if(is.character(y)) parms <- sort(unique(y[y!="" & y!=" "]))
      else parms <- as.character(sort(unique(y[!is.na(y)])))
								}
if(!is.factor(y))x <- factor(y, levels=parms) else x <- y
if((is.character(y) && any(y!="" & y!=" " & is.na(x))) ||
   (is.numeric(y) & any(!is.na(y) & is.na(x))))
      stop(paste(nam,"has non-allowable values"))
if(all(is.na(x)))stop(paste(nam,"has no non-missing observations"))
lp <- length(parms)
if(lp<2)stop(paste(nam,"has <2 category levels"))
attributes(x) <- list(levels=parms,class=c("factor","Design"),
	name=nam,label=z$label,assume="category",assume.code=5,
	parms=parms,nonlinear=rep(FALSE,lp-1),
	colnames=paste(nam,"=",parms[-1],sep=""))
x
									}

# Scored expansion   parms=unique values

scored <- function(...)							{

cal <- sys.call()
xx <- list(...)
z <- des.args(xx,TRUE,cal)
parms <- z$parms
nam <- z$name
x <- xx[[1]]
if(is.category(x))			{
	levx <- as.single(levels(x))
	if(any(is.na(levx))) stop(paste("levels for",nam,"not numeric"))
	if(is.null(parms)) parms <- levx
#	.Options$warn <- -1   #suppress warning about NAs  14Sep00
    oldopt <- options(warn=-1)
    on.exit(options(oldopt))
	x <- levx[x]			}
if(!is.numeric(x))stop(paste(nam,"is not a numeric variable"))

y <- sort(unique(x[!is.na(x)]))
if(is.null(parms)) parms <- y

parms <- sort(parms)
n.unique <- length(parms)
if(n.unique<3) stop("scored specified with <3 levels")
lp <- length(parms)-1

#Form contrast matrix of the form linear | dummy | dummy ...

xd <- matrix(single(1),nrow=length(y),ncol=lp)
xd[,1] <- y
name <- character(lp)
name[1] <- nam
i <- 1
for(k in parms[3:length(parms)])	{
	i <- i+1
	name[i] <- paste(nam,"=",k,sep="")
	xd[,i] <- y==k			}

dimnames(xd) <- list(NULL, name)

x <- ordered(x)
if(!.SV4.) oldClass(x) <- c("ordered","factor","Design")  # 17Jul01
attributes(x) <- c(attributes(x),
	list(name=nam,label=z$label,assume="scored",assume.code=7,
	parms=parms,
	nonlinear=c(FALSE,rep(TRUE,lp-1)),colnames=name,contrasts=xd))
x
									}

# strat	parms=value labels

strat <- function(...)							{

cal <- sys.call()
xx <- list(...)
y <- xx[[1]]
z <- des.args(xx,TRUE,cal)
if(is.category(y) & !is.factor(y)) oldClass(y) <- "factor"
parms <- z$parms
if(is.null(parms)) parms <- levels(y)
if(is.null(parms))				         {
   if(is.character(y)) parms <- sort(unique(y[y!="" & y!=" "]))
   else parms <- as.character(sort(unique(y[!is.na(y)])))}
nam <- z$name
if(!is.factor(y)) x <- factor(y,levels=parms) else x <- y
if((is.character(y) & any(y!="" & y!=" " & is.na(x))) ||
   (is.numeric(y) & any(!is.na(y) & is.na(x))))
      stop(paste(nam," has a non-allowable value"))
name <- nam
attributes(x) <- list(levels=parms,class=c("factor","Design"),
	name=nam, label=z$label, assume="strata", assume.code=8,
	parms=parms, nonlinear=FALSE)
x
									}

#Function to subscript a variable, keeping attributes
#Is similar to [.smooth, but does not keep attribute NAs

"[.Design" <- function(x, ..., drop = FALSE)		{
	ats <- attributes(x)
	ats$dimnames <- NULL
	ats$dim <- NULL
	ats$names <- NULL
	oldClass(x) <- NULL
	y <- x[..., drop = drop]
	attributes(y) <- c(attributes(y), ats)
	y						}

#Function to get parms of factor in fit or design object "fit" with name
#given by second argument (without quotes)

gparms <- function(fit,...)	{
  name <- as.character(sys.call())[3]
  atr <- fit$Design
  if(!length(atr)) atr <- getOldDesign(fit)   # 13Apr01
  atr$parms[[name]]
}

#value.chk - if x=NA, returns list of possible values of factor i defined
#	in object f's attributes.  For continuous factors, returns n values
#	in default prediction range.  Use n=0 to return trio of effect
#	limits.  Use n<0 to return pretty(plotting range,nint=-n).
#       If type.range="full" uses the full range instead of default plot rng.
# If x is not NA, checks that list to see that each value is allowable
#	for the factor type, and returns x
# Last argument is object returned from Getlim (see Design.Misc)
# First argument is Design list

value.chk <- function(f, i, x, n, limval, type.range="plot")		{

as     <- f$assume.code[i]
name   <- f$name[i]
parms  <- f$parms[[name]]
isna   <- length(x)==1 && is.na(x)
values <- limval$values[[name]]
charval <- !is.null(values) && is.character(values)
if(isna & as!=7) {
   if(is.null(limval) || match(name, dimnames(limval$limits)[[2]], 0)==0 ||
	is.na(limval$limits["Adjust to",name]))
	stop(paste("variable",name,"does not have limits defined by datadist"))
   limits <- limval$limits[,name]
   lim    <- if(type.range=="full") limits[6:7] else limits[4:5]
}

if(as<5 | as==6) {
  if(isna) {
    if(is.null(values))	{
      if(n==0) x <- limits[1:3]	else {
        if(n>0) x <- seq(oldUnclass(lim[1]), #handles chron
		oldUnclass(lim[2]),length=n)
        else x <- pretty(oldUnclass(lim[1:2]), n=-n)
        oldClass(x) <- oldClass(lim)
      }
    } else x <- values
  } else {
    if(is.character(x) && !charval)
	stop(paste("character value not allowed for variable",
	name))   #Allow any numeric value
    if(charval) {
      j <- match(x, values, 0)
      if(any(j==0)) 
	stop(paste("illegal values for categorical variable:",
		paste(x[j==0],collapse=" "),"\nPossible levels:",
		paste(values,collapse=" ")))
      }	
    }
  }
else if(as==5|as==8) {
# if(isna) x <- lim[1]:lim[2]	else
  if(isna) x <- parms	else {
    j <- match(x, parms, 0)  #match converts x to char if needed
    if(any(j==0)) 
	stop(paste("illegal levels for categorical variable:",
	 paste(x[j==0],collapse=" "),"\nPossible levels:",
	 paste(parms,collapse=" ")))
    x
  }
}

else if(as==7) {
  if(isna) x <- parms
  else if(is.character(x))stop(paste("character value not allowed for",
		"variable",name))
  else {
    j <- match(x, parms, 0)
    if(any(j==0)) 
      stop(paste("illegal levels for categorical variable:",
	paste(x[j==0],collapse=" "),"\n","Possible levels:",
	paste(parms,collapse=" ")))
  }
}

invisible(x)
}


# This is a stripped-down version of terms.inner 
# Returns the inner-most variable expression
# Moved to Hmisc 2Apr01
if(FALSE) var.inner <- function(formula)
{
	if(!inherits(formula,"formula")) formula <- attr(formula,"formula")
	if(length(formula) > 2)
		formula[[2]] <- NULL
	maxch <- 100
	z <- .C("all_names",
		list(formula),
		as.integer(FALSE),
		labels = character(maxch),
		n = as.integer(maxch),
		expr = character(maxch),
		as.logical(TRUE),
		NAOK = TRUE)
	z$labels[1:z$n]
}


#ia.operator.s - restricted interaction operators for use with Design
#F. Harrell  8 Nov 91

#Set up proper attributes for a restricted interaction for a model
#such as y ~ rcs(x1) + rcs(x2) + x1 %ia% x2 or x1 %ia% rcs(x2)
#or rcs(x1) %ia% x2

"%ia%" <- function(x1, x2)						{
a1 <- attributes(x1)
a2 <- attributes(x2)
nam <- as.character(sys.call())[-1]

redo <- function(x,nam)	{
   if(is.null(attr(x,"assume.code")))		{
      if(!is.null(oldClass(x)) && oldClass(x)[1]=="ordered")
         x <- scored(x, name=nam)
      else if(is.character(x) | is.category(x)) x <- catg(x, name=nam)
      else if(is.matrix(x)) x <- matrx(x, name=nam)
      else x <- asis(x, name=nam)		}
   ass <- attr(x,"assume.code")
   nam <- attr(x,"name")

   if(ass==5)				{
      colnames <- attr(x,"colnames")
      len <- length(attr(x,"parms"))-1	}
   else if(ass==8)			{
      prm <- attr(x,"parms")
      colnames <- paste(nam,"=",prm[-1],sep="")
      len <- length(prm)-1		}
   else if(ass==7)			{
      prm <- attr(x,"parms")
      colnames <- c(nam,paste(nam,"=",prm[-(1:2)],sep=""))
      len <- length(prm)-1		}
   else							{
      if(is.null(ncol(x)))		{
         len <- 1
         colnames <- nam		} else 	{
         colnames <- dimnames(x)[[2]]
         len <- ncol(x)				}
							}

   attr(x,"colnames") <- colnames
   attr(x,"len") <- len
   if(ass==8) attr(x,"nonlinear") <- rep(FALSE, len)
   x							}

x1 <- redo(x1,nam[1])
x2 <- redo(x2,nam[2])
a1 <- attributes(x1)
a2 <- attributes(x2)
n1 <- a1$colnames
n2 <- a2$colnames
nl1 <- a1$nonlinear
nl2 <- a2$nonlinear
as1 <- a1$assume.code
as2 <- a2$assume.code

l1 <- a1$len
l2 <- a2$len
if(any(nl1) & any(nl2))	nc <- l1+l2-1   
else nc <- l1*l2
if(is.matrix(x1)) nr <- nrow(x1) 
else nr <- length(x1)
x <- matrix(single(1),nrow=nr,ncol=nc)
name <- character(nc)
parms <- matrix(integer(1),nrow=2,ncol=nc+1)
nonlinear <- logical(nc)

k <- 0
if(!is.factor(x1)) x1 <- as.matrix(x1)
if(!is.factor(x2)) x2 <- as.matrix(x2)
for(i in 1:l1)					{
   if(as1==5 | as1==8) x1i <- oldUnclass(x1)==(i+1)
   else x1i <- x1[,i]
   for(j in 1:l2)			{
	#Remove doubly nonlinear terms
	if(nl1[i] & nl2[j]) break
		k <- k + 1
		if(as2==5 | as2==8) x2j <- oldUnclass(x2)==(j+1)
		else x2j <- x2[,j]
		x[,k] <- x1i * x2j
		name[k] <- paste(n1[i],"*",n2[j])
		parms[,k+1] <- c(nl1[i],nl2[j])
		nonlinear[k] <- nl1[i] | nl2[j]
					}	}

dimnames(x) <- list(NULL, name)
attr(x,"ia") <- c(a1$name, a2$name)
attr(x,"parms") <- parms
attr(x,"nonlinear") <- nonlinear
attr(x,"assume.code") <- 9
attr(x,"name") <- paste(a1$name,"*",a2$name)
attr(x,"label") <- attr(x,"name")
attr(x,"colnames") <- name
attr(x,"class") <- "Design"
storage.mode(x) <- "single"
x
}


Function.Design <- function(fit, intercept, digits=max(8,.Options$digits))
{

# .Options$digits <- digits  14Sep00
oldopt <- options(digits=digits)
on.exit(options(oldopt))
at <- fit$Design
if(!length(at)) at <- getOldDesign(fit)
name <- at$name
ac <- at$assume.code
p <- length(name)
nrp <- num.intercepts(fit)
name.main <- name[ac!=9]  #non-intercepts
pm <- length(name.main)
adj.to <- Getlim(at, allow.null=TRUE, need.all=TRUE)$limits['Adjust to',]


chr <- function(y, digits) if(is.category(y)||is.character(y)) 
  paste('"',as.character(y),'"',sep='') else format.sep(y, digits)

adj.to <- unlist(lapply(adj.to,chr,digits=digits))
z <- paste('function(',paste(name.main,'=',adj.to,collapse=','), ') {', sep='')


#f$term.labels does not include strat
TL <- attr(terms(fit),"term.labels")
#Get inner transformations
#from <- c("asis","pol","lsp","rcs","catg","scored","strat","matrx","I")
#from <- paste(from,"(\\(.*\\))",sep="")
from <- c('asis(*)','pol(*)','lsp(*)','rcs(*)','catg(*)','scored(*)',
  'strat(*)','matrx(*)','I(*)')
to   <- rep('*',9)

#trans <- paste("h(",translate(TL[ac!=9], from, "\\1"),")",sep="")  
trans <- paste("h(",sedit(TL[ac!=9], from, to),")",sep="")
#change wrapping function to h()
h <- function(x,...) deparse(substitute(x))
for(i in (1:pm)) trans[i] <- eval(parse(text=trans[i]))
j <- trans!=name.main
if(any(j)) z <- paste(z, paste(name.main[j],'<-',trans[j],collapse=';'),
					  ';',sep='')

interaction <- at$interactions
if(length(interaction)==0) interaction <- 0

parms <- at$parms

Two.Way <- function(prm,Nam,nam.coef,cof,coef,f,varnames,at,digits) {
  i1 <- prm[1,1]; i2 <- prm[2,1]
  num.nl <- any(prm[1,-1])+any(prm[2,-1])
 #If single factor with nonlinear terms, get it as second factor
 #Otherwise, put factor with most # terms as second factor
  rev <- FALSE
  if((num.nl==1 & any(prm[1,-1])) || (length(Nam[[i1]])>length(Nam[[i2]])))
	{ i1 <- i2; i2 <- prm[1,1]; rev <- TRUE }
  N1 <- Nam[[i1]]; N2 <- Nam[[i2]]
  n1 <- nam.coef[[i1]]; n2 <- nam.coef[[i2]]
  v <- ""
  for(j1 in 1:length(N1)) {
   nam1 <- nam.coef[[i1]][j1]
   lN2 <- length(N2)
   cnam <- if(rev) paste(nam.coef[[i2]],"*",nam1) else
		paste(nam1, "*", nam.coef[[i2]])
   mnam <- match(cnam, names(cof), nomatch=0)
   act <- mnam[mnam>0]
   lN2.act <- length(act)
  #Check if restricted interaction between a rcs and another nonlinear
  #var, i.e. >1 2nd term possible, only 1 (linear) there, and at first 
  #nonlinear term of rcs
   if(lN2.act==1 & lN2>1 & at$assume.code[i1]==4 & j1==2) {
     v <- paste(v,"+",N2[1],"*(",sep="")
     cnam <- paste(nam.coef[[if(rev)i2 else i1]][1], "*",
                   nam.coef[[if(rev)i1 else i2]][-1])
     ## rev and re-order -1 1 4Dec00
     vv <- attr(rcspline.restate(at$parms[[at$name[i1]]], c(0, coef[cnam]), 
		x=varnames[i1], digits=digits),'function.text')
	 v <- paste(v, vv, ')', sep='')
     break
  }
   else if(lN2.act==1) {
     vv <- paste(cof[act],"*",N1[j1],"*", N2[mnam>0], sep="")
	 v <- paste(v, vv, sep='')
   }
   else	if(lN2.act>0) {
    vv <- paste("+",N1[j1],"*(",sep="")
	v <- paste(v, vv, sep='')

    if(at$assume.code[i2]==4 & !any(mnam==0)) {
   #rcspline, interaction not restricted
     vv <- attr(rcspline.restate(at$parms[[at$name[i2]]], coef[act], 
                                 x=varnames[i2], digits=digits),
                'function.text')
	 v <- paste(v, vv, ')', sep='')
     }

    else {
      for(j2 in 1:lN2) {
       l <- mnam[j2]
       if(l>0) {	    #not a restricted-out nonlinear term
        if(j2==1 && substring(cof[l],1,1)=="+") cof[l] <- substring(cof[l],2)
        vv <- paste(cof[l],"*",N2[j2],sep="")
        v <- paste(v, vv, sep='')
	   }
  				 }
     v <- paste(v, ")", sep='')
					}
							}	}
 v
									}

Three.Way <- function(prm,Nam,nam.coef,cof,coef,f,at,digits){
  i1 <- prm[1,1]; i2 <- prm[2,1]; i3 <- prm[3,1]
  N1 <- Nam[[i1]]; N2 <- Nam[[i2]]; N3 <- Nam[[i3]]
  v <- ""; l <- 0
  for(j3 in 1:length(N3)) {
    for(j2 in 1:length(N2)) {
	  for(j1 in 1:length(N1)) {
		l <- l+1
		v <- paste(v,cof[l], "*", N1[j1], "*", N2[j2],
				   "*", N3[j3], sep="")
	  }
	}
  }
  v
}


Coef <- fit$coef
if(nrp==1 | !missing(intercept))	{
  cof <- if(missing(intercept))format.sep(Coef[1],digits) else 
          format.sep(intercept,digits)
  z <- paste(z, cof, sep='')
}

Nam <- list();  nam.coef <- list()
assig <- fit$assign  #DesignAssign(at, f$non.slopes, formula(f))  ## 10Apr02

for(i in (1:p)) {
   ass <- ac[i]
   nam <- name[i]
   prm <- at$parms[[nam]]
   if(any(ass==c(5,7,8))) prm <- chr(at$parms[[nam]],digits=digits)

      k <- assig[[TL[i]]]   ## was f$assign 10Apr02
      coef <- Coef[k]
      nam.coef[[i]] <- names(coef)
      cof <- format.sep(coef,digits)
      cof <- ifelse(coef<=0, cof, paste("+", cof, sep=""))

      switch(ass,
         {  nam <- name[i]; Nam[[i]] <- nam
            q <- paste(cof, '*', nam, sep="")
         },

         {     q <- ""; pow <- 1:prm
               nams <- ifelse(pow==1,nam,paste(nam,"^",pow,"",sep=""))
               Nam[[i]] <- nams
               for(j in pow) q <- paste(q, cof[j], "*", nams[j], sep="")
         },

         {  
            q <- paste(cof[1], "*", nam, sep="")
            nams <- nam
            kn <- format.sep(-prm,digits)
            for(j in 1:length(prm))	{
               zz <- paste("pmax(", nam, if(prm[j]<0) "+" else NULL, 
						   if(prm[j]!=0) kn[j] else NULL, 
						  ",0)", sep="")
               nams <- c(nams, zz)
               q <- paste(q, cof[j+1], "*", zz, sep="")
						}
            Nam[[i]] <- nams
		 },

         {  q <- attr(rcspline.restate(prm, coef, x=nam, digits=digits),
					  'function.text')
	    if(coef[1]>=0) q <- paste('+',q,sep='')
            nn <- nam
            for(j in 1:(length(prm)-2))	{
              nam <- paste(nam, "'", sep=""); nn <- c(nn, nam)
            }
            Nam[[i]] <- nn       #Two.Way only needs first name
								 #for 2nd-order ia with 1 d.f. (restr ia)
                                 #Three.Way needs original design matrix
         } ,
         {  nn <- paste('(',nam,'==',prm[-1],')',sep='')
            ## was prm[-1],, - R barked
		    Nam[[i]] <- nn
			q <- ''
            for(j in 1:(length(prm)-1))	{
               vv <- paste(cof[j], nn[j], sep="*")
               q <- paste(q, vv, sep="")
            }
		 },

		 q <- '',
         
         {
            q <- paste(cof[1], "*", nam, sep="")
            nams <- nam
            for(j in 3:length(prm))	{
               zz <- prm[j]
               vv <- paste(cof[j-1], "*(", nam, "==", zz, ")", sep="")
               nams <- c(nams, zz)
               q <- paste(q, vv, sep="")
               }
             Nam[[i]] <- nams 
		 },
         #Strat factor doesn't exist as main effect, but keep variable
         #names and their lengths if they will appear in interactions later
         {  if(!length(Nam[[i]]) && any(interaction==i)) {
            nam.coef[[i]] <- paste(name[i], "=", prm[-1], sep="")
            ## prm[-1] was oprm[-1] 26Mar01; thanks to
            ## Geskus, Ronald <RGeskus@gggd.amsterdam.nl>
            Nam[[i]] <- prm[-1]
		 }
	    q <- "" },

	 {  
	     if(prm[3,1]==0) 
              q <- Two.Way(prm,Nam,nam.coef,cof,coef,fit,	name, at, digits)
              else q <- Three.Way(prm,Nam,nam.coef,cof,coef,fit,at, digits)
            
	 }, 
         {  nam <- names(coef)
            q <- ""
            nam <- paste("(", nam, ")", sep="")
            Nam[[i]] <- nam
            for(j in 1:length(prm))	{
               vv <- paste(cof[j], '*', nam[j], sep="")
               q <- paste(q, vv, sep="")
               }
		 }) 
     z <- paste(z, q, sep='')
}
z <- paste(z, '}')
eval(parse(text=z))
									}

Function.cph <-
  function(fit, intercept=-fit$center, ...)
  Function.Design(fit, intercept=intercept, ...)


sascode <- function(object, file="", append=FALSE) {

chr <- function(y) if(is.category(y)||is.character(y))
  paste('"',as.character(y),'"',sep='') else as.character(y)

n <- names(object)[names(object)!='']
# 14Sep00: In S+ 6.0, can't append to tty
for(i in n) if(file=='') cat(i,'=',chr(object[[i]]),';\n') else
  cat(i,'=',chr(object[[i]]),';\n',file=file, append=append|i>1)

if(.R.) {
  tf <- tempfile()
  dput(object, file=tf)
  object <- scan(file=tf, what='', sep='\n', quiet=TRUE)
  object <- paste(paste(object[3:(length(object)-1)],collapse='\n'),';',sep='')
} else object <- paste(as.character(object[[length(object)]]), ';')

#com <- 'sed -e "s/pmax/max/g" -e "s/pmin/min/g" -e "s/==/=/g" 
#-e "s/<-/=/g" -e "s/\\^/\*\*/g"'
#w <- sys(com, w)
object <- sedit(object, c('pmax','pmin','==','<-','^'),c('max','min','=','=','**'),
           wild.literal=TRUE)
if(file=='') cat(object, sep='\n') else cat(object, sep="\n", file=file, append=TRUE)
invisible()
}

#SCCS 7/10/92 @(#)Surv.s	4.10
# Package up surivival type data as a structure
#
Surv <- function(time, time2, event,
	      type=c('right', 'left', 'interval', 'counting', 'interval2'),
		       origin=0) {
    nn <- length(time)
    ng <- nargs()
    nam <- as.character(sys.call())[-1]
    if (missing(type)) {
      if (ng==1 || ng==2) type <- 'right'
      else if (ng==3)     type <- 'counting'
      else stop("Invalid number of arguments")
	}
    else {
      type <- match.arg(type)
      ng <- ng-1
      if (ng!=3 && (type=='interval' || type =='counting'))
		stop("Wrong number of args for this type of survival data")
      if (ng!=2 && (type=='right' || type=='left' ||  type=='interval2'))
        stop("Wrong number of args for this type of survival data")
	}
    who <- !is.na(time)
    if (ng==1) {
      if (!is.numeric(time)) stop ("Time variable is not numeric")
      ss <- cbind(time, 1)
      dimnames(ss) <- list(NULL, c("time", "status"))
      tvar <- time; svar <- NULL; nam <- nam[1]	#FEH
	}
    else if (type=='right' || type=='left') {
      if(missing(time2) && !missing(event))	time2 <- event
      if (!is.numeric(time)) stop ("Time variable is not numeric")
      if (length(time2) != nn) stop ("Time and status are different lengths")
      if (is.logical(time2)) status <- 1*time2
      else  if (is.numeric(time2)) {
		who2 <- !is.na(time2)
        if(max(time2[who2]) == 2)
          status <- time2 - 1
        else status <- time2
        if(any(status[who2] != 0 & status[who2] != 1))
				stop("Invalid status value")
      }
      else stop ("Invalid status value")
      ss <- cbind(time, status)
      dimnames(ss) <- list(NULL, c("time", "status"))
      tvar <- time; svar <- time2; nam <- nam[1:2]	#FEH
	}
    else if(type == 'counting') {
      if (length(time2) !=nn) stop ("Start and stop are different lengths")
      if (length(event)!=nn) stop ("Start and event are different lengths")
      if (!is.numeric(time))stop("Start time is not numeric")
      if (!is.numeric(time2)) stop("Stop time is not numeric")
      who3 <- who & !is.na(time2)
      if (any (time[who3]>= time2[who3]))stop("Stop time must be > start time")
      tvar <- time2; svar <- event; nam <- nam[2:3]	#FEH
      if (is.logical(event)) status <- 1*event
      else  if (is.numeric(event)) {
          who2 <- !is.na(event)
          status <- event - min(event[who2])
          if(all(status == 0)) status <- status + 1
          if(any(status[who2] != 0 & status[who2] != 1))
            stop("Invalid status value")
		}
		else stop("Invalid status value")
      ss <- cbind(time-origin, time2-origin, status)
    }

    else {
      ##interval censored data
      if(type == "interval2") {
        event <- ifelse(is.na(time), 2, 
                        ifelse(is.na(time2), 0, 
                               ifelse(time == time2, 1,	3)))
        if(any(time[event == 3] > time2[event == 3]))
          stop("Invalid interval: start > stop")
        time <- ifelse(event != 2, time, time2)
        type <- "interval"
      }
      else {
        temp <- event[!is.na(event)]
        if(!is.numeric(temp))
          stop("Status indicator must be numeric")
        if(length(temp) > 0 && any(temp !=
                                   floor(temp) | temp < 0 | 
                                   temp > 3))
          stop("Status indicator must be 0, 1, 2 or 3")
      }
      status <- event
      ss <- cbind(time, ifelse(!is.na(event) & event == 3, time2, 1), status)
      tvar <- time2; svar <- event; nam <- nam[2:3]  ## needs checking FEH
	}

    oldClass(ss) <- "Surv"
    attr(ss, "type")  <- type
    ## Below is all FEH
    uni <- attr(tvar,"units")
    if(is.null(uni)) {
      uni <- "Day"
##	warning('Time variable has no units() attribute. Assuming Day.\nFor cph with surv=T may need to specify time.inc.')
    }
    tlab <- attr(tvar,"label")
    if(is.null(tlab)) tlab <- nam[1]
    elab <- attr(svar,"label")
    if(is.null(elab) & length(nam)>1) elab <- nam[2]
    attr(ss,"units") <- uni
    attr(ss,"time.label") <- tlab
    attr(ss,"event.label") <- elab
    ss
  }

if(FALSE && !.SV4.) as.character.Surv <- function(xx) {
    attr(xx,'class') <- NULL
    type <- attr(xx, 'type')
    if (type=='right') {
	temp <- xx[,2]
	temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "+"," "))
	paste(format(xx[,1]), temp, sep='')
	}
    else if (type=='counting') {
	temp <- xx[,3]
	temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "+"," "))
	paste('(', format(xx[,1]), ',', xx[,2], temp, ']', sep='')
	}
    else if (type=='left') {
	temp <- xx[,2]
	temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "<"," "))
	paste(temp, format(xx[,1]), sep='')
	}
    else {   #interval type
	stat <- xx[,3]
	temp <- c("+", "", "-", "]")[stat+1]
	temp2 <- ifelse(stat==3,
		paste("[", format(xx[,1]), ", ",format(xx[,2]), sep=''),
			 format(xx[,1]))
	ifelse(is.na(stat), "NA", paste(temp2, temp, sep=''))
	}
    }

"[.Surv" <- function(x, ..., drop=FALSE) {  # was i,j 13Nov00
    atr <- attributes(x)
    atr$dim <- NULL; atr$dimnames <- NULL
    if (missing(..2)) {
      cl <- oldClass(x)
      oldClass(x) <- NULL
      x <- NextMethod('[')
      oldClass(x) <- cl
      attributes(x) <- c(attributes(x), atr)
	x
	}
    else {
      oldClass(x) <- NULL
      NextMethod("[")
    }
  }

if(FALSE && !.SV4.) {
  is.na.Surv <- function(x) {
    oldClass(x) <- NULL
    as.vector( (1* is.na(x))%*% rep(1, ncol(x)) >0)
    }

Math.Surv <- function(...)  stop("Invalid operation on a survival time")
Ops.Surv  <- function(...)  stop("Invalid operation on a survival time")
Summary.Surv<-function(...) stop("Invalid operation on a survival time")
is.Surv <- function(x) inherits(x, 'Surv')
}
#main.effect=F to suppress printing main effects when the factor in
#question is involved in any interaction.

anova.Design <- function(object,...,main.effect=FALSE, tol=1e-9, 
						 test=c('F','Chisq'), ss=TRUE)	{

ava <- function(idx,coef,cov,tol) {
	chisq <- coef[idx] %*% solvet(cov[idx,idx], coef[idx], tol=tol)
	c(chisq, length(idx))
  }

obj.name <- as.character(sys.call())[2]
itype <- 1	#Wald stats. Later sense score stats from object$est
misstest <- missing(test)   ## R updates missing  8Apr02
test <- match.arg(test)
is.ols <- inherits(object,'ols') ||
 (length(object$fitFunction) && any(object$fitFunction=='ols')) ##14Nov00 22May01
if(misstest) test <- if(is.ols) 'F' else 'Chisq'
if(!is.ols && test=='F') stop('F-test not allowed for this type of model')
if(!is.ols) ss <- FALSE

at <- object$Design
if(!length(at)) at <- getOldDesign(object)
assign <- object$assign
name <- at$name
nama <- names(assign)[1]
asso <- 1*(nama=="(Intercept)" | nama=="Intercept")
names(assign)[-asso] <- name

ia <- at$interactions
if(!length(ia))nia <- 0 else nia <- ncol(ia)
assume <- at$assume.code
#if(is.null(assume))stop("fit does not have Design information")
parms <- at$parms
f <- length(assume)
dotlist <- if(!.SV4. && !.R.) (((sys.frame())[["..."]])[[1]])[-1] else
{  ## 12Nov00
  ncall <- names(sys.call())[-(1:2)]
  other.arg <- as.character(sys.call())[-(1:2)]
  if(length(other.arg) && length(ncall))
    other.arg <- other.arg[ncall=='']
  other.arg
}
#  (as.character(sys.call()[-1])[-1])[if(length(ncall))ncall=='' else TRUE]

if(length(dotlist)==0) which <- 1:f	else		{
  if(!.SV4. && !.R.) {
	alist <- NULL
	for(i in 1:length(dotlist))
		alist <- c(alist, deparse(dotlist[[i]]))
  } else alist <- dotlist ## 12Nov00
	jw <- charmatch(alist,name,0)
	if(any(jw==0))	stop(paste("factor names not in design: ",
			paste(alist[jw==0],collapse=" ")))
	which <- jw	}

if(length(object$est) && !length(object$u))
	stop("est in fit indicates score statistics but no u in fit")

if(itype==1)	{
	if(!length(object$coefficients))
	   stop("estimates not available for Wald statistics")
	coef <- object$coefficients	}	else
	{
	if(!length(object$u)) stop("score statistics not available")
#	if(pr)cat("\n                     Score Statistics\n\n")
	coef <- object$u	}
np <- length(coef)

#Compute # intercepts to skip in testing
nrp <- num.intercepts(object)
if(itype==2 & nrp!=0)stop("fit score statistics and x are incompatible")

nc <- length(coef)

cov <- Varcov(object, regcoef.only=TRUE)  #Omit row/col for scale parameters

stats <- NULL
lab <- NULL
W <- list()
s <- 0
all.slopes <- rep(FALSE, nc)
all.ia <- rep(FALSE, nc)
all.nonlin <- rep(FALSE, nc)
num.ia <- 0
num.nonlin <- 0
issue.warn <- FALSE

for(i in which)	{
	j <- assume[i]
	parmi <- parms[[name[i]]]
	if(j!=9) low.fact <- i
	else low.fact <- (parmi[,1])[parmi[,1]>0]
	if(!length(names(at$nonlinear))) nl <- at$nonlinear[[i]]
	else nl <- at$nonlinear[[name[i]]]
	if(!length(nl)) nl <- rep(FALSE,length(assign[[name[i]]]))
#Factor no. according to model matrix is 1 + number of non-strata factors
#before this factor
       if(j!=8)	{				#ignore strata
	if(i==1) jfact <- 1
	else jfact <- 1 + sum(assume[1:(i-1)]!=8)
	main.index <- assign[[jfact+asso]]
	nonlin.ia.index <- NULL	#Should not have to be here. Bug in S?
	all.slopes[main.index] <- TRUE
	if(nia==0)ni <- 0 else ni <- sum(ia==i)
	if(nia==0)ni <- 0
	else for(k in 1:ncol(ia)) ni <- ni + !any(is.na(match(low.fact,ia[,k])))
	if(ni==0 | main.effect)			{
		w <- ava(main.index,coef,cov,tol=tol)
        s <- s+1; W[[s]] <- main.index
		stats <- rbind(stats,w)
		lab <- c(lab, name[i])		}
	#If term is involved in any higher order effect, get pooled test
	#by adding in all high-order effects containing this term
	#For 2nd order interaction, look for 3rd order interactions
	#containing both factors
#	nonlin.ia.index <- NULL	#Used to be here.  Bug in S?
	if(ni>0)	{
		ia.index <- NULL
		mm <- (1:f)[assume==9]
		mm <- mm[mm!=i]
		for(k in mm)		{

			parmk <- parms[[name[k]]]
			hi.fact <- parmk[,1]
			m <- match(low.fact, hi.fact)
			if(!any(is.na(m)))		{
				if(k==1)kfact <- 1
				else kfact <- 1 + sum(assume[1:(k-1)]!=8)
				idx <- assign[[kfact+asso]]
				ia.index <- c(ia.index,idx)
				if(ncol(parmk)>1)for(jj in 1:length(m))	{
				  nonlin.ia.index <- c(nonlin.ia.index,
				     idx[parmk[m[jj],-1]==1])		}
				  nonlin.ia.index <- if(length(nonlin.ia.index))
					unique(nonlin.ia.index) else NULL
				#Highest order can be counted twice
#				c(nonlin.ia.index, added 17 Sep 91
							}
									 }
			idx <- c(main.index,ia.index)
			all.slopes[idx] <- TRUE
			w <- ava(idx,coef,cov,tol=tol)
            s <- s+1; W[[s]] <- idx
			stats <- rbind(stats,w)
			lab <- c(lab, paste(name[i], 
			   " (Factor+Higher Order Factors)"))
			#If factor i in >1 interaction, print summary
			#Otherwise, will be printed later
			if(j!=9 & ni>1)	{
				w <- ava(ia.index,coef,cov,tol=tol)
                s <- s+1; W[[s]] <- ia.index
				stats<-rbind(stats,w)
				lab <- c(lab, " All Interactions")
									}
					    				 }
#	if((any(nl) & j!=9) | (j==9 && parmi[3,i]==1))	{
	if(any(nl))					{
		# Tests of adequacy of linear relationship
		idx <- c(main.index[nl], nonlin.ia.index)
		num.nonlin <- num.nonlin+1
		all.nonlin[idx] <- TRUE
		w <- ava(idx,coef,cov,tol=tol)
        s <- s+1; W[[s]] <- idx
		stats <- rbind(stats,w)
		lab <- c(lab, if(!length(nonlin.ia.index))" Nonlinear"
			else " Nonlinear (Factor+Higher Order Factors)")	
							} 
		#If interaction factor involves a non-linear term from an
		#expanded polynomial, lspline, rcspline, or scored factor,
		#do tests to see if a simplification (linear interaction) is
		#adequate.  Do for second order only.
		if(j==9)	{
			num.ia <- num.ia+1
			all.ia[main.index] <- TRUE
			if(parmi[3,1]>0) issue.warn <- TRUE
			if(parmi[3,1]==0 && ncol(parmi)>1)		{
			nonlin.x <- as.logical(parmi[1,2:ncol(parmi)])
			nonlin.y <- as.logical(parmi[2,2:ncol(parmi)])
			nonlin.xy <- nonlin.x | nonlin.y
			nonlin.xandy <- nonlin.x & nonlin.y
			idx <- main.index[nonlin.xy]
			li <- length(idx)
			if(li>0)	{
			  num.nonlin <- num.nonlin+1
			  all.nonlin[idx] <- TRUE
			  w <- ava(idx,coef,cov,tol=tol)
              s <- s+1; W[[s]] <- idx
			  stats<-rbind(stats,w)
			  lab<-c(lab," Nonlinear Interaction : f(A,B) vs. AB")
			  idx <- main.index[nonlin.xandy]
			  li <- length(idx)
			  if(li>0)	{
			  w <- ava(idx,coef,cov,tol=tol)
              s <- s+1; W[[s]] <- idx
			  stats<-rbind(stats,w)
			  lab<-c(lab," f(A,B) vs. Af(B) + Bg(A)")	}
			  idx <- main.index[nonlin.x]
			  li <- length(idx)
			  if(li>0 & any(nonlin.x!=nonlin.xy)) {
			    w <- ava(idx,coef,cov,tol=tol)
                s <- s+1; W[[s]] <- idx
			    stats<-rbind(stats,w)
			    lab<-c(lab,paste(" Nonlinear Interaction in",
			         name[parmi[1,1]],"vs. Af(B)"))	}
			    idx <- main.index[nonlin.y]
			    li <- length(idx)
			  if(li>0 & any(nonlin.y!=nonlin.xy)) {
			    w <- ava(idx,coef,cov,tol=tol)
                s <- s+1; W[[s]] <- idx
			    stats<-rbind(stats,w)
			    lab<-c(lab,paste(" Nonlinear Interaction in",
				 name[parmi[2,1]],"vs. Bg(A)"))	}
									}}  }
									  } }

#If >1 test of adequacy, print pooled test of all nonlinear effects
if(num.nonlin>1)					{
	idx <- (1:nc)[all.nonlin]
	li <- length(idx)
	w <- ava(idx,coef,cov,tol=tol)
    s <- s+1; W[[s]] <- idx
	stats <- rbind(stats,w)
	lab <- c(lab, "TOTAL NONLINEAR")			}
#If >1 test of interaction, print pooled test of all interactions in list
if(num.ia>1)						{
	idx <- (1:nc)[all.ia]
	li <- length(idx)
	w <- ava(idx,coef,cov,tol=tol)
    s <- s+1; W[[s]] <- idx
	stats <- rbind(stats,w)
	lab <- c(lab,"TOTAL INTERACTION")			}
#If >0 test of adequacy and >0 test of interaction, print pooled test of
#all nonlinear and interaction terms
if(num.nonlin>0 & num.ia>0)	{
	idx <- (1:nc)[all.nonlin | all.ia]
	li <- length(idx)
	w <- ava(idx,coef,cov,tol=tol)
    s <- s+1; W[[s]] <- idx
	stats <- rbind(stats,w)
	lab <- c(lab,"TOTAL NONLINEAR + INTERACTION")
				}
#Get total test for all factors listed
	idx <- (1:nc)[all.slopes | all.ia]
	w <- ava(idx,coef,cov,tol=tol)
s <- s+1; W[[s]] <- idx
stats <- rbind(stats,w)
lab <- c(lab,"TOTAL")

statnam <- c('Chi-Square','d.f.')
if(is.ols) {
  sigma2 <- object$stats['Sigma']^2
  dfe <- object$df.residual
}

if(ss) {
  stats <- cbind(stats[,2], stats[,1]*sigma2, stats[,1]*sigma2/stats[,2], 
				 stats[,1])
  statnam <- c('d.f.','Partial SS','MS','Chi-Square')
  stats <- rbind(stats, Error=c(dfe, sigma2*dfe, sigma2, NA))
  s <- s+1; W[[s]] <- NA
  lab <- c(lab, 'ERROR')
}

j <- statnam=='Chi-Square'
dfreg <- stats[,statnam=='d.f.']
if(test=='F') {
  stats[,j] <- stats[,j] / dfreg
  statnam[j] <- 'F'
  stats <- cbind(stats, P=1-pf(stats[,j], dfreg, dfe))
  attr(stats,'df.residual') <- dfe
} else stats <- cbind(stats,1-pchisq(stats[,j], dfreg))

statnam <- c(statnam, 'P')
dimnames(stats) <- list(lab, statnam)
## attr(stats,"formula") <- formula(object$terms) 30may02
attr(stats,'formula') <- formula(object)
## was attr(object$terms,"formula") 17Apr02
attr(stats,"obj.name") <- obj.name
attr(stats,"class") <- if(.SV4.)'anova.Design' else c("anova.Design","matrix")
names(W) <- lab
attr(stats,"which") <- W
attr(stats,"coef.names") <- names(coef)
attr(stats,"non.slopes") <- nrp
if(issue.warn) 
 warning("tests of nonlinear interaction with respect to single component \nvariables ignore 3-way interactions")
stats
}


print.anova.Design <-
  function(x, which=c('none','subscripts','names','dots'), ...) {
  stats <- x
  digits <- c('Chi-Square'=2, F=2, 'd.f.'=0, 'Partial SS'=15, MS=15, P=4)
  cstats <- matrix('', nrow=nrow(stats), ncol=ncol(stats), 
                   dimnames=dimnames(stats))

which <- match.arg(which)

do.which <- which!='none' && length(W <- attr(stats,'which'))
if(do.which) {
  if(which=='subscripts') simplifyr <- function(x) {
    x <- sort(unique(x))
    n <- length(x)
    ranges <- character(n)
    m <- 0
    s <- x
    while(length(s) > 0) {
      j <- s == s[1] + (1:length(s))-1
      m <- m+1
      ranges[m] <- if(sum(j)>1) paste(range(s[j]),collapse='-') else s[1]
      s <- s[!j]
    }
    ranges[1:m]
  }
  
  k <- length(W)
  w <- character(k)
  coef.names <- attr(stats,'coef.names')
  nrp <- attr(stats,'non.slopes')
  for(i in 1:k) {
    z <- W[[i]]
    if(all(is.na(z))) w[i] <- '' else {
      z <- sort(z)
      w[i] <- switch(which,
                     subscripts=paste(simplifyr(z - nrp), collapse=','),
                     names=paste(coef.names[z],collapse=','),
                     dots={
                       dots <- rep(' ',length(coef.names)-nrp)
                       dots[z - nrp] <- '.'
                       paste(dots,collapse='')} )
    }
  }
}
sn <- dimnames(cstats)[[2]]
for(j in 1:ncol(cstats)) cstats[,j] <- format(round(stats[,j], digits[sn[j]]))

cstats[is.na(stats)] <- ''
j <- sn=='P'
cstats[stats[,j] < 0.00005,j] <- '<.0001'
cstats <- cbind(dimnames(stats)[[1]], cstats)

#cstats<-cbind(dimnames(stats)[[1]],format(round(stats[,1],2)),
#	format(stats[,2]),format(round(stats[,3],4)))
dimnames(cstats) <- list(rep("",nrow(stats)),
	c("Factor    ",dimnames(stats)[[2]]))

heading <- paste("                ",
				 if(any(dimnames(stats)[[2]]=='F'))"Analysis of Variance" else
				 "Wald Statistics", "          Response: ", 
		as.character(attr(stats, "formula")[2]), sep = "")
cat(heading,"\n\n")
if(any(sn=='MS')) cstats[cstats[,1]=='TOTAL',1] <- 'REGRESSION'
if(do.which) cstats <- cbind(cstats, Tested=w)
print(cstats,quote=FALSE)
if(do.which && which!='names') {
  cat('\nSubscripts correspond to:\n')
  print(if(nrp > 0)coef.names[-(1:nrp)] else coef.names, quote=FALSE)
}
if(!any(sn=='MS') && length(dfe <- attr(stats,'df.residual'))) 
  cat('\nError d.f.:', dfe, '\n')
invisible()
}

latex.anova.Design <- function(object,
  title=if(under.unix) paste('anova',attr(object,'obj.name'),sep='.') else
   paste("ano",substring(first.word(attr(object,"obj.name")),
                         1,5),sep=""), 
  psmall=TRUE, dec.chisq=2, dec.F=2, dec.ss=NA, dec.ms=NA, dec.P=4, ...) {

## expr in first.word 18Nov00 removed 25May01

rowl <- dimnames(object)[[1]]
#Translate interaction symbol (*) to times symbol
#rowl <- translate(rowl, "*", "$\\\\times$")
rowl <- sedit(rowl, "*", "$\\times$", wild.literal=TRUE)
#Put TOTAL rows in boldface
rowl <- ifelse(substring(rowl,1,5) %in% c("TOTAL","ERROR"), paste("{\\bf",rowl,"}"),rowl)
rowl <- ifelse(substring(rowl,1,1)==" ",
	paste("~~{\\it ",substring(rowl,2),"}",sep=""), rowl) # preserve leading blank
P <- object[,3]

dstats <- as.data.frame(object)
attr(dstats, 'row.names') <- rowl

## 4may03
if(psmall) {
  psml <- !is.na(dstats$P) & dstats$P < 0.00005
  if(any(psml)) dstats$P <- ifelse(is.na(dstats$P),'',ifelse(psml, 
#if(psmall && any(dstats$P <0.00005)) dstats$P <- ifelse(dstats$P <0.00005,
				"$<0.0001$",
				paste("~",format(round(dstats$P,4)),sep="")))
}

digits <- c('Chi-Square'=dec.chisq, F=dec.F, 'd.f.'=0,
            'Partial SS'=dec.ss, MS=dec.ms, P=dec.P)

sn <- dimnames(object)[[2]]
dig <- digits[sn]
sn[sn=='Chi-Square'] <- '\\chi^2'
names(dstats) <- paste('$',sn,'$',sep='')

#dstats <- structure(list("$\\chi^2$"=stats[,1],"$d.f.$"=stats[,2],
#		"$P$"=P), row.names=rowl, class="data.frame")
#Make LaTeX preserve spaces in heading
head <- paste(if(any(sn=='F'))"Analysis of Variance" else "Wald Statistics", "for {\\tt",
	as.character(attr(object,"formula")[2]),"}")
latex(dstats, cdec=dig, title=title, caption=head, rowlabel="",
	col.just=rep('r',length(sn)), ...)
}

text.anova.Design <- function(x, at, cex=.5, font=2, ...)	{

#Note: a bug in text() prevents writing long character strings
ltext <- function(z, line, label, cex = 0.5, font=2, adj = 0) {
  zz <- z
  zz$y <- z$y - ((line - 1) * 1.2 * cex * par("csi") * (
	par("usr")[4] - par("usr")[3]))/(par("fin")[2])
  text(zz, label, cex = cex, adj = adj, font=font)
}

fi <- tempfile()
sink(fi)
print.anova.Design(x)
sink()

k <- if(.R.) scan(fi, list(z=""), sep="\n", quiet=TRUE)$z else
             scan(fi, list(z=""), sep="\n")$z
if(!.R. && existsFunction('unlink')) unlink(fi)

for(l in 1:length(k)) ltext(at, l, k[l], font=font, cex=cex)

invisible(k)
}

plot.anova.Design <- function(x,
    what=c("chisqminusdf","chisq","aic","P","partial R2","remaining R2",
      "proportion R2"),
	xlab=switch(what, chisq=if(.R.)expression(chi^2) else "Chi-square", 
      chisqminusdf=if(.R.)expression(chi^2~-~df) else
       "Chi-Square Minus Degrees of Freedom", 
      aic="Akaike Information Criterion",
      P="P-value",
      "partial R2"=if(.R.)expression(paste("Partial",~R^2)) else "Partial R^2",
      "remaining R2"=if(.R.)expression(paste("Remaining~",R^2,
          "~After Removing Variable")) else
      "Remaining R^2 After Removing Variable",
      "proportion R2"=
      if(.R.)expression(paste("Proportion of Overall",~R^2))
      else "Proportion of Overall R^2"),
	pch=if(FALSE)183 else 16, rm.totals=TRUE, rm.ia=FALSE, rm.other=NULL, newnames,
	sort=c("descending","ascending","none"), pl=TRUE, ...) {

what <- match.arg(what)
sort <- match.arg(sort)

if(.SV4.) x <- matrix(oldUnclass(x), nrow=nrow(x),
                          dimnames=dimnames(x))  ##14Nov00

rm <- c(if(rm.totals) c("TOTAL NONLINEAR","TOTAL NONLINEAR + INTERACTION",
	"TOTAL INTERACTION","TOTAL"), 
        " Nonlinear"," All Interactions", "ERROR", rm.other)
rn <- dimnames(x)[[1]]
rm <- c(rm, rn[substring(rn,2,10)=="Nonlinear"])
k <- !(rn %in% rm)
if(rm.ia) k[grep("\\*", rn)] <- FALSE
an <- x[k,,drop=FALSE]

dof <- an[,'d.f.']
P <- an[,'P']
chisq <- if(any(dimnames(an)[[2]]=='F')) an[,'F']*dof else an[,'Chi-Square']

if(what %in% c("partial R2","remaining R2","proportion R2")) {
  if("Partial SS" %nin% dimnames(x)[[2]])
    stop('to plot R2 you must have an ols model and must not have specified ss=F to anova')
  sse <- x['ERROR','Partial SS']
  ssr <- x['TOTAL','Partial SS']
  sst <- sse + ssr
}

an <- switch(what,
             chisq=chisq,
             chisqminusdf=chisq-dof,
             aic=chisq-2*dof,
             P=P,
             "partial R2" = an[,"Partial SS"]/sst,
             "remaining R2" = (ssr - an[,"Partial SS"]) / sst,
             "proportion R2" = an[,"Partial SS"] / ssr)

if(missing(newnames)) newnames <- sedit(names(an),
	"  (Factor+Higher Order Factors)", "")
names(an) <- newnames
an <- switch(sort, descending=-sort(-an), ascending=sort(an), none=an)

if(pl) dotchart2(an, xlab=xlab, pch=pch, ...)
invisible(an)
}

bj <- function(formula=formula(data), data,
               subset, na.action=na.delete, 
               link="log",	control=NULL,
               method='fit', x=FALSE, y=FALSE, time.inc) {

    call <- match.call()
    m <- match.call(expand=FALSE)
    if(.R.) library(survival)
    m$x <- m$y <- m$control <- m$method <- m$link <- m$time.inc <- NULL
    m$na.action <- na.action
    if(.R.) m$drop.unused.levels <- TRUE  ## 31jul02
    m[[1]] <- as.name("model.frame")
      if(.R.) {
        dul <- .Options$drop.unused.levels
        if(!length(dul) || dul) {
          on.exit(options(drop.unused.levels=dul))
          options(drop.unused.levels=FALSE)
        }
      }
    X <- Design(eval(m, sys.parent()))    # 24Apr01
	if(method=='model.frame') return(X)
    atrx <- attributes(X)
    nact <- atrx$na.action
    Terms <- atrx$terms
    atr <- atrx$Design              # 24Apr01
    if(!is.null(nact$nmiss))
      names(nact$nmiss) <- 
        c(as.character(formula[2]), atr$name[atr$assume.code!=9])

    lnames <- if(.R.) c("logit","probit","cloglog","identity","log","sqrt",
      "1/mu^2","inverse") else dimnames(glm.links)[[2]]

    link <- pmatch(link, lnames, 0)
    if(link==0) stop("invalid link function")
    link <- lnames[link]
    Y <- model.extract(X, "response")
    atY <- attributes(Y)
    ncy <- ncol(Y)
    maxtime <- max(Y[,-ncy])
    nnn <- c(nrow(Y),sum(Y[,ncy]))
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")

    type <- attr(Y, "type")

    linkfun <- if(.R.) make.link(link)$linkfun else
                       glm.links["link", link][[1]]

    if (type != 'right') stop ("Surv type must by 'right' censored")
    Y <- cbind(linkfun(Y[,1]), Y[,2])

    X <- model.matrix(Terms, X)
##	assgn <- attr(X,'assign') 20may02
    assgn <- DesignAssign(atr, 1, Terms)

	if(method=='model.matrix') return(X)

    time.units <- attr(Y, "units")
    if(is.null(time.units)) time.units <- "Day"
    if(missing(time.inc))	{
      time.inc <- switch(time.units,Day=30,Month=1,Year=1,maxtime/10)
      if(time.inc>=maxtime | maxtime/time.inc>25) 
		time.inc <- max(pretty(c(0,maxtime)))/10
	}
    rnam <- dimnames(Y)[[1]]
    dimnames(X) <- list(rnam, c("(Intercept)",atr$colnames))

    n <- nrow(X)
    nvar <- ncol(X)

    fit <- bj.fit(X, Y, control=control)

	if(fit$fail) {
	  cat("Failure in bj.fit\n")
	  return(fit)
	}

    fit$linear.predictors <- matxv(X, fit$coefficients)

    if (length(nact)) fit$na.action <- nact

    fit <- c(fit, list(maxtime=maxtime, units=time.units,
	time.inc=time.inc,non.slopes=1,assign=assgn,fitFunction='bj'))
    oldClass(fit) <-  if(.SV4.)'Design' else c("bj", "Design") ##13Nov00
    fit$terms <- Terms
    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$call <- call
    fit$Design <- atr    # 24Apr01
    if (x)     fit$x <- X
    if (y) {
      oldClass(Y) <- 'Surv'
      attr(Y,'type') <- atY$type
      fit$y <- Y
    }
	scale.pred <- if(link=="log") c("log(T)","Survival Time Ratio") else "T"
	fit$scale.pred <- scale.pred
	fit$link       <- link
    fit
    }

bj.fit <- function(x, y, control = NULL) {
  if(ncol(y) != 2)
	stop("y is not a right-censored Surv object")
  status <- y[, 2]
  yy <- y[, 1]
  iter.max <- control$iter.max
  eps <- control$eps
  trace <- control$trace
  tol <- control$tol
  max.cycle <- control$max.cycle
  if(length(iter.max) == 0)
	iter.max <- 20
  if(length(eps) == 0)
	eps <- 0.001
  if(length(trace) == 0)
	trace <- FALSE
  if(length(tol) == 0)
	tol <- 1e-007
  if(length(max.cycle) == 0)
	max.cycle <- 30
  x <- as.matrix(x)
  if(all(x[, 1] == 1))
	x <- x[, -1, drop = FALSE]
  d <- dim(x)
  nvar <- d[2]
  if(length(nvar) == 0)
	nvar <- 0
  N <- length(yy)
  if(nvar > 0) {
	xbar <- apply(x, 2, mean)
	xm <- x - rep(xbar, if(.R.)rep(N, nvar) else rep.int(N, nvar))
  }
	else xm <- 0
  timeorig <- yy
  order.orig <- 1:N
  dummystrat <- factor(rep(1, N))
  betahat <- rep(0, max(nvar, 1))
  betamatrix <- NULL
  sse <- 0
  n <- 0	##
  ## new stuff
  nonconv <- FALSE	##
  repeat {
	oldbeta <- betahat
	oldsse <- sse
	if(nvar == 0)
	  ypred <- 0
	  else {
		betahat <- solvet(t(xm) %*% xm, t(xm) %*% yy, tol = tol)
		ypred <- x %*% betahat
	  }
	alphahat <- mean(yy - ypred)
	sse <- sum((yy - ypred)^2)
	razlika <- oldsse/sse
	if(trace)
	  cat("iteration = ", n, "   sse ratio = ", format(razlika), "\n")
	n <- n + 1
	if(trace)
	  cat("  alpha = ", format(alphahat), "   beta = ", format(betahat), "\n\n")	##
	ehat <- timeorig - ypred
	if(!nonconv) {
	  if(abs(razlika - 1) <= eps)
		break
		else if(n > iter.max) {
		  cyclesse <- NULL
		  cycleperiod <- 0
		  nonconv <- TRUE
		  firstsse <- sse
		}
	}
	  else {
		betamatrix <- cbind(betamatrix, c(alphahat, betahat))
		cyclesse <- c(cyclesse, sse)
		cycleperiod <- cycleperiod + 1
		if(any(abs(firstsse - cyclesse) < 1e-007)) {
		  cat("\nCycle period = ", cycleperiod, "\n")
		  meanbeta <- apply(betamatrix, 1, mean)
		  alphahat <- meanbeta[1]
		  betahat <- meanbeta[2:length(meanbeta)]
		  ypred <- x %*% betahat
		  ehat <- timeorig - ypred
		  break
		}
		  else if(cycleperiod >= max.cycle)
			break
	  }
	state <- status
	state[ehat == max(ehat)] <- 1
	S <- structure(cbind(ehat, state), class = "Surv", type = "right")
	KM.ehat <- survfit.km(dummystrat, S, conf.type = "none", se.fit = FALSE)
	n.risk <- KM.ehat$n.risk
	surv <- KM.ehat$surv
	repeats <- c(diff( - n.risk), n.risk[length(n.risk)])
	surv <- rep(surv, repeats)
	w <-  - diff(c(1, surv))
	m <- order(ehat,  - status)
	bla <- cumsum((w * ehat[m]))
	bla <- (bla[length(bla)] - bla)/(surv + state[m])	## Put bla back into original order
	bl <- bla
	bl[(1:N)[m]] <- bla
		yhat <- if(nvar == 0) bl else x %*% betahat + bl
	yy[state == 0] <- yhat[state == 0]
  }
  n <- n - 1
  if(nonconv) {
	if(cycleperiod < max.cycle)
	  cat("\nNo convergence in", n, "steps, but cycle found - average beta returned\n")
	  else {
		cat("\nNo convergence in", n, "steps\n")
		return(list(fail = TRUE))
	  }
  }
  f <- list(fail = FALSE, iter = n)
  cof <- if(nvar == 0) alphahat else c(alphahat, betahat)
  dx <- dimnames(x)[[2]]
  if(length(dx) == 0 && nvar > 0)
	dx <- paste("x", 1:nvar, sep = "")
  names(cof) <- c("Intercept", dx)
  f$coefficients <- cof
  ehat.u <- ehat[status == 1]
  edf <- sum(status) - nvar - 1
  sigma <- sqrt(sum((ehat.u - mean(ehat.u))^2)/edf)
  if(nvar > 0) {
	x <- cbind(Intercept = 1, x)[status == 1,  , drop = FALSE]
	f$var <- solvet(t(x) %*% x, tol = tol) * sigma * sigma
  }
	else f$var <- (sigma * sigma)/N
  stats <- c(N, sum(status), nvar, edf, sigma)
  names(stats) <- c("Obs", "Events", "d.f.", "error d.f.", "sigma")
  f$stats <- stats
  if(any(status == 0))
	yy <- structure(yy, class = "impute", imputed = (1:N)[status == 0])
  f$y.imputed <- yy
  f
}

bjplot <- function(fit, which=1:dim(X)[[2]]) {
  if(!all(c('x','y') %in% names(fit)))
	stop('must specify x=T,y=T to bj to use bjplot')
  X <- (fit$x)[,-1,drop=FALSE]
  Y <- fit$y
#  predicted <- fit$linear.predictors
  xnam <- dimnames(X)[[2]]
  yy <- fit$y.imputed
  imp <- is.imputed(yy)
  trans <- if(fit$link=='identity') '' else fit$link

## Do Hillis plot first
  N <- length(fit$y[, 1])
  dummystrat <- factor(rep(1, N))
  S <- resid(fit)
  S[S[, 1] == max(S[, 1]), 2] <- 1
  m <- order(fit$y[, 1],  - fit$y[, 2])
  resd <- S[m, 1]
  cens <- S[m, 2]
  KM.ehat <- survfit.km(dummystrat, S, 
						conf.type = "none", se.fit = FALSE)
  repeats <- c(diff( - KM.ehat$n.risk), KM.ehat$n.risk[length(KM.ehat$n.risk)])
  if(length(KM.ehat$time) != N) {
	time <- rep(KM.ehat$time, repeats)
	surv <- rep(KM.ehat$surv, repeats)
  }
	else {
	  time <- KM.ehat$time
	  surv <- KM.ehat$surv
	}
  u <- runif(N-1, 0, surv[1:(N - 1)])
  w <- approx(surv, time, xout=u, method='constant', f=0)
  t.i <- c(w$y, max(time))
  surv.i <- c(w$x, min(surv))
  residnew <- resd
  residnew[cens == 0] <- t.i[cens == 0]
  retlist <- list(predictor = fit$linear.predictor[m], 
				  x = fit$x[m,  ], res.cens = resd, hillis = residnew, 
				  cens = cens)
  predictor <- fit$linear.predictor[m]
  plot(predictor, resd, type = "n", 
	   xlab = "Linear Predictor", ylab = "Residuals")
  points(predictor[cens == 0], resd[cens == 0], pch = 1)
  points(predictor[cens == 1], resd[cens == 1], pch = 16)
  plot(predictor, residnew, type = 	"n", xlab = "Linear Predictor", 
	   ylab = "Residuals")
  points(predictor[cens == 0], residnew[cens == 0], pch = 1)
  points(predictor[cens == 1], residnew[cens == 1], pch = 16)
 
  for(i in which) {
	xi <- X[,i]
  ##	plot(predicted, residual, xlab = "", ylab
  ##		 = "", type = "n")
  ##	title(xlab = "predicted ", ylab = 
  ##		"residual")	##
  ##points(predicted[status == 0], residual[status == 0], pch = 4)
  ##	points(predicted[status == 1], residual[
  ##		status == 1], pch = 16)
  ##	abline(h = 0)	##

	ry <- range(yy,Y)

	plot(xi, Y[,1], xlab=xnam[i], ylab=paste('Observed',trans,'Time'),
		 type='n', ylim=ry)
	points(xi[!imp], Y[!imp,1], pch=16)
	if(any(imp)) {
	  points(xi[imp],  Y[imp,1],  pch=1)

	  plot(xi, yy, xlab = xnam[i], ylab=paste('Imputed',trans,'Time'), 
		   type = "n", ylim=ry)
	  points(xi[imp],   yy[imp],  pch = 1)
##	points(xi[imp],   Y[imp,1], pch = 1,cex=.48)
	  segments(xi[imp], Y[imp,1], xi[imp], yy[imp])
	points(xi[!imp],  yy[!imp], pch = 16)

	plot(xi, yy, xlab=xnam[i], ylab=paste('Observed or Imputed',trans,'Time'),
		 type='n', ylim=ry)
	points(xi[!imp], yy[!imp], pch=16)
	points(xi[imp],  yy[imp],  pch=1)
	##	plot(as.matrix(x)[, i], residual, type = "n")
	##	points(as.matrix(x)[, i][status == 0], residual[status == 0], pch = 16)
	##	points(as.matrix(x)[, i][status == 1], residual[status == 1], pch = 4)
	}
  }
  invisible(retlist)
}

print.bj <- function(x, digits=4, long=FALSE, ...) {
  
  if(x$fail) {
	  warning(" bj failed, no summary provided\n")
	}
		
#	.Options$digits <- digits 14Sep00
  old <- options(digits=digits)
  on.exit(options(old))
		
	cat("Buckley-James Censored Data Regression\n\n")
	dput(x$call)			; cat('\n')

	if(length(z <- x$na.action)) naprint(z)
	stats <- x$stats
	if(.R.) print(format.sep(stats),quote=FALSE) else print(stats); cat("\n")

	cof <- x$coefficients
	cnames <- names(cof)

    cof <- matrix(rep(cof, 4), ncol = 4)
    dimnames(cof) <- list(cnames, c("Value", "Std. Error", "Z", "Pr(>|Z|)"))
    stds <- sqrt(diag(x$var))
    cof[, 2] <- stds
    cof[, 3] <- cof[, 1]/stds
    cof[, 4] <- 2*pnorm(-abs(cof[,3]))
	print(cof)

	p <- length(cof)
    if(long && p>1) {
	  ss <- diag(1/stds)
	  correl <- ss %*% x$var %*% ss
	  dimnames(correl) <- list(cnames, cnames)
	  cat('\n\nCorrelation Matrix for Parameter Estimates\n\n')
	  ll <- lower.tri(correl)
	  correl[ll] <- format(round(correl[ll], digits=max(digits-2,2)))
	  correl[!ll] <- ""
	  print(correl[-1,  - p, drop = FALSE], quote = FALSE)
	}
  invisible()
}

predict.bj <- 
  function(object, newdata,
           type=c("lp","x","data.frame","terms","adjto","adjto.data.frame",
             "model.frame"),
           se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
           incl.non.slopes, non.slopes, kint=1,
           na.action=na.keep, expand.na=TRUE, center.terms=TRUE, ...)
  predictDesign(object, newdata, type, se.fit, conf.int, conf.type,
                incl.non.slopes, non.slopes, kint,
                na.action, expand.na, center.terms, ...)


residuals.bj <- function(object, 
						 type = c("censored","censored.normalized"), ...) {
    type <- match.arg(type)
    
    y <- object$y
    aty <- attributes(y)
    if('y' %nin% names(object)) stop('did not use y=T with fit')
    ## 13Nov00 - handled not finding y.imputed in place of y
    ncy <- ncol(y)
    r <- y[,-ncy,drop=FALSE] - object$linear.predictors
	if(type=='censored.normalized') r <- r/object$stats['sigma']
    r <- cbind(r, y[,ncy])
    attr(r,'type') <- aty$type
    attr(r,'units') <- ' '
    attr(r,'time.label') <- if(type=='censored') 
	  'Residual' else 'Normalized Residual'
    attr(r,'event.label') <- aty$event.label
    oldClass(r) <- if(.SV4.)'Surv' else
						 c('residuals.bj','Surv') ##13Nov00
    if (!is.null(object$na.action))
	 naresid(object$na.action, r)
    else r
    }


validate.bj <- function(fit,method="boot",B=40,
		   bw=FALSE,rule="aic",type="residual",sls=.05,aics=0,pr=FALSE,
		   dxy=TRUE, tol=1e-7, rel.tolerance=1e-3, maxiter=15, ...) {

  if(!(length(fit$x) && length(fit$y)))
    stop('you must specify x=T and y=T to bj')
xb <- fit$linear.predictors
ny <- dim(fit$y)
nevents <- sum(fit$y[,ny[2]])

#Note: fit$y already has been transformed by the link function by psm

distance <- function(x,y,fit,iter,evalfit=FALSE,fit.orig,dxy=TRUE,
	maxiter=15, tol=1e-7, rel.tolerance=1e-3, ...){
	#Assumes y is matrix with 1st col=time, 2nd=event indicator

  Dxy <- rcorr.cens(x,y)["Dxy"]; z <- Dxy; nam <- "Dxy"
  names(z) <- nam
  z
}

predab.resample(fit, method=method,
		fit=bj.fit2, measure=distance,
		pr=pr, B=B, bw=bw, rule=rule, type=type,  
		dxy=dxy,
		sls=sls, maxiter=maxiter, tol=tol,
		rel.tolerance=rel.tolerance, ...)
								}


bj.fit2 <- function(x,y,iter=0,maxiter=15, 
	init=NULL, rel.tolerance=1e-3, tol=1e-7, ...) {
	e <- y[,2]
	if(sum(e)<1)return(list(fail=TRUE))
	x <- x	#Get around lazy evaluation creating complex expression
	f <- bj.fit(as.matrix(x),y,
		control=list(iter.max=maxiter,eps=rel.tolerance,tol=tol))
	if(f$fail) warning('bj.fit failed')
	f
}
  

latex.bj <- function(...) latexDesign(...)
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



#Resampling optimism of reliability of a Cox survival model
#For predicting survival at a fixed time u, getting grouped K-M estimates
#with avg. of m subjects in a group, or using cutpoints cuts if present,
#e.g. cuts=c(0,.2,.4,.6,.8,1).
#B: # reps  method=see predab.resample
#bw=T to incorporate backward stepdown (using fastbw) with params rule,type,sls
#pr=T to print results of each rep
#what="observed" to get optimism in observed (Kaplan-Meier) survival for
#groups by predicted survival
#what="observed-predicted" to get optimism in KM - Cox - more suitable if
#distributions of predicted survival vary greatly withing quantile groups
#defined from original sample's predicted survival

calibrate.cph <- function(fit,method="boot",u,m=150,cuts,B=40,
		bw=FALSE,rule="aic",
		type="residual",sls=.05,aics=0,
		pr=FALSE,what="observed-predicted",tol=1e-12, ...)	{

call <- match.call()
#.Options$digits <- 3  14Sep00
oldopt <- options(digits=3)
on.exit(options(oldopt))
unit <- fit$units
if(unit=="") unit <- "Day"
ssum <- fit$surv.summary  ##14Nov00
if(!length(ssum)) stop('did not use surv=T for cph( )')
cat("Using Cox survival estimates at ", dimnames(ssum)[[1]][2],
	" ", unit,"s\n",sep="")
surv.by.strata <- ssum[2,,1] #2nd time= at u, all strata
xb <- fit$linear.predictors
if(length(stra <- attr(xb,"strata"))) 
   surv.by.strata <- surv.by.strata[stra]
survival <- surv.by.strata^exp(xb)
if(missing(cuts)) {
  g <- max(1,floor(length(xb)/m))
  cuts <- quantile(c(0,1,survival), seq(0,1,length=g+1),na.rm=TRUE)
}


distance <- function(x,y,fit,iter,u,fit.orig,what="observed",
                     orig.cuts, ...) {
  ##Assumes y is matrix with 1st col=time, 2nd=event indicator, 3rd=strata

  if(sum(y[,2])<5)return(NA)
  surv.by.strata <- fit$surv.summary[2,,1]
  ##2 means to use estimate at first time past t=0 (i.e., at u)
  surv.by.strata <- surv.by.strata[y[,3]] #Get for each stratum in data
  cox <- surv.by.strata^exp(x - fit$center)
  ##Assumes x really= x * beta
  pred.obs <- groupkm(cox,Surv(y[,1],y[,2]),u=u,cuts=orig.cuts)
  if(what=="observed") dist <- pred.obs[,"KM"]	else
  dist <- pred.obs[,"KM"] - pred.obs[,"x"]
  if(iter==0) {
    print(pred.obs)
    ##		assign("pred.obs", pred.obs, where=1)   18Apr01
    storeTemp(pred.obs, "pred.obs")  #Store externally for plotting
  }
  
  dist
}

coxfit <- function(x,y,u,iter=0, ...) {
  etime <- y[,1]
  e <- y[,2]
  stra <- y[,3]
  ##	attr(x,"strata") <- strata
  if(sum(e)<5)return(list(fail=TRUE))
  ##	f <- coxph(x,etime,e,surv="summary",time.inc=u,pr=F)
  x <- x	#Get around lazy evaluation creating complex expression
  f <- if(length(x))
    cph(Surv(etime,e) ~ x + strat(stra), surv=TRUE, time.inc=u) else
  cph(Surv(etime,e) ~ strat(stra), surv=TRUE, time.inc=u)
  ## length(x)==0 case 25apr03
  ##Get predicted survival at times 0, u, 2u, 3u, ...
  attr(f$terms, "Design") <- NULL
  ##	f$non.slopes <- f$assume <- f$assume.code <- f$assign <- f$name <- NULL
  ##Don't fool fastbw called from predab.resample
  f
}

b <- min(10,B)
overall.reps <- max(1,round(B/b)) #Bug in S prevents>10 loops in predab.resample
cat("\nAveraging ", overall.reps," repetitions of B=",b,"\n\n")
rel <- 0
opt <- 0
nrel <- 0
B <- 0

for(i in 1:overall.reps) {

  reliability <-
    predab.resample(fit, method=method,
                    fit=coxfit,measure=distance,
                    pr=pr, B=b, bw=bw, rule=rule, type=type,  
                    u=u, m=m, what=what, sls=sls, aics=aics, strata=TRUE,
                    orig.cuts=cuts, tol=tol, ...)
  n <- reliability[,"n"]
  rel <- rel + n * reliability[,"index.corrected"]
  opt <- opt + n * reliability[,"optimism"]
  nrel <- nrel + n
  B <- B + max(n)	
  ##	cat("Memory used after ",i," overall reps:",memory.size(),"\n")
}

mean.corrected <- rel/nrel
mean.opt <- opt/nrel
rel <- cbind(mean.optimism=mean.opt,mean.corrected=mean.corrected,n=nrel)
cat("\nMean over ",overall.reps," overall replications\n\n")
print(rel)

pred <- pred.obs[,"x"]
KM <- pred.obs[,"KM"]
obs.corrected <- KM - mean.opt

e <- fit$y[,2]

structure(cbind(reliability[,c("index.orig","training","test"),drop=FALSE],
	rel,mean.predicted=pred,KM=KM,
	KM.corrected=obs.corrected,std.err=pred.obs[,"std.err",drop=FALSE]),
	class="calibrate", u=u, units=unit, n=length(e), d=sum(e),
	p=length(fit$coef), m=m, B=B, what=what, call=call)
}

calibrate.default <- function(fit, predy, 
			      method=c("boot","crossvalidation",".632","randomization"),
			      B=40, bw=FALSE, rule=c("aic","p"),
			      type=c("residual","individual"),
			      sls=.05, pr=FALSE, kint, smoother="lowess", ...) {
call <- match.call()
method <- match.arg(method)
rule <- match.arg(rule)
type <- match.arg(type)

ns <- num.intercepts(fit)
if(missing(kint)) kint <- floor((ns+1)/2)
clas <- attr(fit,"class")
model <- if(any(clas=="lrm"))"lr" else if(any(clas=="ols"))"ol" else
  stop("fit must be from lrm or ols")
lev.name <- NULL
yvar.name <- as.character(formula(fit))[2]
y <- fit$y
n <- length(y)
if(length(y)==0) stop("fit did not use x=T,y=T")
if(model=="lr") {
  y <- factor(y); lev.name <- levels(y)[kint+1]; fit$y <- as.integer(y)-1
  ## was category(y)   y-1  11Apr02
}

predicted <- if(model=="lr") 
  1/(1+exp(-(fit$linear.predictors-fit$coef[1]+fit$coef[kint])))  else
fit$linear.predictors
 
if(missing(predy)) {
 if(n<11) stop("must have n>10 if do not specify predy")
  p <- sort(predicted)
  predy <- seq(p[5],p[n-4],length=50)
  p <- NULL
}

penalty.matrix <- fit$penalty.matrix

cal.error <- function(x, y, iter, smoother, predy, kint, model, ...) {
  if(model=="lr") {
    x <- 1/(1+exp(-x))
    y <- y>=kint
  }
  smo <- if(is.function(smoother)) smoother(x,y) else lowess(x,y,iter=0)
  cal <- if(.R.) approx(smo, xout=predy, ties=function(x)x[1])$y else
                 approx(smo, xout=predy)$y
  ## 11Apr01  .R. lowess has duplicates
#  if(iter==0) assign(".orig.cal",cal,where=1)   17Apr01
  if(iter==0) storeTemp(cal,".orig.cal")
  cal-predy
}

fitit <- function(x, y, model, penalty.matrix=NULL, xcol=NULL, ...) {
  if(length(penalty.matrix) && length(xcol)) {
    if(model=='ol') xcol <- xcol[-1]-1   # take off intercept position
    penalty.matrix <- penalty.matrix[xcol,xcol,drop=FALSE]
  }
  switch(model,
	 lr=lrm.fit(x, y, penalty.matrix=penalty.matrix,tol=1e-13),
	 ol=c(if(length(penalty.matrix)==0) lm.fit.qr.bare(x, y, intercept=FALSE) else 
                             lm.pfit(x, y, 
                                     penalty.matrix=penalty.matrix),fail=FALSE))
  ## Was lm.fit.qr 14Sep00
}

z <- predab.resample(fit, method=method, fit=fitit, measure=cal.error,
		     pr=pr, B=B, bw=bw, rule=rule, type=type, sls=sls,
		     non.slopes.in.x=model=="ol",
		     smoother=smoother, predy=predy, model=model, kint=kint,
		     penalty.matrix=penalty.matrix, ...)

z <- cbind(predy, calibrated.orig=.orig.cal,
	          calibrated.corrected=.orig.cal-z[,"optimism"],
	   z)
structure(z, class="calibrate.default", call=call, kint=kint, model=model,
		  lev.name=lev.name, yvar.name=yvar.name, n=n, freq=fit$freq,
		  non.slopes=ns, B=B, method=method, 
		  predicted=as.single(predicted), smoother=smoother)
}

print.calibrate.default <- function(x, ...) {
  at <- attributes(x)
  cat("\nEstimates of Calibration Accuracy by ",at$method," (B=",at$B,")\n\n",
      sep="")
  dput(at$call)
  if(at$model=="lr") {
    lab <- paste("Pr{",at$yvar.name,sep="")
    if(at$non.slopes==1) lab <- paste(lab,"=",at$lev.name,"}",sep="") else
    lab <- paste(lab,">=",at$lev.name,"}",sep="")
  } else lab <- at$yvar.name
  cat("\nPrediction of",lab,"\n\n")
  predicted <- at$predicted
  if(length(predicted)) {  ## for downward compatibility
	s <- !is.na(x[,'predy'] + x[,'calibrated.corrected'])
	err <- predicted - approx(x[s,'predy'],x[s,'calibrated.corrected'], 
							  xout=predicted)$y
	cat('\nn=',length(err),    '   Mean absolute error=',
		format(mean(abs(err),na.rm=TRUE)),'   Mean squared error=',
		format(mean(err^2,na.rm=TRUE)),   '\n0.9 Quantile of absolute error=',
		format(quantile(abs(err),.9,na.rm=TRUE)),	   '\n\n',sep='')
  }
  print.matrix(x)
  invisible()
}

plot.calibrate.default <- function(x, xlab, ylab, xlim, ylim, legend=TRUE, 
                                   subtitles=TRUE, ...){
  at <- attributes(x)
  if(missing(ylab)) ylab <- if(at$model=="lr") "Actual Probability" else
      paste("Observed",at$yvar.name)
  if(missing(xlab)) {
    if(at$model=="lr") {
      xlab <- paste("Predicted Pr{",at$yvar.name,sep="")
      if(at$non.slopes==1) {
        xlab <- if(at$lev.name=="TRUE") paste(xlab,"}",sep="") else
                paste(xlab,"=",at$lev.name,"}",sep="")
      } else
      xlab <- paste(xlab,">=",at$lev.name,"}",sep="")
    } else xlab <- paste("Predicted",at$yvar.name)
  }
  p <- x[,"predy"]
  p.app <- x[,"calibrated.orig"]
  p.cal <- x[,"calibrated.corrected"]
  if(missing(xlim) & missing(ylim)) xlim <- ylim <- range(c(p,p.app,p.cal),
							  na.rm=TRUE) else {
    if(missing(xlim)) xlim <- range(p)
    if(missing(ylim)) ylim <- range(c(p.app,p.cal,na.rm=TRUE))
  }
  plot(p, p.app, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, type="n", ...)
  predicted <- at$predicted
  if(length(predicted)) {  ## for downward compatibility
	s <- !is.na(p + p.cal)
	err <- predicted - approx(p[s],p.cal[s],xout=predicted)$y
	cat('\nn=',n <- length(err),    '   Mean absolute error=',
		format(mae <- mean(abs(err),na.rm=TRUE)),'   Mean squared error=',
		format(mean(err^2,na.rm=TRUE)),   '\n0.9 Quantile of absolute error=',
		format(quantile(abs(err),.9,na.rm=TRUE)),	   '\n\n',sep='')
	if(subtitles) title(sub=paste('Mean absolute error=',format(mae),
						  ' n=',n,sep=''), cex=.65, adj=1)
	scat1d(predicted)
  }
  lines(p, p.app, lty=3)
  lines(p, p.cal, lty=1)
  abline(a=0,b=1,lty=2)
  if(subtitles) title(sub=paste("B=",at$B,"repetitions,",at$method),adj=0)
  if(!(is.logical(legend) && !legend)) {
    if(is.logical(legend)) legend <- list(x=xlim[1]+.55*diff(xlim), #was .57
					    y=ylim[1]+.32*diff(ylim))
    legend(legend, c("Apparent","Bias-corrected","Ideal"),
	   lty=c(3,1,2), bty="n")
  }
  invisible()
}


calibrate.psm <- function(fit,method="boot",u,m=150,cuts,B=40,
		bw=FALSE,rule="aic",
		type="residual",sls=.05,aics=0,
		pr=FALSE,what="observed-predicted",tol=1e-12, maxiter=15, 
	    rel.tolerance=1e-5, ...) {

call <- match.call()
if(!length(fit$y)) stop("fit did not store y")
# .Options$digits <- 3  14Sep00
oldopt <- options(digits=3)
on.exit(options(oldopt))
unit <- fit$units
if(unit=="") unit <- "Day"
ny <- dim(fit$y)
nevents <- sum(fit$y[,ny[2]])

#Note: fit$y already has been transformed by the link function by psm

if(missing(cuts))	{
  g <- max(1,floor(ny[1]/m))
  survival <- survest.psm(fit,times=u,conf.int=FALSE)$surv
  cuts <- quantile(c(0,1,survival), seq(0,1,length=g+1),na.rm=TRUE)
}

if(.newSurvival.) {
  dist <- fit$dist
  if(!.R.) survreg.distributions <- survReg.distributions
  inverse <- survreg.distributions[[dist]]$itrans
  if(!length(inverse)) inverse <- function(x) x
  parms <- fit$parms
} else {
  link <- fit$family["link"]
  inverse <- glm.links["inverse",link][[1]]
  family <- fit$family
  fixed <- fit$fixed
  if(length(fixed)==1 && is.logical(fixed) && !fixed) fixed <- list()
  dist <- fit$family["name"]
}


distance <- function(x,y,fit,iter,u,fit.orig,what="observed",inverse,
                     orig.cuts, ...)						{
  ##Assumes y is matrix with 1st col=time, 2nd=event indicator
  if(sum(y[,2])<5)return(NA)
  oldClass(fit) <- 'psm'   # for survest.psm which uses Survival.psm
  if(.newSurvival.) fit$dist <- fit.orig$dist
  psurv <- survest.psm(fit, linear.predictors=x, times=u, conf.int=FALSE)$surv
  ##Assumes x really= x * beta
  pred.obs <- 
    groupkm(psurv,Surv(inverse(y[,1]),y[,2]),u=u,cuts=orig.cuts)
  if(what=="observed") dist <- pred.obs[,"KM"]	else
  dist <- pred.obs[,"KM"] - pred.obs[,"x"]
  if(iter==0) 	{
    print(pred.obs)
    storeTemp(pred.obs)
    ##Store externally for plotting
  }
  dist
}

if(FALSE) {
sfit <- function(x,y,u,iter=0,dist,fixed,family,tol=1e-12,
				 rel.tolerance=1e-5,maxiter=15,...) 	{
  e <- y[,2]
  if(sum(e)<5)return(list(fail=TRUE))
  x <- x	#Get around lazy evaluation creating complex expression
  f <- survreg.fit(as.matrix(x),y,dist=dist,
                   fixed=fixed,offset=rep(0,length(e)),init=NULL,
                   controlvals=
                   survreg.control(failure=2,maxiter=maxiter,
                                   rel.tolerance=rel.tolerance)) #1Jul96
  if(is.character(f)) {warning(f); return(list(fail=TRUE)) }
  
  if(!length(f$coefficients)) {   ## 14Aug01
    f$coefficients <- f$coef
    f$coef <- NULL
  }

  p <- length(f$coefficients)            ## was coef 14Aug01
  f$coefficients <- f$coefficients[-p]   ## " "
  f$family <- family
  f$fail <- FALSE
  f$var <- solvet(f$imat[-p,-p,drop=FALSE], tol=tol)
  f
}
NULL
}

b <- min(10,B)
overall.reps <- max(1,round(B/b)) #Bug in S prevents>10 loops in predab.resample
cat("\nAveraging ", overall.reps," repetitions of B=",b,"\n\n")
rel <- 0
opt <- 0
nrel <- 0
B <- 0

for(i in 1:overall.reps)		{

  reliability <- predab.resample(fit, method=method,
                                 fit=survreg.fit2,measure=distance,
                                 pr=pr, B=b, bw=bw, rule=rule, type=type,  
                                 u=u, m=m, what=what, 
                                 dist=dist, inverse=inverse, parms=parms,
                                 fixed=fixed, family=family,
                                 sls=sls, aics=aics, strata=FALSE,
                                 tol=tol, orig.cuts=cuts, maxiter=maxiter,
                                 rel.tolerance=rel.tolerance, ...)
  n <- reliability[,"n"]
  rel <- rel + n * reliability[,"index.corrected"]
  opt <- opt + n * reliability[,"optimism"]
  nrel <- nrel + n
  B <- B + max(n)	
  print(reliability)
  ##	cat("Memory used after ",i," overall reps:",memory.size(),"\n")
					}

mean.corrected <- rel/nrel
mean.opt <- opt/nrel
rel <- cbind(mean.optimism=mean.opt,mean.corrected=mean.corrected,n=nrel)
cat("\nMean over ",overall.reps," overall replications\n\n")
print(rel)

pred <- pred.obs[,"x"]
KM <- pred.obs[,"KM"]
se <- pred.obs[,"std.err"]
obs.corrected <- KM - mean.opt

structure(cbind(reliability[,c("index.orig","training","test"),drop=FALSE],
	rel,mean.predicted=pred,KM=KM,
	KM.corrected=obs.corrected,std.err=se), 
	class="calibrate", u=u, units=unit, n=ny[1], d=nevents, 
	p=length(fit$coef)-1, m=m, B=B, what=what, call=call)
}







calibrate <- function(fit, ...) 
   UseMethod("calibrate")

print.calibrate <- function(x, ...) {
 at <- attributes(x)
 dput(at$call); cat("\n")
 print.matrix(x)
 invisible()
}

plot.calibrate <- function(x, xlab, ylab, subtitles=TRUE,
                           conf.int=TRUE, ...) { 
   at <- attributes(x)
   u <- at$u
   units <- at$units
   pred <- x[,"mean.predicted"]
   KM <- x[,"KM"]
   obs.corrected <- x[,"KM.corrected"]
   se <- x[,"std.err"]
   if(missing(xlab)) xlab <- paste("Predicted ",format(u),units,"Survival")
   if(missing(ylab)) ylab <- paste("Fraction Surviving ",format(u)," ",units,
				"s",sep="")

#Remember that groupkm stored the se of the log(-log) survival
if(conf.int) errbar(pred, KM,
	ifelse(KM==0 | KM==1, NA, exp(-exp(logb(-logb(KM))-1.96*se))),
	ifelse(KM==0 | KM==1, NA, exp(-exp(logb(-logb(KM))+1.96*se))),
	xlab=xlab, ylab=ylab, type="b", ...)
else plot(pred, KM, xlab=xlab, ylab=ylab, type="b", ...)

if(subtitles)	{
	title(sub=paste("n=",at$n," d=",at$d," p=",at$p,
		", ",at$m," subjects per group",sep=""),adj=0,cex=1)
	title(sub=paste("X - resampling optimism added, B=",at$B,
		"\nBased on ",at$what,sep=""),
		adj=1,cex=1)
		}
abline(0,1,lty=2)	#ideal
points(pred, obs.corrected, pch=4)

invisible()
}

contrast <- function(fit, ...) UseMethod("contrast")

contrast.Design <- function(fit, a, b, cnames=NULL,
                            type=c('individual','average'),
                            weights='equal', conf.int=0.95, ...) {
type <- match.arg(type)

zcrit <- if(length(idf <- fit$df.residual)) qt((1+conf.int)/2, idf) else
         qnorm((1+conf.int)/2)

da <- do.call('gendata', list(fit, factors=a))
xa <- predict(fit, da, type='x')
ma <- nrow(xa)

if(missing(b)) {
  xb <- 0*xa
  db <- da
} else {
  db <- do.call('gendata', list(fit, factors=b))
  xb <- predict(fit, db, type='x')
}
mb <- nrow(xb)

vary <- NULL
if(type=='individual' && !length(cnames)) {
  ## If two lists have same length, label contrasts by any variable
  ## that has the same length and values in both lists
  if(ma==mb) {
    if(ncol(da) != ncol(db)) stop('program logic error')
    if(any(sort(names(da)) != sort(names(db)))) stop('program logic error')
    k <- integer(0)
    nam <- names(da)
    for(j in 1:length(da)) {
      if(all(as.character(da[[nam[j]]])==as.character(db[[nam[j]]])))
        k <- c(k, j)
    }
    if(length(k)) vary <- da[k]
  } else if(max(ma,mb)>1) {
  ## Label contrasts by values of longest variable in list if
  ## it has the same length as the expanded design matrix
    d <- if(ma>1) a else b
    l <- sapply(d, length)
    vary <- if(sum(l==max(ma,mb))==1)d[l==max(ma,mb)]
  }
}

if(max(ma,mb)>1 && min(ma,mb)==1) {
  if(ma==1) xa <- matrix(xa, nrow=mb, ncol=ncol(xb), byrow=TRUE) else
            xb <- matrix(xb, nrow=ma, ncol=ncol(xa), byrow=TRUE)
} else if(mb != ma)
  stop('number of rows must be the same for observations generated\nby a and b unless one has one observation')

X <- xa - xb
p <- ncol(X)
m <- nrow(X)

if(is.character(weights)) {
  if(weights!='equal') stop('weights must be "equal" or a numeric vector')
  weights <- rep(1, m)
} else if(length(weights) > 1 && type=='individual')
    stop('can specify more than one weight only for type="average"')
  else if(length(weights) != m) stop(paste('there must be',m,'weights'))
weights <- as.vector(weights)  # 28Oct99
if(m > 1 && type=='average')
  X <- matrix(apply(weights*X, 2, sum) / sum(weights), nrow=1,
              dimnames=list(NULL,dimnames(X)[[2]]))

est <- drop(X %*% coef(fit))
v <- drop(X %*% Varcov(fit, regcoef.only=FALSE) %*% t(X))   # 28Apr99 regcoef
ndf <- if(is.matrix(v))nrow(v) else 1
se <- if(ndf==1) sqrt(v) else sqrt(diag(v))
Z <- est/se
P <- if(length(idf)) 2*(1-pt(abs(Z), idf)) else 2*(1-pnorm(abs(Z)))
res <- list(Contrast=est, SE=se,
            Lower=est - zcrit*se, Upper=est + zcrit*se,
            Z=Z, Pvalue=P, 
            var=v, df.residual=idf,
            X=X, 
            cnames=if(type=='average')NULL else cnames, nvary=length(vary))
if(type=='individual') res <- c(vary, res)
structure(res, class='contrast.Design')
}

print.contrast.Design <- function(x, X=FALSE, fun=function(u)u, ...) {
  edf <- x$df.residual
  sn <- if(length(edf))'t' else 'Z'
  pn <- if(length(edf))'Pr(>|t|)' else 'Pr(>|z|)'
  w <- x[1:(x$nvary + 6)]
  w$Z <- round(w$Z, 2)
  w$Pvalue <- round(w$Pvalue, 4)
  no <- names(w)
  no[no=='SE'] <- 'S.E.'
  no[no=='Z'] <- sn
  no[no=='Pvalue'] <- pn
  names(w) <- no
  
  cnames <- x$cnames
  if(!length(cnames)) cnames <- if(x$nvary)rep('',length(x[[1]])) else
    as.character(1:length(x[[1]]))
  attr(w,'row.names') <- cnames
  attr(w,'class') <- 'data.frame'
  w$Contrast <- fun(w$Contrast)
  w$SE       <- fun(w$SE)
  w$Lower    <- fun(w$Lower)
  w$Upper    <- fun(w$Upper)
  print(as.matrix(w),quote=FALSE)  ## was print(w) 20may02
  if(length(edf))cat('\nError d.f.=',edf,'\n')
  if(X) {
    cat('\nDesign Matrix for Contrasts\n\n')
    if(is.matrix(x$X)) dimnames(x$X) <- list(cnames, dimnames(x$X)[[2]])
    print(x$X)
  }
  invisible()
}

cph <- function(formula=formula(data),
                data=if(.R.)parent.frame() else sys.parent(),
                weights, subset, na.action=na.delete, 
                method=c("efron","breslow","exact",
                  "model.frame", "model.matrix"),
                singular.ok=FALSE, robust=FALSE,
                model=FALSE, x=FALSE, y=FALSE, se.fit=FALSE,
                eps=.0001, init, iter.max=10, tol=1e-9,
                surv=FALSE, time.inc, type, vartype, conf.type, ...) {

  if(.R.) require('survival')
  method <- match.arg(method)
  call <- match.call()
  m <- match.call(expand=FALSE)
  m$na.action <- na.action    ## 30apr02
  m$method <- m$model <- m$x <- m$y <- m$... <- m$se.fit <-
    m$type <- m$vartype <-
    m$surv <- m$time.inc <- m$eps <- m$init <- m$iter.max <- m$tol <-
      m$weights <- m$singular.ok <- m$robust <- NULL
  m$na.action <- na.action
  if(.R.) m$drop.unused.levels <- TRUE  ## 31jul02
  m[[1]] <- as.name("model.frame")

  if (!inherits(formula,"formula")) {
	## I allow a formula with no right hand side
	## The dummy function stops an annoying warning message "Looking for
	##  'formula' of mode function, ignored one of mode ..."
	if (inherits(formula,"Surv")) {
      xx <- function(x) formula(x)
      formula <- xx(paste(deparse(substitute(formula)), 1, sep="~"))
    }
	else stop("Invalid formula")
  }
  m$formula <- formula

  nstrata <- 0; Strata <- NULL

  if(!missing(data) || (length(z <- attr(terms(formula),"term.labels"))>0 &&
                        any(z!="."))) { #X's present
    if(.R.) {
      dul <- .Options$drop.unused.levels
      if(!length(dul) || dul) {
        on.exit(options(drop.unused.levels=dul))
        options(drop.unused.levels=FALSE)
      }
    }
    X <- Design(eval(m, if(.R.) parent.frame() else sys.parent()))
    atrx <- attributes(X)
    atr  <- atrx$Design
    nact <- atrx$na.action
    if(method=="model.frame") return(X)

#    atr <- attr(attr(X,'terms'),'Design')
    Terms <- if(missing(data)) terms(formula, 
                                     specials=c("strat","cluster")) else
    terms(formula, 
          specials=c("strat","cluster"), data=data)

    asm   <- atr$assume.code
    name  <- atr$name

    cluster <- attr(Terms, "specials")$cluster
    stra <- attr(Terms, "specials")$strat
    Terms.ns <- Terms
    if(length(cluster)) {
      if(missing(robust)) robust <- TRUE
      tempc <- untangle.specials(Terms.ns, "cluster", 1:10)
      ord <- attr(Terms, 'order')[tempc$terms]
      if(any(ord>1)) stop("Cluster can not be used in an interaction")
      cluster <- strata(X[,tempc$vars], shortlabel=TRUE)  #allow multiples
      Terms.ns <- Terms.ns[-tempc$terms]
    }
    if(length(stra)) {
      temp <- untangle.specials(Terms.ns, "strat", 1)
      Terms.ns <- Terms.ns[-temp$terms]	#uses [.terms function
      ##  Set all factors=2 (-> interaction effect not appearing in main effect
      ##  that was deleted strata effect
      if(!.R.) {
        tfac <- attr(Terms,'factors')
        ## For some reason attr(...) <- pmin(attr(...)) changed a detail
        ## in factors attribute in R but doesn't seem to be needed in
        ## R or SV4 anyway 30apr02
        if(length(tfac) && any(tfac > 1))
          attr(Terms,'factors') <- pmin(tfac, 1)
        tfac <- attr(Terms.ns,'factors')
        if(length(tfac) && any(tfac > 1))
          attr(Terms.ns,'factors') <-  pmin(tfac, 1)
      }

      ##added 21Apr94, if( ) 28Apr02
      Strata <- list()
      for(i in (1:length(asm))[asm==8])	{
        nstrata <- nstrata+1
        xi <- X[[i+1]]
        levels(xi) <- paste(name[i],"=",levels(xi),sep="")
        Strata[[nstrata]] <- xi			}
      names(Strata) <- paste("S",1:nstrata,sep="")
      Strata <- interaction(as.data.frame(Strata),drop=TRUE)
    }
    offs <- offset<- attr(Terms, "offset")  ## offs 23nov02 moved up 6dec02
    ## 23nov02
    if(length(nact$nmiss)) {
      jia <- grep('*%ia%*',names(nact$nmiss))  ## 8feb03
      if(length(jia)) nact$nmiss <- nact$nmiss[-jia]
      s <- if(length(offs)) names(nact$nmiss) !=  atrx$names[offs] else TRUE
      names(nact$nmiss)[s] <- 
        c(as.character(formula[2]), atr$name[atr$assume.code!=9])
    }
    ## [s] 23nov02
    xpres <- length(asm) && any(asm!=8)
    Y <- model.extract(X, 'response')
    if(!inherits(Y,"Surv"))stop("response variable should be a Surv object")
    weights <- model.extract(X, 'weights')
    tt <- length(offset)
    offset <- if(tt == 0) rep(0, nrow(Y))
    else if(tt == 1)
      X[[offset]]
    else {
      ff <- X[[offset[1]]]
      for(i in 2:tt)   # for case with multiple offset terms
        ff <- ff + X[[offset[i]]]
      ff
    }
    
    if(model) m <- X

    ##No mf if only strata factors

    if(!xpres)
      {
        X <- if(.R.)matrix(nrow=0,ncol=0) else NULL
        assign <- NULL
      }
    else
      {
        X <- model.matrix(Terms.ns, X)[,-1,drop=FALSE]
        assign <- attr(X, "assign")
        assign[[1]] <- NULL  # remove intercept position, renumber
      }
    nullmod <- FALSE
  }

  else	# model with no right-hand side
    {
      X <- NULL
      Terms <- terms(formula)
      yy <- attr(terms(formula),"variables")[1]
      Y <- eval(yy,data=data)     #sys.parent(2))
      if(!inherits(Y,"Surv"))stop("response variable should be a Surv object")
      Y <- Y[!is.na(Y)]
      assign <- NULL
      xpres <- FALSE
      nullmod <- TRUE
      nact <- NULL
    }
  ny <- ncol(Y)
  time.units <- attr(Y, "units")
  maxtime <- max(Y[,ny-1])

  rnam <- dimnames(Y)[[1]]
  if(xpres) dimnames(X) <- list(rnam, atr$colnames)

  if(method=="model.matrix") return(X)
  if(!length(time.units)) time.units <- "Day"
  if(missing(time.inc)) {
    time.inc <- switch(time.units,Day=30,Month=1,Year=1,maxtime/10)
    if(time.inc>=maxtime | maxtime/time.inc>25) time.inc <- max(pretty(c(0,maxtime)))/10
  }
  if(nullmod) f <- NULL
  else
    {
      ytype <- attr(Y, "type")
      if( method=="breslow" || method =="efron") {
        if (ytype== 'right')  fitter <- get("coxph.fit")
        else if (ytype=='counting') fitter <- get("agreg.fit")
        else stop(paste("Cox model doesn't support \"", ytype,
                        "\" survival data", sep=''))
      }
      else if (method=='exact') fitter <- get("agexact.fit")
      else stop(paste ("Unknown method", method))
      if (missing(init)) init <- NULL

      ## S-Plus 6 has control parameter.  S-Plus 5 has toler.chol.
      ## Previous to 5 has neither.  R has names(fitter)=NULL
      nf <- names(fitter)
      if(any(nf=='toler.chol'))
        f <- fitter(X, Y, strata=Strata, offset=offset, iter.max=iter.max,
                    eps=eps, weights=weights, init=init,
                    method=method, rownames=rnam, toler.chol=tol) else
      if(.R. || any(nf=='control'))
        f <- fitter(X, Y, strata=Strata, offset=offset,
                    weights=weights, init=init,
                    method=method, rownames=rnam,
                    control=coxph.control(eps=eps, toler.chol=tol,
                      toler.inf=1, iter.max=iter.max)) else
      f <- fitter(X, Y, strata=Strata, offset=offset, iter.max=iter.max,
                  eps=eps, weights=weights, init=init,
                  method=method, rownames=rnam)
    }
  if (is.character(f)) {
    cat("Failure in cph:\n",f,"\n")
    if(.SV4.)return(structure(list(fail=TRUE,fitFunction='cph'),
                              class='Design')) else 
    return(structure(list(fail=TRUE),class="cph"))
  }   else {
    if(length(f$coef) && any(is.na(f$coef))) {  ## length 1may02
      vars <- names(f$coef)[is.na(f$coef)]
      msg <- paste("X matrix deemed to be singular; variable",
                   paste(vars, collapse=" "))
      if(singular.ok) warning(msg)
	  else {
		cat(msg,"\n")
        if(.SV4.)return(structure(list(fail=TRUE,fitFunction='cph'),
                                  class='Design')) else 
        return(structure(list(fail=TRUE),class="cph")) ##13Nov00
	  }
    }
  }
#  attr(Terms, 'Design') <- atr
  f$terms <- Terms

  if(robust) {
    f$naive.var <- f$var
    ## Terry gets a little tricky here, calling resid before adding
    ## na.action method to avoid re-inserting NAs.  Also makes sure
    ## X and Y are there
    if(missing(cluster)) cluster <- FALSE
    fit2 <- c(f, list(x=X, y=Y, method=method))
    if(length(stra)) fit2$strata <- Strata
    temp <- residuals.coxph(fit2, type='dfbeta', collapse=cluster)
    f$var <- t(temp) %*% temp
  }
  
  if(length(weights) && any(weights!=1)) f$weights <- weights

  nvar <-  length(f$coefficients)

  temp <- factor(Y[,ny], levels=0:1, labels=c("No Event","Event"))
  ## was category 30apr02
  n.table <- if(.R.) {
    if(!length(Strata)) table(temp,dnn='Status') else
    table(Strata, temp, dnn=c('Stratum','Status'))
  } else {
    if(!length(Strata)) table(temp) else table(Strata, temp)
  }
  f$n <- n.table
  nnn <- nrow(Y)
  nevent <- sum(Y[,ny])
  if(xpres)	{
    logtest <- -2 * (f$loglik[1] - f$loglik[2])
    R2.max <- 1 - exp(2*f$loglik[1]/nnn)
    R2 <- (1 - exp(-logtest/nnn))/R2.max
    P <- 1-pchisq(logtest,nvar)
    stats <- c(nnn, nevent, logtest, nvar, P, f$score, 
               1-pchisq(f$score,nvar), R2)
    names(stats) <- c("Obs", "Events", "Model L.R.", "d.f.", "P", 
                      "Score", "Score P","R2")
  }
  else { stats <- c(nnn, nevent); names(stats) <- c("Obs","Events") }
  f$method <- NULL
  if(xpres) dimnames(f$var) <- list(atr$colnames, atr$colnames)  #2Apr95
  f <- c(f, list(call=call, Design=atr, 
                 assign=DesignAssign(atr, 0, atrx$terms), # was =assign 1may02
                 na.action=nact,  #terms=Terms
                 fail = FALSE, non.slopes = 0, stats = stats, method=method,
                 maxtime = maxtime, time.inc = time.inc,
                 units = time.units, fitFunction=c('cph','coxph')))

  if(xpres) {
    f$center <- sum(f$means*f$coef)
    f$scale.pred <- c("log Relative Hazard","Hazard Ratio")
    attr(f$linear.predictors,"strata") <- Strata
    names(f$linear.predictors) <- rnam	#23Feb93
    if(se.fit) {
      XX <- X - rep(f$means,rep.int(nnn,nvar))   # see scale() function
      ##  XX <- sweep(X, 2, f$means)	# center   (slower)
      se.fit <- drop(((XX %*% f$var) * XX) %*% rep(1,ncol(XX)))^.5
      if(!.R.) storage.mode(se.fit) <- "single"
      names(se.fit) <- rnam
      f$se.fit <- se.fit  	
    }
  }
  if(model) f$model <- m
  if(nstrata > 0) {
    attr(X, "strata") <- attr(Y, "strata") <- Strata
    f$strata <- levels(Strata)
  }
  if(x) f$x <- X
  if(y) f$y <- Y

  if(is.character(surv) || surv) {
    if(!length(Strata)) Strata <- rep(1, nnn)
    else Strata <- oldUnclass(Strata)
    nstr <- max(Strata, na.rm=TRUE)
    srv <- NULL
    tim <- NULL
    s.e. <- NULL
    timepts <- seq(0, maxtime, by=time.inc)
    s.sum <- array(if(.R.)double(1) else single(1),c(length(timepts),nstr,3),
                   list(format(timepts),paste("Stratum",1:nstr),
                        c("Survival","n.risk","std.err")))
    
    if(xpres) {
      g <- f; if(!x) g$x <- X; if(!y) g$y <- Y
      fname <- 'survfit.cph'
    } else {
      g <- f; if(!y) g$y <- Y
      fname <- 'survfit.cph.null'
      }
    g <- list(g)
    if(!missing(type))      g$type      <- type
    if(!missing(vartype))   g$vartype   <- vartype
    if(!missing(conf.type)) g$conf.type <- conf.type
    g <- do.call(fname, g)

    if(nstr==1) stemp <- rep(1, length(g$time))
    else stemp <- rep(1:nstr,g$strata)
    i <- 0
    for(k in 1:nstr) {
      j <- stemp==k; i <- i+1
      yy <- Y[Strata==i,ny-1]
      maxt <- max(yy)
      ##n.risk from surv.fit does not have usual meaning if not Kaplan-Meier
      
      tt <- c(0,g$time[j])
      su <- c(1,g$surv[j])
      se <- c(NA,-g$std.err[j]/logb(g$surv[j]))
      if(!.R.) {
        storage.mode(tt) <- 'single'
        storage.mode(su) <- 'single'
        storage.mode(se) <- 'single'
      }
      if(maxt>tt[length(tt)]) {
        tt <- c(tt, maxt); su <- c(su, su[length(su)]); se <- c(se, NA)
      }
      kk <- 0
      for(tp in timepts) {
        kk <- kk + 1
        ##	  t.choice <- max((1:length(tt))[max(tt[tt<=tp])==tt])
        t.choice <- max((1:length(tt))[tt<=tp+1e-6])
        if(tp > max(tt)+1e-6 & su[length(su)]>0) {
          Su <- NA
          Se <- NA			}	else {
            Su <- su[t.choice]
            Se <- se[t.choice]		 }
        n.risk <- sum(yy>=tp)
        s.sum[kk,i,1:3] <- c(Su, n.risk, Se)		}
      if(!is.character(surv)) {
        if(nstr==1)	{
          tim <- tt
          srv <- su
          s.e. <- se		} else {
            tim <- c(tim, list(tt))
            srv <- c(srv, list(su))
            s.e. <- c(s.e., list(se))		}	}
    }
    if(is.character(surv)) f$surv.summary <- s.sum
    else {
      attr(srv, "type") <- if(missing(type)) method else type #7Jun99
      if(nstr>1) { names(srv) <- names(tim) <- names(s.e.) <- f$strata }
      f <- c(f, list(time=tim, surv=srv,
                     std.err=s.e., surv.summary=s.sum))		
    }
  }
  oldClass(f) <- if(.SV4.) 'Design' else c("cph", "Design", "coxph") ##13Nov00
  f
}

## Define a version of coxph.fit that works in S-Plus post version
## 4.5 as well as in earlier versions, as toler.chol argument was
## added in survival5 for S-Plus 2000 and Unix S-Plus 5.0
## coxphFit also handles the case where toler.chol and 4 other
## arguments were forgotten.  In S-Plus 6 control= is used.

coxphFit <- if(.R. || .SV4.)
  function(..., strata=NULL, rownames=NULL, offset=NULL,
           init=NULL, toler.chol=1e-9, eps=.0001, iter.max=10) {

    res <- coxph.fit(..., strata=strata, rownames=rownames,
                     offset=offset, init=init,
                     control=coxph.control(toler.chol=toler.chol, toler.inf=1,
                       eps=eps, iter.max=iter.max))
    if(is.character(res)) return(list(fail=TRUE))
    if(iter.max > 1 && res$iter >= iter.max) return(list(fail=TRUE))
    res$fail <- FALSE
    res
} else  function(..., strata=NULL, rownames=NULL, offset=NULL,
                     init=NULL, toler.chol=1e-9, eps=.0001, iter.max=10) {

    nf <- names(coxph.fit)  ##13Nov00
    res <-
      if(any(nf=='control'))
        coxph.fit(..., strata=strata, rownames=rownames,
                  offset=offset, init=init,
                  control=coxph.control(toler.chol=toler.chol, toler.inf=1,
                    eps=eps, iter.max=iter.max))
      else if(all(c('toler.chol','eps','iter.max') %in% nf))
        coxph.fit(..., strata=strata, rownames=rownames,
                  offset=offset, init=init, toler.chol=toler.chol,
                  eps=eps, iter.max=iter.max)
      else if(all(c('iter.max','eps') %in% nf))
        coxph.fit(..., strata=strata, rownames=rownames,
                  offset=offset, init=init, eps=eps, iter.max=iter.max)
      else coxph.fit(..., strata=strata, rownames=rownames,
                     offset=offset, init=init)
    if(is.character(res)) return(list(fail=TRUE))
    if(length(res$iter) && iter.max > 1 && res$iter >= iter.max)
      return(list(fail=TRUE))
    res$fail <- FALSE
    res
  }

Survival.cph <- function(fit, ...) {
if(!length(fit$time) || !length(fit$surv))
  stop("did not specify surv=T with cph")
f <- function(times, lp=0, stratum=1, type=c("step","polygon"),
              time, surv) {
  type <- match.arg(type)
  if(length(stratum)>1) stop("does not handle vector stratum")
  if(length(times)==0) {
    if(length(lp)>1) stop("lp must be of length 1 if times=NULL")
    return(surv[[stratum]]^exp(lp))
  }
  s <- matrix(NA, nrow=length(lp), ncol=length(times),
              dimnames=list(names(lp), format(times)))
  if(is.list(time)) {time <- time[[stratum]]; surv <- surv[[stratum]]}
  if(type=="polygon") {
    if(length(lp)>1 && length(times)>1)
      stop('may not have length(lp)>1 & length(times>1) when type="polygon"')
    su <- approx(time, surv, times)$y
    return(su ^ exp(lp))
  }
  for(i in 1:length(times)) {
    tm <- max((1:length(time))[time <= times[i]+1e-6])
    su <- surv[tm]
    if(times[i] > max(time)+1e-6) su <- NA
    s[,i] <- su^exp(lp)
  }
  drop(s)
}
formals(f) <- list(times=NULL, lp=0, stratum=1,
                   type=c("step","polygon"),
                   time=fit$time, surv=fit$surv)
f
}

Quantile.cph <- function(fit, ...) {
if(!length(fit$time) || !length(fit$surv))
  stop("did not specify surv=T with cph")
f <- function(q=.5, lp=0, stratum=1, type=c("step","polygon"), time, surv) {
  type <- match.arg(type)
  if(length(stratum)>1) stop("does not handle vector stratum")
  if(is.list(time)) {time <- time[[stratum]]; surv <- surv[[stratum]]}
  Q <- matrix(NA, nrow=length(lp), ncol=length(q),
              dimnames=list(names(lp), format(q)))
  for(j in 1:length(lp)) {
    s <- surv^exp(lp[j])
    if(type=="polygon") Q[j,] <- approx(s, time, q)$y
    else for(i in 1:length(q))
      if(any(s <= q[i])) Q[j,i] <- min(time[s<=q[i]])  #is NA if none
    ## if(any()) 20may02
  }
  drop(Q)
}
formals(f) <- list(q=.5, lp=0, stratum=1,
                   type=c('step','polygon'),
                   time=fit$time, surv=fit$surv)
f
}


Mean.cph <- function(fit, method=c("exact","approximate"),
type=c("step","polygon"), n=75, tmax, ...) {
method <- match.arg(method)
type   <- match.arg(type)

if(!length(fit$time) || !length(fit$surv))
  stop("did not specify surv=T with cph")

if(method=="exact") {
  f <- function(lp=0, stratum=1, type=c("step","polygon"),
                tmax=NULL, time, surv) {
    type <- match.arg(type)
    if(length(stratum)>1) stop("does not handle vector stratum")
    if(is.list(time)) {time <- time[[stratum]]; surv <- surv[[stratum]]}
    Q <- lp
    if(!length(tmax)) {
      if(min(surv)>1e-3)
        warning(paste("Computing mean when survival curve only defined down to",
                      format(min(surv)),"\n Mean is only a lower limit"))
      k <- rep(TRUE,length(time))
    } else {
      if(tmax>max(time)) stop(paste("tmax=",format(tmax),
                                    "> max follow-up time=",format(max(time))))
      k <- (1:length(time))[time<=tmax]
    }
    for(j in 1:length(lp)) {
      s <- surv^exp(lp[j])
      Q[j] <- if(type=="step") sum(c(diff(time[k]),0) * s[k]) else 
      trap.rule(time[k], s[k])
    }
    Q
  }
  formals(f) <- list(lp=0, stratum=1,
                     type=if(!missing(type))type else c("step","polygon"),
                     tmax=tmax,
                     time=fit$time, surv=fit$surv)
##    if(!missing(tmax)) f$tmax <- tmax ??
} else {
  lp     <- fit$linear.predictors
  lp.seq <- if(length(lp)) lp.seq <- seq(min(lp), max(lp), length=n) else 0
  
  time   <- fit$time
  surv   <- fit$surv
  nstrat <- if(is.list(time)) length(time) else 1
  areas  <- list()

  for(is in 1:nstrat) {
    tim <- if(nstrat==1) time else time[[is]]
    srv <- if(nstrat==1) surv else surv[[is]]
    if(!length(tmax)) {
      if(min(srv)>1e-3)
        warning(paste("Computing mean when survival curve only defined down to",
                      format(min(srv)),"\n Mean is only a lower limit"))
      k <- rep(TRUE,length(tim))
    } else {
      if(tmax>max(tim)) stop(paste("tmax=",format(tmax),
                                   "> max follow-up time=",format(max(tim))))
      k <- (1:length(tim))[tim<=tmax]
    }
    ymean <- lp.seq
    for(j in 1:length(lp.seq)) {
      s <- srv^exp(lp.seq[j])
      ymean[j] <- if(type=="step") sum(c(diff(tim[k]),0) * s[k]) else 
      trap.rule(tim[k], s[k])
    }
    if(!.R.) storage.mode(ymean) <- "single"
    areas[[is]] <- ymean
  }
  if(nstrat>1) names(areas) <- names(time)

  f <- function(lp=0, stratum=1, lp.seq, areas) {

    if(length(stratum)>1) stop("does not handle vector stratum")
    area <- areas[[stratum]]
    if(length(lp.seq)==1 && all(lp==lp.seq)) ymean <- rep(area,length(lp)) else
    ymean <- approx(lp.seq, area, xout=lp)$y
    if(any(is.na(ymean)))
      warning("means requested for linear predictor values outside range of linear\npredictor values in original fit")
    names(ymean) <- names(lp)
    ymean
  }
if(!.R.) storage.mode(lp.seq) <- "single"
formals(f) <- list(lp=0, stratum=1, lp.seq=lp.seq, areas=areas)
}
f
}

## cox.zph demands that the fit object inherit 'coxph'  14Nov00
## The following slightly changed cox.zph also explicitly invokes
## residuals.cph
cox.zph <- function(fit, transform = "km", global = TRUE)
{
	call <- match.call()
    clas <- c(oldClass(fit), fit$fitFunction)  ##FEH
	if(!any(c('coxph','cph') %in% clas))       ##FEH
		stop("Argument must be the result of coxph or cph")
      if('coxph.null' %in% clas)               ##FEH
		stop("The are no score residuals for a Null model")
	sresid <- if(any(clas=='cph')) residuals.cph(fit, "schoenfeld") else
     residuals.coxph(fit, 'schoenfeld')        ##FEH
	varnames <- names(fit$coef)
	nvar <- length(varnames)
	ndead <- length(sresid)/nvar
	if(nvar == 1)
		times <- as.numeric(names(sresid))
	else times <- as.numeric(dimnames(sresid)[[1]])
	if(missing(transform) && attr(fit$y, "type") != "right")
		transform <- "identity"
	if(is.character(transform)) {
		tname <- transform
		ttimes <- switch(transform,
			identity = times,
			rank = rank(times),
			log = log(times),
			km = {
				temp <- survfit.km(factor(rep(1, nrow(fit$
					y))), fit$y, se.fit = FALSE)
				# A nuisance to do left cont KM
				t1 <- temp$surv[temp$n.event > 0]
				t2 <- temp$n.event[temp$n.event > 0]
				km <- rep(c(1, t1), c(t2, 0))
				if(!length(attr(sresid, "strata")))
					1 - km
				else (1 - km[sort.list(sort.list(times))])
			}
			,
			stop("Unrecognized transform"))
	}
	else {
		tname <- deparse(substitute(transform))
		ttimes <- transform(times)
	}
	xx <- ttimes - mean(ttimes)
	r2 <- sresid %*% fit$var * ndead
	test <- xx %*% r2
	# time weighted col sums
	corel <- c(cor(xx, r2))
	z <- c(test^2/(diag(fit$var) * ndead * sum(xx^2)))
	Z.ph <- cbind(corel, z, 1 - pchisq(z, 1))
	if(global && nvar > 1) {
		test <- c(xx %*% sresid)
		z <- (c(test %*% fit$var %*% test) * ndead)/sum(xx^2)
		Z.ph <- rbind(Z.ph, c(NA, z, 1 - pchisq(z, ncol(sresid))))
		dimnames(Z.ph) <- list(c(varnames, "GLOBAL"), c("rho", "chisq",
			"p"))
	}
	else dimnames(Z.ph) <- list(varnames, c("rho", "chisq", "p"))
	dimnames(r2) <- list(times, names(fit$coef))
	temp <- list(table = Z.ph, x = ttimes, y = r2 + outer(rep(1, ndead),
		fit$coef), var = fit$var, call = call, transform = tname)
	oldClass(temp) <- "cox.zph"
	temp
}

predict.cph <- 
  function(object, newdata,
           type=c("lp","x","data.frame","terms","adjto","adjto.data.frame",
             "model.frame"),
           se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
           incl.non.slopes, non.slopes, kint=1,
           na.action=na.keep, expand.na=TRUE, center.terms=TRUE, ...)
  predictDesign(object, newdata, type, se.fit, conf.int, conf.type,
                incl.non.slopes, non.slopes, kint,
                na.action, expand.na, center.terms, ...)

cr.setup <- function(y) {

yname <- as.character(substitute(y))

if(!is.category(y)) y <- factor(y, exclude=NA)  # was category 20may02
y       <- oldUnclass(y)   # in case is.factor
ylevels <- levels(y)
kint    <- length(ylevels) - 1
y       <- as.integer(y-1)

reps   <- ifelse(is.na(y), 1, ifelse(y < kint-1, y+1, kint))
subs   <- rep(1:length(y), reps)

cuts   <- vector('list',kint+2)
cuts[[1]] <- NA
for(j in 0:kint) cuts[[j+2]] <- if(j < kint-1) 0:j else 0:(kint-1)

cuts   <- unlist(cuts[ifelse(is.na(y),1,y+2)])
y      <- rep(y, reps)
Y      <- 1*(y==cuts)
labels <- c('all', paste(yname,'>=',ylevels[2:kint],sep=''))
cohort <- factor(cuts, levels=0:(kint-1), labels=labels)

list(y=Y, cohort=cohort, subs=subs, reps=reps)
}

datadist <- function(..., data, q.display, q.effect=c(.25,.75),
                     adjto.cat=c('mode','first'), n.unique=10) {
adjto.cat <- match.arg(adjto.cat)
X <- list(...)

argnames <- as.character(sys.call())[-1]

if(inherits(x <- X[[1]],"datadist")) {
  Limits <- x$limits
  Values <- x$values
  X[[1]] <- NULL
  argnames <- argnames[-1]
} else {
  Limits <- list()
  Values <- list()
}

if(is.data.frame(X[[1]])) {
  if(length(X) > 1) stop('when the first argument is a data frame, no other variables may be specified')  ## 2Apr99
  X <- X[[1]]
}
  
else if(length(Terms <- X[[1]]$terms) && length(D <- attr(Terms,"Design"))){
   n <- D$name[D$assume!="interaction"]
   X <- list()
   if(missing(data)) for(nm in n) X[[nm]] <- eval(as.name(nm),local=FALSE) #27May99
   else if(length(names(data)))
   {
     j <- match(n, names(data), 0)
     if(any(j==0)) stop(paste("variable(s)",
              paste(n[j==0],collapse=" "),
              "in model not found on data=, \nwhich has variables",
              paste(names(data),collapse=" ")))
     for(nm in n) X[[nm]] <- data[[nm]]
   }
   else for(nm in n) X[[nm]] <- get(nm, data)  
 }
else {

  if(length(X) & !length(names(X))) names(X) <- argnames[1:length(X)]

  if(!missing(data))
    {	# This duplicative code is for efficiency for large data frames
      if(length(X))
        {
          if(is.numeric(data)) X <- c(X,database.object(data))
          else X <- c(X, data)
        }
      else
        {
          if(is.numeric(data)) X <- database.object(data)
          else X <- data
        }
    }
}
nam <- names(X)
p <- length(nam)
if(p==0) stop("you must specify individual variables or a data frame")

maxl <- 0
for(i in 1:p) {
  values <- NULL
  x <- X[[i]]
  if(is.character(x)) x <- as.factor(x)
  lx <- length(x)
  lev <- levels(x)
  ll <- length(lev)
  limits <- rep(NA, 5)
  if(is.matrix(x) | (i>1 && lx!=maxl)) {
    warning(paste(nam[i],"is a matrix or has incorrect length; ignored"))
  }
  else {
    if(ll && (ll<length(x))) values <- lev   # if # levels=length(x) is ID variable
    ## First look for ordered variable with numeric levels (scored() var)
#    if(is.ordered(x) && !any(is.na(levx <- sort(as.numeric(lev))))) {
## 6may03 next 2 lines
    if(is.ordered(x) && all.is.numeric(lev)) {
      levx <- sort(as.numeric(lev))
      limits <- c(levx[1],levx[(ll+1)/2],levx[ll],levx[1],levx[ll],
                  levx[1],levx[ll])
      values <- levx
    }
      
    else if(ll) {
      adjto <- if(adjto.cat=='first') lev[1] else {
        tab <- table(x)
        (names(tab)[tab==max(tab)])[1]
      }
      limits <- factor(c(NA,adjto,NA,lev[1],lev[ll],lev[1],lev[ll]),
                       levels=lev) 
      ## was c("",lev[1],"",lev[1],lev[ll],lev[1],lev[ll]) 6Feb95
      ## non-ordered categorical
    }
    else {	 	# regular numeric variable
      clx <- oldClass(x)
      y <- x[!is.na(x)]
      n <- length(y)
      if(n<2)stop(paste("fewer than 2 non-missing observations for",nam[i]))
      values <- sort(unique(y))
      names(values) <- NULL
      nunique <- length(values)
      if(nunique < 2) {
        warning(paste(nam[i],"is constant"))
        limits <- rep(y[1], 7)
      }
      else	{
        r <- range(values)
        limits[6:7] <- r
        if(nunique<4) q <- r else {
          if(missing(q.display))	{
            q.display <- 10/max(n,200)
            q.display <- c(q.display,1-q.display)	}
          q <- quantile(oldUnclass(y),q.display)	}  #chron obj. not work here
        limits[4] <- q[1]; limits[5] <- q[2]
        ## check for very poorly distributed categorical numeric variable
        if(limits[4]==limits[5]) limits[4:5] <- r

        ## Use low category if binary var, middle if 3-level, median otherwise
        if(nunique < 3) limits[2] <- values[1] else
        if(nunique==3) limits[2] <- values[2] else
        limits[2] <- median(oldUnclass(y))

        if(nunique < 4) q <- r else
        q <- quantile(oldUnclass(y), q.effect)
        limits[1] <- q[1]; limits[3] <- q[2]
        if(limits[1]==limits[3]) limits[c(1,3)] <- r
        if(nunique > n.unique) values <- NULL
        oldClass(limits) <- clx
      }
    }
    Limits[[nam[i]]] <- limits
    if(length(values)) Values[[nam[i]]] <- values
    maxl <- max(maxl, lx)
  }
}

Limits <- structure(Limits, class="data.frame", 
                    row.names=c("Low:effect","Adjust to",
                      "High:effect","Low:prediction",
                      "High:prediction","Low","High"))
##data.frame(Limits) gives error with chron objects

d <- list(limits=Limits, values=Values)
oldClass(d) <- "datadist"
d
}

print.datadist <- function(x, ...) {
  lim <- x$limits
  for(n in names(lim))	{
    z <- lim[[n]]
    if(inherits(z,"dates") | inherits(z,"times"))
      lim[[n]] <- factor(format(z))
  }
  if(length(lim)) print(lim)
  ##print.data.frame doesn't print chron objects correctly
  if(length(V <- x$values)) {
    cat("\nValues:\n\n")
    wid <- .Options$width
    for(n in names(V)) {
      v <- V[[n]]
      if(length(v)==0) next  # for gendata
      if(is.character(v) && length(v)>80) 
        v <- c(v[1:20],paste("+",length(v),"others"))
      w <- if(is.character(v)) v else format(v)
      nc <- nchar(paste(w,collapse=" "))
      if(nc+nchar(n)+4>wid) {cat(n,":\n"); print(v, quote=FALSE)}
      else cat(n,":",w,"\n")
	}
  }
  invisible()
}

# Fast backward elimination using a slow but numerically stable version
# of the Lawless-Singhal method (Biometrics 1978), used in the SAS
# PHGLM and LOGIST procedures
# Uses function solvet, a slightly edited version of solve that passes
# the tol argument to qr.
# Modified 12Oct92 - if scale parameter present, ignore last row and col of cov
# Modified 22Sep93 - new storage format for design attributes
# Modified 1Mar94 - add k.aic
# Modified 4Mar96 - use S commands instead of avia if not under UNIX
#
# F. Harrell 18Jan91

fastbw <- function(fit, rule="aic", type="residual", sls=.05, aics=0, 
			eps=1e-9, k.aic=2) {

ns <- num.intercepts(fit)
L <- if(ns==0) NULL else 1:ns

pt <- length(fit$coef)
p <- pt-ns
atr <- fit$Design
if(!length(atr)) atr <- getOldDesign(fit)

assume <- atr$assume.code
if(!length(assume))stop("fit does not have design information")
assign <- fit$assign
nama <- names(assign)[1]
asso <- 1*(nama=="(Intercept)" | nama=="Intercept")

f <- sum(assume!=8)
strt <- integer(f)
len <- strt
j <- 0
for(i in 1:length(assume)) {
   if(assume[i]!=8)	{
      j <- j+1; aj <- assign[[j+asso]]
      strt[j] <- min(aj)
      len[j]  <- length(aj)
    }
 }
name <- atr$name[assume!=8]
ed <- as.integer(strt+len-1)

rule <- charmatch(rule,c("aic","p"),0)
if(rule==0)stop("rule must be aic or p for Akaike's info criterion or p-value")
type <- charmatch(type,c("residual","individual","total"),0)
if(type==0)stop("type must be residual or individual")
if(type==3) type <- 1

factors.in <- 1:f
parms.in <- 1:pt

##library.dynam(section="local",file="mlmats.o")

##Delete this block of code if using solve() instead of avia
## Allocate work areas for avia
if(under.unix || .R.) {
  s1 <- double(pt)
  s2 <- s1
  s3 <- double(2*pt)    #28Jul95
  s4 <- s1
  vsub <- double(pt*pt)
  pivot <- integer(pt)
}

factors.del <- integer(f)
chisq.del <- single(f)
df.del <- integer(f)
resid.del <- single(f)
df.resid <- integer(f)
beta <- fit$coef
Cov  <- Varcov(fit, regcoef.only=TRUE)  #Ignore scale parameters
cov <- Cov
Coef <- matrix(NA, nrow=f, ncol=pt, dimnames=list(NULL,names(beta)))
d <- 0

##14Sep99
fcl <- oldClass(fit)  # added this and length(fcl) below 8Dec99
## added fitFunction 11Apr02
dor2 <- length(fcl) && (any(fcl=='ols') || (length(fit$fitFunction) &&
                              any(fit$fitFunction=='ols'))) && 
  (length(fit$y) || (length(fit$fitted.values) &&
                     length(fit$residuals)))
if(dor2) {
  ## X <- fit$x
  Y <- if(length(fit$y))fit$y else fit$fitted.values + fit$residuals
  r2 <- single(f)
  sst <- sum((Y-mean(Y))^2)
  sigma2 <- fit$stats['Sigma']^2
  ## Get X'Y using b=(X'X)^-1 X'Y, X'X^-1 = var matrix / sigma2
  xpy <- matrix(solve(Cov, beta)*sigma2, ncol=1)
  ypy <- sum(Y^2)
}

for(i in 1:f) {
	fi <- length(factors.in)
	ln <- len[factors.in]
	st <- as.integer(ns+c(1,1+cumsum(ln[-fi]))[1:fi])
	en <- as.integer(st+ln-1)
	crit.min <- 1e10
	chisq.crit.min <- 1e10
	jmin <- 0
	dfmin <- 0
	k <- 0
	factors.in.loop <- factors.in  #indirect reference prob in S 3.1
	for(j in factors.in.loop) {
	k <- k+1
	if(under.unix && !.R.) {  # can't get this to work in R - CHECK
		z <- if(.R.)
          .Fortran("avia",beta,cov,chisq=double(1),length(beta),
                   st[k]:en[k],
                   ln[k],df=integer(1),eps,vsub,s1,s2,s3,s4,pivot,NAOK=TRUE,
                   PACKAGE="Design") else
          .Fortran("avia",beta,cov,chisq=double(1),length(beta),
                   st[k]:en[k],
                   ln[k],df=integer(1),eps,vsub,s1,s2,s3,s4,pivot,NAOK=TRUE)
		chisq <- z$chisq
		df <- z$df
	} else {
##replace previous 5 statements with following 3 to use slow method
	   q <- st[k]:en[k]
	   chisq <- beta[q] %*% solve(cov[q,q], beta[q])
	   df <- length(q)
     }
		switch(rule, crit <- chisq-k.aic*df, crit <- pchisq(chisq,df))
		if(crit < crit.min)	{
			jmin <- j
			crit.min <- crit
			chisq.crit.min <- chisq
			df.min <- df
          }	
  }
##	kk <- factors.in[factors.in!=jmin]
##	factors.in <- NULL
##	factors.in <- kk		#funny bug in S 3.1 !
	factors.in <- factors.in[factors.in!=jmin]
	parms.in <- parms.in[parms.in<strt[jmin] | parms.in>ed[jmin]]
	if(length(parms.in)==0) q <- 1:pt else q <- (1:pt)[-parms.in]
	if(under.unix && !.R.) {
	z <- if(.R.)
      .Fortran("avia",fit$coef,Cov,chisq=double(1),
               pt,q,as.integer(pt-length(parms.in)),
               df=integer(1),eps,vsub,s1,s2,s3,s4,pivot,NAOK=TRUE,
               PACKAGE="Design") else
      .Fortran("avia",fit$coef,Cov,chisq=double(1),
               pt,q,as.integer(pt-length(parms.in)),
               df=integer(1),eps,vsub,s1,s2,s3,s4,pivot,NAOK=TRUE)
	resid <- z$chisq
	resid.df <- z$df
	} else {
##replace previous 5 statements with following 2 to use slow method
	resid <- fit$coef[q] %*% solve(Cov[q,q], fit$coef[q])
	resid.df <- length(q)
  }
	switch(type,
           switch(rule, del <- resid-k.aic*resid.df <= aics,
                  del <- 1-pchisq(resid,resid.df)>sls),
           switch(rule, del <- crit.min <= aics,
                  del <- 1-crit.min > sls)	)
	if(del) {
		d <- d+1
		factors.del[d] <- jmin
		chisq.del[d] <- chisq.crit.min
		df.del[d] <- df.min
		resid.del[d] <- resid
		df.resid[d] <- resid.df
		if(length(parms.in)) {
		   cov.rm.inv <- solvet(Cov[-parms.in,-parms.in], tol=eps)
		   cov.cross <- Cov[parms.in,-parms.in,drop=FALSE]
		   w <- cov.cross %*% cov.rm.inv
		   beta <- fit$coef[parms.in] - w %*% fit$coef[-parms.in]
		   cov <- Cov[parms.in,parms.in] - w %*% t(cov.cross)
           cof <- rep(0, pt)
           cof[parms.in] <- beta
           Coef[d,] <- cof
           if(dor2) {
             ## yhat <- matxv(X[,parms.in,drop=F], beta)
             ## r2[d] <- 1 - sum((yhat-Y)^2)/sst
             ## sse = Y'(I - H)Y, where H = X*inv(X'X)*X'
             ##     = Y'Y - Y'X*inv(X'X)*X'Y
             ##     = Y'Y - Y'Xb
             sse <- ypy - t(xpy[parms.in,,drop=FALSE])%*%beta
             r2[d] <- 1 - sse/sst
           }
         }
		else {
          beta <- NULL; cov <- NULL
          if(dor2) r2[d] <- 0
        }
      }	else break
  }

if(d>0) {
  fd <- factors.del[1:d]
  if(dor2) {
    r2 <- r2[1:d]
    Coef <- Coef[1:d,,drop=FALSE]
  }
  res <- cbind(chisq.del[1:d],df.del[1:d],
               1-pchisq(chisq.del[1:d],df.del[1:d]),
               resid.del[1:d],df.resid[1:d],
               1-pchisq(resid.del[1:d],df.resid[1:d]),resid.del[1:d]-k.aic*
               df.resid[1:d])
  labs <- c("Chi-Sq","d.f.","P","Residual","d.f.","P","AIC")
  dimnames(res) <- list(name[fd],labs)
  if(length(fd)==f) fk <- NULL else fk <- (1:f)[-fd]
}	else {
  fd <- NULL
  res <- NULL
  fk <- 1:f
}

nf <- name[fk]

pd <- NULL
##if(d>0) for(i in 1:d) pd <- c(pd, (strt[fd[i]]:ed[fd[i]])-ns)
if(d>0) for(i in 1:d) pd <- c(pd, (strt[fd[i]]:ed[fd[i]]))

if(length(fd)==f) fk <- NULL else if(d==0) fk <- 1:f else fk <- (1:f)[-fd]
##if(length(pd)==p) pk <- NULL else if(d==0) pk <- 1:p else pk <- (1:p)[-pd]
if(length(pd)==p) pk <- L else if(d==0) pk <- 1:pt else pk <- (1:pt)[-pd]

if(length(pd)!=p) {
  beta <- as.vector(beta)
  names(beta) <- names(fit$coef)[pk]
  dimnames(cov) <- list(names(beta),names(beta))
}
##if(ns>0)	{	#removed 11Jan94
##   if(length(pk)>ns)pk <- pk[-(1:ns)]-ns else pk <- NULL
##   pd <- pd-ns	}

if(dor2) res <- cbind(res, R2=r2)
r <- list(result=res,names.kept=nf,factors.kept=fk,
          factors.deleted=fd,
          parms.kept=pk,parms.deleted=pd, coefficients=beta, var=cov,
          Coefficients=Coef)
oldClass(r) <- "fastbw"
r
}

		
print.fastbw <- function(x, digits=4, ...) {

res <- x$result
fd <- x$factors.deleted
if(length(fd)) {
cres <- cbind(dimnames(res)[[1]],format(round(res[,1],2)),format(res[,2]),
	format(round(res[,3],4)),format(round(res[,4],2)),
	format(res[,5]),format(round(res[,6],4)),
	format(round(res[,7],2)),if(ncol(res)>7)format(round(res[,8],3)))
dimnames(cres) <- list(rep("",nrow(cres)), c("Deleted", dimnames(res)[[2]]))
cat("\n")
print(cres, quote=FALSE)
if(length(x$coef)) {
    cat("\nApproximate Estimates after Deleting Factors\n\n")
    cof <- coef(x)
    vv <- if(length(cof)>1) diag(x$var) else x$var
    z <- cof/sqrt(vv)
    stats <- cbind(cof, sqrt(vv), z, 1-pchisq(z^2,1))
    dimnames(stats) <- list(names(cof), c("Coef","S.E.","Wald Z","P"))
    print(stats, digits=digits)
  }
}
else cat("\nNo Factors Deleted\n")
cat("\nFactors in Final Model\n\n")
nk <- x$names.kept
if(length(nk))print(nk, quote=FALSE)
else cat("None\n")
}
gendata <- function(fit, ...) UseMethod("gendata")

gendata.default <- function(fit, ...) gendata.Design(obj, ...)

gendata.Design <- function(fit, nobs, viewvals=FALSE,
	editor=.Options$editor, ..., factors) {
 
at <- fit$Design
if(!length(at)) at <- getOldDesign(fit)

nam <- at$name[at$assume!="interaction"]

if(!length(editor) && exists('using.X') && using.X()) editor <- "xedit"

if(!missing(nobs) && !is.logical(nobs)) {
  df <- predictDesign(fit, type="adjto.data.frame")
  df[1:nobs,] <- df
  cat("Edit the list of variables you would like to vary.\nVariables not listed will be set to reference values.\n")
  if(editor=="xedit") cat("To delete an individual variable, type Cntl-k\nTo delete blocks of variables, highlight the block by holding down the left\nmouse button, then type Cntl-w.\n")
  nam.sub <- if(.R.)edit(nam, editor=editor) else ed(nam, editor=editor)
  if(!all(nam.sub %in% nam)) stop("misspelled a variable name")
  df.sub <- as.data.frame(df[,nam.sub])  #df[,] was returning list (?)
  cat("Edit the predictor settings to use.\n")
  if(viewvals && 
    length(val <- Getlim(at, allow.null=TRUE, need.all=FALSE)$values[nam.sub])) {
    cat("A window is being opened to list the valid values of discrete variables.\n")
    sink(tf <- tempfile())
    print.datadist(list(values=val))
    sink()
    if(.R.)file.show(tf) else page(filename=tf)
  }
  if(existsFunction('Edit.data')) {
    stop('use of S-PLUS 4.x GUI not yet implemented for gendata')
    assign('.df.sub.', df.sub, where=1)
    Edit.data(.df.sub., '.df.sub.')
    df.sub <- get('.df.sub.', where=1)
    remove('.df.sub.', where=1)
  }
  else if(existsFunction('data.ed')) {
#    if(!(exists('using.X') && using.X()))
#      stop("must be using X-windows to use interactive data.ed")
    df.sub <- data.ed(df.sub)
  }
  else if(existsFunction('data.entry')) df.sub <- data.entry(df.sub)
  df[nam.sub] <- df.sub
  return(structure(df, names.subset=nam.sub))
}

factors <- if(missing(factors)) list(...) else factors
fnam <- names(factors)
nf <- length(factors)
if(nf==0) return(predictDesign(fit, type="adjto.data.frame"))
which <- charmatch(fnam, nam, 0)
if(any(which==0)) stop(paste("factor(s) not in design:",
	paste(names(factors)[which==0],collapse=" ")))
settings <- if(nf<length(nam)) predictDesign(fit, type="adjto.data.frame") else
	list()
settings <- oldUnclass(settings)
if(nf>0) for(i in 1:nf) settings[[fnam[i]]] <- factors[[i]]
if(nf==0) return(as.data.frame(settings))
expand.grid(settings)

}
glmD <- if(.R.)
  function(formula, family = gaussian, data = list(), weights = NULL,
    subset = NULL, na.action = na.fail, start = NULL, offset = NULL,
    control = glm.control(...), model = TRUE, method = "glm.fit",
    x = FALSE, y = TRUE, contrasts = NULL, ...)
{
    call <- match.call()
    if (is.character(family))
        family <- get(family)
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("`family' not recognized")
    }
    mt <- terms(formula, data = data)
    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    mf$family <- mf$start <- mf$control <- mf$maxit <- NULL
    mf$model <- mf$method <- mf$x <- mf$y <- mf$contrasts <- NULL
    mf$... <- NULL
    mf$drop.unused.levels <- TRUE   # FEH 31jul02
    mf[[1]] <- as.name("model.frame")

    dul <- .Options$drop.unused.levels  # FEH 31jul02
    if(!length(dul) || dul) {
      on.exit(options(drop.unused.levels=dul))
      options(drop.unused.levels=FALSE)
    }

    mf <- Design(eval(mf, parent.frame()))   # FEH 13Apr01
    desatr <- attr(mf,'Design')
    attr(mf,'Design') <- NULL
    
    switch(method, model.frame = return(mf), glm.fit = 1, glm.fit.null = 1,
        stop(paste("invalid `method':", method)))
    xvars <- as.character(attr(mt, "variables"))[-1]
    if ((yvar <- attr(mt, "response")) > 0)
        xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(mf[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts)
    Y <- model.response(mf, "numeric")
    weights <- model.weights(mf)
    offset <- model.offset(mf)
    if (!is.null(weights) && any(weights < 0))
        stop("Negative wts not allowed")
    if (!is.null(offset) && length(offset) != NROW(Y))
        stop(paste("Number of offsets is", length(offset), ", should equal",
            NROW(Y), "(number of observations)"))
    fit <- (if (is.empty.model(mt))
        glm.fit.null
    else glm.fit)(x = X, y = Y, weights = weights, start = start,
        offset = offset, family = family, control = control,
        intercept = attr(mt, "intercept") > 0)
    if (any(offset) && attr(mt, "intercept") > 0) {
        fit$null.deviance <- if (is.empty.model(mt))
            fit$deviance
        else glm.fit(x = X[, "(Intercept)", drop = FALSE], y = Y,
            weights = weights, start = start, offset = offset,
            family = family, control = control, intercept = TRUE)$deviance
    }
    if (model)
        fit$model <- mf
    if (x)
        fit$x <- X
    if (!y)
        fit$y <- NULL
    fit <- c(fit, list(call = call, formula = formula, terms = mt,
        data = data, offset = offset, control = control, method = method,
        contrasts = attr(X, "contrasts"), xlevels = xlev,
                       Design=desatr,
                       assign=DesignAssign(desatr,1,mt)))
    ##FEH 13Apr01 24nov02 above
    class(fit) <- c("Design", 'glmD',
                    if (is.empty.model(mt)) "glm.null", "glm",
        "lm") # FEH 13Apr01  glmD 26nov02
    fit
} else function(formula = formula(data), family = gaussian,
                data = sys.parent(),
        weights, subset, na.action, start = eta, control = glm.control(...),
        method = "glm.fit", model = FALSE, x = FALSE, y = TRUE, contrasts = NULL, ...)
{
        call <- match.call()
        m <- match.call(expand = FALSE)
        m$family <- m$method <- m$model <- m$x <- m$y <- m$control <- m$
                contrasts <- m$... <- NULL
        m$drop.unused.levels <- TRUE
        m[[1]] <- as.name("model.frame")

        m <- Design(eval(m, sys.parent()))  # FEH 13Apr01
        desatr <- attr(m,'Design')
        attr(m,'Design') <- NULL
        
        Terms <- attr(m, "terms")
        if(method == "model.frame")
                return(m)
        xvars <- as.character(attr(Terms, "variables"))
        if(length(xvars) > 0) {
                xlevels <- lapply(m[xvars], levels)
                xlevels <- xlevels[!sapply(xlevels, is.null)]
                if(length(xlevels) == 0)
                        xlevels <- NULL
        }
        else xlevels <- NULL
        a <- attributes(m)
        Y <- model.extract(m, response)
        X <- model.matrix(Terms, m, contrasts)
        w <- model.extract(m, weights)
        if(!length(w))
                w <- rep(1, nrow(m))
        else if(any(w < 0))
                stop("negative weights not allowed")
        start <- model.extract(m, start)
        offset <- model.extract(m, offset)
        family <- as.family(family)
        if(missing(method))
                method <- attr(family, "method")
        if(!is.null(method)) {
                if(!existsFunction(method))
                        stop(paste("unimplemented method:", method))
        }
        else method <- "glm.fit"
        glm.fitter <- get(method)
        fit <- glm.fitter(x = X, y = Y, w = w, start = start, offset = offset,
                family = family, maxit = control$maxit, epsilon = control$
                epsilon, trace = control$trace, null.dev = TRUE, ...)
        #
        # If an offset and intercept is present, iterations are needed to
        # compute the Null deviance; these are done here, unless the model
        # is NULL, in which case the computations have been done already
        #
        if(any(offset) && attr(Terms, "intercept")) {
                null.deviance <- if(length(Terms)) glm.fitter(X[, "(Intercept)",                                drop = FALSE], Y, w, offset = offset, family =
                                family, maxit = control$maxit, epsilon =
                                control$epsilon, null.dev = NULL)$deviance
                         else fit$deviance
                fit$null.deviance <- null.deviance
        }
        oldClass(fit) <- if(.SV4.) 'Design' else
            c("Design","glmD","glm","lm")  # FEH 13Apr01 16aug02
        ## glmD 2dec02 8p
        if(!is.null(xlevels))
                attr(fit, "xlevels") <- xlevels
        fit$terms <- Terms
        fit$formula <- as.vector(attr(Terms, "formula"))
        fit$call <- call
        fit$Design <- desatr   # FEH 13Apr01
        fit$assign <- DesignAssign(desatr,1,Terms) ## 24nov02
        if(model)
                fit$model <- m
        if(x)
                fit$x <- X
        if(!y)
                fit$y <- NULL
        fit$control <- control
        if(!is.null(attr(m, "na.action")))
                fit$na.action <- attr(m, "na.action")
        fit$fitFunction <- c('glmD','glm','lm')  ## glmD 26nov02
        fit
}

## 26nov02
print.glmD <- function(x, digits=4, ...) {
  cat('General Linear Model\n\n')
  dput(x$call); cat('\n\n')
  cof <- coef(x)
  lr <- x$null.deviance - x$deviance
  names(cof) <- ifelse(names(cof)=='(Intercept)','Intercept',names(cof))
  dof <- x$rank - (names(cof)[1]=='Intercept')
  pval <- 1 - pchisq(lr, dof)
  print(c('Model L.R.'=format(lr,digits=2), 'd.f.'=format(dof),
          'P'=format(pval,digits=4)), quote=FALSE)
  cat('\n')
  se <- sqrt(diag(Varcov(x)))
  z <- cof/se
  p <- 1 - pchisq(z^2, 1)
  w <- cbind(format(cof, digits=digits),
             format(se,  digits=digits),
             format(z,   digits=2),
             format(p,   digits=4))
  dimnames(w) <- list(names(cof), c('Coef','S.E.','Wald Z','P'))
  print(w, quote=FALSE)
  invisible()
}

## 26nov02
summary.glmD <- function(...) summary.Design(...)
## 2dec02
Varcov.glmD <- function(object, regcoef.only=FALSE, ...)
  Varcov.glm(object, regcoef.only, ...)
## 6dec02
predict.glmD <- 
  function(object, newdata,
           type=c("lp","x","data.frame","terms","adjto","adjto.data.frame",
             "model.frame"),
           se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
           incl.non.slopes, non.slopes, kint=1,
           na.action=na.keep, expand.na=TRUE, center.terms=TRUE, ...)
  predictDesign(object, newdata, type, se.fit, conf.int, conf.type,
                incl.non.slopes, non.slopes, kint,
                na.action, expand.na, center.terms, ...)


latex.glmD <- function(...) latexDesign(...)
# B= pr= opmeth= added FEH 29mar03
glsD <-
  function (model, data = sys.frame(sys.parent()), correlation = NULL, 
    weights = NULL, subset, method = c("REML", "ML"), na.action = na.fail, 
    control = list(), verbose = FALSE, B=0, dupCluster=FALSE,
            pr=FALSE, opmeth=c('optimize','optim'))
            
{
    if(.R.) require('nlme') else stop('not implemented for S-Plus')
    # In R 1.7 nlme namespace does not export glsEstimate
    if(.R. && !existsFunction('glsEstimate'))
      glsEstimate <- getFromNamespace('glsEstimate','nlme')
    
    Call <- match.call()
    opmeth <- match.arg(opmeth)
    controlvals <- glsControl()
    if (!missing(control)) {
        if (!is.null(control$nlmStepMax) && control$nlmStepMax < 0) {
            warning("Negative control$nlmStepMax - using default value")
            control$nlmStepMax <- NULL
        }
        controlvals[names(control)] <- control
    }
    if (!inherits(model, "formula") || length(model) != 3) {
        stop("\nModel must be a formula of the form \"resp ~ pred\"")
    }
    method <- match.arg(method)
    REML <- method == "REML"
    if (!is.null(correlation)) {
        groups <- getGroupsFormula(correlation)
    }
    else groups <- NULL
    glsSt <- glsStruct(corStruct = correlation, varStruct = varFunc(weights))
    mfArgs <- list(formula = asOneFormula(formula(glsSt), model, 
        groups), data = data, na.action = na.action)
    if (!missing(subset)) {
        mfArgs[["subset"]] <- asOneSidedFormula(Call[["subset"]])[[2]]
    }
    mfArgs$drop.unused.levels <- TRUE
    dataMod <- do.call("model.frame", mfArgs)
    rn <- origOrder <- row.names(dataMod)  ## rn FEH 6apr03
    if (!is.null(groups)) {
        groups <- eval(parse(text = paste("~1", deparse(groups[[2]]), 
            sep = "|")))
        grps <- getGroups(dataMod, groups,
                          level = length(getGroupsFormula(groups, 
                            asList = TRUE)))
        ord <- order(grps)
        grps <- grps[ord]
        dataMod <- dataMod[ord, , drop = FALSE]
        rn <- rn[ord]  ## FEH 4apr03
        revOrder <- match(origOrder, rn)  ## was row.names(dataMod)) FEH
    }
    else grps <- NULL

    X <- model.frame(model, dataMod)
    dul <- .Options$drop.unused.levels  # FEH 29mar03
    if(!length(dul) || dul) {
      on.exit(options(drop.unused.levels=dul))
      options(drop.unused.levels=FALSE)
    }
    X <- Design(X)
    atrx <- attributes(X)
    desatr <- atrx$Design
    mt <- atrx$terms
    attr(X,'Design') <- NULL

    contr <- lapply(X, function(el) if (inherits(el, "factor")) 
        contrasts(el))
    contr <- contr[!unlist(lapply(contr, is.null))]
    X <- model.matrix(model, X)
    dimnames(X)[[2]] <- cn <- c('Intercept',desatr$colnames)  ## FEH 3apr03
    y <- eval(model[[2]], dataMod)
    N <- nrow(X)
    p <- ncol(X)
    parAssign <- attr(X, "assign")
    fTerms <- terms(as.formula(model))
    namTerms <- attr(fTerms, "term.labels")
    if (attr(fTerms, "intercept") > 0) {
      namTerms <- c("Intercept", namTerms)  # FEH 4apr03
    }
    namTerms <- factor(parAssign, labels = namTerms)
    parAssign <- split(order(parAssign), namTerms)
    ## Start FEH 4apr03
    if(B > 0) {
      bootcoef <- matrix(NA, nrow=B, ncol=p, dimnames=list(NULL,cn))
      bootcorr <- numeric(B)
      Nboot    <- integer(B)
      if(length(grps)) {
        obsno    <- split(1:N,grps)
        levg     <- levels(grps)
        ng       <- length(levg)
        if(!length(levg)) stop('program logic error')
      } else {
        obsno <- 1:N
        levg  <- NULL
        ng    <- N
      }
    }
    for(j in 0:B) {
      if(j == 0) s <- 1:N else {
        if(ng == N) s <- sample(1:N, N, replace=TRUE) else {
          grps.sampled <- sample(levg, ng, replace=TRUE)
          s <- unlist(obsno[grps.sampled])
          dataMods <- dataMod[s,]
          if(!dupCluster) {
            grp.freqs <- table(grps)
            newgrps <- factor(rep(paste('C',1:ng,sep=''),
                                  table(grps)[grps.sampled]))
            dataMods$id <- newgrps
          }
        }
        Nboot[j] <- Nb <- length(s)
        if(pr) cat(j,'')
      }
    attr(glsSt, "conLin") <-
      if(j==0)
        list(Xy = array(c(X, y), c(N, p + 1),
               list(rn, c(cn, deparse(model[[2]])))), 
             dims = list(N = N, p = p, REML = as.integer(REML)),
             logLik = 0) else
      list(Xy = array(c(X[s,,drop=FALSE], y[s]), c(Nb, p + 1),
               list(rn[s], c(cn, deparse(model[[2]])))), 
             dims = list(N = Nb, p = p, REML = as.integer(REML)),
             logLik = 0)
        ## FEH colnames(X) -> cn, ncol(X) -> p, j>0 case above 4apr03
      glsEstControl <- controlvals[c("singular.ok", "qrTol")]
      glsSt <- Initialize(glsSt, if(j==0) dataMod else dataMods,
                          glsEstControl)
      parMap <- attr(glsSt, "pmap")
      numIter <- numIter0 <- 0
      repeat {
        oldPars <- c(attr(glsSt, "glsFit")[["beta"]], coef(glsSt))
        if (length(coef(glsSt))) {
          co <- c(coef(glsSt))  ## FEH
          if(opmeth == 'optimize' && co > 1) opmeth <- 'optim'
          best <- if(opmeth=='optimize')
            optimize(function(z)
                     -logLik(glsSt,z), lower=-12, upper=12)$minimum else
          optim(fn = function(glsPars)
                -logLik(glsSt, glsPars),
                par = co,  ## FEH
                method = "BFGS", 
                control = list(trace = controlvals$msVerbose, 
                  reltol = if (numIter == 0)
                  controlvals$msTol else 100 * 
                  .Machine$double.eps,
                  maxit = controlvals$msMaxIter))$par
          coef(glsSt) <- best   ## FEH
        }
        attr(glsSt, "glsFit") <- glsEstimate(glsSt, control = glsEstControl)
        if (!needUpdate(glsSt)) 
          break
        numIter <- numIter + 1
        glsSt <- update(glsSt, if(j==0) dataMod else dataMods)  ## FEH
        aConv <- c(attr(glsSt, "glsFit")[["beta"]], coef(glsSt))
        conv <- abs((oldPars - aConv)/ifelse(aConv == 0, 1, aConv))
        aConv <- c(beta = max(conv[1:p]))
        conv <- conv[-(1:p)]
        for (i in names(glsSt)) {
          if (any(parMap[, i])) {
            aConv <- c(aConv, max(conv[parMap[, i]]))
            names(aConv)[length(aConv)] <- i
          }
        }
        if (verbose) {
          cat("\nIteration:", numIter)
          cat("\nObjective:", format(aNlm$value), "\n")
          print(glsSt)
          cat("\nConvergence:\n")
          print(aConv)
        }
        if (max(aConv) <= controlvals$tolerance) {
          break
        }
        if (numIter > controlvals$maxIter) {
          stop("Maximum number of iterations reached without convergence.")
        }
      }
      if(j > 0) {
        bootcoef[j,] <- attr(glsSt, "glsFit")[["beta"]]
        bootcorr[j]  <- coef(glsSt$corStruct, unconstrained=FALSE)
      }
      if(j==0) glsSt0 <- glsSt  ## FEH 4apr03
    }  ## end bootstrap reps
    if(pr && B > 0) cat('\n')
    glsSt <- glsSt0   ## FEH
        
    glsFit <- attr(glsSt, "glsFit")
    namBeta <- names(glsFit$beta)
    attr(parAssign, "varBetaFact") <- varBeta <- glsFit$sigma * 
      glsFit$varBeta * sqrt((N - REML * p)/(N - p))
    varBeta <- crossprod(varBeta)
    dimnames(varBeta) <- list(namBeta, namBeta)
    Fitted <- fitted(glsSt)
    if (!is.null(grps)) {
      grps <- grps[revOrder]
      Fitted <- Fitted[revOrder]
      Resid <- y[revOrder] - Fitted
      attr(Resid, "std") <- glsFit$sigma/(varWeights(glsSt)[revOrder])
    }
    else {
      Resid <- y - Fitted
      attr(Resid, "std") <- glsFit$sigma/(varWeights(glsSt))
    }
    if (controlvals$apVar && FALSE) ## FEH 3apr03
      apVar <-
        glsApVar(glsSt, glsFit$sigma, .relStep = controlvals[[".relStep"]],  
                 minAbsPar = controlvals[["minAbsParApVar"]],
                 natural = controlvals[["natural"]])
    else {
      apVar <- "Approximate variance-covariance matrix not available"
    }
    dims <- attr(glsSt, "conLin")[["dims"]]
    dims[["p"]] <- p
    attr(glsSt, "conLin") <- NULL
    attr(glsSt, "glsFit") <- NULL
    estOut <- list(modelStruct = glsSt, dims = dims, contrasts = contr, 
                   coefficients = glsFit[["beta"]], varBeta = varBeta,
                   sigma = glsFit$sigma, 
                   apVar = apVar, logLik = glsFit$logLik,
                   numIter = if(needUpdate(glsSt)) numIter else numIter0,  
                   groups = grps, call = Call, method = method,
                   fitted = Fitted,  
                   residuals = Resid, parAssign = parAssign,
                   Design=desatr,assign=DesignAssign(desatr,1,mt),
                   formula=model, opmeth=opmeth,
                   B=B, boot.Coef=if(B > 0) bootcoef,
                   boot.Corr=if(B > 0) bootcorr,
                   Nboot=if(B > 0) Nboot,
                   var=if(B > 0) var(bootcoef))
    ## Last 2 lines FEH 29mar03
    if (inherits(data, "groupedData")) {
      attr(estOut, "units") <- attr(data, "units")
      attr(estOut, "labels") <- attr(data, "labels")
    }
    attr(estOut, "namBetaFull") <- colnames(X)
    class(estOut) <- c('glsD','Design','gls')  ## FEH 29mar03
    estOut
  }

#summary.glsD <- function(...) summary.Design(...)
#predict.glsD <- function(...) predict.Design(...)

print.glsD <- function(x, digits=4, ...) {
## Following taken from print.gls with changes marked FEH
  ## In R 1.7 nlme namespace does not export glsEstimate
  if(.R. && !existsFunction('summary.gls'))   # FEH
      summary.gls <- getFromNamespace('summary.gls','nlme')    # FEH

    dd <- x$dims
    mCall <- x$call
    if (inherits(x, "gnls")) {
        cat("Generalized nonlinear least squares fit\n")
    }
    else {
        cat("Generalized least squares fit by ")
        cat(ifelse(x$method == "REML", "REML\n", "maximum likelihood\n"))
    }
    cat("  Model:", deparse(mCall$model), "\n")
    cat("  Data:", deparse(mCall$data), "\n")
    if (!is.null(mCall$subset)) {
        cat("  Subset:", deparse(asOneSidedFormula(mCall$subset)[[2]]), 
            "\n")
    }
    if (inherits(x, "gnls")) {
        cat("  Log-likelihood: ", format(x$logLik), "\n", sep = "")
    }
    else {
        cat("  Log-", ifelse(x$method == "REML", "restricted-", 
            ""), "likelihood: ", format(x$logLik), "\n", sep = "")
    }
##    cat("\nCoefficients:\n")  FEH
##    print(coef(x))            FEH and following 9 lines
    cat('\n')
    if(any(names(x)=='var') && length(x$var)) {
      cat('Using bootstrap variance estimates\n\n')
      se <- sqrt(diag(x$var))
      beta <- coef(x)
      zTable <- cbind(Coef=format(beta,digits=digits),
                     'S.E'=format(se, digits=digits),
                      Z   =format(beta/se, digits=2),
                      'Pr(>|Z|)'=format.pval(2*pnorm(-abs(beta/se)),digits=4))
      print(zTable, quote=FALSE)
    } else print(summary.gls(x)$tTable)
    
    cat("\n")
    if (length(x$modelStruct) > 0) {
        print(summary(x$modelStruct))
    }
    cat("Degrees of freedom:", dd[["N"]], "total;", dd[["N"]] - 
        dd[["p"]], "residual\n")
    cat("Residual standard error:", format(x$sigma), "\n")

    cat('Clusters:',length(unique(x$groups)),'\n')
    if(x$B > 0) {
      cat('Bootstrap repetitions:',x$B,'\n')
      tn <- table(x$Nboot)
      if(length(tn) > 1) {
        cat('Table of Sample Sizes used in Bootstraps\n')
        print(tn)
      } else cat('Bootstraps were all balanced with respect to clusters\n')
      dr <- diag(x$varBeta)/diag(x$var)
      cat('Ratio of Original Variances to Bootstrap Variances\n')
      print(round(dr,2))
      cat('Bootstrap Nonparametric 0.95 Confidence Limits for Correlation Parameter\n')
      r <- round(quantile(x$boot.Corr, c(.025,.975)),3)
      names(r) <- c('Lower','Upper')
      print(r)
    }
    invisible()
  }

Varcov.glsD <- function(object, regcoef.only=FALSE, ...)
  if(any(names(object)=='var') && length(object$var))
  object$var else object$varBeta

predict.glsD <- 
  function(object, newdata,
           type=c("lp","x","data.frame","terms","adjto","adjto.data.frame",
             "model.frame"),
           se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
           incl.non.slopes, non.slopes, kint=1,
           na.action=na.keep, expand.na=TRUE, center.terms=TRUE, ...)
  predictDesign(object, newdata, type, se.fit, conf.int, conf.type,
                incl.non.slopes, non.slopes, kint,
                na.action, expand.na, center.terms, ...)


latex.glsD <- function(...) latexDesign(...)
#Function to divide x (e.g. Cox predicted survical at time y created by
#survest) into g quantile groups, get Kaplan-Meier estimates at time u
#(a scaler), and to return a matrix with columns x=mean x in
# quantile, n=#subjects, events=#events, and KM=K-M survival at time u,
# std.err = s.e. of log-log K-M
#Failure time=y   Censoring indicator=e
#Instead of supplying g, the user can supply the number of subjects to have
#in the quantile group on the average; then g will be computed.  The default
#m is 50, so the default g is n/50.
#If cuts is given (e.g. cuts=c(0,.1,.2,...,.9,.1)), it overrides m and g.
#Set pl=T to plot results.  If pl=T, units attribute of y applies.  
#Default is "Day".
#xlab and ... are passed to plot() if pl=T.  Default xlab is label(x)
#if it is defined, otherwise the name of the calling argument for x.
#
#Author: Frank Harrell   8 May 91
#Modified               18 May 91 - made x general
#			24 Jun 91 - added loglog argument
#			 4 Aug 91 - use cut2 for grouping
#			27 Mar 92 - use Surv()
#			17 Apr 92 - add lines with err bars, allow conf.int=F,
#					added lty, add
#			3 Jun 92  - add cex.subtitle
#		       23 Jul 92  - add ylab
#		       05 Sep 92  - added dimension check
groupkm <- function(x, Srv, m=50, g, 
                    cuts, u, pl=FALSE, loglog=FALSE, conf.int=.95, xlab, ylab,
                    lty=1, add=FALSE,
                    cex.subtitle=.7, ...) {
  if(.R.) library(survival)
  if(missing(u))stop("u (time point) must be given")
  if(missing(xlab)) xlab <- label(x)
  if(xlab=="") xlab <- as.character(sys.call())[2]
  s <- !(is.na(x)|is.na(Srv[,1])|is.na(Srv[,2]))
  x <- x[s]; Srv <- Srv[s,]
  x[abs(x)<1e-10] <- 0 #cut2 doesn't work with tiny x
  e <- Srv[,2]
  if(nrow(Srv)!=length(x))stop("lengths of x and Srv must match")
  unit <- attr(Srv,"units")
  if(is.null(unit) || unit=="") unit <- "Day"
  if(!missing(cuts)) q <- cut2(x, cuts)
  else if(!missing(g)) q <- cut2(x, g=g)
  else q <- cut2(x, m=m)
  if(any(table(q)) < 2) warning('one interval had < 2 observations')
  q <- oldUnclass(q)  #need integer
  g <- length(levels(q))

  km <- single(g)
  pred <- km
  std.err <- km
  events <- integer(g)
  numobs <- events

#f <- survfit.km(q, Srv, conf.int=conf.int, conf.type="log-log")
#if(is.null(f$strata)) {nstrat <- 1; stemp <- rep(1, length(f$time))}
#else { nstrat <- length(f$strata); stemp <- rep(1:nstrat,f$strata)}
#This is more efficient but doesn't handle empty strata

  for(i in 1:g) {
	s <- q==i
	nobs <- sum(s); ne <- sum(e[s])
	if(nobs < 2) {   ## was ==0 25apr03
      numobs[i] <- 0
      events[i] <- 0
      pred[i] <- if(nobs==1) mean(x[s], na.rm=TRUE) else NA
      km[i] <- NA
      std.err[i] <- NA	} else {
        pred[i] <- mean(x[s], na.rm=TRUE)
        ##	f <- surv.fit(y[s], e[s])
        dummystrat <- rep(1, nobs)
        attributes(dummystrat) <- list(class="factor",levels="1")
        f <- survfit.km(dummystrat,Srv[s,], conf.type="log-log") 
        ##doesn't need conf.int since only need s.e.
        tt <- c(0, f$time)
        ss <- c(1, f$surv)
        se <- c(0, -f$std.err/logb(f$surv))
        ##	tm <- max((1:length(tt))[max(tt[tt<=u+1e-6])==tt])
        tm <- max((1:length(tt))[tt<=u+1e-6])
        km[i] <- ss[tm]
        std.err[i] <- se[tm]
        numobs[i] <- nobs
        events[i] <- ne
        n <- length(tt)
        if(u > tt[n]+1e-6 & ss[n]>0) {
          km[i] <- NA
          std.err[i] <- NA
        }
      }				}
  z <- cbind(x=pred, n=numobs, events=events, KM=km, 
             std.err=std.err)

  if(pl) {
	y <- km
    if(conf.int) {
      zcrit <- qnorm((conf.int+1)/2)
      low <- km^exp(zcrit*std.err); hi <- km^exp(-zcrit*std.err)
    }
	if(missing(ylab))
      ylab <- paste("Kaplan-Meier ",format(u),"-",unit," Survival",sep="")
	if(loglog) {
      y <- logb(-logb(y))
      if(conf.int) {
		low <- logb(-logb(low))
		hi <- logb(-logb(hi))
      }
      if(missing(ylab))
		ylab <- paste("log(-log Kaplan-Meier ",format(u),unit,
                      " Survival",sep="")
    }
	if(!add)plot(pred, y, xlab=xlab, ylab=ylab, type="n", ...)
	lines(pred, y, lty=lty)
	if(conf.int)errbar(pred, y, hi, low, add=TRUE, ...)
	if(!is.logical(cex.subtitle)) {
      nn <- sum(numobs,na.rm=TRUE)
      mm <- round(nn/g)
      title(sub=paste("n=",nn," d=",sum(events,na.rm=TRUE),
              ", avg. ",mm," patients per group",sep=""),
            adj=0,cex=cex.subtitle)
    }
  }
z
}
hazard.ratio.plot<-function(x,Srv,which,
                            times,e=30,subset,conf.int=.95,legendloc=NULL,
                            smooth=TRUE,pr=FALSE,pl=TRUE,add=FALSE,ylim,cex=.5,xlab="t",ylab,
                            antilog=FALSE, ...)						{

  if(missing(ylab)) ylab <- if(antilog)"Hazard Ratio" else "Log Hazard Ratio"

  trans <- if(antilog) function(x) exp(x) else function(x) x

  if(is.matrix(x)) 		{
	nam <- dimnames(x)[[2]]
	if(!length(nam)) nam <- paste("x[",1:ncol(x),"]",sep="")
  }	else	{
	nam <- label(x)
	x <- as.matrix(oldUnclass(x))
	if(!length(nam)) nam <- ""		}

  y <- Srv[,1];   event <- Srv[,2]
  if(length(y)!=nrow(x))stop("number of rows in x must be length of y")
  nx <- ncol(x)
  if(missing(which)) which <- 1:nx

  labele<-attr(Srv, "event.label")
  if(!length(labele)) labele <- ""

  isna <- is.na(matxv(x,rep(1,nx)) + y + event)

  if(!missing(subset))isna <- isna | (!subset)
  x <- x[!isna,,drop=FALSE]
  if(length(dimnames(x)[[2]])==0) dimnames(x) <- list(NULL,paste("x",1:nx,sep=""))
  y <- y[!isna]
                                        #Srv <- Srv[!isna,]
  event <- event[!isna]

  if(!missing(times))uft<-c(0,sort(times),1000000) else	{
	nblock<-max(round(sum(event)/e),2)
	uft<-c(0,quantile(y[event==1],seq(0,1,length=nblock+1))[2:nblock]
           ,1000000)
	uft <- unique(uft)
  }
  thr<-NULL
  lhr<-NULL
  se<-NULL
  storage.mode(x) <- "double"
  for(i in seq(length(uft)-1))	{
	s<-y>=uft[i]
	tt<-pmin(y[s],uft[i+1])
	ev<-event[s] & (y[s]<=uft[i+1])
	if(sum(ev)>nx)			{
      cox <- coxphFit(x[s,,drop=FALSE], cbind(tt,ev),
                       iter.max=10, eps=.0001, method="efron")
      if(!is.character(cox))		{
        if(pr)	{
          r <- range(tt)
          cat(paste("Time interval:",format(r[1]),"-",
                    format(r[2]),"  At Risk:",sum(s),
                    "  Events:",sum(ev),"\n"))
          k <- cbind(cox$coefficients,sqrt(diag(cox$var)))
          dimnames(k) <- list(names(cox$coefficients),
                              c("Coef","S.E."))
          print(k)	}
        tmid<-mean(y[y>=uft[i] & y<=uft[i+1]])
        thr<-c(thr,tmid)
        lhr<-cbind(lhr,cox$coef)
        se<-cbind(se,sqrt(diag(cox$var)))
      }	}
  }

  if(!pl) return(list(time=thr,log.hazard.ratio=lhr,se=se))

  zcrit<-qnorm((1+conf.int)/2)
  for(j in which)					{
	lhrj <- lhr[j,]
	sej <- se[j,]
	labelx <- nam[j]
	if(missing(ylim)) ylim <- trans(range(c(lhrj+zcrit*sej,lhrj-zcrit*sej)))
	if(!add)	{
      oldpar <- par(c('err','mar'))
      on.exit(par(oldpar))
      oldmar <- oldpar$mar
      if(labelx!="" & labele!="")oldmar[1]<-oldmar[1]+1
      par(err=-1,mar=oldmar)
      plot(thr,trans(lhrj),xlab=xlab,ylim=ylim,ylab=ylab,...)
    }
    ##Next line had ...
    else points(thr,trans(lhrj))
    lines(thr,trans(lhrj))
    lines(thr,trans(lhrj+zcrit*sej),lty=2)
    lines(thr,trans(lhrj-zcrit*sej),lty=2)
    leg<-c("Subset Estimate",paste(format(conf.int),"C.L."))
    ltype<-1:2
    if(smooth & length(thr)>3)	{
      ##Next line did have ...
      lines(supsmu(thr,trans(lhrj)),lty=3)
      leg<-c(leg,"Smoothed")
      ltype<-c(ltype,3)	}

    if(!add)	{
      
      labels<-""
      if(labelx!="")labels<-paste("Predictor:",labelx,"\n",sep="")
      if(labele!="")labels<-paste(labels,"Event:",labele,sep="")
      title(sub=labels,adj=1,cex=cex)

      if(!interactive() && !length(legendloc))legendloc<-"ll"
      if(!length(legendloc))	{
        cat("Click left mouse button at upper left corner for legend\n")
        z<-locator(1)
        legendloc<-"l"	}
      else if(legendloc[1]!="none")	{
        if(legendloc[1]=="ll")z<-list(x=par("usr")[1],y=par("usr")[3])
        else z<-list(x=legendloc[1],y=legendloc[2]) 		  }	
      if(legendloc[1]!="none")legend(z,leg,lty=ltype,cex=cex,bty="n")
    }  }
  list(time=thr,log.hazard.ratio=lhr,se=se)
}

#ia.operator.s - restricted interaction operators for use with Design
#F. Harrell  8 Nov 91

#Set up proper attributes for a restricted interaction for a model
#such as y ~ rcs(x1) + rcs(x2) + x1 %ia% x2 or x1 %ia% rcs(x2)
#or rcs(x1) %ia% x2

"%ia%" <- function(x1, x2)						{
a1 <- attributes(x1)
a2 <- attributes(x2)
nam <- as.character(sys.call())[-1]

redo <- function(x,nam)					{
   if(is.null(attr(x,"assume.code")))		{
      if(!is.null(oldClass(x)) && oldClass(x)[1]=="ordered")
         x <- scored(x, name=nam)
      else if(is.character(x) | is.category(x)) x <- catg(x, name=nam)
      else if(is.matrix(x)) x <- matrx(x, name=nam)
      else x <- asis(x, name=nam)		}
   ass <- attr(x,"assume.code")
   nam <- attr(x,"name")

   if(ass==5)				{
      colnames <- attr(x,"colnames")
      len <- length(attr(x,"parms"))-1	}
   else if(ass==8)			{
      prm <- attr(x,"parms")
      colnames <- paste(nam,"=",prm[-1],sep="")
      len <- length(prm)-1		}
   else if(ass==7)			{
      prm <- attr(x,"parms")
      colnames <- c(nam,paste(nam,"=",prm[-(1:2)],sep=""))
      len <- length(prm)-1		}
   else							{
      if(is.null(ncol(x)))		{
         len <- 1
         colnames <- nam		} else 	{
         colnames <- dimnames(x)[[2]]
         len <- ncol(x)				}
							}

   attr(x,"colnames") <- colnames
   attr(x,"len") <- len
   if(ass==8) attr(x,"nonlinear") <- rep(FALSE, len)
   x							}

x1 <- redo(x1,nam[1])
x2 <- redo(x2,nam[2])
a1 <- attributes(x1)
a2 <- attributes(x2)
n1 <- a1$colnames
n2 <- a2$colnames
nl1 <- a1$nonlinear
nl2 <- a2$nonlinear
as1 <- a1$assume.code
as2 <- a2$assume.code

l1 <- a1$len
l2 <- a2$len
if(any(nl1) & any(nl2))	nc <- l1+l2-1   
else nc <- l1*l2
if(is.matrix(x1)) nr <- nrow(x1) 
else nr <- length(x1)
x <- matrix(single(1),nrow=nr,ncol=nc)
name <- character(nc)
parms <- matrix(integer(1),nrow=2,ncol=nc+1)
nonlinear <- logical(nc)

k <- 0
if(!is.factor(x1)) x1 <- as.matrix(x1)
if(!is.factor(x2)) x2 <- as.matrix(x2)
for(i in 1:l1)					{
   if(as1==5 | as1==8) x1i <- oldUnclass(x1)==(i+1)
   else x1i <- x1[,i]
   for(j in 1:l2)			{
	#Remove doubly nonlinear terms
	if(nl1[i] & nl2[j]) break
		k <- k + 1
		if(as2==5 | as2==8) x2j <- oldUnclass(x2)==(j+1)
		else x2j <- x2[,j]
		x[,k] <- x1i * x2j
		name[k] <- paste(n1[i],"*",n2[j])
		parms[,k+1] <- c(nl1[i],nl2[j])
		nonlinear[k] <- nl1[i] | nl2[j]
					}	}

dimnames(x) <- list(NULL, name)
attr(x,"ia") <- c(a1$name, a2$name)
attr(x,"parms") <- parms
attr(x,"nonlinear") <- nonlinear
attr(x,"assume.code") <- 9
attr(x,"name") <- paste(a1$name,"*",a2$name)
attr(x,"label") <- attr(x,"name")
attr(x,"colnames") <- name
attr(x,"class") <- "Design"
storage.mode(x) <- "single"
x									}
ie.setup <- function(failure.time, event, ie.time, break.ties=FALSE) {

s <- !is.na(ie.time)
if(all(s)) warning('every subject had an intervening event')
if(!any(s)) stop('no intervening events')
if(any(ie.time[s] > failure.time[s])) 
  stop('an ie.time was after a failure.time')
if(break.ties) {
	mindif <-
      min(diff(sort(unique(failure.time[!is.na(failure.time)]))))
    ## 8Nov01 Thanks: Josh Betcher
	k <- s & (ie.time==failure.time)
	if(sum(k)==0) warning('break.times=T but no ties found')
	ie.time[k] <- ie.time[k] - runif(sum(k),0,mindif)
}
 
if(any(ie.time[s]==failure.time[s])) 
  stop('ie.times not allowed to equal failure.times')

n <- length(failure.time)
reps <- ifelse(is.na(ie.time), 1, 2)
subs <- rep(1:n, reps)

start <- end <- ev <- ie.status <- vector('list', n)
start[]     <- 0
end[]       <- failure.time
ev[]        <- event
ie.status[] <- 0
for(i in seq(along=s)[s]) {
  start[[i]]  <- c(0, ie.time[i])
  end[[i]]    <- c(ie.time[i], failure.time[i])
  ev[[i]]     <- c(0, event[i])
  ie.status[[i]] <- c(0, 1)
}

start     <- unlist(start)
end       <- unlist(end)
ev        <- unlist(ev)
ie.status <- unlist(ie.status)

u <- units(failure.time)
units(end) <- if(u=='')'Day' else u

s <- !is.na(start+end) & (end <= start)
if(any(s)) {
  cat('stop time <= start time:\n')
  print(cbind(start=start[s], end=end[s]))
  stop()
}

S <- Surv(start, end, ev)

list(S=S, ie.status=ie.status, subs=subs, reps=reps)
}

latexDesign <- function(object,
		file=paste(first.word(deparse(substitute(f))),".tex",sep=""),
		append=FALSE, which=1:p, varnames, columns=65, prefix=NULL, 
		inline=FALSE, before=if(inline)"" else "& &", intercept, 
		pretrans=TRUE, digits=.Options$digits)
{
f <- object	
at <- f$Design
if(!length(at)) at <- getOldDesign(f)
name <- at$name
ac <- at$assume.code
p <- length(name)
nrp <- num.intercepts(f)
#f$term.labels does not include strat
TL <- attr(terms(f),"term.labels")
tl <- TL
#Get inner transformations
#from <- c("asis","pol","lsp","rcs","catg","scored","strat","matrx","I")
#from <- paste(from,"(\\(.*\\))",sep="")

from <- c('asis(*)','pol(*)','lsp(*)','rcs(*)','catg(*)','scored(*)',
  'strat(*)','matrx(*)','I(*)')
to   <- rep('*',9)

TLi <- paste("h(",sedit(TL, from, to),")",sep="")
#TLi <- paste("h(",translate(TL, from, "\\1"),")",sep="")  
#change wrapping function to h()
h <- function(x,...) deparse(substitute(x))
for(i in (1:p)[ac!=9]) TLi[i] <- eval(parse(text=TLi[i]))
TLi <- ifelse(TLi==name | ac==1 | ac==9, "", TLi)
anytr <- any(TLi!="")
if(!missing(varnames)) {
   if(length(varnames)!=sum(ac!=9)) stop("varnames is wrong length")
   vn <- name; vn[ac!=9] <- varnames; varnames <- vn
   tl <- sedit(tl, name, varnames)  # was translate
   if(anytr) TLi <- sedit(TLi, name, varnames) # was translate
			}
else varnames <- name
lnam <- nchar(varnames)
# digits at end of name -> subscript, change font
#vnames  <- sys('sed -e "s/\\([0-9]\\)$/_{\\\\mit \\1}/g"', varnames)
#vnames2 <- sys('sed -e "s/\\([0-9]\\)$/_{\\\\\\\\mit \\1}/g"',
#			varnames)
vnames <- sedit(varnames, '*$', '_{\\mit *}', test=all.digits)

if(is.character(which)) {
   wh <- charmatch(which, name, 0)
   if(any(wh==0))stop(paste("variable name not in model:",
	paste(which[wh==0], collapse=" ")))
}

interaction <- at$interactions
if(length(interaction)==0) interaction <- 0

parms <- at$parms

if(TRUE)					{
#If any interactions to be printed, make sure all main effects are included
ia <- ac[which]==9
if(length(which)<p & any(ia))	{
  for(i in which[ia]) which <- c(which,parms[[name[i]]][,1])
  which <- which[which>0]
  which <- sort(unique(which))	} 
					} 

#trancom <- 'sed -e "s/sqrt(\\(.*\\))/\\\\sqrt{\\1}/" -e "s/log(/\\\\log(/" -e "s/I(\\(.*\\))/[\\1]/" -e "s/\\[1\\/(\\(.*\\))\\]/[\\1]^{-1}/" -e "s/pmin(/\\\\min( /" -e "s/pmax(/\\\\max( /"'
from <- c('sqrt(*)',  'log(',  'I(*)', '1/(*)',   'pmin(', 'pmax(')
to   <- c('\\sqrt{*}','\\log(','[*]',  '(*)^{-1}','\\min(','\\max(')
#tl <- sys(trancom, tl)
tl  <- sedit(tl, from, to)
#tl <- ifelse(TL==varnames,vnames,tl)
tl <- sedit(tl, varnames, vnames)  # was translate.  was vnames2 29Apr96
#varnames <-  paste("{\\rm ", vnames, "}", sep="")
ltl <- nchar(tl)
tl <- paste("{\\rm ", tl, "}", sep="")
if(anytr)	{
#  TLi <- sys(trancom, TLi)
  TLi <- sedit(TLi, from, to)
  TLi <- sedit(TLi, varnames, vnames)  # was translate  was vnames2 29Apr96
  TLi <- ifelse(TLi=="","",paste("{\\rm ", TLi, "}", sep=""))
		}

varnames <- paste("{\\rm ", vnames, "}", sep="")

Two.Way <- function(prm,Nam,nam.coef,lNam,cof,coef,f,columns,lcof,varnames,
			lnam, at, digits=digits)					{
  i1 <- prm[1,1]; i2 <- prm[2,1]
  num.nl <- any(prm[1,-1])+any(prm[2,-1])
 #If single factor with nonlinear terms, get it as second factor
 #Otherwise, put factor with most # terms as second factor
  rev <- FALSE
  if((num.nl==1 & any(prm[1,-1])) || (length(Nam[[i1]])>length(Nam[[i2]])))
	{ i1 <- i2; i2 <- prm[1,1]; rev <- TRUE }
  N1 <- Nam[[i1]]; N2 <- Nam[[i2]]
  n1 <- nam.coef[[i1]]; n2 <- nam.coef[[i2]]
  q <- NULL; cur <- ""; m <- 0
  for(j1 in 1:length(N1))					{
   nam1 <- nam.coef[[i1]][j1]
   l1 <- lNam[[i1]][j1]
   lN2 <- length(N2)
   cnam <- if(rev) paste(nam.coef[[i2]],"*",nam1) else
		paste(nam1, "*", nam.coef[[i2]])
   mnam <- match(cnam, names(cof), nomatch=0)
   act <- mnam[mnam>0]
   lN2.act <- length(act)
  #Check if restricted interaction between a rcs and another nonlinear
  #var, i.e. >1 2nd term possible, only 1 (linear) there, and at first 
  #nonlinear term of rcs
   if(lN2.act==1 & lN2>1 & at$assume.code[i1]==4 & j1==2)	{
     if(cur!="") { q <- c(q, cur); m <- 0; cur <- ""}
     v <- paste("+",N2[1],"[",sep=""); n <- lNam[[i2]][1]
     if(m + n > columns) { q <- c(q, cur); cur <- ""; m <- 0}
     cur <- paste(cur, v, sep="")
     m <- m+n
     cnam <- paste(nam.coef[[if(rev)i2 else i1]][1], "*",
                   nam.coef[[if(rev)i1 else i2]][-1])  # rev 4Dec00
     v <- rcspline.restate(at$parms[[at$name[i1]]], c(0, coef[cnam]), 
		x=varnames[i1],
		lx=lnam[i1], columns=columns, before="", after="",
		begin=cur, nbegin=m, digits=digits)
     m <- attr(v, "columns.used")+1   #+1 for "]"
     v <- attr(v, "latex")
     j <- length(v)
     if(j>1) q <- c(q, v[-j])
     cur <- paste(v[j], "]")
     break
								}
   else if(lN2.act==1)		{
     v <- paste(cof[act],"\\:",N1[j1],"\\:\\times\\:", N2[mnam>0], sep="")
     n <- l1+lNam[[i2]][mnam>0]+2
     if(m + n > columns) { q <- c(q, cur); cur <- ""; m <- 0}
     cur <- paste(cur, v, sep="")
     m <- m + n			}
   else	if(lN2.act>0)					{
    if(cur!="") { q <- c(q, cur); m <- 0; cur <- ""}
    v <- paste("+",N1[j1],"[",sep=""); n <- l1 + 1
    if(m + n > columns) { q <- c(q, cur); cur <- ""; m <- 0}
    cur <- paste(cur, v, sep="")
    m <- m + n

    if(at$assume.code[i2]==4 & !any(mnam==0))	{
   #rcspline, interaction not restricted
     v <- rcspline.restate(at$parms[[at$name[i2]]], coef[act], x=varnames[i2],
 		lx=lnam[i2],
		columns=columns, before="", after="", 
		begin=cur, nbegin=m, digits=digits)
     m <- attr(v, "columns.used") + 1   #1 for "]"
     v <- attr(v, "latex")
     j <- length(v)
     if(j>1) q <- c(q, v[-j])
     cur <- paste(v[j],"]")			}

    else 				{
      for(j2 in 1:lN2)		{
       l <- mnam[j2]
       if(l>0)	    	{	#not a restricted-out nonlinear term
        if(j2==1 && substring(cof[l],1,1)=="+") cof[l] <- substring(cof[l],2)
        v <- paste(cof[l],"\\:",N2[j2],sep="")
        n <- lcof[l]+lNam[[i2]][j2]
        if(m + n > columns) { q <- c(q, cur); cur <- ""; m <- 0}
        cur <- paste(cur, v, sep="")
        m <- m + n
			}
  				 }
     cur <- paste(cur, "]")
					}
							}	}
 if(cur!="") q <- c(q, cur)
 attr(q, "columns.used") <- m
 q
									}

Three.Way <- function(prm,Nam,nam.coef,lNam,cof,coef,f,columns,lcof,at){
  i1 <- prm[1,1]; i2 <- prm[2,1]; i3 <- prm[3,1]
  N1 <- Nam[[i1]]; N2 <- Nam[[i2]]; N3 <- Nam[[i3]]
  q <- NULL; cur <- ""; m <- 0; l <- 0
   for(j3 in 1:length(N3))		{
    for(j2 in 1:length(N2))		 {
     for(j1 in 1:length(N1))		  {
      l <- l + 1
      v <- paste(cof[l], "\\:", N1[j1], "\\:\\times\\:", N2[j2],
                 "\\:\\times\\:", N3[j3], sep="")
      n <- lcof[l] + lNam[[i1]][j1]+lNam[[i2]][j2]+lNam[[i3]][j3] + 3
      if(m + n > columns)	{
       q <- c(q, cur)
       cur <- ""
       m <- 0			}
      cur <- paste(cur, v, sep="")
      m <- m + n			}}}
  q <- c(q, cur)
  attr(q, "columns.used") <- m
  q
}

if(!inline)	{
  tex <- "\\begin{eqnarray*}"
  if(length(prefix)) tex <- c(tex,
      paste("\\lefteqn{",prefix,"=}\\\\",sep=""))
		} else tex <- NULL

cur <- ""
cols <- 0
Coef <- f$coef
if((length(which)==p)&& (nrp==1 | !missing(intercept)))	{
   cof <- if(missing(intercept))format(Coef[1]) else format(intercept)
   cur <- cof
   cols <- nchar(cof)					}

anybrace <- FALSE; anyplus <- FALSE
Nam <- list();  lNam <- list();  nam.coef <- list()

for(i in (1:p)[which])							{
   ass <- ac[i]
   nam <- varnames[i]
   prm <- at$parms[[at$name[i]]]
   if(any(ass==c(5,7,8)))	{
     if(ass==7) prm <- format(prm)
     oprm <- prm
     lprm <- nchar(prm)
     z <- substring(prm,1,1)=="["
     u <- !z & ass==7
#     prm <- sys('sed -e "s/ /\\\\ /g"',prm)
#     prm <- sys('sed -e "s/&/\\\\&/g"',prm)
     prm <- sedit(prm, c(' ','&'), c('\\ ','\\&'))
     prm <- ifelse(z | u, prm, paste("{\\rm ", prm, "}", sep=""))
     prm <- ifelse(z,paste(nam,"\\in ",prm),prm)
     prm <- ifelse(u,paste(nam,"=",prm),prm)
     lprm <- lprm + (z | u)*(lnam[i]+1)
     prm <- paste("\\{", prm, "\\}", sep="")
     anybrace <- TRUE
				}
 #						{
    if(ass != 8) {
      k <- f$assign[[TL[i]]]
      coef <- Coef[k]
      nam.coef[[i]] <- names(coef)
      cof <- format(coef)
      lcof <- nchar(cof)
#cof <- sys('sed -e "s/e+00//" -e "s/e-0\\(.\\)/\\\\!\\\\times\\\\!10^{-\\1}/" -e "s/e-\\(..\\)/\\\\!\\\\times\\\\!10^{-\\1}/" -e "s/e+0\\(.\\)/\\\\!\\\\times\\\\!10^{\\1}/" -e "s/e+\\(..\\)/\\\\!\\\\times\\\\!10^{\\1}/"', cof)
cof <- latexSN(cof)  # 23jun03
#      cof <- sedit(cof, c('e+00','e-0*',                'e-*',
#                    'e+0*',               'e+*'),
#                  c('',    '\\\!\\times\\\!10^{-*}','\\\!\\times\\\!10^{-*}',
#                    '\\\!\\times\\\!10^{*}','\\\!\\times\\\!10^{*}'))

      cof <- ifelse(coef<=0, cof, paste("+", cof, sep=""))
      cof.sp <- cof
      if(ass==2 | ass==10)	{
        r <- grep("times",cof)
        r <- if(length(r)==0) 1:length(cof) else -r
        cof.sp[r] <- paste(cof.sp[r],"\\:",sep="")    
				}
      else if(length(grep("time",cof[1]))==0)
              cof.sp[1] <- paste(cof[1],"\\:",sep="")
      # medium space between constant and variable names if constant
      # does not end in 10^x
    }
      newline <- FALSE
      switch(ass,
         {  nam <- tl[i]; Nam[[i]] <- nam; lNam[[i]] <- ltl[i]
            q <- paste(cof.sp, nam, sep="")
            m <- ltl[i]+lcof},

         {     q <- ""; m <- 0; pow <- 1:prm
               nams <- ifelse(pow==1,nam,paste(nam,"^{",pow,"}",sep=""))
               Nam[[i]] <- nams; lNam[[i]] <- rep(lnam[i],prm)
               for(j in pow) q <- paste(q,cof.sp[j], nams[j], sep="")
            m <- prm*lnam[i]+sum(lcof) },

         {  if(cols>0) {tex <- c(tex, cur); cur <-""; cols <- 0}
            anyplus <- TRUE
            q <- paste(cof.sp[1], nam, sep="")
            m <- lcof[1]+lnam[i]
            nams <- nam; lnams <- lnam[i]
            kn <- format(-prm)
            lkn <- nchar(kn)
            for(j in 1:length(prm))		{
               z <- paste("(", nam, if(prm[j]<0) "+" else NULL, 
			if(prm[j]!=0) kn[j] else NULL, ")_{+}",
                          sep="")
               nams <- c(nams, z); u <- lnam[i]+lkn[j]+2; lnams <- c(lnams,u)
               v <- paste(cof[j+1], z, sep="")
               n <- lcof[j+1]+u
               if(m + n > columns)	{ 
                  cur <- paste(cur, q)
                  tex <- c(tex, cur)
                  cur <- ""
                  cols <- 0
                  q <- ""
                  m <- 0                }
               q <- paste(q, v, sep="")
               m <- m + n
						}
            Nam[[i]] <- nams; lNam[[i]] <- lnams},

         {  q <- rcspline.restate(prm, coef, x=nam, lx=lnam[i], columns=columns,
                                  before="",after="",digits=digits)
            anyplus <- TRUE
            m <- attr(q, "columns.used")
            nn <- nam; ln <- lnam[i]
            for(j in 1:(length(prm)-2))	{
              nam <- paste(nam, "'", sep=""); nn <- c(nn, nam)
              ln <- c(ln, lnam[i]+j)	}
            Nam[[i]] <- nn       #Two.Way only needs first name
            lNam[[i]] <- ln      #for 2nd-order ia with 1 d.f. (restr ia)
                                 #Three.Way needs original design matrix
            q <- attr(q, "latex")
            if(substring(sedit(q[1]," ",""),1,1)!="-")   # was translate
              q[1] <- paste("+", q[1], sep="")
            j <- length(q)
            if(cur!="") {tex <- c(tex,cur); cur <- ""; cols <- 0}
            if(j>1) { tex <- c(tex, q[-j]); q <- q[j]}} ,
         {  Nam[[i]] <- prm[-1]; lNam[[i]] <- lprm[-1]
            if(cols>0) {tex <- c(tex,cur); cur <- ""; cols <- 0}
            q <- ""
            m <- 0
            for(j in 2:length(prm))		{
               v <- paste(cof[j-1], prm[j], sep="")
               n <- lcof[j-1]+lprm[j]
               if(m + n > columns)	{
                  cur <- paste(cur,q)
                  tex <- c(tex, cur)
                  cur <- ""
                  cols <- 0
                  q <- ""
                  m <- 0		}
               q <- paste(q, v, sep="")
               m <- m + n			}},
	q <- "",
         
         {  if(cols>0) {tex <- c(tex,cur); cur <- ""; cols <- 0}
            q <- paste(cof.sp[1], nam, sep="")
            m <- nchar(q)
            nams <- nam; lnams <- lnam[i]
            for(j in 3:length(prm))		{
               z <- prm[j]
               v <- paste(cof[j-1], z, sep="")
               u <- lprm[j]+lnam[i]+3
               n <- lcof[j-1]+u
               nams <- c(nams, z); lnams <- c(lnams,u)
               if(m + n > columns)	{
                  cur <- paste(cur, q)
                  tex <- c(tex, cur)
                  cur <- ""
                  cols <- 0
                  q <- ""
                  m <- 0		}
               q <- paste(q, v, sep="")
               m <- m + n			}
             Nam[[i]] <- nams; lNam[[i]] <- lnams },
         #Strat factor doesn't exist as main effect, but keep variable
         #names and their lengths if they will appear in interactions later
         {  if(length(Nam[[i]])==0 && any(interaction==i))	{
            nam.coef[[i]] <- paste(name[i], "=", oprm[-1], sep="")
            Nam[[i]] <- prm[-1]; lNam[[i]] <- lprm[-1]
								}
	    q <- "" },

	 {  if(prm[3,1]==0) 
              q <- Two.Way(prm,Nam,nam.coef,lNam,cof,coef,f,columns,lcof,
			varnames,lnam,at,digits=digits)
              else q <- Three.Way(prm,Nam,nam.coef,lNam,cof,coef,f,columns,lcof,at)
            m <- attr(q, "columns.used")
            j <- length(q)
            if(cur!="") {tex <- c(tex,cur); cur <- ""; cols <- 0}
            if(j>1) {tex <- c(tex,q[-j]); q <- q[j]}   }, 
         {  nam <- names(coef)
            if(cols>0) {tex <- c(tex,cur); cur <- ""; cols <- 0}
            q <- ""
            m <- 0
            lnam <- nchar(nam)
            nam <- paste("{\\rm ", nam, "}", sep="")
            Nam[[i]] <- nam; lNam[[i]] <- lnam
            for(j in 1:length(prm))		{
               v <- paste(cof.sp[j], nam[j], sep="")
               n <- lcof[j]+lnam[j]
               if(m + n > columns)	{
                  cur <- paste(cur, q)
                  tex <- c(tex, cur)
                  cur <- ""
                  cols <- 0
                  q <- ""
                  m <- 0		}
               q <- paste(q, v, sep="")
               m <- m + n			}}) 

     
     if(length(q) && q!="")			{
       if(cols+m > columns)	{
         tex <- c(tex, cur)
         cur <- ""
         cols <- 0		}

      cur <- paste(cur, q)
      cols <- cols + m				}
								}	#}

if(cur!="") tex <- c(tex, cur)

if(inline)	{
  cat(tex, sep="\n", file=file, append=append)
#  return(structure(fi, class=c("latex","file")))   17Jul01
  return(structure(list(file=file,style=NULL), class='latex'))
}

tex <- c(tex,"\\end{eqnarray*}")
tex <- ifelse(substring(tex,1,1)=="\\",tex,paste(before,tex,"\\\\"))

if(anybrace | anyplus)	{
  s <- if(length(which)==p) "and $" else "where $"
  if(anybrace) s <- paste(s,"\\{c\\}=1 {\\rm\\ if\\ subject\\ is\\ in\\ group\\ } c, \\ 0 {\\rm\\ otherwise}")
  if(anybrace & anyplus) s <- paste(s, ";\\ ")  # ; was , 2Mar01
  if(anyplus) s <- paste(s, 
                    "(x)_{+}=x {\\rm\\ if\\ } x>0, \\ 0 {\\rm\\ otherwise}")
  s <- paste(s, "$.")
#  if(anytr) s <- paste(s, "\\\\")
  tex <- c(tex, s)	}

if(anytr & pretrans)		{
  i <- TLi!=""
  if(sum(i)==1) tr <- paste("$",varnames[i],
	"$ is pre--transformed as $",TLi[i],"$.",sep="")
  else		{
    tr <- c("\\vspace{0.5ex}\\begin{center}{\\bf Pre--Transformations}\\\\",
     "\\vspace{1.5ex}\\begin{tabular}{|l|l|} \\hline",
     "\\multicolumn{1}{|c|}{Variable} & \\multicolumn{1}{c|}{Transformation} \\\\ \\hline",
     paste("$",varnames[i],"$ & $",TLi[i],"$ \\\\",sep=""),
     "\\hline", "\\end{tabular}\\end{center}")
		}
  tex <- c(tex, tr) 		}

cat(tex, sep="\n", file=file, append=append)
#if(.SV4.) new("latexFile", file=file, style=NULL) else
 structure(list(file=file, style=NULL),class='latex')
}

latex.cph <- function(object, title, 
		file=paste(first.word(deparse(substitute(object))),".tex",sep=""),
		append=FALSE, surv=TRUE, maxt=FALSE, which, varnames, 
		columns=65, inline=FALSE, 
		before=if(inline)"" else "& &", dec=3, pretrans=TRUE,
		caption, ...) {

  f <- object
  
  atr <- f$Design
  if(!length(atr)) atr <- getOldDesign(f)  ## 30may02
  
lev <- names(f$freq)
Intercept <- -f$center
strata <- f$strata
w <- if(length(caption)) paste('\\begin{center} \\bf',caption,'\\end{center}')
if(missing(which) & !inline)			{
if(length(strata)==0)	{
  w <- c(w,paste("\\[{\\rm Prob}\\{T\\geq t\\} = S_{0}(t)^{{\\textstyle e}^{X\\beta}}, {\\rm \\ \\ where} \\\\ \\]",sep=""))
			}
else			{
  sname <- atr$name[atr$assume.code==8]
  strata.sub <- letters[8+(1:length(sname))]
  s <- paste("{\\rm ",sname,"}=",strata.sub,sep="")
  s <- paste(s, collapse=",")
  w <- c(w,paste("\\[{\\rm Prob}\\{T\\geq t\\ |\\ ",s,"\\}=S_{",
	paste(strata.sub,collapse=""),
	"}(t)^{{\\textstyle e}^{X\\beta}}, {\\rm \\ \\ where} \\\\ \\]", sep=""))
			}
						}
if(missing(which)) which <- 1:length(atr$name)
if(missing(varnames)) varnames <- atr$name[atr$assume.code!=9]
cat(w, sep=if(length(w))"\n" else "", file=file, append=append)
z <- latexDesign(f, file=file, append=TRUE, which=which, varnames=varnames, 
                 columns=columns, 
                 before=before,
                 prefix=if(missing(which))"X\\hat{\\beta}" else NULL, 
                 intercept=Intercept, inline=inline, pretrans=pretrans)

if(inline) return(z)   ## 4Dec00

ss <- f$surv.summary
if(surv && length(ss))			{
  fs <- f$strata
  nstrat <- 0; if(length(fs)) nstrat <- length(fs)
  times <- as.numeric(dimnames(ss)[[1]])
  maxtime <- f$maxtime
  if(max(times)>=maxtime) maxt <- FALSE
  if(nstrat==0)		{
    s <- matrix(ss[,,1],ncol=1)
    
    if(maxt) 	{
      s <- cbind(s, f$surv[L <- length(f$surv)])
      times <- c(times, f$time[L]) 
		}
    dimnames(s) <- list(format(times), "$S_{0}(t)$")
    latex.default(s, file=file, append=TRUE, rowlabel="$t$", rowlabel.just="r",
		dec=dec, table.env=FALSE)
			}

  else					{
#    com <- paste(paste("-e 's/",sname,"=\\(.*\\),/",
#	"\\1, /' ",sep=""),collapse="")
#   # Adding \\: spacer in subscript caused LaTeX to barf  (19May95)
#    n <- sys(paste('sed -e "s/[.]/, /g"',com, # "-e 's/ /\\\\\\\\:/g'"
#		), 
#		paste(fs,",",sep=""))
    # Change . to ,blank
    n <- sedit(paste(fs,',',sep=''), '.', ', ')
    # Change sname=*, to *,
    n <- sedit(n, paste(sname,'=*,',sep=''), rep('*, ', length(sname)))
    n <- substring(n, 1, nchar(n)-sum(atr$assume.code==8)-1)
   #was -3*sum()-1   19May95
    s <- ss[,,1]
    if(maxt)		{
      smax <- rep(NA,nstrat)
      for(i in 1:nstrat) 	{
        smax[i] <- f$surv[[i]][abs(f$time[[i]]-maxtime)<.001]
				}
      s <- rbind(s, smax)
      times <- c(times, maxtime)
			}    
    
    dimnames(s) <- list(format(times),
			paste("$S_{", n, "}(t)$", sep=""))
    latex.default(s, file=file, append=TRUE, rowlabel="$t$", rowlabel.just="r",
		dec=dec, table.env=FALSE)
  }
}
z
}


latex.lrm <- function(object, title, 
   file=paste(first.word(deparse(substitute(object))),".tex",sep=""),
   append=FALSE, which, varnames, columns=65, inline=FALSE, 
   before=if(inline)"" else "& &",pretrans=TRUE,caption=NULL,...) {

  f <- object
  
if(missing(which) & !inline)				{
Y <- paste("{\\rm ",as.character(attr(f$terms,"formula")[2]),"}",sep="")
lev <- names(f$freq)
nrp <- f$non.slopes

w <- '\\['

j <- if(lev[2]=="TRUE") "" else paste("=",lev[2],sep="")
if(nrp==1) w <- paste(w,"{\\rm Prob}\\{",Y, j,
   "\\} = \\frac{1}{1+\\exp(-X\\beta)}", sep="")

else w <- paste(w,"{\\rm Prob}\\{", Y, 
   "\\geq j\\} = \\frac{1}{1+\\exp(-\\alpha_{j}-X\\beta)}", sep="")

w <- paste(w, ", {\\rm \\ \\ where} \\\\ \\]", sep="")

if(length(caption)) w <- c(paste('\\begin{center} \\bf',caption,
                                 '\\end{center}'), w)

if(nrp>1)			{
  w <- c(w,"\\begin{eqnarray*}")
  cof <- format(f$coef[1:nrp])
  for(i in 1:nrp)
    w <- c(w, paste("\\hat{\\alpha}_{\\rm ",lev[i+1],"} &=&",cof[i],"\\\\",sep=""))
  w <- c(w,"\\end{eqnarray*}",sep="")
				}
							}
else w <- NULL
if(missing(which) | missing(varnames)) { # 17Jul01
  at <- f$Design
  if(!length(at)) at <- getOldDesign(f)
}
if(missing(which)) which <- 1:length(at$name)
if(missing(varnames)) varnames <- at$name[at$assume.code!=9]
cat(w, file=file, append=append, sep=if(length(w))"\n" else "")
latexDesign(f, file=file, append=TRUE, which=which, varnames=varnames, 
	columns=columns, 
	before=before, prefix="X\\hat{\\beta}", inline=inline,
	pretrans=pretrans)  ## 4Dec00
}


latex.ols <- function(object, title,
   file=paste(first.word(deparse(substitute(object))),".tex",sep=""),
   append=FALSE, which, varnames, columns=65, inline=FALSE, 
   before=if(inline)"" else "& &", pretrans=TRUE, caption=NULL, ...)	{

  f <- object
  
w <- if(length(caption)) paste('\\begin{center} \\bf',caption,'\\end{center}')

if(missing(which) & !inline)			{
Y <- paste("{\\rm ",as.character(attr(f$terms,"formula")[2]),"}",sep="")

w <- c(w, paste("\\[{\\rm E(",Y,
           "}) = X\\beta, {\\rm \\ \\ where} \\\\ \\]", sep=""))
						}
at <- f$Design
if(!length(at)) at <- getOldDesign(f)

if(missing(which)) which <- 1:length(at$name)

if(missing(varnames)) varnames <- at$name[at$assume.code!=9]
cat(w, file=file, sep=if(length(w)) "\n" else "", append=append)
latexDesign(f, file=file, append=TRUE, which=which, varnames=varnames, 
	columns=columns, 
	before=before, prefix="X\\hat{\\beta}", inline=inline, 
	pretrans=pretrans)  ## 4Dec00
}


latex.pphsm <- function(object, title,
    file=paste(first.word(deparse(substitute(object))),".tex",sep=""),
    append=FALSE, which, varnames, 
    columns=65, inline=FALSE, 
    before=if(inline)"" else "& &",pretrans=TRUE, caption=NULL, ...) {

w <- if(length(caption)) paste('\\begin{center} \\bf',caption,'\\end{center}')

sc <- exp(f$parms)
at <- f$Design
if(!length(at)) at <- getOldDesign(f)

if(missing(which) & !inline)			{
  dist <- paste("\\exp\\{-t^{",format(1/sc),"} \\exp(X\\hat{\\beta})\\}")
  w <- c(w,paste("\\[{\\rm Prob}\\{T\\geq t\\} = ",dist,
	"{\\rm \\ \\ where} \\\\ \\]",sep=""))
						}				
if(missing(which)) which <- 1:length(at$name)
if(missing(varnames)) varnames <- at$name[at$assume.code!=9]
cat(w, file=file, sep=if(length(w))"\n" else "", append=append)
latexDesign(f, file=file, append=TRUE, which=which, varnames=varnames, 
	columns=columns, 
	before=before, prefix=if(missing(which))"X\\hat{\\beta}" else NULL, 
	inline=inline,pretrans=pretrans)  ## 4Dec00
}


latex.psm <- function(object,  title,
   file=paste(first.word(deparse(substitute(object))),".tex",sep=""),
   append=FALSE, which, varnames, 
   columns=65, inline=FALSE, 
   before=if(inline)"" else "& &",pretrans=TRUE, caption=NULL, ...) {

  f <- object
  
w <- if(length(caption)) paste('\\begin{center} \\bf',caption,'\\end{center}')

if(missing(which) & !inline)			{
  if(.SV4. || .R.) {
    dist <- f$dist
    w <- c(w, paste("\\[{\\rm Prob}\\{T\\geq t\\} = ",
	survreg.auxinfo[[dist]]$latex(f$scale),
	"{\\rm \\ \\ where} \\\\ \\]",sep=""))
  } else {
  fam <- f$family[1:2]
  dist <- fam[1]
  transform <- fam[2]
  w <- c(w,paste("\\[{\\rm Prob}\\{T\\geq t\\} = ",
	survreg.auxinfo[[dist]]$latex(f$parms, transform),
	"{\\rm \\ \\ where} \\\\ \\]",sep=""))
}
}
atr <- f$Design
if(!length(atr)) atr <- getOldDesign(f)

if(missing(which)) which <- 1:length(atr$name)
if(missing(varnames)) varnames <- atr$name[atr$assume.code!=9]

cat(w, sep=if(length(w)) "\n" else "", file=file, append=append)
latexDesign(f, file=file, append=TRUE, which=which,
            varnames=varnames, columns=columns, 
            before=before,
            prefix=if(missing(which))"X\\hat{\\beta}" else NULL, 
            inline=inline,pretrans=pretrans)  ## 4Dec00
}


if(FALSE) latex <-  function(x, ...)   # duplicates what's in print.display
{
  if(is.null(oldClass(x)))
    oldClass(x) <- data.class(x)
  UseMethod("latex")
}
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
               penmat, weights)
## 17jan03
##      .Fortran("lrmfit",coef=initial,as.integer(0),0,x,y,offset,
##               u=double(kint),
##               double(kint*(kint+1)/2),loglik=double(1),n,as.integer(0),
##               sumw,kint,
##               v=double(kint*kint),double(kint),double(kint),
##               double(kint),pivot=integer(kint),opts=opts,ftable,
##               penmat)
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
               pivot=integer(nvi),opts=opts,ftable,penmat,weights)
    ## 17jan03
##      .Fortran("lrmfit",coef=initial,nxin,est,x,y,offset,
##               u=double(nvi),
##               double(nvi*(nvi+1)/2),loglik=double(1),n,nx,sumw,nvi,
##               v=double(nvi*nvi),double(nvi),double(2*nvi),double(nvi),
##               pivot=integer(nvi),opts=opts,ftable,penmat)   #2*nvi 28Jul95
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
           double(nvi),integer(nvi),opts=opts,ftable,penmat,weights)
  ## 17jan03
##    .Fortran("lrmfit",coef=initial,nx,1:nx,x,y,offset,
##           u=double(nvi),double(nvi*(nvi+1)),double(1),n,nx,
##           sumw,nvi,v=double(nvi*nvi),double(nvi),double(nvi),
##           double(nvi),integer(nvi),opts=opts,ftable,penmat)
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
           penmat, weights)
  ## 17jan03
##    .Fortran("lrmfit",coef=initial,as.integer(0),0,x,y,offset,
##           u=double(kint),
##           double(kint*(kint+1)/2),loglik=double(1),n,as.integer(0),
##           numy,kint,
##           v=double(kint*kint),double(kint),double(kint),
##           double(kint),pivot=integer(kint),opts=opts,ftable,
##           penmat)
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
               penmat)
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
  m$formula <- formula  ## 16Dec97
  if(.R.) m$drop.unused.levels <- TRUE  ## 31jul02
  m[[1]] <- as.name("model.frame")
  nact <- NULL
  tform <- terms(formula, specials='strat')  ##specials= 16Dec97
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
    X <- Design(eval(m, sys.parent()))   # 24Apr01
    atrx <- attributes(X)
    nact <- atrx$na.action
    if(method=="model.frame") return(X)
    Terms <- atrx$terms             # 16Dec97
    attr(Terms, "formula") <- formula
    atr <- atrx$Design              # 24Apr01
    ## 22nov02
    if(length(nact$nmiss)) {
      jia <- grep('*%ia%*',names(nact$nmiss))  ## 8feb03
      if(length(jia)) nact$nmiss <- nact$nmiss[-jia]
      s <- if(length(offs)) names(nact$nmiss) !=  atrx$names[offs] else TRUE
      names(nact$nmiss)[s] <-
        c(as.character(formula[2]), atr$name[atr$assume.code!=9])
    }
    ## added [s] 22nov02

    Y <- model.extract(X, response)
    weights <- wt <- model.extract(X, 'weights')  ##6jun02 18jan03
    if(length(weights))
      warning('currently weights are ignored in model validation and bootstrapping lrm fits')
    ##   offset <- attr(X,"offset")
    ofpres <- length(offs) > 0
    if(ofpres) offs <- X[[offs]]  else offs <- 0
    if(!.R.)storage.mode(offs) <- "single"
    if(model) m <- X
    stra <- attr(tform,'specials')$strat ## 16Dec97 ----
    Strata <- NULL
    Terms.ns <- Terms
    if(length(stra)) {
      temp <- untangle.specials(Terms.ns, 'strat', 1)
      Terms.ns <- Terms.ns[-temp$terms]
      attr(Terms,   "factors") <- pmin(attr(Terms,"factors"),1)
      attr(Terms.ns,"factors") <- pmin(attr(Terms.ns,"factors"),1)
      Strata <- X[[stra]]
      nstrata <- length(levels(Strata))
    }                                    ## 16Dec97 ----
    X <- model.matrix(Terms.ns, X)  ## 8Apr02
    ##      ass <- attr(X, 'assign')   ## 9Apr02
    ##	  X <- model.matrix(Terms.ns, X)[,-1,drop=FALSE]  ## .ns 16Dec97
    X <- X[,-1,drop=FALSE]   ## 8Apr02  R drops assign with []
    if(!.R.)storage.mode(X) <- "single"
    dimnames(X)[[2]] <- atr$colnames
    ##	  oldClass(X) <- c("model.matrix", "Design")   14Sep00
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
	if(var.penalty=='sandwich') f$var.from.info.matrix <- v  # 25Mar00
	f.nopenalty <- if(ofpres) fitter(X,Y,offset=offs,initial=f$coef, maxit=1, tol=tol) else
	fitter(X,Y,initial=f$coef, maxit=1, tol=tol)
	##  info.matrix.unpenalized <- solvet(f.nopenalty$var, tol=tol)  ##6May96
	info.matrix.unpenalized <- f.nopenalty$info.matrix
	dag <- diag(info.matrix.unpenalized %*% v)
	f$effective.df.diagonal <- dag
	f$var <- if(var.penalty=='simple')v else
             v %*% info.matrix.unpenalized %*% v   # 25Mar00
	df <- sum(dag[-(1:nrp)])   ## 6May96
	lr <- f.nopenalty$stats["Model L.R."]
	pval <- 1-pchisq(lr, df)
	f$stats[c('d.f.','Model L.R.','P')] <- c(df, lr, pval)  
}
  }
    ass <- if(xpres) DesignAssign(atr, nrp, Terms) else list()
    ## 8Apr02
    ## ass[[1]] <- 1:nrp
    ## if(length(ass) > 1) for(i in 2:length(ass)) ass[[i]] <- ass[[i]]+nrp
    ##[,-1 after model.matrix had subtract intercept position, but was
	## only 1 d.f.

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
  f <- c(f, list(call=call, Design=if(xpres)atr,  ## 4Apr02
                 scale.pred=c("log odds","Odds Ratio"),
				 terms=Terms, assign=ass, na.action=nact, fail=FALSE,
                 nstrata=nstrata, fitFunction=c('lrm','glm')))
  
  oldClass(f) <- if(.SV4.)'Design' else c("lrm","Design","glm") ##13Nov00
  f
}
#Uses matinv Fortran function, which uses ginv and sweep
#Returns matrix inverse with attributes rank (integer rank of x)
# and swept (logical - whether or not ith variable has been swept)
#Input matrix should set its swept attribute before the first invocation
# of matinv for that matrix.  If swept isn't set, it defaults to all F.
#
#Inverse is with respect to diagonal elements which[1],which[2],...
#For collinearities, the appropriate rows and columns of inv are set to 0
#Caller must negate matrix when finished with all partial inversions if
# negate is false.  The default is to automatically negate the which
# portion of the inverse, i.e., to assume that no further operations are
# to be done on the matrix 
#
#Eps is singularity criterion, like 1-Rsquare
#
#F. Harrell 1 Aug 90

matinv <- function(a,which,negate=TRUE,eps=1E-12)			{

	swept <- attr(a,"swept")
	if(!is.matrix(a)) a <- as.matrix(a)
	storage.mode(a) <- "double"
	m<-nrow(a)
	if(missing(which))which <- 1:m	else{
		rw <- range(which)
		if(rw[1]<1 | rw[2]>m)stop("illegal elements to invert") }
	storage.mode(which) <- "integer"
	if(!length(swept))swept <- rep(FALSE, m)
	if(m!=ncol(a))stop("matrix must be square")

#	library.dynam(section="local",file="mlmats.o") 
	y <- if(.R.)
      .Fortran("matinv",x = a, as.integer(m), 
               as.integer(length(which)),which,
               swept=swept, logical(m), double(m*(m+1)/2), 
               double(m), rank = integer(1), as.double(eps),
               as.logical(negate), PACKAGE="Design") else
    .Fortran("matinv",x = a, as.integer(m), 
             as.integer(length(which)),which,
             swept=swept, logical(m), double(m*(m+1)/2), 
             double(m), rank = integer(1), as.double(eps),
             as.logical(negate))
	x <- y$x
	attr(x,"rank") <- y$rank
	attr(x,"swept") <- y$swept
	dimnames(x) <- dimnames(a)
	x
 
}
nomogram <- function(fit, ...) UseMethod("nomogram")

nomogram.Design <- function(fit, ..., adj.to, 
			lp=TRUE, lp.at, lplabel="Linear Predictor",
			fun, fun.at, fun.lp.at, funlabel="Predicted Value",
			fun.side,
			interact=NULL, intercept=1,
			conf.int=FALSE, 
			col.conf=c(1,if(under.unix).3 else 12),
			conf.space=c(.08,.2), 
			conf.lp=c("representative", "all", "none"),
			est.all=TRUE, abbrev=FALSE, minlength=4, maxscale=100, nint=10, 
			label.every=1, force.label=FALSE, 
			xfrac=.35, cex.axis=.85, cex.var=1,
			col.grid=FALSE, vnames=c("labels","names"),
			varname.label=TRUE, varname.label.sep="=", ia.space=.7, 
			tck=-.009, lmgp=.4, omit=NULL, naxes,
			points.label='Points', total.points.label='Total Points',
			total.sep.page=FALSE, total.fun, verbose=FALSE) {	

conf.lp <- match.arg(conf.lp)
vnames  <- match.arg(vnames)

abb <- (is.logical(abbrev) && abbrev) || is.character(abbrev)
if(is.logical(conf.int) && conf.int) conf.int <- c(.7,.9)
if(!is.logical(conf.int) && (length(conf.int)!=length(col.conf)))
  stop("conf.int and col.conf must have same length")
if(is.logical(col.grid) && col.grid) col.grid <- if(under.unix).2 else 5

oldpar <- oPar()  # in Hmisc Misc.s
mgp <- oldpar$mgp
mar <- oldpar$mar
par(mgp=c(mgp[1],lmgp,mgp[3]),mar=c(mar[1],1.1,mar[3],mar[4]))
on.exit(setParNro(oldpar))  ## was par(oldpar) 11Apr02
tck2 <- tck/2

se <- FALSE
if(any(conf.int>0)) {
  se <- TRUE
  zcrit <- qnorm((conf.int+1)/2)
  bar <- function(x, y, zcrit, se, col.conf, nlev=4) {
    y <- rep(seq(y[1], y[2], length=nlev), length.out=length(x))
    for(j in 1:length(x)) {
      xj <- x[j]; yj <- y[j]
      W <- c(0,zcrit)*se[j]
      for(i in 1:length(zcrit)) {
        segments(xj-W[i+1], yj, xj-W[i], yj, col=col.conf[i], lwd=1)
        segments(xj+W[i+1], yj, xj+W[i], yj, col=col.conf[i], lwd=1)
      }
    }
  }
}

nfun <- if(missing(fun)) 0 else if(is.list(fun)) length(fun) else 1
if(nfun>1 && length(funlabel)==1) funlabel <- rep(funlabel, nfun)
if(nfun>0 && is.list(fun) && length(names(fun))) funlabel <- names(fun)

if(!missing(fun.at)) {
  if(!is.list(fun.at)) fun.at <- rep(list(fun.at),nfun)
}
if(!missing(fun.lp.at)) {
  if(!is.list(fun.lp.at)) fun.lp.at <- rep(list(fun.lp.at),nfun)
}
if(!missing(fun.side)) {
  if(!is.list(fun.side)) fun.side <- rep(list(fun.side),nfun)
  if(any(!(unlist(fun.side) %in% c(1,3))))
    stop('fun.side must contain only the numbers 1 and 3')
}
at <- fit$Design
if(!length(at)) at <- getOldDesign(fit)
assume <- at$assume.code
if(any(assume==10))
 warning("does not currently work with matrix factors in model")
name  <- at$name
names(assume) <- name
parms <- at$parms
label <- if(vnames=="labels") at$label else name
if(any(d <- duplicated(name))) stop(paste("duplicated variable names:",
	paste(name[d],collapse=" ")))
label <- name
if(vnames=="labels") {
  label <- at$label
  if(any(d <- duplicated(label))) stop(paste("duplicated variable labels:",
	paste(label[d],collapse=" ")))
}

ia    <- at$interactions

factors <- list(...)
nf <- length(factors)

which <- if(est.all) (1:length(assume))[assume!=8] else 
	(1:length(assume))[assume!=8 & assume!=9]
if(nf>0) {
  jw <- charmatch(names(factors),name,0)
  if(any(jw==0))stop(paste("factor name(s) not in the design:",
		paste(names(factors)[jw==0],collapse=" ")))
  if(!est.all) which <- jw
}

Limval <- Getlim(at, allow.null=TRUE, need.all=FALSE)
values <- Limval$values
lims <- Limval$limits[c(6,2,7),,drop=FALSE]

#Next 4 lines 27Nov99 - keep character variables intact
lims <- oldUnclass(lims)
for(i in 1:length(lims))
   if(is.factor(lims[[i]]))lims[[i]] <- as.character(lims[[i]])
attr(lims, 'class') <- 'data.frame'  # so can subscript later

#Find underlying categorical variables
ucat <- rep(FALSE, length(assume))
names(ucat) <- name
for(i in (1:length(assume))[assume!=5 & assume<8]) {
   ucat[i] <- !is.null(V <- values[[name[i]]]) # did add && is.character(V)
   if(ucat[i]) parms[[name[i]]] <- V
}

discrete <- assume==5 | assume==8 | ucat
names(discrete) <- name

#Number of non-slopes:
nrp <- num.intercepts(fit)
Intercept <- if(nrp>0) fit$coefficients[intercept] else 
             if(!is.null(fit$center)) -fit$center else 0
#non.slopes <- rep(0, nrp)   23Jun95
#non.slopes[intercept] <- 1

intercept.offset <- if(nrp<2) 0 else
   fit$coefficients[intercept]-fit$coefficients[1]
#linear.predictors stored in fit always used first intercept

settings <- list()
for(i in which[assume[which]<9]) {
  ni <- name[i]
  z <- factors[[ni]]
  lz <- length(z)
  if(lz < 2) settings[[ni]] <- value.chk(at, i, NA, -nint, Limval, 
				type.range="full") else
  if(lz > 0 && any(is.na(z))) stop("may not specify NA as a variable value")
  if(lz==1) lims[2,i] <- z else if(lz>1) {
    settings[[ni]] <- z
    if(is.null(lims[[ni]]) || is.na(lims[2,ni])) {
      lims[[ni]] <- c(NA,z[1],NA)
      warning(paste("adjustment values for ",ni,
        " not defined in datadist; taken to be first value specified (", 
        z[1],")" ,sep=""))
    }
  }
}

adj <- lims[2,,drop=FALSE]
if(!missing(adj.to)) for(nn in names(adj.to)) adj[[nn]] <- adj.to[[nn]]
isna <- sapply(adj, is.na)
if(any(isna)) stop(
   paste("adjustment values not defined here or with datadist for",
	paste(name[assume!=9][isna],collapse=" ")))

num.lines <- 0
space.used <- 0
entities <- 0
set <- list()
iset <- 0
start <- len <- NULL
end <- 0

#Sort to do continuous factors first if any interactions present
main.effects <- which[assume[which]<8]  # this logic not handle strata w/intera.
if(any(assume==9)) main.effects <- main.effects[order(10*discrete[main.effects]+
	(name[main.effects] %in% names(interact)))]

#For each predictor, get vector of predictor numbers directly or
#indirectly associated with it
rel <- related.predictors(at)   # Function in Design.Misc.s

already.done <- structure(rep(FALSE,length(name)), names=name)
for(i in main.effects) {
  nam <- name[i]
  if(already.done[nam] || (nam %in% omit)) next
  r <- sort(rel[[i]])
  if(length(r)==0) { #main effect not contained in any interactions
    num.lines <- num.lines + 1
    space.used <- space.used + 1
    entities <- entities+1
    x <- list()
    x[[nam]] <- settings[[nam]]
    iset <- iset+1
    attr(x,'info') <- list(predictor=nam, effect.name=nam,
                           type='main')
    set[[label[i]]] <- x
#    set[[label[i]]] <- structure(x, predictor=nam, effect.name=nam,
#                                 type="main")  24Sep00
    start <- c(start, end+1)
    n <- length(settings[[nam]])
    len <- c(len, n)
    end <- end+n
    NULL	#handles bug in S
  }  else {
    namo <- name[r]
    s <- !(name[r] %in% names(interact))
    if(any(s)) {
      if(is.null(interact)) interact <- list()
      for(j in r[s]) {
        nj <- name[j]
        if(discrete[j]) interact[[nj]] <-  parms[[nj]]
	NULL
      }
      s <- !(name[r] %in% names(interact))
    }
    if(any(s)) stop(paste("factors not defined in interact=list(...):",
                          paste(name[r[s]],collapse=",")))
    combo <- expand.grid(interact[namo]) #list[vector] gets sublist
    oldClass(combo) <- NULL   # so combo[[n]] <- as.character will really work
    acombo <- combo
    if(abb) for(n in if(is.character(abbrev))abbrev else names(acombo)) {
      if(discrete[n]) { 
       acombo[[n]] <- abbreviate(parms[[n]],
         minlength=if(minlength==1)4 else minlength)[combo[[n]]]  
       #lucky that abbreviate function names its result
      }
    }
    for(n in names(combo)) if(is.factor(combo[[n]])) {
      combo[[n]] <- as.character(combo[[n]])   #so row insertion will work xadj
      acombo[[n]] <- as.character(acombo[[n]]) #so format() will work
      NULL
    }
    entities <- entities+1
    already.done[namo] <- TRUE
    for(k in 1:length(combo[[1]])) {		# was nrow(combo))
      num.lines <- num.lines+1
      space.used <- space.used + if(k==1) 1 else ia.space
      x <- list()
      x[[nam]] <- settings[[nam]]   #store fastest first
      for(nm in namo) x[[nm]] <- combo[[nm]][k]   #was combo[k,nm] 2Dec94
      iset <- iset+1
      set.name <- paste(nam, " (",sep="")
      for(j in 1:length(acombo)) {		# was ncol
        set.name <- paste(set.name, 
          if(varname.label) paste(namo[j],varname.label.sep,sep="") else "",
          format(acombo[[j]][k]),sep="")   # was acombo[k,j]
        if(j<length(acombo)) set.name <- paste(set.name," ",sep="")
      }
      set.name <- paste(set.name,")",sep="")
     #Make list of all terms needing inclusion in calculation
     #Include interation term names  - interactions.containing in Design.Misc.s
      ia.names <- NULL
      for(j in r) ia.names <- c(ia.names, name[interactions.containing(at, j)])
      ia.names <- unique(ia.names)
      attr(x,'info') <- list(predictor=nam,
                             effect.name=c(nam,namo[assume[namo]!=8],ia.names), 
                             type=if(k==1) "first" else "continuation")
      set[[set.name]] <- x
#      set[[set.name]] <- list(x, predictor=nam,
#                          effect.name=c(nam,namo[assume[namo]!=8],ia.names), 
#                          type=if(k==1) "first" else "continuation") #24Sep00
			##Don't include strata main effects
      start <- c(start, end+1)
      n <- length(settings[[nam]])
      len <- c(len, n)
      end <- end+n
      NULL
    }
    NULL
  }    
}
# xadj <- Design.levels(adj, at)  # in Design.Misc.s  20jul03 + next
# xadj[1:sum(len),] <- adj        # for R 1.8
xadj <- oldUnclass(Design.levels(adj, at))
for(k in 1:length(xadj)) xadj[[k]] <- rep(xadj[[k]], sum(len))

j <- 0
for(S in set) {
  j <- j+1
  ns <- names(S)
  nam <- names(S)
  for(k in 1:length(nam)) {
    xadj[[nam[k]]][start[j]:(start[j]+len[j]-1)] <- S[[k]]
    NULL
  }
}
#xadj <- data.frame(xadj)  27Nov99 to preserve character variables
xadj <- structure(xadj, class='data.frame',
                  row.names=as.character(1:sum(len)))   # 27Nov99

xx <- predictDesign(fit, newdata=xadj, type="terms", center.terms=FALSE, se.fit=FALSE,
	kint=intercept)   #non.slopes=non.slopes) 23Jun95
##20Nov00:
if(any(is.infinite(xx)))
  stop("variable limits and transformations are such that an infinite axis value has resulted.\nRe-run specifying your own limits to variables.")

if(se) xse <- predictDesign(fit, newdata=xadj, se.fit=TRUE, 
	kint=intercept)   #non.slopes=non.slopes) 23Jun95

R <- matrix(NA, nrow=2, ncol=length(main.effects),
            dimnames=list(NULL,name[main.effects]))
R[1,] <- 1e30
R[2,] <- -1e30
# R <- apply(xx, 2, range)  - does not work since some effects are for
# variable combinations that were never used in constructing axes

for(i in 1:num.lines) {
  is <- start[i]; ie <- is+len[i]-1
  s <- set[[i]]
#  nam <- attr(s, "effect.name") 24Sep00
  setinfo <- attr(s,'info')
  nam <- setinfo$effect.name
  xt <- xx[is:ie, nam]
  if(length(nam)>1) xt <- apply(xt, 1, sum)  # add all terms involved
  set[[i]]$Xbeta <- xt
  r <- range(xt)
#  pname <- attr(s, "predictor") 24Sep00
  pname <- setinfo$predictor
  R[1,pname] <- min(R[1,pname], r[1])
  R[2,pname] <- max(R[2,pname], r[2])
  if(se) {
    set[[i]]$Xbeta.whole <- xse$linear.predictors[is:ie] #note-has right interc.
    set[[i]]$se.fit      <- xse$se.fit[is:ie]
    NULL
  }
  NULL
}


R <- R[,R[1,]<1e30,drop=FALSE]
sc <- maxscale/max(R[2,]-R[1,])
Intercept <- Intercept + sum(R[1,])

xl <- -xfrac*maxscale
if(missing(naxes)) naxes <- 
  if(total.sep.page) max(space.used + 1, nfun + lp + 1) else
                     space.used + 1 + nfun + lp + 1

Format <- function(x) { # like format but does individually
  f <- character(l <- length(x))
  for(i in 1:l) f[i] <- format(x[i])
  f
}

newpage <- function(naxes,xl,maxscale,cex.var,nint,space.used,col.grid,
		cex.axis,tck,tck2,label.every,force.label,
		points=TRUE, points.label='Points',usr) {
  y <- naxes-1
  plot(0,0,xlim=c(xl,maxscale),ylim=c(0,y),
	   type="n",axes=FALSE,xlab="",ylab="")
  if(!missing(usr)) par(usr=usr)
  if(!points) return(y + 1)

  ax <- c(0,maxscale)
  text(xl, y, points.label, adj=0, cex=cex.var)
  x <- pretty(ax, n=nint)
  x2 <- seq((x[1]+x[2])/2, max(x), by=x[2]-x[1])
  if(col.grid>0) {
	segments(x ,y,x, y-space.used,col=col.grid,lwd=1)
	segments(x2,y,x2,y-space.used,col=if(col.grid==1)1 else col.grid/2,lwd=1)
  }  
  axisf(3, at=x,  pos=y, cex=cex.axis, tck=tck, label.every=label.every, 
	   force.label=force.label)
  axisf(3, at=x2, labels=FALSE, pos=y, tck=tck2, cex=cex.axis)
  y
}

y <- newpage(naxes,xl,maxscale,cex.var,nint,space.used,col.grid,
			cex.axis,tck,tck2,label.every=label.every,
			force.label=force.label,points.label=points.label)

i <- 0
ns <- names(set)
Abbrev <- list()
for(S in set) {
  i <- i+1
#  type <- attr(S,"type")  24Sep00
  setinfo <- attr(S,'info')
  type <- setinfo$type
  y <- y - (if(type=="continuation") ia.space else 1)
  if(y < -.05) {
	y <- newpage(naxes,xl,maxscale,cex.var,nint,space.used,col.grid,
			cex.axis,tck,tck2,
			label.every=label.every,force.label=force.label,
			points.label=points.label) -
			  (if(type=="continuation") ia.space else 1)
  }

  text(xl, y, ns[[i]], adj=0, cex=cex.var)
  x <- S[[1]]
  nam <- names(S)[1]  #stored with fastest first
  fx <- if(is.character(x)) x else 
          sedit(Format(x)," ","") #axis not like bl   - was translate()
  if(abb && discrete[nam] && (is.logical(abbrev) || nam %in% abbrev)) {
    old.text <- fx
    fx <- if(abb && minlength==1)letters[1:length(fx)] else 
          abbreviate(fx, minlength=minlength)
    Abbrev[[nam]] <- list(abbrev=fx, full=old.text)
  }

  j <- match(nam, name, 0)
  if(any(j==0)) stop("program logic error 1")
#  if(!discrete[nam] && label.every>1) {
#    sq <- seq(along=x, by=label.every)
#    fx[-sq] <- ""}
  is <- start[i]; ie <- is+len[i]-1
  xt <- (S$Xbeta - R[1,nam])*sc
  set[[i]]$points <- xt
 #Find flat pieces and combine their labels
  r <- rle(xt)
  if(any(r$length>1)) {
    is <- 1
    for(j in r$length) {
      ie <- is+j-1
      if(j>1) {
        fx[ie] <- if(discrete[nam] || ie < length(xt))
		  paste(fx[is], "-", fx[ie],sep="") else
		paste(fx[is], '+', sep='')

        fx[is:(ie-1)] <- ""
	xt[is:(ie-1)] <- NA
      }
      is <- ie+1
    }
  fx <- fx[!is.na(xt)]
  xt <- xt[!is.na(xt)]
  }
 #Find direction changes
  ch <- if(length(xt)>2) c(FALSE,FALSE,diff(diff(xt)>0)!=0) else rep(FALSE, length(xt))
  if(discrete[nam] && length(xt)>1) { # categorical - alternate adjacent levels
    j <- order(xt)
    lines(range(xt),rep(y,2))   # make sure all ticks are connected
    for(k in 1:2) {
      is <- j[seq(k, length(j), by=2)]
	  ##If using 1-letter abbreviations, force all to print
	  ##      if(abb && minlength==1 && (is.logical(abbrev) || nam %in% abbrev)) 
	  ##        for(jj in is) axis(1+2*(k==2), at=xt[jj], label=fx[jj], pos=y,
	  ##                           cex=cex.axis, tck=tck) else
	  ##      axis(1+2*(k==2), at=xt[is], label=fx[is], pos=y, cex=cex.axis, tck=tck)
	  axisf(1+2*(k==2), at=xt[is], labels=fx[is], pos=y, cex=cex.axis,
		   tck=tck, force.label=force.label || (
		   abb && minlength==1 && (is.logical(abbrev) || nam %in% abbrev)),
		   disc=TRUE)
      if(se) bar(xt[is], if(k==1) y-conf.space-.32 else y+conf.space+.32, 
				 zcrit, sc*S$se.fit[is], col.conf)
    }
  } else if(!any(ch)) {
    axisf(1, at=xt, labels=fx, pos=y, cex=cex.axis, tck=tck,
		 label.every=label.every, force.label=force.label, disc=discrete[nam])
    if(se)bar(xt, y+conf.space, zcrit, sc*S$se.fit, col.conf)
  }
  else {
    lines(range(xt), rep(y,2))  # make sure all ticks are connected
    j <- (1:length(ch))[ch]
    if(max(j)<length(ch)) j <- c(j, length(ch)+1)
    side <- 1
    is <- 1
    for(k in j) {
      ie <- k-1
      axisf(side, at=xt[is:ie], labels=fx[is:ie], pos=y, cex=cex.axis, tck=tck,
		   label.every=label.every, force.label=force.label, 
		   disc=discrete[nam])
      if(se) bar(xt[is:ie], if(side==1) y-conf.space-.32 else y+conf.space+.32, 
                 zcrit, sc*S$se.fit[is:ie], col.conf)
      side <- if(side==3)1 else 3
      is <- ie+1
    }
  }
}

if(missing(lp.at)) {
  xb <- fit$linear.predictors
  if(is.null(xb)) xb <- fit$fitted
  if(is.null(xb)) stop("lp.at not given and fit did not store linear.predictors or fitted.values")
  if(nrp>1) xb <- xb + intercept.offset
  lp.at <- pretty(range(xb), n=nint)
}

sum.max <- if(entities==1) maxscale else max(maxscale,sc*max(lp.at-Intercept))
x <- pretty(c(0, sum.max), n=nint)

new.max <- max(x)
xl.old <- xl
xl <- -xfrac*new.max
u <- par()$usr
if(!missing(total.fun)) total.fun()
usr <- c(xl*u[1]/xl.old, new.max*u[2]/maxscale, u[3:4])
par(usr=usr)

x.double <- seq(x[1], new.max, by=(x[2]-x[1])/2)

y <- y-1

if(y < -.05 || total.sep.page) 
  y <- newpage(naxes,xl,maxscale,cex.var,nint,space.used,col.grid,
			cex.axis,tck,tck2,
			label.every=label.every,force.label=force.label,
			points=FALSE,usr=usr) - 1

text(xl, y, total.points.label, adj=0, cex=cex.var)
axisf(1, at=x, pos=y, cex=cex.axis, tck=tck, label.every=label.every,
	 force.label=force.label)
axisf(1, at=x.double, labels=FALSE, pos=y, tck=tck2, cex=cex.axis)
set$total.points <- list(x=x, y=y)

if(lp) {
  x2 <- seq(lp.at[1], max(lp.at), by=(lp.at[2]-lp.at[1])/2)
  scaled.x <- (lp.at-Intercept)*sc
  scaled.x2 <- (x2-Intercept)*sc
  y <- y-1
  if(y < -.05) 
	y <- newpage(naxes,xl,maxscale,cex.var,nint,space.used,col.grid,
			cex.axis,tck,tck2,
			label.every=label.every,force.label=force.label,
			points=FALSE,usr=usr) - 1

  text(xl, y, lplabel, adj=0, cex=cex.var)
  axisf(1, at=scaled.x,  labels=Format(lp.at), pos=y, cex=cex.axis, tck=tck,
	   label.every=label.every, force.label=force.label)
  axisf(1, at=scaled.x2, labels=FALSE, tck=tck2,   pos=y, cex=cex.axis)
  set$lp <- list(x=scaled.x, y=y, x.real=lp.at)
  if(se && conf.lp!="none") {
    xxb <- NULL
    xse <- NULL
    for(S in set) { xxb <- c(xxb, S$Xbeta.whole); xse <- c(xse, S$se.fit) }
    i <- order(xxb)
    if(length(xxb)<16 | conf.lp=="representative") 
	{nlev <- 4; w <- 1} else {nlev <- 8; w <- 2}
    if(conf.lp=="representative") {
      deciles <- cut2(xxb[i], g=10)
      mean.xxb <- tapply(xxb[i], deciles, mean)
      median.se <- tapply(xse[i], deciles, median)
      bar((mean.xxb-Intercept)*sc, 
	  y+c(conf.space[1],conf.space[1]+w*diff(conf.space)),
          zcrit, sc*median.se, col.conf, nlev=nlev)
    } else
    bar((xxb[i]-Intercept)*sc, y+c(conf.space[1],
	conf.space[1]+w*diff(conf.space)), 
        zcrit, sc*xse[i], col.conf, nlev=nlev)
  }
}

if(nfun>0) {
  if(!is.list(fun)) fun <- list(fun)
  i <- 0
  for(func in fun) {
    i <- i+1
   #Now get good approximation to inverse of fun evaluated at fat
   #Unless inverse function given explicitly
    if(!missing(fun.lp.at)) {
      xseq <- fun.lp.at[[i]]
      fat <- func(xseq)
      w <- xseq
    } else {
      if(missing(fun.at)) fat <- pretty(func(range(lp.at)), n=nint)
      else fat <- fun.at[[i]]
      if(verbose) {
        cat('Function',i,'values at which to place tick marks:\n')
        print(fat)
      }
      xseq <- seq(min(lp.at),max(lp.at),length=1000)
      fu <- func(xseq)
      s <- !is.na(fu)
      w <- approx(fu[s], xseq[s], fat)$y
      if(verbose) {
        cat('Estimated inverse function values (lp):\n')
        print(w)
      }
    }
    s <- !(is.na(w)|is.na(fat))
    w <- w[s]
    fat <- fat[s]
    fat.orig <- fat
    fat <- if(is.category(fat)) as.character(fat) else Format(fat)
    scaled <- (w-Intercept)*sc
    y <- y-1
	if(y < -.05) 
	  y <- newpage(naxes,xl,maxscale,cex.var,nint,space.used,col.grid,
			cex.axis,tck,tck2,
			label.every=label.every,force.label=force.label,
			points=FALSE,usr=usr) - 1

    text(xl, y, funlabel[i], adj=0, cex=cex.var)
    sides <- if(missing(fun.side)) rep(1, length(fat)) else (fun.side[[i]])[s]
    if(length(sides)!=length(fat)) 
      stop('fun.side vector not same length as fun.at or fun.lp.at')
    for(jj in 1:length(fat)) axis(sides[jj], at=scaled[jj], label=fat[jj],
                                  pos=y, cex=cex.axis, tck=tck)
    lines(range(scaled),rep(y,2))  #make sure all ticks are connected
    set[[funlabel[i]]] <- list(x=scaled, y=y, x.real=fat.orig)
  }
}
set$abbrev <- Abbrev
oldClass(set) <- "nomogram"
invisible(set)
}

print.nomogram <- function(x, dec=0, ...) {

  obj <- x
  
  w <- diff(range(obj$lp$x))/diff(range(obj$lp$x.real))
  cat('Points per unit of linear predictor:',format(w),
	  '\nLinear predictor units per point   :',format(1/w),'\n\n')

fun <- FALSE
for(x in names(obj)) {

  k <- x=='total.points' || x=='lp' || x=='abbrev'
  if(k) { fun <- TRUE; next }
  y <- obj[[x]]
  if(fun) {
    z <- cbind(round(y[[1]],dec), y$x.real)
    dimnames(z) <- list(rep('',nrow(z)), c('Total Points',x))
  } else {
    if(.R.) {
      z <- cbind(format(y[[1]]), format(round(y$points,dec)))
      dimnames(z)  <- list(rep('',length(y$points)), c(x, 'Points'))
    } else {
      z <- list(x=y[[1]], Points=round(y$points,dec))
      names(z) <- c(x,'Points')
      attr(z,'row.names') <- rep('',length(y$points))
      attr(z,'class') <- 'data.frame'
    }
    ## didn't use data.frame since wanted blank row names
  }
  cat('\n')
  if(.R.) print(z, quote=FALSE) else print(z)
  cat('\n')
}
invisible()
}

legend.nomabbrev <- function(object, which, x, y=NULL, ncol=3, ...) {
abb <- object$abbrev[[which]]
if(length(abb)==0) stop(paste('no abbreviation information for',which))
if(max(nchar(abb$abbrev))==1) if(length(y)) legend(x, y, abb$full, ncol=ncol,
                                     pch=paste(abb$abbrev,collapse=''), ...)
	else legend(x, abb$full, ncol=ncol, pch=paste(abb$abbrev,collapse=''),
		 ...)
else if(length(y)) legend(x, y, paste(format(abb$abbrev),':',abb$full,sep=''), 
		ncol=ncol, ...) else
	legend(x, paste(format(abb$abbrev),':',abb$full,sep=''), ncol=ncol, 
            ...)
invisible()
}


##Version of axis allowing tick mark labels to be forced, or to
##label every 'label.every' tick marks

axisf <- function(side, at, labels=TRUE, pos, cex, tck, 
				  label.every=1, force.label=FALSE, disc=FALSE) {

  if(.R.) {   # R doesn't pass cex to axis   16may03
    ocex <- par('cex')
    par(cex=cex)
    on.exit(par(cex=ocex))
  }
  
  axis(side, at, labels=FALSE, pos=pos, cex=cex, tck=tck)

  if(is.logical(labels) && !labels) return(invisible())

  if(label.every>1 && !disc) {
	sq <- seq(along=at, by=label.every)
	at[-sq] <- NA
  }
  if(is.logical(labels)) labels <- format(at)

  if(force.label) {
    for(i in 1:length(labels))
	if(!is.na(at[i])) axis(side, at[i], labels[i], pos=pos, cex=cex, tck=tck)
    }
  else axis(side, at[!is.na(at)], labels[!is.na(at)], 
            pos=pos, cex=cex, tck=tck)

  invisible()
}



if(FALSE) structure <- function(.Data = NULL, ...)
{
        s <- .Data
        a <- list(...)
        nn <- sort(names(a))
        #To allow sourcing of dumped time series
        for(i in nn) {
          prn(i)
                attr(s, i) <- a[[i]]
          prn(names(attributes(s)))
            }
        s
}                                                                               
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
    if(.R.) m$drop.unused.levels <- TRUE  ## 31jul02
	m[[1]] <- as.name("model.frame")
    ##	if(!missing(data) || any(attr(terms(formula),"term.labels")!=".")){
    ##X's present) 
	if(length(attr(terms(formula),"term.labels"))) {
      ## R's model.frame.default gives wrong model frame if [.factor
      ## removes unused factor levels  31jul02
      if(.R.) {
        dul <- .Options$drop.unused.levels
        if(!length(dul) || dul) {
          on.exit(options(drop.unused.levels=dul))
          options(drop.unused.levels=FALSE)
        }
      }
      X <- Design(eval(m, sys.parent()))   # 17Apr01 + next
      if(.R.) options(drop.unused.levels=dul)
    atrx <- attributes(X)
    atr <- atrx$Design
	nact <- atrx$na.action
	Terms <- atrx$terms
    assig <- DesignAssign(atr, 1, Terms)  ## 11Apr02
    
	penpres <- !(missing(penalty) && missing(penalty.matrix))
    if(penpres && missing(var.penalty))
      warning('default for var.penalty has changed to "simple"')

	if(method == "model.frame") return(X)
	scale <- as.character(formula[2])
	attr(Terms, "formula") <- formula
	if(length(nact$nmiss)) {
      jia <- grep('*%ia%*',names(nact$nmiss))  ## 8feb03
      if(length(jia)) nact$nmiss <- nact$nmiss[-jia]
      names(nact$nmiss) <- 
		c(scale,atr$name[atr$assume.code!=9])
    }
	weights <- model.extract(X, weights)
	if(length(weights) && penpres)
	  stop('may not specify penalty with weights')

	Y <- model.extract(X, response)
	n <- length(Y)
	if(model) m <- X
	X <- model.matrix(Terms, X)
	if(!.R.) storage.mode(X) <- "single"
	if(length(atr$colnames)) 
	   dimnames(X)[[2]] <- c("Intercept",atr$colnames)
	else dimnames(X)[[2]] <- c("Intercept",dimnames(X)[[2]][-1])
#	oldClass(X) <- c("model.matrix","Design")   14Sep00
	if(method=="model.matrix") return(X)				   }

#Model with no covariables:

else {
  assig <- NULL
  yy <- attr(terms(formula),"variables")[1]
  Y <- eval(yy,sys.parent(2))
  nmiss <- sum(is.na(Y))
  if(nmiss==0) nmiss <- NULL else names(nmiss) <- as.character(yy)
  Y <- Y[!is.na(Y)]
  ##Note: this logic doesn't work in presence of weights
  yest <- mean(Y)
  coef <- yest
  n <- length(Y)
  if(!length(sigma)) sigma <- sqrt(sum((Y-yest)^2)/(n-1))
  cov <- matrix(sigma*sigma/n, nrow=1, ncol=1,
                dimnames=list("Intercept","Intercept"))
  fit <- list(coefficients=coef, var=cov,
              non.slopes=1, fail=FALSE, residuals=Y-yest,
              df.residual=n-1, intercept=TRUE)
  if(linear.predictors) {fit$linear.predictors <- rep(yest,n); 
                         names(fit$linear.predictors) <- names(Y)}
  if(model) fit$model <- m
  if(x) fit$x <- matrix(1, ncol=1, nrow=n, 
                        dimnames=list(NULL,"Intercept"))
  if(y) fit$y <- Y
  fit$fitFunction <- c('ols','lm')  ## 13Nov00
  oldClass(fit) <- if(.SV4.)"Design" else c("ols","Design","lm")
  ## 13Nov00  17Apr01
  return(fit)
}

	if(!penpres) {
	  fit <- if(length(weights))
        lm.wfit(X, Y, weights, method=method,  ...) else 
	    lm.fit(X, Y, method=method, ...)  ## added method=x2 8Apr02
      if(.R.) cov.unscaled <- chol2inv(fit$qr$qr) else {
        rinv <- solve(fit$R, diag(length(fit$coefficients)))
        cov.unscaled <- rinv %*% t(rinv)
      }
	  sse <- sum(fit$residuals^2)
	  if(!length(sigma)) sigma <- sqrt(sse/fit$df.residual)
	  fit$var <- sigma*sigma*cov.unscaled
	  cnam <- dimnames(X)[[2]]
	  dimnames(fit$var) <- list(cnam, cnam)
	  r2 <- 1-sse/sum((Y-mean(Y))^2)
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
                     var.penalty=var.penalty)  #25Mar00
	  fit$penalty <- penalty
	}

	if(model)
		fit$model <- m
	if(linear.predictors) 	{
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
    fit$assign <- assig   ## 11Apr02
	oldClass(fit) <- if(.SV4.)'Design' else c("ols","Design","lm") ##13Nov00
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
if(regcoef.only) return(list(coefficients=coef))  ## 16Oct97
res  <- drop(Y - X %*% coef)
n <- length(Y)
sse <- sum(res^2)
s2 <- drop( (sse + t(coef) %*% pm %*% coef) / n )
var <- if(var.penalty=='simple') s2 * Z else s2 * Z %*% xpx %*% Z  #25Mar00
cnam <- dimnames(X)[[2]]
dimnames(var) <- list(cnam, cnam)
sst <- sum((Y-mean(Y))^2)
lr <- n*(1+logb(sst/n))-n*logb(s2)-sse/s2
s2.unpen <- sse/n
dag <- diag((xpx / s2.unpen) %*% (s2 * Z))
df <- sum(dag) - 1
stats <- c(n=n, 'Model L.R.'=lr, 'd.f.'=df, R2=1-sse/sst, Sigma=sqrt(s2))

## was assign=attr(X,'assign') 11Apr01
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
pentrace <- function(fit, penalty, penalty.matrix,
					 method=c('grid','optimize'),
					 which=c('aic.c','aic','bic'), target.df=NULL,
					 fitter, pr=FALSE,
                     tol=1e-7, keep.coef=FALSE, complex.more=TRUE, verbose=FALSE,
                     maxit=12, subset)  {

# Need to check Strata for cph

method <- match.arg(method)
which  <- match.arg(which)
tdf    <- length(target.df)
if(tdf) method <- 'optimize'

if(!length(X <- fit$x) | !length(Y <- as.matrix(fit$y)))
  stop("you did not specify x=T and y=T in the fit")
fit$x <- fit$y <- NULL

if(length(pn <- fit$penalty)>0 && max(unlist(pn))!=0) 
 warning('you should not have specified penalty= in fit so that unpenalized model can be a candidate for the best model')

sc.pres <- match("parms",names(fit),0)>0
## 20Apr02
if(!.newSurvival.) {
  fixed <- fit$fixed   #psm only
  fixed <- if(length(fixed)==1 && is.logical(fixed) && !fixed) list() else
           list(scale=TRUE)
}

if(.R.) {
  fixed <- NULL
  dist  <- fit$dist
  parms <- fit$parms
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
  storeTemp(NULL, 'Strata')
}

clas <- fit$fitFunction[1]  ## 14Nov00
if(!length(clas)) clas <- oldClass(fit)[1]

isols <- clas=='ols'

if(!(isols | clas=='lrm')) 
  stop("at present pentrace only works for lrm or ols")

if(missing(fitter)) fitter <- switch(clas,
				ols=function(x,y,maxit,...)lm.pfit(x,y,...), 
				lrm=function(x,y,maxit=12,...)
                                 lrm.fit(x,y,maxit=maxit,...), 
				cph=function(x,y,maxit=12,...)coxphFit(x,y,
				 strata=Strata, iter.max=maxit, 
				 eps=.0001, method="efron", toler.chol=tol),
				psm=function(x,y,maxit=12,...) survreg.fit2(x, y,
                        dist=dist, parms=parms, fixed=fixed, offset=NULL,
                        init=NULL, maxiter=maxit))

## psm= was (20Apr02) the following.  survreg.fit2 is in validate.psm.s
##psm=function(x,y,maxit=12,...) survreg.fit(x, y, 
## dist=dist, fixed=fixed, offset=NULL, init=NULL,
## controlvals=survreg.control(maxiter=maxit,rel.tol=.0001))
## survreg.fit2 ignores unneeded parms vs. fixed

if(!length(fitter))stop("fitter not valid")

str.pres <- FALSE
if(clas=="cph")	{  ##14Nov00 was inherits
   str <- attr(X, "strata")
   str.pres <- length(str) > 0
}

if(!missing(subset)) {
  Y <- Y[subset,,drop=FALSE]
  X <- X[subset,,drop=FALSE]
  if(str.pres) str <- str[subset,drop=FALSE]
}
n <- nrow(Y)
#if(str.pres) assign("Strata", str, 1)  # 28Jul98  17Apr01
if(str.pres) {  # 17Jul01
  if(.R.) Strata <- str else storeTemp(str, "Strata")
}


atr <- fit$Design
if(!length(atr)) atr <- getOldDesign(fit)

if(missing(penalty.matrix)) 
  penalty.matrix <- Penalty.matrix(atr, X)

obj.best <- -1e10

ns <- num.intercepts(fit)


islist <- is.list(penalty)

if(islist) {
  penalty <- expand.grid(penalty)
  if(complex.more && ncol(penalty) > 1 && nrow(penalty) > 1) {
    ikeep <- NULL
    for(i in 1:nrow(penalty)) {
      ok <- TRUE
      peni <- penalty[i,]
      for(j in 2:length(peni)) if(peni[[j]] < peni[[j-1]]) ok <- FALSE
      if(ok) ikeep <- c(ikeep, i)
    }
  penalty <- penalty[ikeep,,drop=FALSE]
  }
  np <- nrow(penalty)
} else {
  if(method=='grid') penalty <- c(0, penalty[penalty>0])
  np <- length(penalty)
}

if(method=='optimize') {
  if(.R.) stop('method="optimize" not yet implemented in R')
  
  if((islist && nrow(penalty)>1) || (!islist && length(penalty)>1)) 
    stop('may not specify multiple potential penalties when method="optimize"')

  objective <- function(pen, X, Y, z) { 

	##Problem with sending so many auxiliary parameters to nlminb -
	##nlminb's internal parameters got shifted somehow
	n <- z$n; penalty.matrix <- z$penalty.matrix; pennames <- z$pennames
	isols <- z$isols; islist <- z$islist; tol <- z$tol; maxit <- z$maxit
	ns <- z$ns; fitter <- z$fitter; pr <- z$pr; atr <- z$atr;
	tdf <- length(z$target.df)

	if(length(pen) > 1) {
	  pen <- structure(as.list(pen), names=pennames)
	  penfact <- Penalty.setup(atr, pen)$multiplier
	} else penfact <- pen
	
	if(length(penfact)==1 || !islist) pm <- penfact*penalty.matrix
	  else {
		a <- diag(sqrt(penfact))
		pm <- a %*% penalty.matrix %*% a
	  }
    f <- fitter(X, Y, penalty.matrix=pm, tol=tol, maxit=maxit)
    if(length(f$fail) && f$fail) 
      stop('fitter failed.  Try changing maxit or tol')

	if(isols) {
	  ## ols (from lm.pfit) already stored correct LR chisq and effective df
	  stats <- f$stats
	  df <- stats['d.f.']
	  lr <- stats['Model L.R.']
	  dag <- f$effective.df.diagonal
	} else {
	  v <- f$var   #Later: Varcov(f, regcoef.only=T)
	  f.nopenalty <- fitter(X, Y, initial=f$coef, maxit=1, tol=tol)
	  if(length(f.nopenalty$fail) && f.nopenalty$fail)
		stop('fitter failed.  Try changing tol')
	  info.matrix.unpenalized <- 
		if(length(f.nopenalty$info.matrix)) f.nopenalty$info.matrix else
	      solvet(f.nopenalty$var, tol=tol) # -> Varcov
	  dag <- diag(info.matrix.unpenalized %*% v)
	  df <- if(ns==0)sum(dag) else sum(dag[-(1:ns)])
	  lr <- f.nopenalty$stats["Model L.R."]
	}
    obj <- switch(z$which,
				  aic.c <- lr - 2*df*(1 + (df+1)/(n-df-1)),
				  aic   <- lr - 2*df,
				  bic   <- lr - df*logb(n))
    if(tdf) obj <- abs(df - z$target.df)
	if(pr) {
	  w <- if(tdf)df else obj
	  names(w) <- NULL
	  pp <- if(islist) unlist(pen) else c(Penalty=pen)
	  print(c(pp, Objective=w))
	}
	if(!tdf) obj <- -obj else attr(obj,'df') <- df
	obj
  }
  res <- nlminb(unlist(penalty), objective, lower=0, X=X, Y=Y, z=list(n=n,
		 penalty.matrix=penalty.matrix, pennames=names(penalty), 
		 isols=isols, islist=islist, tol=tol, maxit=maxit, ns=ns, 
		 fitter=fitter, atr=atr, pr=pr, which=which, target.df=target.df), 
		 control=nlminb.control(abs.tol=.00001,
		   rel.tol=if(tdf)1e-6 else .00001))
  return(list(penalty=if(islist)	
			  structure(as.list(res$parameters),names=names(penalty)) 
	           else res$parameters, 
			  objective=if(tdf)res$aux$df else -res$objective))
}

if(.R.) df <- aic <- bic <- aic.c <- 
  if(islist) double(length(penalty[[1]])) else double(length(penalty)) else
df <- aic <- bic <- aic.c <- 
  if(islist) single(length(penalty[[1]])) else single(length(penalty))

for(i in 1:np) {
  if(islist) {
    pen     <- penalty[i,]
    penfact <- Penalty.setup(atr, pen)$multiplier
  } else {
    pen     <- penalty[i]
    penfact <- pen
  }
  unpenalized <- all(penfact==0)

  if(i==1) Coef <- if(keep.coef) matrix(NA,ncol=length(f$coef),nrow=np)
      else NULL

  if(unpenalized) f <- fit else{ 
    if(length(penfact)==1 || !islist) pm <- penfact*penalty.matrix
    else {
      a <- diag(sqrt(penfact))
      pm <- a %*% penalty.matrix %*% a
    }
    f <- fitter(X, Y, penalty.matrix=pm, tol=tol, maxit=maxit)
    if(length(f$fail) && f$fail) 
      stop('fitter failed.  Try changing maxit or tol')
  }

  if(keep.coef) Coef[i,] <- f$coef

  if(unpenalized || isols) {
    # ols (from lm.pfit) already stored correct LR chisq and effective df
    stats <- f$stats
    df[i] <- stats['d.f.']
    lr    <- stats['Model L.R.']
    dag <- if(unpenalized) rep(1, length(df[i])) else f$effective.df.diagonal
  } else {
    v <- f$var   #Later: Varcov(f, regcoef.only=T)
    f.nopenalty <- fitter(X, Y, initial=f$coef, maxit=1, tol=tol)
    if(length(f.nopenalty$fail) && f.nopenalty$fail)
      stop('fitter failed.  Try changing tol')
    info.matrix.unpenalized <- 
      if(length(f.nopenalty$info.matrix)) f.nopenalty$info.matrix else
      solvet(f.nopenalty$var, tol=tol) # -> Varcov
    dag <- diag(info.matrix.unpenalized %*% v)
    df[i] <- if(ns==0)sum(dag) else sum(dag[-(1:ns)])
    lr <- f.nopenalty$stats["Model L.R."]
    if(verbose) {
      cat('non slopes',ns,'\neffective.df.diagonal:\n')
      print(dag)
    }
  }
  aic[i]   <- lr - 2*df[i]
  bic[i]   <- lr - df[i]*logb(n)
  aic.c[i] <- lr - 2*df[i]*(1 + (df[i]+1)/(n-df[i]-1))
  obj <- switch(which, aic.c=aic.c[i], aic=aic[i], bic=bic[i])

  if(obj > obj.best) {
    pen.best <- pen
    df.best <- df[i]
    obj.best <- obj
    f.best <- f
    var.adj.best <- if(unpenalized || isols) f$var else v %*% info.matrix.unpenalized %*% v
    diag.best <- dag
  }
  if(pr) {
    d <- if(islist)as.data.frame(pen, row.names='') else
      data.frame(penalty=pen, row.names='')
    d$df <- df[i]
    d$aic <- aic[i]
    d$bic <- bic[i]
    d$aic.c <- aic.c[i]
    cat('\n'); print(d)
  }
}
mat <- if(islist)as.data.frame(penalty) else
  data.frame(penalty=penalty)
mat$df <- df
mat$aic <- aic
mat$bic <- bic
mat$aic.c <- aic.c

structure(list(penalty=pen.best, df=df.best, objective=obj.best, 
          fit=f.best, var.adj=var.adj.best, diag=diag.best,
          results.all=mat, Coefficients=Coef), class="pentrace")
}


plot.pentrace <- function(x, method=c('points','image'),
						  which=c('effective.df','aic','aic.c','bic'), 
						  pch=2, add=FALSE, ylim, ...) {

method <- match.arg(method)

x     <- x$results.all

penalty      <- x[[1]]
effective.df <- x$df
aic          <- x$aic
bic          <- x$bic
aic.c        <- x$aic.c

if(length(x) == 5) {  ## only one variable given to penalty=

##  domfrow <- all(par('mfrow')==1) && ('effective.df' %in% which)
##  if(domfrow) par(mfrow=c(1,2))

  if('effective.df' %in% which)
	if(!add) plot(penalty, effective.df, xlab="Penalty", ylab="Effective d.f.",
	   type="l", ...)
  if(!add) plot(penalty, aic, 
				ylim=if(missing(ylim))range(c(aic,bic)) else ylim,
				xlab="Penalty", ylab="Information Criterion", 
				type=if('aic' %in% which)"l" else "n", lty=1, ...)
  else if('aic' %in% which) lines(penalty, aic, lty=2, ...)
  if('bic' %in% which) lines(penalty, bic,  lty=3, ...)
  if('aic.c' %in% which) lines(penalty, aic.c, lty=1, ...)
  if(!add)title(sub=paste(if('aic.c' %in% which) "Solid: AIC_c",
				  if('aic' %in% which) "Dotted: AIC",
				  if('bic' %in% which) "Dashed: BIC",sep='  '),
				adj=0,cex=.75)

##  if(domfrow) par(mfrow=c(1,1))
  return(invisible())
}

## At least two penalty factors
if(add) stop('add=T not implemented for >=2 penalty factors')

X1 <- x[[1]]
X2 <- x[[2]]
nam <- names(x)
x1 <- sort(unique(X1))
x2 <- sort(unique(X2))
n1 <- length(x1)
n2 <- length(x2)

aic.r <- rank(aic); aic.r <- aic.r/max(aic.r)

if(method=='points') {
  plot(0, 0, xlim=c(1,n1), ylim=c(1,n2), xlab=nam[1], ylab=nam[2], 
	   type='n', axes=FALSE, ...)
  mgp.axis(1, at=1:n1, labels=format(x1))
  mgp.axis(2, at=1:n2, labels=format(x2))
  ix <- match(X1, x1)
  iy <- match(X2, x2)
  for(i in 1:length(aic)) points(ix[i], iy[i], pch=pch, cex=(.1+aic.r[i])*3)
  return(invisible(aic.r))
}

z <- matrix(NA,nrow=n1,ncol=n2)
for(i in 1:n1) 
  for(j in 1:n2) z[i,j] <- aic.r[X1==x1[i] & X2==x2[j]]
image(x1, x2, z, xlab=nam[1], ylab=nam[2], zlim=range(aic.r), ...)
invisible(aic.r)
}

print.pentrace <- function(x, ...) {
cat('\nBest penalty:\n\n')
pen <- if(is.list(x$penalty)) as.data.frame(x$penalty,row.names='')
  else data.frame(penalty=x$penalty, row.names='')
pen$df <- x$df
pen$aic <- x$aic
print(pen)
cat('\n')
if(is.data.frame(x$results.all)) print(x$results.all) else
print(as.data.frame(x$results.all, 
  row.names=rep('',length(x$results.all[[1]]))))
## added if(is.data.frame()) 28apr02
invisible()
}

effective.df <- function(fit, object) {

atr <- fit$Design
if(!length(atr)) atr <- getOldDesign(fit)

dag <- if(missing(object)) fit$effective.df.diagonal else object$diag
if(length(dag)==0) stop('object not given or fit was not penalized')

ia.or.nonlin <- param.order(atr, 2)
nonlin       <- param.order(atr, 3)
ia           <- param.order(atr, 4)
ia.nonlin    <- param.order(atr, 5)

ns <- num.intercepts(fit)
if(ns>0) dag <- dag[-(1:ns)]

z <- rbind(c(length(dag),        sum(dag)),
           c(sum(!ia.or.nonlin), sum(dag[!ia.or.nonlin])),
           c(sum(ia.or.nonlin),  sum(dag[ia.or.nonlin])),
           c(sum(nonlin),        sum(dag[nonlin])),
           c(sum(ia),            sum(dag[ia])),
           c(sum(ia.nonlin),     sum(dag[ia.nonlin])))

dimnames(z) <- list(c('All','Simple Terms','Interaction or Nonlinear',
  'Nonlinear', 'Interaction','Nonlinear Interaction'),
  c('Original','Penalized'))

cat('\nOriginal and Effective Degrees of Freedom\n\n')
print(round(z,2))
invisible(z)
}
## First variable in ... cannot be named x - S methods try to call plot.default

plot.Design <- function(x, ..., xlim, ylim, fun, xlab, ylab,
						conf.int=.95, conf.type=c('mean','individual'),
						add=FALSE,  label.curves=TRUE,  eye,
						lty, col=1, lwd=par("lwd"), lwd.conf=1, pch=1,
						adj.zero=FALSE, ref.zero=FALSE,
						adj.subtitle, cex.adj,
						non.slopes, time=NULL, loglog=FALSE,
                        val.lev=FALSE, digits=4, 
						log="", perim, 
						method=c("persp","contour","image",
						  "dotchart","default"),
						sortdot=c('neither','ascending','descending'),
						nlevels=10, name,
						zlim=range(zmat,na.rm=TRUE),
                        vnames=c('labels','names'), abbrev=FALSE) {

  fit <- x
  conf.type <- match.arg(conf.type)
  vnames <- match.arg(vnames)
  dotlist <- list(...)
  fname  <- if(missing(name)) '' else name
  at     <- fit$Design
  if(!length(at)) at <- getOldDesign(fit)
  assume <- at$assume.code


#  if(length(assume)==0)stop("fit does not have design information")
  name <- at$name	##interactions are placed at end by design

  if(any(name == 'time')) {  ## 7apr03
    dotlist$time <- time
    time <- NULL
  }

  if(length(fname)>1 || (length(dotlist)==0 && fname=='')) {
	m <- match.call(expand=FALSE)
	m[[1]] <- as.name('plot.Design')
	for(nam in if(length(fname)>1)fname else name[assume!=9]) {
	  m$name <- nam
      if(vnames=='names') m$xlab <- nam  # 10Mar01
	  lv <- eval(m)
	}
	return(invisible(lv))
  }

#  .Options$digits <- digits  14Sep00
  oldopt <- options(digits=digits)
  on.exit(options(oldopt))

  method <- match.arg(method)
  sortdot <- match.arg(sortdot)

  cex <- par('cex')
  if(missing(cex.adj))   cex.adj   <- .75*cex

#  abbrev <- NULL   10Mar01

  Pretty <- function(x){ ##handles chron objects
	if(inherits(x,"dates") | inherits(x,"times"))
	  structure(pretty(oldUnclass(x)), class=oldClass(x))
	  else pretty(x)
  }
  

  f <- sum(assume!=9)	##limit list to main effects factors
  parms <- at$parms
  label <- at$label
  values <- at$values
  yunits <- fit$units   ## was units  8oct02
  units  <- at$units    ## 8oct02
  scale <- fit$scale.pred
  if(!length(scale)) scale <- "X * Beta"
  Center <- fit$center
  if(length(Center)==0) Center <- 0

  if(missing(ylab))	 {
	if(!missing(fun)) ylab <- ""
	  else	if(length(time)) {
		if(loglog) ylab <- paste("log[-log S(",format(time),")]",sep="")
		  else ylab <- paste(format(time),yunits,"Survival Probability")
	  }
		else ylab <- scale[1]
  }
  trlab <- if(.R. && ylab=='X * Beta') function(k) expression(X*hat(beta))
	else function(k) k   ## 24nov02 and everywhere trlab is used

  if(ref.zero & length(time))
	stop("ref.zero=T does not make sense with time given")

  labelc <- is.list(label.curves) || label.curves

  if(fname=='') factors <- dotlist else {
	factors <- list(NA)
	names(factors) <- fname
  }

  nf <- length(factors)

  if(nf<1)stop("must specify 1 or 2 factors to plot")
  which <- charmatch(names(factors), name, 0)
  if(any(which==0))stop(paste("factor(s) not in design:",
		   paste(names(factors)[which==0],collapse=" ")))
  if(any(assume[which]==9))stop("cannot plot interaction terms")

  ##Number of non-slopes
  nrp <- num.intercepts(fit)
  if(missing(non.slopes)) {
	non.slopes <- rep(0, nrp)
	if(!adj.zero) non.slopes[1] <- 1
  }
  if(nrp>0 && length(non.slopes)!=nrp)stop("wrong # values in non.slopes")

  if(is.logical(conf.int)) {
	if(conf.int) conf.int <- .95	else conf.int <- 0
  }

  if(conf.int) {
    vconstant <- 0
    if(conf.type=='individual') {
      vconstant <- fit$stats['Sigma']^2
      if(is.na(vconstant))
        stop('conf.type="individual" requires that fit be from ols')
    }
    zcrit <- if(length(idf <- fit$df.residual)) qt((1+conf.int)/2, idf) else
      qnorm((1+conf.int)/2)
  }

  ix <- which[1]
  ixt <- assume[ix]
  Limval <- Getlim(at, allow.null=TRUE, need.all=FALSE)
  parmx <- parms[[name[ix]]]
  if(ixt!=7 && is.numeric(parmx)) parmx <- NULL
  ##if(ixt!=5 && ixt<7 && !is.null(V <- Limval$values[[name[ix]]])) parmx <- V
  ## Was ixt<8 above   14Jun97
  ## Replaced above line with following 4 3apr03 - was using partial
  ## matching of names
  if(ixt != 5 && ixt < 7) {
    Lv <- Limval$values
    if(any(names(Lv)==name[ix])) parmx <- Lv[[name[ix]]]
  }
  if(ixt>8)stop("matrix or interaction factor may not be displayed")
#  val.lev <- val.lev & (ixt==5 | ixt==8)  commented out 20Jun99

  xadj <- Limval$limits[2,,drop=FALSE]
  ##for(i in 1:length(xadj))    14Feb95
  ##   if(is.factor(xadj[[i]])) xadj[[i]] <- as.character(xadj[[i]])

  if(adj.zero) for(i in 1:length(xadj)) {
	if(assume[i]==5 | assume[i]==8) xadj[[i]] <- parms[[name[i]]][1]
	  else if(!is.null(V <- Limval$values[[name[i]]]) & is.character(V))
		xadj[[i]] <- V[1]
		else xadj[[i]] <- 0
  }

  ##Use default adjusted to, replace some later.  Will be NA if 
  ## datadist doesn't have the variable

  nv <- 1
  xseq <- factors[[1]]

  if(nf>1) {
	y <- factors[[2]]
	ny <- length(y)
	if(ny>1 || (ny==1 && is.na(y))) nv <- 2	
  }

  if(missing(adj.subtitle)) {
	if(add)adj.subtitle <- FALSE	else {
	  adj.subtitle <- f-nv <= 6	
	}
  }
  
  jf <- nv
  if(nf>nv) for(i in which[(nv+1):nf])	 {
	jf <- jf+1
	z <- factors[[jf]]
	if(!is.na(z)) z <- value.chk(at, i, z, 0, Limval)
	if(length(z)!=1)stop("must specify single value for adjustment factors")
	if(!is.na(z)) xadj[[name[i]]] <- z
  }
  beta <- fit$coef
  if(length(beta) & conf.int>0) cov <- Varcov(fit, regcoef.only=TRUE)

  plot.type <- "curves"
  curve.labels <- NULL
  xval <- parmx
  if(nv>1) {
	iy <- which[2]
	isna <- sapply(xadj,is.na); isna[c(ix,iy)] <- FALSE
	if(any(isna))
	  stop(paste("variables not set to values or defined with datadist:",
				 paste(name[isna],collapse=" ")))
	iyt <- assume[iy]
	parmy <- parms[[name[iy]]]
	if(iyt!=5 && iyt<8 && !is.null(V <- Limval$values[[name[iy]]]))
	  parmy <- V

	if(iyt>8)stop("matrix or interaction factor may not be displayed")
	if(((!is.null(xval)) | ixt==5 | ixt==7 | ixt==8) &
	   (iyt!=5 & iyt!=7 & iyt!=8))
	  warning("plot more effective if continuous factor on x-axis")
	y <- value.chk(at, iy, y, 40, Limval)
	ny <- length(y)
	if(!(ixt==5 | ixt==7 | ixt==8 | iyt==5 | iyt==7 | iyt==8 | 
		 ny<=10)) plot.type <- "3d" }

  ## Use 40 x 40 grid if two continuous factors (for perspective or contour plot)

  if(plot.type=="3d") conf.int <- 0
  xseq <- value.chk(at, ix, xseq, 
					if(plot.type=="curves")100 else 40, Limval)
  if(is.factor(xseq)) xseq <- as.character(xseq)

  if(missing(xlab))
    xlab <- labelPlotmath(label[ix],units[ix],
                          plotmath=!(plot.type=='3d' && method=='persp'))
  ## was label[ix] 8oct02; persp( ) not support expressions in R

  ##No intercept in model -> expand factors at adjustment values for later
  ##subtraction in variance.  Also needed for ref.zero option.

  if(ref.zero | (nrp==0 & conf.int>0)) {
	if(any(sapply(xadj,is.na))) 
	  stop("with conf.int>0 for ref.zero=T or Cox model, all variables must have been defined with datadist")
	xadjdf <- structure(xadj, class="data.frame", row.names="1")
	xsub <- predictDesign(fit,xadjdf,type="x",non.slopes=non.slopes)
  }
  if(ref.zero) ycenter <- matxv(xsub, beta) - Center else ycenter <- 0
  ## - Center added 4 Oct 96

#  xseqn <- if(length(parmx) && is.character(parmx)) match(xseq, parmx)
#	else xseq   5Sep97
  xseqn <- if(length(parmx) && is.character(parmx)) 1:length(xseq) else xseq
  if(val.lev) xseqn <- as.numeric(xseq)  ## was as.single 21apr03

  if(nv==1) {
	nna <- sapply(xadj, is.na); nna[ix] <- FALSE  ##ignore this one
	if(any(nna))
	  stop(paste("variables not set to values here or defined with datadist:"
				 , paste(name[nna],collapse=" ")))
	if(missing(lty)) lty <- 1
	##Expand xadj to ##rows=length(xseq), replace col. ix with xseq
	adj <- oldUnclass(xadj)  ##so can expand one column
	adj[[name[ix]]] <- xseq   ## was adj[[ix]] 16jul02
	adj <- expand.grid(adj)
	if(!length(time)) {
	  xx <- predictDesign(fit, adj, type="x",non.slopes=non.slopes)
	  if(!is.null(attr(xx,"strata")) && any(is.na(attr(xx,"strata"))))
		warning("Computed stratum NA.  Requested stratum may not\nexist or reference values may be illegal strata combination\n")

	  if(length(xx)==0)
		stop("model has no covariables and survival not plotted")
	  xb <- matxv(xx, beta) - Center    ## Center 31May94

	  if(conf.int>0) {
		##Subtract from rows if need to center for variance
		if(nrp==0) xxx <- sweep(xx,2,xsub) else xxx <- xx
		var <- ((xxx %*% cov) * xxx) %*% rep(1,ncol(xxx))
		lower <- xb - zcrit*sqrt(vconstant + var)
		upper <- xb + zcrit*sqrt(vconstant + var)
	  }
	  if(ref.zero) {
		xb <- xb-ycenter
		if(conf.int>0) {lower <- lower-ycenter;upper <- upper-ycenter}
	  }
	}
	  else {
		xb <- survest(fit, adj, times=time,loglog=loglog, 
					  conf.int=conf.int)
		lower <- xb$lower; upper <- xb$upper; xb <- xb$surv
	  }

	if(!missing(fun)) {
	  xb <- fun(xb)
	  if(conf.int>0) {
		lower <- fun(lower)
		upper <- fun(upper)
	  }	
	}
	if(missing(ylim)) {
	  if(conf.int>0) ylim <- range(c(lower,upper),na.rm=TRUE)	else
	  ylim <- range(xb,na.rm=TRUE)
	}
	if((ixt==5 | ixt==8 | (!is.null(xval))) & !val.lev)	{
	  if(method=='dotchart') {  ## 21jul02 next line not always format()
		w <- xb; names(w) <- if(is.numeric(xseq)) format(xseq) else xseq
		isd <- switch(sortdot, ascending=order(w), descending=order(-w),
					  neither=1:length(w))
	  }
	  if(add) {
		if(method=='dotchart') 
		  dotchart2(w[isd], add=TRUE, pch=pch, reset.par=FALSE) else
		points(xseqn,xb,col=col[1],pch=pch)
	  } else {
		labs <- format(xseq)
		if(is.character(xseq) || ((ixt==5 | ixt==8) || 
								  (length(xseq)==length(xval) &&
								   all(abs(xseq-as.numeric(xval))<.00001))))
		  {att <- 1; at <- xseqn}
		  else {att <- 2; at <- Pretty(xseqn)}
		xlm <- if(missing(xlim)) range(at) else xlim
		if(method=='dotchart') dotchart2(w[isd], pch=pch, reset.par=FALSE,
			 xlim=ylim, xlab=trlab(ylab), ylab='')  # was ylab=xlab 21jul02
		  else {
			plot(xseqn,xb,xlab=xlab, xlim=xlm, ##Handle NAs in Y
				 axes=FALSE, ylim=ylim, ylab=trlab(ylab), log=log, type='n') #26Oct96
			points(xseqn,xb,col=col[1],pch=pch)
			if(att==1) mgp.axis(1,at=at,
                 labels=if(abbrev && ixt %in% c(5,8))
                 abbreviate(labs) else labs)  #7Feb98,2Jun99,10Mar01
			  else mgp.axis(1,at=at)
			mgp.axis(2,at=pretty(ylim))
		  }
	  }
	  if(conf.int>0) {
		if(method=='dotchart') {
		  dotchart2(lower[isd], add=TRUE, pch='[', reset.par=FALSE)
		  dotchart2(upper[isd], add=TRUE, pch=']', reset.par=FALSE)
		} else {
		  points(xseqn,lower,pch="-",col=col[1])
		  points(xseqn,upper,pch="-",col=col[1])
		}
	  }
	} else {
	  if(!add) {
		xlm <- if(missing(xlim)) range(Pretty(xseqn)) else xlim
		plot(xseqn,xb,xlab=xlab, xlim=xlm, ylim=ylim,
             ylab=trlab(ylab), log=log,
			 type='n', axes=FALSE) # 26Oct96, 7Feb98
		mgp.axis(1, at=pretty(xlm))  #2Jun99
		mgp.axis(2, at=pretty(ylim))
	  }
	  lines(xseqn,xb,lty=lty[1],lwd=lwd[1],col=col[1])
	  if(conf.int>0) {
		lines(xseqn,lower,lty=2,col=col[1],lwd=lwd.conf)
		lines(xseqn,upper,lty=2,col=col[1],lwd=lwd.conf)
	  }
	}
	xx <- list(); xx[[name[ix]]] <- xseq
	xx <- structure(xx,class="data.frame",
					row.names=as.character(1:length(xseq)))	
  }  ## end nv=1
	else {
	  ##Expand xadj into length(xseq)*ny rows, replace columns
	  ##ix and iy with xseq and y
	  adj <- oldUnclass(xadj)   ##so can expand a column
	  adj[[name[ix]]] <- NULL	##Guarantee y moves fastest (expand.grid moves
	  ##first factor first
	  adj[[name[iy]]] <- y
	  adj[[name[ix]]] <- xseq
	  adj <- expand.grid(adj)
	  xx <- predictDesign(fit,adj,type="x",non.slopes=non.slopes)
	  if(!is.null(attr(xx,"strata")) && any(is.na(attr(xx,"strata"))))
		warning("Computed stratum NA.  Requested stratum may not\nexist or reference values may be illegal strata combination\n")
	  if(length(xx)==0) {
		xb <- 0
		if(!length(time))
		  stop("model has no covariables and survival not plotted") 
	  }
		else xb <- matxv(xx, beta) - Center   ## Center 31May94
	  if(length(time)) {
		xb <- survest(fit, adj, times=time, loglog=loglog, 
					  conf.int=conf.int)
		lower <- xb$lower; upper <- xb$upper; xb <- xb$surv
	  }	##was conf.int=F
		else {
		  if(conf.int>0) {
			if(nrp==0) xxx <- sweep(xx,2,xsub) else xxx <- xx
			var <- ((xxx %*% cov) * xxx) %*% rep(1,ncol(xxx))
			lower <- xb - zcrit*sqrt(vconstant + var)
			upper <- xb + zcrit*sqrt(vconstant + var)			
		  }
		  if(ref.zero) {
			xb <- xb-ycenter
			if(conf.int>0){lower <- lower-ycenter;upper <- upper-ycenter}	
		  }
		}
	  if(!missing(fun))	{
		xb <- fun(xb)
		if(conf.int>0) {
		  lower <- fun(lower)
		  upper <- fun(upper)
		}
	  }

	  if(!missing(perim) && plot.type=="3d") {
        if(.SV4.) perim <- matrix(oldUnclass(perim), nrow=nrow(perim),
                                  dimnames=dimnames(perim))
		Ylo <- approx(perim[,1],perim[,2],adj[[name[ix]]])$y
		Yhi <- approx(perim[,1],perim[,3],adj[[name[ix]]])$y
		Ylo[is.na(Ylo)] <-  1e30
		Yhi[is.na(Yhi)] <- -1e30
		xb[adj[[name[iy]]] < Ylo] <- NA
		xb[adj[[name[iy]]] > Yhi] <- NA
	  }

	  if(missing(ylim))	 {
		if(conf.int>0)ylim <- range(c(lower,upper),na.rm=TRUE)	else
		ylim <- range(xb,na.rm=TRUE)	
	  }

	  xx <- adj[,c(name[ix],name[iy])]

	  if(plot.type=="3d") {
		zmat <- matrix(pmin(ylim[2],pmax(ylim[1],xb)),
					   nrow=length(xseqn),
					   ncol=ny, byrow=TRUE)

        laby <- labelPlotmath(label[iy],units[iy],
                              plotmath=method!='persp')  ## 8oct02
		if(method=="contour") {
		  contour(xseqn, y, zmat,
				  nlevels=nlevels, xlab=xlab, ylab=laby)
          ## was  ylab=label[iy]) 8oct02; + changes in persp image()
		} else if(method=="persp") {
		  if(missing(eye))
			persp(xseqn, y, zmat,zlim=zlim,
				  xlab=xlab,ylab=laby,zlab=trlab(ylab),box=TRUE)
			else persp(xseqn, y, zmat,zlim=zlim,
					   xlab=xlab,ylab=laby,zlab=trlab(ylab),eye=eye,box=TRUE)
		} else image(xseqn, y, zmat, xlab=xlab, ylab=laby)
	  } else {

		## One curve for each value of y, excl style used for C.L.

		lty <- if(missing(lty)) seq(ny+1)[-2] else
		           rep(lty, length=ny)
		i <- 0
		if(labelc) curves <- vector('list',ny)
		col <- rep(col, length=ny)
		lwd <- rep(lwd, length=ny)
		for(ay in y) {
		  i <- i+1
		  index.y <- (1:nrow(xx))[xx[,2]==ay]
		  xseq <- xx[index.y, name[ix]]
		  if(is.factor(xseq)) xseq <- as.character(xseq)
		  if(val.lev) xl <- as.single(xseq)
			else if(is.character(xseq)) xl <- match(xseq, parmx)
			  else xl <- xseq
		  curve.labels <-c(curve.labels, format(ay))
          if(!missing(perim)) { # 26Nov00
            if(!is.function(perim))stop('perim must be a function')
            show.pt <- perim(xl, ay)
            xb[index.y][!show.pt] <- NA
            if(conf.int) {
              lower[index.y][!show.pt] <- NA
              upper[index.y][!show.pt] <- NA
            }
          }
		  if(labelc) curves[[i]] <- list(xl, xb[index.y])
		  if(i>1 | add) lines(xl,xb[index.y],
							  lty=lty[i],col=col[i],lwd=lwd[i]) else {
								if((ixt==5 | ixt==8 | (!is.null(xval))) &
								   !val.lev)	 {
								  labs <- format(xseq)
								  if(is.character(xseq) || 
									 ((ixt==5 | ixt==8) || (length(xseq)==
												  length(xval) &&
												  all(abs(xseq-as.numeric(xval))<.00001))))
									{att <- 1; at <- xl}
									else {att <- 2; at <- Pretty(xl)}
								  xlm <- if(missing(xlim)) range(at) else xlim
								  plot(xl,xb[index.y],log=log,
									   xlab=xlab,ylab=trlab(ylab),
									   xlim=xlm, ylim=ylim,
									   type='n', axes=FALSE) #26Oct96,7Feb98
								  lines(xl,xb[index.y],lty=lty[i],lwd=lwd[i],
										col=col[i]) # 26Oct96
								  if(att==1) mgp.axis(1,at=at,
                                       label=if(abbrev && ixt %in% c(5,8))
                                       abbreviate(labs) else labs) #2Jun99,10Mar01
									else mgp.axis(1,at=at)
								  mgp.axis(2,at=pretty(ylim))
								} else
								{
								  xlm <- if(missing(xlim)) range(Pretty(xl)) else xlim
								  plot(xl,xb[index.y],xlab=xlab,
                                       ylab=trlab(ylab),
									   xlim=xlm,ylim=ylim,log=log,type='n',
									   axes=FALSE)  #7Feb98
								  mgp.axis(1,at=pretty(xlm)) #2Jun99
								  mgp.axis(2,at=pretty(ylim))
								  lines(xl,xb[index.y],col=col[i],
										lwd=lwd[i],lty=lty[i]) #26Oct96
								}
							  }
		  if(conf.int>0) {
			lines(xl,lower[index.y],
				  lty=2,col=col[i],lwd=lwd.conf)
			lines(xl,upper[index.y],
				  lty=2,col=col[i],lwd=lwd.conf) 
		  }
		}
	  if(labelc) labcurve(curves, curve.labels, 
						  opts=label.curves, lty=lty, lwd=lwd, col=col)
	  }
	}

  xx[[if(ylab=="")"Z" else ylab]] <- xb

  if(conf.int>0) {xx$lower <- lower; xx$upper <- upper}

  adjust <- ""
  if(f>nv) for(i in 1:f) {
	if(!any(i==which[1:nv]))
	  adjust <- paste(adjust, name[i], "=", 
					  if(is.factor(xadj[[i]])) as.character(xadj[[i]]) else
					  format(xadj[[i]])," ",sep="")
  }

  if(adjust!="" & adj.subtitle)
	title(sub=paste("Adjusted to:",adjust),adj=0,cex=cex.adj)

  R <- list(x.xbeta=xx, adjust=adjust, curve.labels=curve.labels, 
			plot.type=plot.type, method=method, col=col, 
			lty=if(plot.type=='curves') lty, lwd=lwd) #, abbrev=abbrev)
  oldClass(R) <- "plot.Design"
  invisible(R)
}

print.plot.Design <- function(x, ...) {
  print(x$x.xbeta)
  if(x$adjust!="") cat("Adjust to:",x$adjust,"\n")
  if(!is.null(cl <- x$curve.labels)) cat("Curves:",cl,"\n")
  invisible()
}

Legend <- function(object, ...) UseMethod("Legend")
#Legend.default <- legend

Legend.plot.Design <- function(object, x, y, size = c(1, 1), 
                               horizontal = TRUE, nint = 50, fun, at, 
                               zlab,  ...) {
  if(missing(x)) 
	if(object$method=="image" && missing(size)) {
	  cat("Using function \"locator(2)\" to place opposite corners of image.legend\n"
		  )
	  x <- locator(2)
	}
	  else {
		cat("Using function \"locator(1)\" to place upper left corner of legend\n"
			)
		x <- locator(1)
	  }
  if(!missing(y)) x <- list(x=x,y=y)   ## 17Sep96
  xb <- object$x.xbeta
##  if(obj$plot.type=="curves") {
	##  legend(x, obj$curve.labels, lty=obj$lty,
	##        col=rep(obj$col,length=length(obj$lty)), ...)
##	if(!exists('key'))
##	  stop('must type "library(trellis)" to have access to the key() function for making legend')
##	if(length(obj$abbrev)) {
##	  if(length(unique(obj$lty)) < 2) 
##		key(x$x, x$y, text=list(paste(obj$abbrev,obj$curve.labels,sep='  ')))
##		else key(x$x, x$y, lines=list(lty=obj$lty),
##				 text=list(paste(obj$abbrev,obj$curve.labels,sep='  ')))
##	} else
##	key(x$x, x$y,
##		lines=list(lty=obj$lty, col=rep(obj$col,length=length(obj$lty))),
##		text=list(obj$curve.labels),
##		...)
	##  return(if(missing(y))x else list(x=x,y=y))
##	return(invisible(x))  ## 17Sep96
##  }
  if(object$method!="image") stop('expecting to use results from method="image"')
  z <- xb[,3]
  irgz <- seq(min(z,na.rm=TRUE), max(z,na.rm=TRUE), length = nint)
  lirgz <- length(irgz)
  if(horizontal) {
	f <- expression({
	  if(length(list(...))) par(...) ##axis() does not respect mgp
	  image(x = irgz, y = 1:lirgz, z = matrix(irgz, lirgz, lirgz), yaxt
			= "n", xaxt=if(missing(fun))"s" else "n")
	  if(!missing(fun)) mgp.axis(1, if(missing(at)) pretty(irgz) else at,
							 labels=format(fun(if(missing(at)) pretty(irgz) else at)))
      ##2Jun99
	  title(xlab=if(missing(zlab)) names(xb)[3] else zlab)
	})
	subplot(x = x, y = y, size = size, fun = f, hadj=0, vadj=1)
  } else {
	f <- expression({
	  if(length(list(...))) par(...)
	  image(x = 1:lirgz, y = irgz, z = matrix(irgz, lirgz, lirgz,byrow=TRUE), 
			xaxt = "n", yaxt=if(missing(fun))"s" else "n")
	  if(!missing(fun)) mgp.axis(2, if(missing(at)) pretty(irgz) else at,
							 labels=format(fun(if(missing(at)) pretty(irgz) else at))) #2Jun99
	  title(ylab=if(missing(zlab)) names(xb)[3] else zlab)
	})
	subplot(x = x, y = y, size = size, fun = f, hadj=0, vadj=1)
  }
  invisible(if(missing(y))x else list(x=x,y=y))
}

perimeter <- function(x, y, xinc=diff(range(x))/10, n=10, lowess.=TRUE) {

  s <- !is.na(x+y)
  x <- x[s]
  y <- y[s]
  m <- length(x)
  if(m<n) stop("number of non-NA x must be >= n")
  i <- order(x)
  x <- x[i]
  y <- y[i]
  s <- n:(m-n+1)
  x <- x[s]
  y <- y[s]

  x <- round(x/xinc)*xinc

  g <- function(y, n) {
	y <- sort(y)
	m <- length(y)
	if(n>(m-n+1)) c(NA,NA) else c(y[n], y[m-n+1])
  }

  r <- unlist(tapply(y, x, g, n=n))
  i <- seq(1, length(r), by=2)
  rlo <- r[i]
  rhi <- r[-i]
  s <- !is.na(rlo+rhi)
  if(!any(s)) stop("no intervals had sufficient y observations")
  x <- sort(unique(x))[s]
  rlo <- rlo[s]
  rhi <- rhi[s]
  if(lowess.) {
	rlo <- lowess(x, rlo)$y
	rhi <- lowess(x, rhi)$y
  }
  structure(cbind(x, rlo, rhi), dimnames=list(NULL,
								  c("x","ymin","ymax")), class='perimeter')
}

lines.perimeter <- function(x, ...) {
  lines(x[,'x'], x[,'ymin'],...)
  lines(x[,'x'], x[,'ymax'],...)
  invisible()
}

datadensity.plot.Design <- function(object, x1, x2, ...) {
  if(missing(x1)) stop('must specify x1')
  r <- object$x.xbeta
  nam <- names(r)
  x1.name <- deparse(substitute(x1))
  if(x1.name!=nam[1]) warning(paste(x1.name,
       ' is not first variable mentioned in plot() (',nam[1],')',sep=''))

  if(missing(x2)) {
	x <- r[[1]]
	y <- r[[2]]
	x1 <- x1[!is.na(x1)]
	y.x1 <- approx(x, y, xout=x1)$y
	scat1d(x1, y=y.x1, ...)
	return(invisible())
  }

  x2.name <- deparse(substitute(x2))
  if(x2.name!=nam[2]) warning(paste(x2.name,
       ' is not second variable mentioned in plot() (',nam[2],')',sep=''))
  x     <- r[[1]]
  curve <- r[[2]]
  y     <- r[[3]]
  for(s in if(is.factor(curve))levels(curve) else unique(curve)) {
	i <- curve==s
	xs <- x[i]
	ys <- y[i]
	x1s <- x1[x2==s]
	x1s <- x1s[!is.na(x1s)]
	y1s <- approx(xs, ys, xout=x1s)$y
	scat1d(x1s, y=y1s, ...)
  }
  invisible()
}

plot.xmean.ordinaly <- function(x, data, subset, na.action,
	subn=TRUE, cr=FALSE, ...) {

X <- match.call(expand=FALSE)
X$subn <- X$cr <- X$... <- NULL
if(missing(na.action)) X$na.action <- na.keep
Terms <- if(missing(data)) terms(x) else terms(x, data=data)
X$formula <- Terms
X[[1]] <- as.name('model.frame')
X <- eval(X, sys.parent())
resp <- attr(Terms, 'response')
if(resp==0) stop('must have a response variable')

nx <- ncol(X) - 1
Y <- X[[resp]]
nam <- as.character(attr(Terms, 'variables'))
if(.R.) nam <- nam[-1]

for(i in 1:nx) {
  x <- X[[resp+i]]
  if(is.category(x)) stop('categorical predictors not allowed')
  s <- !is.na(oldUnclass(Y)+x)
  y <- Y[s]
  x <- x[s]
  n <- length(x)
  f <- lrm.fit(x, y)
  fy <- f$freq/n

  ##Following is pulled out of predict.lrm
  ns <- length(fy) - 1  # number of intercepts
  k <- ns + 1
  intcept <- f$coef[1:ns]
  xb <- f$linear.predictors - intcept[1]
  xb <- sapply(intcept, '+', xb)
  P <- 1/(1+exp(-xb))
  
  P <- matrix(P, ncol=ns)
  P <- cbind(1, P) - cbind(P, 0)  #one column per prob(Y=j)

  xmean.y <- tapply(x, y, mean)
  xp <- x*P/n
  xmean.y.po <- apply(xp, 2, sum)/fy
  yy <- 1:length(fy)
  rr <- c(xmean.y, xmean.y.po)
  if(cr) {
    u <- cr.setup(y)
    s <- u$subs
    yc <- u$y
    xc <- x[s]
    cohort <- u$cohort
    xcohort <- matrix(0, nrow=length(xc), ncol=length(levels(cohort))-1)
    xcohort[col(xcohort)==oldUnclass(cohort)-1] <- 1  # generate dummies
    cof <- lrm.fit(cbind(xcohort, xc), yc)$coefficients
    cumprob <- rep(1, n)
    for(j in 1:k) {
      P[,j] <- cumprob* (if(j==k) 1 else plogis(cof[1] + (if(j>1) cof[j] else 0) + cof[k]*x))
      cumprob <- cumprob - P[,j]
    }
    xp <- x*P/n
    xmean.y.cr <- apply(xp, 2, sum)/fy
    rr <- c(rr, xmean.y.cr)
  }
  plot(yy, xmean.y, type='b', ylim=range(rr),
       axes=FALSE, xlab=nam[resp], ylab=paste('Mean of',nam[resp+i]), ...)
  mgp.axis(1, at=yy, labels=names(fy))
  mgp.axis(2)
  lines(yy, xmean.y.po, lty=2, ...)
  if(cr) points(yy, xmean.y.cr, pch='C')
  if(subn) title(sub=paste('n=',n,sep=''),adj=0)

}
invisible()
}

pphsm <- function(fit) {

warning("at present, pphsm does not return the correct covariance matrix")

clas <- c(oldClass(fit), fit$fitFunction)
if(!any(c('psm','survreg') %in% clas))
   stop("fit must be created by psm or survreg")  ##14Nov00
if(.newSurvival.) {
  if(fit$dist %nin% c('exponential','weibull'))
    stop("fit must have used dist='weibull' or 'exponential'")
} else {
  if(fit$family[1]!="extreme" | fit$family[2]!="log")
    stop('fit must have used dist="extreme" and link="log"')
}

fit$coefficients <- -fit$coefficients/(if(.newSurvival.)fit$scale else exp(fit$parms))
fit$scale.pred   <- c("log Relative Hazard","Hazard Ratio")
if(.SV4.) fit$fitFunction <- c('pphsm',fit$fitFunction) else
 oldClass(fit) <- c("pphsm",oldClass(fit))   ##14Nov00

fit
}

print.pphsm <- function(x, correlation = TRUE, ...) {
  if (length(f <- x$fail) && f) {
    stop(" Survreg failed.  No summary provided")
    return(invisible(x))
  }

cat("Parametric Survival Model Converted to PH Form\n\n")
print.psm(x, correlation=correlation)
invisible()
}
#Requires fastbw

predab.resample <- function(fit.orig,  fit, measure, 
		method=c("boot","crossvalidation",".632","randomization"),
		bw=FALSE, B=50, pr=FALSE, rule="aic", type="residual",
		sls=.05, aics=0, strata=FALSE, tol=1e-12, 
		non.slopes.in.x=TRUE, kint=1, cluster, subset, group=NULL, ...) {

method <- match.arg(method)
# .Options$digits <- 4  14Sep00
oldopt <- options(digits=4)
on.exit(options(oldopt))

#Following logic prevents having to load a copy of a large x object
if(any(match(c("x","y"),names(fit.orig),0)==0))
   stop("must have specified x=T and y=T on original fit")
fparms <- fit.orig[c("non.slopes","assign","terms","Design")]
if(!length(fparms$Design))fparms$Design <- getOldDesign(fit.orig) #10Jul01

non.slopes <- num.intercepts(fit.orig)
x.index <- if(non.slopes==0 || non.slopes.in.x) function(i,...) i else
	function(i, ns) { if(any(i>ns)) i[i>ns]-ns else NULL }  #23May94

Xb <- function(x, b, non.slopes, non.slopes.in.x, n, kint=1) {
  if(length(x)) {
    if(non.slopes==0 || non.slopes.in.x) x %*% b else
      b[kint] + x %*% b[-(1:non.slopes)]
  } else {
    if(non.slopes==0) rep(0,n) else rep(b[kint],n)
  }
}

nac <- fit.orig$na.action

x <- as.matrix(fit.orig$x)
n <- nrow(x)
attr(x,'class') <- NULL	#Remove model.matrix class for subset operations later

y <- fit.orig$y
y <- as.matrix(if(is.category(y)) oldUnclass(y) else y)  ##25Mar98

multi <- !missing(cluster)   # some subjects have multiple records now

# 19Mar99:
if(length(group)) {
  if(multi || method!='boot')
    stop('group is currently allowed only when method="boot" and cluster is not given')
  if(length(group) > n) {
	## Missing observations were deleted during fit
	if(length(nac)) {
	  j <- !is.na(naresid(nac, y) %*% rep(1,ncol(y)))
	  group <- group[j]
	}
  }
  if(length(group) != n)
	stop('length of group does not match # rows used in fit')
  group.inds <- split(1:n, group)  ## see bootstrap()
  ngroup <- length(group.inds)
} else ngroup <- 0
  
if(multi) {
  if(method!='boot') stop('cluster only implemented for method="boot"')
  if(length(cluster) > n) {
	## Missing observations were deleted during fit
	if(length(nac)) {
	  j <- !is.na(naresid(nac, y) %*% rep(1,ncol(y)))
	  cluster <- cluster[j]
	}
  }
  if(length(cluster) != n)
	stop('length of cluster does not match # rows used in fit')
  if(any(is.na(cluster))) stop('cluster has NAs')
  n.orig <- length(unique(cluster))
  cl.samp <- split(1:n, cluster)
} else n.orig <- n

if(!missing(subset)) {
  if(length(subset) > n && length(nac)) {
    j <- !is.na(naresid(nac, y) %*% rep(1,ncol(y)))
    subset <- subset[j]
  }
  if(length(subset) != n  && all(subset>=0))
    stop('length of subset does not match # rows used in fit')
  if(any(is.na(subset))) stop('subset has NAs')
  if(!is.logical(subset)) {
    subset2 <- rep(FALSE, n)
    subset2[subset] <- TRUE
    subset <- subset2
    subset2 <- NULL
    }
}

if(strata)			{
	stra <- attr(fit.orig$x, "strata")
	if(!length(stra)) stra <- rep(1, nrow(y))
	y <- cbind(y, stra)	}

if(bw)				{
#	fit.orig <- fit(x,y,iter=0,tol=tol,...)
	if(fit.orig$fail) return()
	cat("\n		Backwards Step-down - Original Model\n")
	fbw <- fastbw(fit.orig,rule=rule,type=type,sls=sls,aics=aics,eps=tol)
	print(fbw)
	orig.col.kept <- fbw$parms.kept
	if(!length(orig.col.kept))stop("no variables kept in original model")
	xcol <- x.index(orig.col.kept, non.slopes)
	fit.orig <- fit(x[,xcol,drop=FALSE], y, iter=0, tol=tol, xcol=xcol, ...)

				}	else 
	orig.col.kept <- 1:length(fit.orig$coef)

b <- fit.orig$coef
xcol <- x.index(orig.col.kept, non.slopes)
xb <- Xb(x[,xcol,drop=FALSE], b, non.slopes, non.slopes.in.x, n,
         kint=kint)

index.orig <- if(missing(subset))measure(xb, 
	y, fit=fit.orig,
	iter=0,	evalfit=TRUE, fit.orig=fit.orig, kint=kint, ...)   else
  measure(xb[subset], y[subset,,drop=FALSE], fit=fit.orig,
        iter=0, evalfit=FALSE, fit.orig=fit.orig, kint=kint, ...)

test.stat <- single(length(index.orig))
train.stat <- test.stat
#name <- attr(fparms$terms,"Design")$name   10Jul01
name <- fparms$Design$name
if(bw) 	{
	varin <- matrix("", nrow=B, ncol=length(name))
	nvarin <- rep(NA,B)
	}

j <- 0
num <- 0

if(method=="crossvalidation")		{ 
	per.group <- n/B
	if(per.group<2)stop("B>n/2")
	sb <- sample(n, replace=FALSE)	}
#Cross-val keeps using same random set of indexes, without replacement

ntest <- 0 #Used in getting weighted average for .632 estimator

## if(exists('.Random.seed'))   28jul03
##  cat(".Random.seed:",.Random.seed,"in",find(.Random.seed)[1],"\n")
## exists( ) 8oct02

if(method==".632")
{
   #Must do assignments ahead of time so can weight estimates
   #according to representation in bootstrap samples
   S <- matrix(integer(1), nrow=n, ncol=B)
   W <- matrix(TRUE, nrow=n, ncol=B)
   for(i in 1:B)
   {
	S[,i] <- s <- sample(n, replace=TRUE)
	W[s,i] <- FALSE  #now these obs are NOT omitted
   }
   nomit <- drop(W %*% rep(1,ncol(W)))  #no. boot samples omitting each obs
   if(min(nomit)==0) stop("not every observation omitted at least once in bootstrap samples.\nRe--run with larger B")
   W <- apply(W/nomit, 2, sum)/n
   cat("\n\nWeights for .632 method (ordinary bootstrap weights ",
	format(1/B),")\n",sep="")
   print(summary(W))
}

if(!pr) cat("Iteration:\n")

for(i in 1:B)								{
	if(!pr) { cat(i,""); if(i %% 20 == 0) cat("\n") }
	switch(method,
	crossvalidation=
		{	is <- 1 + round((i-1)*per.group)
			ie <- min(n, round(is+per.group-1))
			test <- sb[is:ie]
			train <- -test	}, #cross-val
    boot=	{
      if(ngroup) {
        train <- integer(n.orig)
        for(si in 1:ngroup) {
          gi <- group.inds[[si]]
          lgi <- length(gi)
          train[gi] <- if(lgi==1) gi else sample(gi, lgi, replace=TRUE)
          ## 6May99: sample behaves differently when first arg is a single integer
        }
      } else {
        train <- sample(n.orig, replace=TRUE)
        if(multi) train <- unlist(cl.samp[train])
      }
			test <- 1:n  },    #boot
	".632"=	{	train <- S[,i]
			test <- -train},   #boot .632
	randomization=	
		{	train <- sample(n, replace=FALSE)
			test <- 1:n   })   #randomization
	xtrain <- if(method=="randomization") 1:n else train
	f <- fit(x[xtrain,,drop=FALSE], y[train,,drop=FALSE], iter=i, tol=tol,...)
	f$assign <- NULL  #Some programs put a NULL assign (e.g. ols.val fit)
 
	fail <- f$fail
	if(!fail)			{
      ## Following if..stop was before f$assign above   28Apr99
      if((ni <- num.intercepts(f)) != non.slopes) 
        stop(paste('\nA training sample has a different number of intercepts (',
                   ni,')\n than the original model fit (',non.slopes,').  \nYou probably fit an ordinal model with sparse cells and a re-sample\ndid not select at least one observation for each value of Y.\nAdd the argument group=y where y is the response variable.\nThis will force balanced sampling on levels of y.',sep=''))
      clf <- attr(f,"class")  # class is removed by c() below
      f[names(fparms)] <- fparms  # 23Dec99
      ##      f <- c(f, fparms)     23Dec99
      attr(f, "class") <- clf
      if(!bw) 				{
        coef <- f$coef  # 14Sep00, coefficients->coef 14Aug01
        col.kept <- 1:length(coef)
      }	else	{
        f <- fastbw(f,rule=rule,type=type,
                    sls=sls,aics=aics,eps=tol)
        if(pr)print(f)
        varin[j+1, f$factors.kept] <- "*"   #did have drop=F
        nvarin[j+1] <- length(f$factors.kept)
        col.kept <- f$parms.kept
        if(!length(col.kept)) f <- fit(NULL, y[train,,drop=FALSE],
                                       iter=i, tol=tol,...)	else     {
                                         xcol <- x.index(col.kept, non.slopes)
                                         f <- fit(x[xtrain,xcol,drop=FALSE], y[train,,drop=FALSE],
                                                  iter=i, tol=tol, xcol=xcol, ...) }
        if(f$fail) fail <- TRUE else coef <- f$coef  #14Sep00 14Aug01
      }	}
	if(!fail)	{
      j <- j+1
      xcol <- x.index(col.kept, non.slopes)
      xb <- Xb(x[,xcol,drop=FALSE], coef, non.slopes, non.slopes.in.x, n,
               kint=kint)
      if(missing(subset)) {
		train.statj <- measure(xb[xtrain], y[train,,drop=FALSE], 
                               fit=f, iter=i,fit.orig=fit.orig,evalfit=TRUE, 
                               kint=kint, ...)
		test.statj <- measure(xb[test], y[test,,drop=FALSE], fit=f, 
                              iter=i,fit.orig=fit.orig, evalfit=FALSE, kint=kint, ...)
      } else {
		ii <- xtrain
		if(any(ii<0)) ii <- (1:n)[ii]
		ii <- ii[subset[ii]]
		train.statj <- measure(xb[ii], y[ii,,drop=FALSE],
                               fit=f, iter=i, fit.orig=fit.orig,evalfit=FALSE,
                               kint=kint, ...)
		ii <- test
		if(any(ii<0)) ii <- (1:n)[ii]
		ii <- ii[subset[ii]]
		test.statj <- measure(xb[ii], y[ii,,drop=FALSE], fit=f,
                              iter=i, fit.orig=fit.orig, evalfit=FALSE, kint=kint, ...)
      }
      na <- is.na(train.statj+test.statj)
      num <- num + !na
      if(pr) print(cbind(training=train.statj, test=test.statj))
      train.statj[na] <- 0
      test.statj[na] <- 0
      if(method==".632") 
		{
          ##wt <- length(xb[test])*(!na)  else wt <- 1
		  wt <- W[i]
		  if(any(na))warning('method=".632" does not properly handle missing summary indexes')
		}
      else wt <- 1
      train.stat <- train.stat + train.statj
      test.stat <- test.stat + test.statj * wt
      ntest <- ntest + 1   #was +wt
    } 
  }
if(!pr)cat("\n\n")
if(j!=B) cat("\nDivergence or singularity in",B-j,"samples\n")
train.stat <- train.stat/num
if(method!=".632") 	{
  test.stat <- test.stat/num
  optimism <- train.stat - test.stat
}	else	{
	optimism <- .632 * (index.orig - test.stat)
  }
res <- cbind(index.orig=index.orig,training=train.stat,test=test.stat,
	optimism=optimism,index.corrected=index.orig-optimism,n=num)

if(bw) {
	varin <- varin[1:j, ,drop=FALSE]
	nvarin <- nvarin[1:j]
#	dimnames(varin) <- list(rep("",j), abbreviate(name,1:2))
dimnames(varin) <- list(rep("",j), name)
	cat("\n		Factors Retained in Backwards Elimination\n\n")
	print(varin, quote=FALSE)
	cat("\n         Frequencies of Numbers of Factors Retained\n\n")
    tvarin <- table(nvarin)
    if(.R.) names(dimnames(tvarin)) <- NULL
	print(tvarin)
  }

res
}

##newdata=data frame, vector,  matrix, or list.  All but first assume data
##need coding, e.g. categorical variables are given as integers
##variable missing for all obs -> use adjust-to value in limits
##(means (parms) for matrx)

## Renamed from predict.Design 6dec02; let predict.Design be
## dispatcher (see Design.Misc.s)

predictDesign <- function(fit, newdata,
   type=c("lp","x","data.frame","terms","adjto","adjto.data.frame",
    "model.frame"),
   se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
   incl.non.slopes, non.slopes, kint=1,
   na.action=na.keep, expand.na=TRUE, center.terms=TRUE, ...)	{

type <- match.arg(type)
conf.type <- match.arg(conf.type)
## R does not preserve missing(x):   31jul02
mnon.slopes <- missing(non.slopes)


at <- fit$Design
if(!length(at)) at <- getOldDesign(fit)
assume <- at$assume.code
Limval <- Getlim(at, allow.null=TRUE, need.all=FALSE)
Values <- Limval$values
non.ia <- assume!=9
non.strat <- assume!=8
f <- sum(non.ia)
nstrata <- sum(assume==8)
somex <- any(non.strat)
rnam <- NULL
cox <- inherits(fit, "cph") ||
 (length(fit$fitFunction) && any(fit$fitFunction=='cph'))  ##14Nov00 22May01
naa <- fit$na.action
if(!expand.na) naresid <- function(a,b) b #don't really call naresid if drop NAs

parms <- at$parms
name <- at$name
coeff <- fit$coefficients
nrp <- num.intercepts(fit)
if(mnon.slopes) {
   non.slopes <- rep(0,nrp)
   non.slopes[kint] <- 1   #13Sep94
}
else if(nrp>0 & length(non.slopes)!=nrp)
stop("length of non.slopes incompatible with fit")

int.pres <- nrp>0  # was !(cox|lrm)
if(somex) cov <- Varcov(fit,regcoef.only=TRUE)    #remove scale params
if(missing(incl.non.slopes)) 
   incl.non.slopes <- !mnon.slopes | (!missing(kint)) | 
                      int.pres | type!="x"
##added 12Feb93   !missing() added 18Feb93, 2nd one 13Sep94
int.pres <- int.pres & incl.non.slopes

assign <- fit$assign

nama <- names(assign)[1]
asso <- 1*(nama=="Intercept" | nama=="(Intercept)")

Center <- if(cox)fit$center else 0

oldopts <- options(contrasts=c(factor="contr.treatment",ordered="contr.poly"),
   Design.attr=at)

##20Nov00   In SV4 options(two lists) causes problems
on.exit({options(contrasts=oldopts$contrasts)
         options(Design.attr=NULL)})

## Terms <- delete.response(terms(attr(fit$terms,"formula"), specials="strat"))

Terms <- if(.R.) delete.response(terms(formula(fit), specials='strat')) else
     delete.response(terms(fit$terms, specials="strat"))  ## 17Apr02  30may02
attr(Terms,"response") <- 0
attr(Terms,"intercept") <- 1    # was int.pres 12Feb93
##Need intercept whenever design matrix is generated to get
##current list of dummy variables for factor variables

stra <- attr(Terms, "specials")$strat

offset <- 0	# used if no newdata

Terms.ns <- Terms
if(length(stra))	{
  Terms.ns <- Terms[-stra]	#Uses [.terms function. 3.0 did not add 1!
  ## was [1-stra], changed 7June94

  ## For some reason attr(...) <- pmin(attr(...)) changed a detail
  ## in factors attribute in R but R and SV4 don't seem to need this
  ## anyway   1may02
  if(!.R.) {
    tfac <- attr(Terms.ns,'factors')
    if(length(tfac) && any(tfac > 1))
      attr(Terms.ns,'factors') <- pmin(tfac, 1)
  }
}

if(conf.int) {
  vconstant <- 0
  if(conf.type=='individual') {
    vconstant <- fit$stats['Sigma']^2
    if(is.na(vconstant))
      stop('conf.type="individual" requires that fit be from ols')
  }
  zcrit <- if(length(idf <- fit$df.residual)) qt((1+conf.int)/2, idf) else
           qnorm((1+conf.int)/2)
}

if(type!="adjto" & type!="adjto.data.frame") {
  X <- NULL
  if(missing(newdata)) {
    if(type=="lp" && length(fit$linear.predictors)) {
      LP <- naresid(naa, fit$linear.predictors)   #changed 8June94
      if(kint>1) LP <- LP-fit$coef[1]+fit$coef[kint]  #added 13Sep94
      if(length(stra <- attr(fit$linear.predictors,"strata")))
        attr(LP, "strata") <- naresid(naa, stra)
      if(!se.fit && !conf.int)return(LP)  ##7Mar99
      else if(length(fit$se.fit)) {
        if(kint>1)
          warning("se.fit is retrieved from the fit but it corresponded to kint=1")
        retlist <- list(linear.predictors=LP, se.fit=naresid(naa,fit$se.fit))
        if(conf.int) {
          plminus <- zcrit*sqrt(retlist$se.fit^2 + vconstant)
          retlist$lower <- LP - plminus
          retlist$upper <- LP + plminus
        }
        return(retlist)
      }
    }
    else if(type=="x") return(structure(naresid(naa,fit$x),
              strata=if(length(stra <- attr(fit$x,"strata")))
              naresid(naa,stra) else NULL))
    X <- fit$x
    rnam <- dimnames(X)[[1]]
    if(!any(names(fit)=="x")) X <- NULL  #fit$x can get fit$xref
    if(!length(X))
      stop("newdata not given and fit did not store x")
  }
  if(!length(X)) {
    if(!is.data.frame(newdata))	{
      if(is.list(newdata)) {
        loc <- if(!length(names(newdata))) 1:f else name[assume!=9]
        new <- matrix(if(.R.)double(1) else single(1),
                      nrow=length(newdata[[1]]),
                      ncol=length(newdata))
        for(j in 1:ncol(new)) new[,j] <- newdata[[loc[j]]]
        newdata <- new
      }
      if(!is.matrix(newdata)) newdata <- matrix(newdata, ncol=f)
      if(ncol(newdata)!=f)stop("# columns in newdata not= # factors in design")
      X <- list()
      k <- 0
      ii <- 0
      for(i in (1:length(assume))[non.ia]) {
        ii <- ii+1
        xi <- newdata[,ii]
        as <- assume[i]
        allna <- all(is.na(xi))
        ##	   if(as!=10 && allna) xi <- at$limits[3,ii]
        if(as==5 | as==8)	{
          xi <- as.integer(xi)
          levels(xi) <- parms[[name[i]]]
          oldClass(xi) <- "factor"
        }
        else if(as==10) {
          if(i==1) ifact <- 1
          else ifact <- 1 + sum(assume[1:(i-1)]!=8)
          ##	Accounts for assign not being output for strata factors
          ncols <- length(assign[[ifact+asso]])
          if(allna) {
            xi <- matrix(if(.R.)double(1) else single(1),
                         nrow=length(xi), ncol=ncols)
            for(j in 1:ncol(xi)) xi[,j] <- parms[[name[i]]][j]	}
          else xi <- matrix(if(.R.)xi else as.single(xi),
                            nrow=length(xi), ncol=ncols) }
        ##	Duplicate single value for all parts of matrix
        k <- k + 1
        X[[k]] <- xi
      }
      names(X) <- name[non.ia]
      attr(X, "row.names") <- as.character(1:nrow(newdata))
      oldClass(X) <- "data.frame"
      newdata <- X
      ##Note: data.frame() converts matrix variables to individual variables
      if(type=="data.frame") return(newdata)
    }
    else {
      ##Need to convert any factors to have all levels in original fit
      ##Otherwise, wrong dummy variables will be generated by model.matrix
      nm <- names(newdata)
      for(i in 1:ncol(newdata))	{
        j <- match(nm[i], name)
        if(!is.na(j)) {
          asj <- assume[j]
          w <- newdata[,i]
          V <- NULL
          if(asj==5 | asj==7 | asj==8 | 
             (length(V <- Values[[name[j]]]) && is.character(V))) {
            if(length(Pa <- parms[[name[j]]])) V <- Pa   #added 8Apr94
            ## if(is.null(V)) V <- parms[[name[j]]]  #subtracted 8Apr94
            newdata[,i] <- factor(w, V)
            ##Handles user specifying numeric values without quotes, that
            ##are levels
            ww <- is.na(newdata[,i]) & !is.na(oldUnclass(w))
            if(any(ww)) 	{
              cat("Error in predictDesign: Values in",names(newdata)[i],
                  "not in",V,":\n")
              print(as.character(w[ww]),quote=FALSE); stop()
            }
          }
        }
      }
    }
  
    newdata <- addOffset4ModelFrame(Terms,newdata)  ## 23nov02
    X <- model.frame(Terms, newdata, na.action=na.action, ...)
    if(type=="model.frame") return(X)
    naa <- attr(X, "na.action")
    rnam <- row.names(X)

    offs <- attr(Terms, "offset")
    if(!length(offs)) offset <- rep(0, length(rnam))
    else offset <- X[[offs]]

    ## if(ncol(X) != sum(non.ia))stop("improperly formed model frame")
    strata <- list()
    nst <- 0
    ii <- 0  ## 23nov02
    for(i in setdiff(1:ncol(X),offs)) {   ## setdiff() was 1:ncol(X) 23nov02
      ii <- ii + 1
      xi <- X[[i]]
      asi <- attr(xi,"assume.code")
      as <- assume[ii]              ## was i 23nov02
      if(!length(asi) && as==7) {
        attr(X[,i],"contrasts") <- 
          attr(scored(xi,name=name[ii]),"contrasts") ## was i 23nov02
        if(length(xi)==1) warning("a bug in model.matrix can produce incorrect results\nwhen only one observation is being predicted for an ordered variable")
      }

      if(as==8) {
        nst <- nst+1
        strata[[nst]] <- paste(name[ii],"=",parms[[name[ii]]][X[,i]],sep="")
        ## was name[i] 23nov02
      }
    }
    if(!somex) X <- NULL
    else if(int.pres && nrp==1) X <- model.matrix(Terms.ns, X) #nrp Jan94
    else X <- model.matrix(Terms.ns, X)[,-1,drop=FALSE]		#12Feb93
  
    if(nstrata>0)	{
      names(strata) <- paste("S",1:nstrata,sep="")
      strata <- factor(interaction(as.data.frame(strata),drop=TRUE),
                       levels=fit$strata)
    }
  }

  else strata <- attr(X,"strata")

  added.col <- if(incl.non.slopes & (nrp>1 | !int.pres)) nrp else 0 #nrp>1 Jan94
  ## & !scale.pres removed from following statement 20Feb93
  if(incl.non.slopes & nrp>0 & somex & added.col>0) {
    xx <- matrix(if(.R.)double(1) else single(1),
                 nrow=nrow(X),ncol=added.col)
    for(j in 1:nrp) xx[,j] <- non.slopes[j]
  }
  else xx <- NULL
}

##For models with multiple intercepts, delete elements of covariance matrix
##containing unused intercepts
elements.to.delete <- 9999
if(somex && nrp>1) {
  i <- (1:nrp)[non.slopes==0]; cov <- cov[-i,-i,drop=FALSE] 
  elements.to.delete <- i
}

if(type=="adjto" | type=="adjto.data.frame" | (center.terms && type=="terms")| 
   (cox & (se.fit | conf.int))) {
  ##Form design matrix for adjust-to values
  adjto <- list()
  ii <- 0
  for(i in (1:length(assume))[non.ia]) {
    ii <- ii+1
    xi <- Getlimi(name[i], Limval, need.all=TRUE)[2] #was =F  5Feb94
    if(assume[i]==5 | assume[i]==8) xi <- factor(xi,parms[[name[i]]])
    else if(assume[i]==7) xi <- scored(xi, name=name[i])
    else if(assume[i]==10) xi <- matrix(parms[[name[i]]],nrow=1) #matrx col medians
    adjto[[ii]] <- xi
  }
  names(adjto) <- name[non.ia]
  ##   adjto <- data.frame(adjto,check.names=FALSE)
  ##   data.frame will take matrix factors and convert into individual vars
  attr(adjto,"row.names") <- "1"
  oldClass(adjto) <- "data.frame"
  if(type=="adjto.data.frame") return(adjto)
  adjto <- addOffset4ModelFrame(Terms, adjto) ## 23nov02
  adjto <- model.frame(Terms, adjto)
  adjto <- if(int.pres) model.matrix(Terms.ns, adjto) else
           model.matrix(Terms.ns,adjto)[,-1,drop=FALSE]   # -1 added 12Feb93
  ## added drop=FALSE 27feb03
  if(type=="adjto")	{
    k <- if(int.pres) 1:length(coeff) else (nrp+1):length(coeff)
    if(is.matrix(adjto))
      dimnames(adjto) <- list(dimnames(adjto)[[1]],names(coeff)[k])
    else names(adjto) <- names(coeff)[k]
    return(adjto)
  }
}

if(length(xx) & type!="terms" & incl.non.slopes)	{
  X <- cbind(xx, X)
  dimnames(X) <- list(rnam, names(coeff))
  if((center.terms && type=="terms") | (cox & (se.fit | conf.int))) 
	adjto <- c(xx[1,], adjto)   #12Feb93
}

else if(somex) dimnames(X) <- 
  ##	list(rnam,names(coeff)[
  ##	 (nrp+1-(int.pres & incl.non.slopes)):length(coeff)])
  list(rnam,names(coeff)[(1+length(coeff)-ncol(X)):length(coeff)]) #22Jun95


storage.mode(X) <- "double"
if(type=="x") return(
     structure(naresid(naa,X), strata=if(nstrata>0) naresid(naa,strata) else NULL,
               offset=if(length(offs)) naresid(naa,offset) else NULL,
               na.action=if(expand.na)NULL else naa)
     )

if(type=="lp") {
  if(somex) {
    ## if( ) 28apr02
    if(elements.to.delete==9999) cof <- coeff else {
      cof <- coeff[-elements.to.delete]
      X <- X[,-elements.to.delete,drop=FALSE]
    }
	xb <- matxv(X, cof)+offset-Center
   	names(xb) <- rnam
   	if(!.R.) storage.mode(xb) <- "single"
  } else {if(!length(offs)) xb <- NULL else xb <- offset}
  xb <- naresid(naa, xb)
  if(nstrata>0)attr(xb,"strata") <- naresid(naa,strata)
  if((se.fit | conf.int) & somex) {
    if(cox) X <- sweep(X,2,adjto) #Center columns
    se <- drop(sqrt(((X %*% cov) * X) %*% rep(1, ncol(X))))
    names(se) <- rnam
    if(!.R.) storage.mode(se) <- "single"
    retlist <- structure(list(linear.predictors=xb, se.fit=naresid(naa,se)),
                         na.action=if(expand.na)NULL else naa)
    if(conf.int) {
      plminus <- zcrit*sqrt(retlist$se.fit^2 + vconstant)
      retlist$lower <- xb - plminus
      retlist$upper <- xb + plminus
    }
    return(retlist)
  }
  else return(structure(xb,na.action=if(expand.na)NULL else naa))
}

if(type=="terms") {
  if(!somex) stop('type="terms" may not be given unless covariables present')
  fitted <- array(0,c(nrow(X),sum(non.strat)),
                  list(rnam,name[non.strat]))
  if(se.fit) se <- fitted
  j <- 0
  if(center.terms) {
    ## 31jul02: lrm and perhaps others put out fit$x without column of
    ## intercepts but model has intercept
    if(ncol(adjto) != ncol(X)) {
      if(dimnames(adjto)[[2]][1] %in% c('Intercept','(Intercept)') &&
         dimnames(X)[[2]][1]    %nin% c('Intercept','(Intercept)'))
        adjto <- adjto[,-1,drop=FALSE]
      if(ncol(adjto) != ncol(X)) stop('program logic error')
    }
    X <- sweep(X, 2, adjto) # center columns
  }
  # PROBLEM: adjto = c(Intercept=1, sexmale=0); no 1s col in f$x
  num.intercepts.not.in.X <- length(coeff)-ncol(X)	#23Jan95
  for(i in (1:length(assume))[non.strat]) {
    j <- j+1
    k <- assign[[j+asso]]   #; m <- k+int.pres
    ko <- k - num.intercepts.not.in.X			#23Jun95
    fitted[,j] <- matxv(X[,ko,drop=FALSE], coeff[k])
    ## was X[,m], coeff[nrp+k]
    if(se.fit) se[,j] <- (((X[,ko,drop=FALSE] %*% cov[ko,ko,drop=FALSE]) * 
                           X[,ko,drop=FALSE]) %*% rep(1,length(ko)))^.5
  }
  if(!.R.) storage.mode(fitted) <- "single"
  fitted <- structure(naresid(naa,fitted), strata=if(nstrata==0) NULL else
                      naresid(naa, strata))
  if(se.fit) {
    if(!.R.) storage.mode(se) <- "single"
    return(structure(list(fitted=fitted, se.fit=naresid(naa,se)),
                     na.action=if(expand.na)NULL else naa)) 	}
  else return(structure(fitted, na.action=if(expand.na)NULL else naa))
}
}   
   
addOffset4ModelFrame <- function(Terms, newdata, offset=0) {
  offs <- attr(Terms,'offset')
  if(!length(offs)) return(newdata)
  offsetVarname <- all.names(attr(Terms,'variables')[offs+1])[1]
  offsetVarname <- offsetVarname[offsetVarname != 'offset']
  if(offsetVarname %nin% names(newdata)) {
    newdata[[offsetVarname]] <- rep(offset, length=nrow(newdata))
    warning(paste('offset variable set to',
                  paste(format(offset),collapse=' ')))
  }
  newdata
}
predict.lrm <- function(object, ..., 
		type=c("lp","fitted","fitted.ind","mean","x","data.frame",
		"terms", "adjto", "adjto.data.frame", "model.frame"),
		se.fit=FALSE, codes=FALSE) {

type <- match.arg(type)
if(!(type=="fitted"|type=="fitted.ind"|type=="mean"))
  return(predictDesign(object,...,type=type, se.fit=se.fit))

xb <- predictDesign(object, ..., type="lp", se.fit=FALSE)
rnam <- names(xb)
ns <- object$non.slopes
cnam <- names(object$coef[1:ns])
if(se.fit)warning('se.fit not supported with type="fitted" or type="mean"')
if(ns==1 & type=="mean")
  stop('type="mean" makes no sense with a binary response')
if(ns==1) return(1/(1+exp(-xb)))
intcept <- object$coef[1:ns]
xb <- xb - intcept[1]
xb <- sapply(intcept, "+", xb)
P <- 1/(1+exp(-xb))
nam <- names(object$freq)
if(is.matrix(P))dimnames(P) <- list(rnam, cnam)
else names(P) <- names(object$coef[1:ns])
if(type=="fitted") return(P)

#type="mean" or "fitted.ind"
vals <- names(object$freq)
k <- ns+1
P <- matrix(P,ncol=ns)
Peq <- cbind(1,P)-cbind(P,0)
if(type=="fitted.ind") 	{
   ynam <- as.character(attr(object$terms,"formula")[2])
   ynam <- paste(ynam,"=", vals, sep="")
   dimnames(Peq) <- list(rnam, ynam)
   return(drop(Peq))	}

#type="mean"
if(codes) vals <- 1:length(object$freq)
else	{
  vals <- as.numeric(vals)
  if(any(is.na(vals)))
    stop('values of response levels must be numeric for type="mean" and codes=F')
	}
m <- drop(Peq %*% vals)
names(m) <- rnam
m

}


#This is Terry Therneau's old print.coxreg with conf.int default to F
#Add Nagelkerke R2 9Jun92
#Remove printing hazard ratios 17Jun92
#Removed stats 23Jun92 since in print.cph now, delete print x$n if 1 stratum
print.cph.fit <- 
function(x, table = TRUE, coef = TRUE, conf.int = FALSE, scale = 1,
         digits = NULL, ...)
{
	if(table && !is.null(x$n) && is.matrix(x$n))
		print(x$n)
	if(is.null(digits))
		digits <- 3
	savedig <- options(digits = digits)
	on.exit(options(savedig))
	beta <- x$coef
	se <- sqrt(diag(x$var))
	if(is.null(beta) | is.null(se))
		stop("Input is not valid")
	if(coef) {
		tmp <- cbind(beta, se, beta/se, 1 - pchisq((beta/
			se)^2, 1))
		dimnames(tmp) <- list(names(beta), c("coef", 
			"se(coef)", "z", "p"))
		cat("\n")
		prmatrix(tmp)
	}
	if(conf.int) {
		z <- qnorm((1 + conf.int)/2, 0, 1)
		beta <- beta * scale
		se <- se * scale
		tmp <- cbind(exp(beta), exp( - beta), exp(beta - z * se), exp(
			beta + z * se))
		dimnames(tmp) <- list(names(beta), c("exp(coef)", "exp(-coef)",
			paste("lower .", round(100 * conf.int, 2), sep = ""),
			paste("upper .", round(100 * conf.int, 2), sep = "")))
		cat("\n")
		prmatrix(tmp)
	}
	invisible(x)
}
print.cph <- function(x, long=FALSE, digits=3, conf.int=FALSE,
                      table=TRUE,  ...) { 

cat("\n")
if(x$fail)	{
	cat("Model Did Not Converge\n")
	return()
		}

cat("Cox Proportional Hazards Model\n\n")
dput(x$call)
cat("\n")
if(!is.null(z <- x$na.action)) naprint(z)
if(!is.null(x$coef))						{
   stats <- x$stats
   stats[3] <- round(stats[3],2)
   stats[5] <- round(stats[5],4)
   stats[6] <- round(stats[6],2)
   stats[7] <- round(stats[7],4)
   stats[8] <- round(stats[8],3)
   if(.R.) print(format.sep(stats), quote=FALSE) else print(stats)
   cat("\n")
   print.cph.fit(x, digits=digits, conf.int=conf.int, table=table, ...)
   if(long)cat("Centering constant:",format(x$center),"\n")
 }
else if(table) print(x$n)
invisible()
									}

print.lrm <- function(x, digits=4, strata.coefs=FALSE, ...) {
sg <- function(x,d) 	{
#  .Options$digits <- d  14Sep00
  oldopt <- options(digits=d)
  on.exit(options(oldopt))
	format(x)	}
rn <- function(x,d) format(round(as.single(x),d))

cat("\n")
if(x$fail)	{
	cat("Model Did Not Converge\n")
	return()
		}

cat("Logistic Regression Model\n\n")
dput(x$call)
cat("\n\nFrequencies of Responses\n")
print(x$freq)
if(length(x$sumwty)) {
  cat('\n\nSum of Weights by Response Category\n')
  print(x$sumwty)
}
cat("\n")
if(!is.null(x$nmiss))	{  #for backward compatibility
   cat("Frequencies of Missing Values Due to Each Variable\n")
   print(x$nmiss)
   cat("\n")		}
else if(!is.null(x$na.action)) naprint(x$na.action)

ns <- x$non.slopes
nstrata <- x$nstrata
if(!length(nstrata)) nstrata <- 1

pm <- x$penalty.matrix
if(length(pm)) {
   psc <- if(length(pm)==1) sqrt(pm) else
	sqrt(diag(pm))
   penalty.scale <- c(rep(0,ns),psc)
   cof <- matrix(x$coef[-(1:ns)], ncol=1)
   cat("Penalty factors:\n\n"); print(as.data.frame(x$penalty, row.names=''))
   cat("\nFinal penalty on -2 log L:",
	rn(t(cof) %*% pm %*% cof,2),"\n\n")
}

#est.exp <- 1:ns
#if(length(f$est)) est.exp <- c(est.exp, ns+f$est[f$est+ns <= length(f$coef)])
vv <- diag(x$var)
cof <- x$coef
if(strata.coefs) {
  cof <- c(cof, x$strata.coef)
  vv  <- c(vv,  x$Varcov(x,which='strata.var.diag'))
  if(length(pm)) penalty.scale <- c(penalty.scale,rep(NA,x$nstrat-1))
}
score.there <- nstrata==1 && (length(x$est) < length(x$coef)-ns)
stats <- x$stats
stats[2] <- signif(stats[2],1)
stats[3] <- round(stats[3],2)
stats[4] <- round(stats[4],2)
stats[5] <- round(stats[5],4)
stats[6] <- round(stats[6],3)
stats[7] <- round(stats[7],3)
if(nstrata==1) { ##17Dec97
  stats[8] <- round(stats[8],3)   ##21Aug97
  stats[9] <- round(stats[9],3)
  stats[10] <- round(stats[10],3)
  if(length(stats)>10) {
    stats[11] <- round(stats[11],3)
    if(length(x$weights)) stats[12] <- round(stats[12],3)
  }
} else stats <- c(stats,Strata=x$nstrat)

if(.R.) {   ## 8Apr02
  nst <- length(stats)
  cstats <- character(nst)
  names(cstats) <- names(stats)
  for(i in 1:nst) cstats[i] <- format(stats[i])
  print(cstats, quote=FALSE)
} else if(!score.there) print(stats)	else	{
	print(stats[1:10])
	cat("\n")
	st <- stats[11:13]
	st[1] <- round(st[1],2)
	st[3] <- round(st[3],4)
	print(st)			}
cat("\n")

##if(length(f$var)==0) vv <- NULL	#doesn't bother with this for x=NULL
z <- cof/sqrt(vv)
stats <- cbind(sg(cof,digits), sg(sqrt(vv),digits), 
	rn(cof/sqrt(vv),2))
stats <- cbind(stats, rn(1-pchisq(z^2,1),4))
dimnames(stats) <- list(names(cof),
	c("Coef","S.E.","Wald Z","P"))
if(length(pm)) stats <- cbind(stats, "Penalty Scale"=sg(penalty.scale,digits))
print(stats,quote=FALSE)
cat("\n")


if(score.there)	{
	q <- (1:length(cof))[-est.exp]
	if(length(q)==1) vv <- x$var[q,q] else vv <- diag(x$var[q,q])
	z <- x$u[q]/sqrt(vv)
	stats <- cbind(rn(z,2), rn(1-pchisq(z^2,1),4))
	dimnames(stats) <- list(names(cof[q]),c("Score Z","P"))
	print(stats,quote=FALSE)
        cat("\n")
}
invisible()

}

print.ols <- function(x, digits=4, long=FALSE, ...)	{

#  .Options$digits <- digits  14Sep00
  oldopt <- options(digits=digits)
  on.exit(options(oldopt))

  cat("\n")

  cat("Linear Regression Model\n\n")
  dput(x$call)
  cat("\n")
  if(!is.null(z <- x$na.action)) naprint(z)
  stats <- x$stats
  if(lst <- length(stats)) {
    if(.R.) {  ## 8Apr02
      cstats <- character(lst)
      names(cstats) <- names(stats)
      for(i in 1:lst) cstats[i] <- format(stats[i])
      print(cstats, quote=FALSE)
    } else print(x$stats); cat('\n')}

  pen <- length(x$penalty.matrix) > 0

#  if(!pen) {    22Dec01
#	x <- summary.lm(f)
#	##The following is part of print.summary.lm
#	resid <- x$residuals
#  } else resid <- f$residuals
  resid <- x$residuals

  n <- length(resid)
  p <- length(x$coef)-(names(x$coef)[1]=="Intercept")
  if(length(x$stats)==0) cat("n=", n,"   p=",p,"\n\n",sep="")
#  if(pen) {  22Dec01
  ndf <- x$stats['d.f.']
  df <- c(ndf, n-ndf-1, ndf)
  r2 <- x$stats['R2']
#  } else {
#  df <- x$df
#  r2 <- x$r.squared
#	##  sigma <- x$sigma
#  }
  sigma <- x$stats['Sigma']
  rdf <- df[2]
  if(rdf > 5) {
	cat("Residuals:\n")
	if(length(dim(resid)) == 2) {
	  rq <- apply(t(resid), 1, quantile)
	  dimnames(rq) <- list(c("Min", "1Q", "Median", "3Q",
							 "Max"), dimnames(resid)[[2]])
	}
	  else {
		rq <- quantile(resid)
		names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
	  }
	print(rq, digits = digits, ...)
  }
	else if(rdf > 0) {
	  cat("Residuals:\n")
	  print(resid, digits = digits, ...)
	}
  if(nsingular <- df[3] - df[1])
	cat("\nCoefficients: (", nsingular, 
		" not defined because of singularities)\n", sep = "")
  
	else cat("\nCoefficients:\n")
##  if(!pen) print(x$coefficients) else {  22Dec01
  se <- sqrt(diag(x$var))
  z <- x$coefficients/se
  P <- 2*(1-pt(abs(z),rdf)) ## was pnorm 8feb03
  co <- cbind(x$coefficients,se,z,P)
  dimnames(co) <- list(names(x$coefficients),
                       c("Value","Std. Error","t","Pr(>|t|)"))
  ## was "Z" "Pr(>|Z|)"
  print(co)
##  } 22Dec01
  if(pen) cat('\n') else
  cat("\nResidual standard error:", format(signif(sigma, digits)),
	  "on", rdf, "degrees of freedom\n")
  rsqa <- 1 - (1 - r2)*(n-1)/rdf
  if(length(x$stats)==0)
	cat("Multiple R-Squared:", format(signif(r2  , digits))," ")
  cat("Adjusted R-Squared:", format(signif(rsqa, digits)), "\n")
  if(!pen) {
#	correl <- x$correlation  22Dec01
    if(long && p > 0) {
      correl <- diag(1/se) %*% x$var %*% diag(1/se)
      dimnames(correl) <- dimnames(x$var)
	  cat("\nCorrelation of Coefficients:\n")
	  ll <- lower.tri(correl)
	  correl[ll] <- format(round(correl[ll], digits), ...)
	  correl[!ll] <- ""
	  print(correl[-1,  - (p+1), drop = FALSE], quote = FALSE, digits = digits,
			...)
    }
  }
  cat("\n")

  invisible()
}

# SCCS @(#)summary.survreg.s	4.5 7/14/92
print.psm <- function(x, correlation = FALSE, ...)
{
    if (length(z <- x$fail) && z) {
	warning(" psm failed.", x$fail, "   No summary provided\n")
	return(invisible(x))
	}

    new <- .newSurvival.
    dist <- if(new) x$dist else x$family["name"]
    if(!.R.) survreg.distributions <- survReg.distributions
    name <- survreg.distributions[[dist]]$name
    cat("Parametric Survival Model:",name,"Distribution\n\n")
    dput(x$call)
    cat("\n")
    if(length(x$nmiss))	{   #backward compatibility
      cat("Frequencies of Missing Values Due to Each Variable\n")
      print(x$nmiss)
      cat("\n")
    } else if(length(z <- x$na.action)) naprint(z)
    stats <- x$stats
    stats[3] <- round(stats[3],2)
    stats[5] <- round(stats[5],4)
    stats[6] <- round(stats[6],2)
    if(.R.) print(format.sep(stats),quote=FALSE) else print(stats)
    cat("\n")

    if(new) {
      if(!x$fail) x$fail <- NULL    # summary.survreg uses NULL for OK
      s <- if(!.R.) summary.survReg(x, correlation=correlation) else
                    summary.survreg(x, correlation=correlation)
      print.summary.survreg2(s, correlation=correlation)
      return(invisible())
    }
    
    wt <- x$weights
    fparms <- x$fixed
    coef <- c(x$coef, x$parms[!fparms])
    resid <- x$residuals
    dresid <- x$dresiduals
    n <- length(resid)
    p <- x$rank
    if(!length(p))
        p <- sum(!is.na(coef))
    if(!p) {
        warning("This model has zero rank --- no summary is provided")
        return(x)
        }
    nsingular <- length(coef) - p
    rdf <- x$df.resid
    if(!length(rdf))
        rdf <- n - p
    R <- x$R   #check for rank deficiencies
    if(p < max(dim(R)))
        R <- R[1:p,     #coded by pivoting
        1:p]
    if(length(wt)) {
        wt <- wt^0.5
        resid <- resid * wt
        excl <- wt == 0
        if(any(excl)) {
            warning(paste(sum(excl), 
                "rows with zero weights not counted"))
            resid <- resid[!excl]
            if(!length(x$df.residual))
                rdf <- rdf - sum(excl)
            }
        }
    famname <- x$family["name"]
    if(!length(famname))
        famname <- "Gaussian"
    scale <- x$fparms
    nas <- is.na(coef)
    cnames <- names(coef[!nas])
    coef <- matrix(rep(coef[!nas], 4), ncol = 4)
    dimnames(coef) <- list(cnames, c("Value", "Std. Error", "z value", "p"))
    stds <- sqrt(diag(x$var[!nas,!nas,drop=FALSE]))
    coef[, 2] <- stds
    coef[, 3] <- coef[, 1]/stds
    coef[, 4] <- 2*pnorm(-abs(coef[,3]))
    if(correlation) {
	if(sum(nas)==1) ss <- 1/stds else ss <- diag(1/stds)
	correl <- ss %*% x$var[!nas, !nas, drop=FALSE] %*% ss
        dimnames(correl) <- list(cnames, cnames)
        }
    else correl <- NULL
    ocall <- x$call
    if(length(form <- x$formula)) {
        if(!length(ocall$formula))
	    ocall <- match.call(get("survreg"), ocall)
        ocall$formula <- form
        }
    dig <- .Options$digits
    print.summary.survreg(
	list(call = ocall, terms = x$terms, coefficients = coef,
	df = c(p, rdf), deviance.resid = dresid,
	var=x$var, correlation = correl, deviance = deviance(x),
	null.deviance = x$null.deviance, loglik=x$loglik,
	iter = x$iter,
	nas = nas))
   options(digits=dig)   #recovers from bug in print.summary.survreg
   invisible()
    }

## Mod of print.summary.survreg from survival5 - suppresses printing a
## few things, added correlation arg
if(.newSurvival.) {
  print.summary.survreg2 <-
    function (x, digits = max(options()$digits - 4, 3),
              correlation=FALSE, ...) 
{
    correl <- x$correl
    n <- x$n
    if (is.null(digits)) 
        digits <- options()$digits
    ##    cat("Call:\n")
    ##   dput(x$call)
    ##    cat('\n')  FEH
    print(x$table, digits = digits)
    if (nrow(x$var) == length(x$coefficients)) 
        cat("\nScale fixed at", format(x$scale, digits = digits), 
            "\n")
    else if (length(x$scale) == 1) 
        cat("\nScale=", format(x$scale, digits = digits), "\n")
    else {
        cat("\nScale:\n")
        print(x$scale, digits = digits, ...)
    }
##    cat("\n", x$parms, "\n", sep = "")
##    df <- sum(x$df) - x$idf
##    cat("Loglik(model)=", format(round(x$loglik[2], 1)), "  Loglik(intercept only)=", 
##        format(round(x$loglik[1], 1)))
##    if (df > 0) 
##        cat("\n\tChisq=", format(round(x$chi, 2)), "on", round(df, 
##            1), "degrees of freedom, p=", format(signif(1 - pchisq(x$chi, 
##            df), 2)), "\n")
##    else cat("\n")
##    cat("Number of Newton-Raphson Iterations:", format(trunc(x$iter)), 
##        "\n")
##    omit <- x$na.action
##    if (length(omit)) 
##        cat("n=", x$n, " (", naprint(omit), ")\n", sep = "")
##    else cat("n=", x$n, "\n")
    if (correlation && !is.null(correl)) {  ## FEH
        p <- dim(correl)[2]
        if (p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            ll <- lower.tri(correl)
            correl[ll] <- format(round(correl[ll], digits = digits))
            correl[!ll] <- ""
            print(correl[-1, -p, drop = FALSE], quote = FALSE)
        }
    }
    cat("\n")
    invisible(NULL)
}
NULL
}
# psm for SV4/R is a modification of Therneau's survreg from survival5
# (survReg in S-Plus 6)   FEH 17Apr02

# SCCS @(#)survreg.s	5.8 07/10/00
#  The newest version of survreg, that accepts penalties and strata
#

# .newSurvival. <- .R. || existsFunction('survReg')

psm <- if(.newSurvival.) 
  function(formula=formula(data),
           data=if(.R.)parent.frame() else sys.parent(),
           weights, subset, na.action=na.delete, dist='weibull', 
           init=NULL,  scale=0,
           control=if(!.R.)survReg.control() else survreg.control(),
           parms=NULL, model=FALSE, x=FALSE, y=TRUE, time.inc, ...) {

    if(.R.) require('survival')
    call <- match.call()
    m <- match.call(expand=FALSE)
    if(dist=='extreme')
warning('Unlike earlier versions of survreg, dist="extreme" does not fit\na Weibull distribution as it uses an identity link.  To fit the Weibull\ndistribution use the default for dist or specify dist="weibull".')
    m$na.action <- na.action  ## FEH
    temp <- c("", "formula", "data", "weights", "subset", "na.action")
    m <- m[ match(temp, names(m), nomatch=0)]
    if(.R.) m$drop.unused.levels <- TRUE  ## 31jul02
    m[[1]] <- as.name("model.frame")
    special <- c("strata", "cluster")
    Terms <- if(missing(data)) terms(formula, special)
             else              terms(formula, special, data=data)
    m$formula <- Terms
    ## Start FEH
    offs <- offset <- attr(Terms, "offset")  ## offs 23nov02 moved 6dec02
    if(!.R.) survreg.distributions <- survReg.distributions
    if(.R.) {
      dul <- .Options$drop.unused.levels
      if(!length(dul) || dul) {
        on.exit(options(drop.unused.levels=dul))
        options(drop.unused.levels=FALSE)
      }
    }
    m <- Design(eval(m, if(.R.)parent.frame() else sys.parent()))
    atrx <- attributes(m)
    nact <- atrx$na.action
    Terms <- atrx$terms
    atr   <- atrx$Design
    if(length(nact$nmiss)) {
      jia <- grep('*%ia%*',names(nact$nmiss))  ## 8feb03
      if(length(jia)) nact$nmiss <- nact$nmiss[-jia]
      s <- if(length(offs)) names(nact$nmiss) !=  atrx$names[offs] else TRUE
      ## 23nov02
      names(nact$nmiss)[s] <- 
        c(as.character(formula[2]), atr$name[atr$assume.code!=9])
    }
    ## End FEH    [s] 23nov02
    
    weights <- model.extract(m, 'weights')
    Y <- model.extract(m, "response")

    ## Start FEH
    atY <- attributes(Y)
    ncy <- ncol(Y)
    maxtime <- max(Y[,-ncy])
    nnn <- c(nrow(Y),sum(Y[,ncy]))
    time.units <- attr(Y,'units')
    if(!length(time.units)) time.units <- "Day"
    if(missing(time.inc))	{
      time.inc <- switch(time.units,Day=30,Month=1,Year=1,maxtime/10)
      if(time.inc>=maxtime | maxtime/time.inc>25)
        time.inc <- max(pretty(c(0,maxtime)))/10
    }
    ## End FEH

    if (!inherits(Y, "Surv")) stop("Response must be a survival object")

    strats <- attr(Terms, "specials")$strata
    cluster<- attr(Terms, "specials")$cluster
    dropx <- NULL
    if (length(cluster)) {
        if (missing(robust)) robust <- TRUE
        tempc <- untangle.specials(Terms, 'cluster', 1:10)
        ord <- attr(Terms, 'order')[tempc$terms]
        if (any(ord>1)) stop ("Cluster can not be used in an interaction")
        cluster <- strata(m[,tempc$vars], shortlabel=TRUE)  #allow multiples
        dropx <- tempc$terms
        }
    if (length(strats)) {
        temp <- untangle.specials(Terms, 'strata', 1)
        dropx <- c(dropx, temp$terms)
        if (length(temp$vars)==1) strata.keep <- m[[temp$vars]]
        else strata.keep <- strata(m[,temp$vars], shortlabel=TRUE)
        strata <- as.numeric(strata.keep)
	nstrata <- max(strata)
        }
    else {
	nstrata <- 1
	strata <- 0
	}

    if (length(dropx)) newTerms<-Terms[-dropx]
    else               newTerms<-Terms
    X <- model.matrix(newTerms,m)

    ## Start FEH
    rnam <- dimnames(Y)[[1]]
    dimnames(X) <- list(rnam, c("(Intercept)",atr$colnames))
    ## End FEH except for 23nov02 and later changes
    
    n <- nrow(X)
    nvar <- ncol(X)

    if (length(offset)) offset <- as.numeric(m[[offset]])  ## length 23nov02
    else                offset <- rep(0, n)

    if (is.character(dist)) {
	dlist <- survreg.distributions[[dist]]
	if (is.null(dlist)) stop(paste(dist, ": distribution not found"))
	}
    else if (is.list(dist)) dlist <- dist
    else stop("Invalid distribution object")
    if (is.null(dlist$dist)) {
	if (is.character(dlist$name) && is.function(dlist$init) &&
	    is.function(dlist$deviance)) {}
	else stop("Invalid distribution object")
	}
    else {
	if (!is.character(dlist$name) || !is.character(dlist$dist) ||
	    !is.function(dlist$trans) || !is.function(dlist$dtrans))
		stop("Invalid distribution object")
	}	

    type <- attr(Y, "type")
    if (type== 'counting') stop ("Invalid survival type")
    
    logcorrect <- 0   #correction to the loglik due to transformations
    if (!is.null(dlist$trans)) {
	tranfun <- dlist$trans
	exactsurv <- Y[,ncol(Y)] ==1
	if (any(exactsurv)) logcorrect <- sum(logb(dlist$dtrans(Y[exactsurv,1])))

	if (type=='interval') {
	    if (any(Y[,3]==3))
		    Y <- cbind(tranfun(Y[,1:2]), Y[,3])
	    else Y <- cbind(tranfun(Y[,1]), Y[,3])
	    }
	else if (type=='left')
	     Y <- cbind(tranfun(Y[,1]), 2-Y[,2])
	else     Y <- cbind(tranfun(Y[,1]), Y[,2])
	if (!all(is.finite(Y))) 
	    stop("Invalid survival times for this distribution")
	}
    else {
	if (type=='left') Y[,2] <- 2- Y[,2]
	else if (type=='interval' && all(Y[,3]<3)) Y <- Y[,c(1,3)]
	}

    if (is.null(dlist$itrans)) itrans <- function(x) x
    else itrans <- dlist$itrans

    if (!is.null(dlist$scale)) {
	if (!missing(scale)) warning(paste(dlist$name, 
			   "has a fixed scale, user specified value ignored"))
	scale <- dlist$scale
	}
    if (!is.null(dlist$dist)) dlist <- survreg.distributions[[dlist$dist]]

    if (missing(control)) control <- if(!.R.)
      survReg.control(...) else survreg.control(...)

    if (scale < 0) stop("Invalid scale value")
    if (scale >0 && nstrata >1) 
	    stop("Cannot have multiple strata with a fixed scale")

    # Check for penalized terms
    pterms <- sapply(m, inherits, 'coxph.penalty')
    if (any(pterms)) {
	pattr <- lapply(m[pterms], attributes)
	# 
	# the 'order' attribute has the same components as 'term.labels'
	#   pterms always has 1 more (response), sometimes 2 (offset)
	# drop the extra parts from pterms
	temp <- c(attr(Terms, 'response'), attr(Terms, 'offset'))
	if (length(dropx)) temp <- c(temp, dropx+1)
	pterms <- pterms[-temp]
	temp <- match((names(pterms))[pterms], attr(Terms, 'term.labels'))
	ord <- attr(Terms, 'order')[temp]
	if (any(ord>1)) stop ('Penalty terms cannot be in an interaction')
	##pcols <- (attr(X, 'assign')[-1])[pterms]
        assign<-attrassign(X,newTerms)
        pcols<-assign[-1][pterms]
  
        fit <- survpenal.fit(X, Y, weights, offset, init=init,
				controlvals = control,
			        dist= dlist, scale=scale,
			        strata=strata, nstrat=nstrata,
				pcols, pattr,assign, parms=parms)
	}
    else fit <- if(!.R.) survReg.fit(X, Y, weights, offset, 
			    init=init, controlvals=control,
			    dist= dlist, scale=scale, nstrat=nstrata, 
			    strata, parms=parms) else
      survreg.fit(X, Y, weights, offset, 
			    init=init, controlvals=control,
			    dist= dlist, scale=scale, nstrat=nstrata, 
			    strata, parms=parms)

    # Next line: FEH added fitFunction='psm'
    if (is.character(fit))
      fit <- list(fail=fit, fitFunction='psm')  #error message
    else {
	if (scale==0) {
	    nvar <- length(fit$coef) - nstrata
	    fit$scale <- exp(fit$coef[-(1:nvar)])
	    if (nstrata==1) names(fit$scale) <- NULL
	    else names(fit$scale) <- levels(strata.keep)
	    fit$coefficients  <- fit$coefficients[1:nvar]
	    fit$idf  <- 1 + nstrata
	    }
	else {
	    fit$scale <- scale
	    fit$idf  <- 1
	    }
	fit$loglik <- fit$loglik + logcorrect
	}

    ## na.action <- attr(m, "na.action")
    ##if (length(na.action)) fit$na.action <- na.action
    if(length(nact)) fit$na.action <- nact  ## FEH
    fit$df.residual <- n - sum(fit$df)
#   fit$fitted.values <- itrans(fit$linear.predictors)
    fit$terms <- Terms
    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$means <- apply(X,2, mean)
    fit$call <- call
    fit$dist <- dist
    fit$df.resid <- n-sum(fit$df) ##used for anova.survreg
    if (model) fit$model <- m
    if (x)     fit$x <- X
##    if (y)     fit$y <- Y   FEH
    if (length(parms)) fit$parms <- parms
    ## Start FEH
    ##    if (any(pterms)) class(fit)<- c('survreg.penal', 'survreg')
    ##    else	     class(fit) <- 'survreg'
    fit$assign <- DesignAssign(atr, 1, Terms)
    fit$formula <- formula
    if(y) {
      oldClass(Y) <- 'Surv'
      attr(Y,'type') <- atY$type
      fit$y <- Y
    }
    if(.newSurvival.) scale.pred <-
      if(dist %in% c('weibull','exponential','lognormal','loglogistic'))
      c('log(T)','Survival Time Ratio') else 'T' else
    scale.pred <- if(substring(dist,1,3)=='log')
      c("log(T)","Survival Time Ratio") else "T"

    logtest <- 2*diff(fit$loglik)
    Nn <- if(length(weights)) sum(weights) else nnn[1]  ## 5jun02
    R2.max <- 1 - exp(2*fit$loglik[1]/Nn)
    R2 <- (1 - exp(-logtest/Nn))/R2.max
    df <- length(fit$coef)-1
    P <- if(df==0) NA else 1-pchisq(logtest,df)
    stats <- c(nnn, logtest, df, P, R2)
    names(stats) <- c("Obs", "Events", "Model L.R.", "d.f.", "P",
                      "R2")
    if(length(weights)) stats <- c(stats, 'Sum of Weights'=sum(weights))
    fit <- c(fit, list(stats=stats, maxtime=maxtime, units=time.units,
                       time.inc=time.inc, scale.pred=scale.pred,
                       non.slopes=1, Design=atr, fail=FALSE,
                       fitFunction=c("psm", "survreg", "glm", "lm")))
    if (any(pterms)) oldClass(fit) <- if(.SV4.)'Design' else
                        c('psm','Design','survreg.penal','survreg')
    else oldClass(fit) <- if(.SV4.)'Design' else c('psm','Design','survreg')
    ## End FEH
                       
    fit
    
  } else function(formula=formula(data), data,
	subset, na.action=na.delete, method="fit",
	link="log",
	dist=c("extreme", "logistic", "gaussian", "exponential","rayleigh","t"),
	init=NULL, fixed=list(), control,
	model=FALSE, x=FALSE, y=FALSE, time.inc, ...) {

    call <- match.call()
    dist <- match.arg(dist)
    
    m <- match.call(expand=FALSE)
    m$dist <- m$link <- m$model <- m$x <- m$y <- m$... <-  NULL
    m$init <- m$fixed <- m$control <- m$time.inc <- NULL
    m$na.action <- na.action
    m[[1]] <- as.name("model.frame")
    X <- Design(eval(m, sys.parent()))    # 24Apr01
    atrx <- attributes(X)                 # 24Apr01 + next 4
    nact <- atrx$na.action
    if(method=="model.frame") return(X)
    Terms <- atrx$terms
    offs <- offset <- attr(Terms, "offset")  ## offs 23nov02 moved 6dec02
    atr   <- atrx$Design
    s <- if(length(offs)) names(nact$nmiss) !=  atrx$names[offs] else TRUE
    ## 23nov02
    if(length(nact$nmiss))
      names(nact$nmiss)[s] <- 
        c(as.character(formula[2]), atr$name[atr$assume.code!=9])
    ## [s] 23nov02
    lnames <- if(.R.) c("logit","probit","cloglog","identity","log","sqrt",
      "1/mu^2","inverse") else dimnames(glm.links)[[2]]
    link <- pmatch(link, lnames, 0)
    if(link==0) stop("invalid link function")
    link <- lnames[link]
    Y <- model.extract(X, "response")
    atY <- attributes(Y)    # 1 Apr 95
    ncy <- ncol(Y)
    maxtime <- max(Y[,-ncy])
    nnn <- c(nrow(Y),sum(Y[,ncy]))
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    if(model) m <- X

    if(length(offset)) offset <- as.numeric(X[[offset]])
    else offset <- rep(0, nrow(X))

    X <- model.matrix(Terms, X)

    time.units <- attr(Y, "units")
    if(!length(time.units)) time.units <- "Day"
    if(missing(time.inc))	{
      time.inc <- switch(time.units,Day=30,Month=1,Year=1,maxtime/10)
      if(time.inc>=maxtime | maxtime/time.inc>25) time.inc <- max(pretty(c(0,maxtime)))/10
      				}
    rnam <- dimnames(Y)[[1]]
    dimnames(X) <- list(rnam, c("(Intercept)",atr$colnames))
    if(method=="model.matrix") return(X)
    n <- nrow(X)
    nvar <- ncol(X)

    type <- attr(Y, "type")
    linkfun <- if(.R.) make.link(link)$linkfun else
                       glm.links["link", link][[1]]
    if (type== 'counting') stop ("Invalid survival type")
    else if (type=='interval') {
	if (any(Y[,3]==3))
	     Y <- cbind(linkfun(Y[,1:2]), Y[,3])
	else Y <- cbind(linkfun(Y[,1]), Y[,3])
	}
    else if (type=='left')
	     Y <- cbind(linkfun(Y[,1]), 2-Y[,2])
    else     Y <- cbind(linkfun(Y[,1]), Y[,2])

    controlvals <- survreg.control()
    if(!missing(control)) controlvals[names(control)] <- control

    if(dist=="exponential") { fixed$scale <- 1; dist <- "extreme" }
    else if (dist=="rayleigh") { fixed$scale <- .5; dist <- "extreme" }

    sd <- survreg.distributions[[dist]]
    if (length(fixed)>0) {
	ifix <- match(names(fixed), names(sd$parms), nomatch=0)
	if (any(ifix==0))
	    stop (paste("Parameter(s)", paste(names(fixed)[ifix==0]),
			"in the fixed list not valid for this dist"))
	}
    if (is.list(init) && length(init)>0) {
	ifix <- match(names(init), c('eta',names(sd$parms)), nomatch=0)
	if (any(ifix==0))
	    stop (paste("Parameter(s)", paste(names(init)[ifix==0]),
			"in the init list not valid for this dist"))
	}

    sfit <- if(.R.) survreg.fit(X, Y, offset=offset, init=init,
                                controlvals=controlvals,
                                dist=dist, parms=fixed) else
    survreg.fit(X, Y, offset, init=init, controlvals=controlvals,
                dist= dist, fixed=fixed)

    if (is.character(sfit))  	{
	cat("Failure in psm:\n",sfit,"\n")
	fit <- list(fail=TRUE, fitFunction='psm')  #error message; 14Nov00
	oldClass(fit) <- if(.SV4.)'Design' else "psm" ##14Nov00
	return(fit)
				}
    else {
	# There may be more clever ways to do this, but ....
	#  In order to make it look like IRLS, and get all the components
	#  that I need for glm inheritance, do one step of weighted least
	#  squares.
	eta <- c(X %*% sfit$coef[1:nvar]) + offset
	wt<- -sfit$deriv[,3]
	fit <- lm.wfit(X, eta + sfit$deriv[,2]/wt, wt, "qr", ...)

    ifun <- if(.R.) make.link(link)$linkinv else
                    glm.links["inverse",link][[1]]
    fit$fitted.values <- ifun(fit$fitted.values)
	fit$family <- c(name=dist, link=link, "")
	fit$linear.predictors <- eta
	fit$iter <- sfit$iter
    fit$parms <- sfit$parms
    fit$df.residual <- fit$df.residual-sum(!sfit$fixed)

	# If singular.ok=T, there may be NA coefs.  The var matrix should
	#   be an inversion of the "non NA" portion.

	var <- 0*sfit$imat
	good <- c(!is.na(fit$coef), rep(TRUE, ncol(var)-nvar))
	var[good,good] <- solve(qr(sfit$imat[good,good], tol=1e-12))
	fit$var <- var
	fit$fixed <- sfit$fixed
	dev <- sd$deviance(Y, fit$parms, sfit$deriv[,1])
	fit$dresiduals <- sign(fit$residuals)*sqrt(dev)
	fit$deviance <- sum(dev)
	fit$null.deviance <- fit$deviance +2*(sfit$loglik[2]- sfit$ndev[2])
	fit$loglik <- c(sfit$ndev[2], sfit$loglik[2])
	}

    if (length(nact)) fit$na.action <- nact

    i <- 1:nvar
    var <- var[i,i,drop=FALSE]	#omit scale row and col.
    fit$se.fit <-  drop(sqrt(((X %*% var) * X) %*% rep(1,nvar)))

    logtest <- fit$null.deviance - fit$deviance
    R2.max <- 1 - exp(-fit$null.deviance/nnn[1])
    R2 <- (1 - exp(-logtest/nnn[1]))/R2.max
    df <- length(fit$coef)-1
    P <- if(df==0) NA else 1-pchisq(logtest,df)
    stats <- c(nnn, logtest, df, P, R2)
    names(stats) <- c("Obs", "Events", "Model L.R.", "d.f.", "P", "R2")
    scale.pred <- if(link=="log") c("log(T)","Survival Time Ratio") else "T"
    fit <- c(fit, list(maxtime=maxtime, units=time.units,
                       time.inc=time.inc,scale.pred=scale.pred,
                       non.slopes=1,
                       fitFunction=c("psm", "survreg", "glm", "lm"))) ##13Nov00
    fit$Design <- atr
    fit$stats  <- stats
    oldClass(fit) <- if(.SV4.)'Design' else
     c("psm", "Design", "survreg", "glm", "lm")  ##14Nov00
    fit$terms <- Terms
    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$call <- call
    fit$fail <- FALSE
    if (model) fit$model <- m
    if (x)     fit$x <- X
    if (y) {
      oldClass(Y) <- 'Surv'    # 1 Apr 95
      attr(Y,'type') <- atY$type
      fit$y <- Y
    }
    fit
    }

Hazard   <- function(fit, ...) UseMethod("Hazard")
Survival <- function(fit, ...) UseMethod("Survival")

Hazard.psm <- if(.newSurvival.) function(fit, ...) {
 dist <- fit$dist
 g <- survreg.auxinfo[[dist]]$hazard
 formals(g) <- list(times=NA, lp=NULL, parms=logb(fit$scale))
 g
} else function(fit) {
  fam <- fit$family
  dist <- fam["name"]
  transform <- fam[2]
  g <- survreg.auxinfo[[dist]]$hazard
  formals(g) <- list(times=NULL, lp=NULL,
                     parms=fit$parms, transform=transform)
  g
}

Survival.psm <- if(.newSurvival.) function(fit, ...) {
 dist <- fit$dist
 g <- survreg.auxinfo[[dist]]$survival
 formals(g) <- list(times=NULL, lp=NULL, parms=logb(fit$scale))
 g
} else function(fit) {
  fam <- fit$family
  dist <- fam["name"]
  transform <- fam[2]
  g <- survreg.auxinfo[[dist]]$survival
  formals(g) <- list(times=NULL, lp=NULL,
                     parms=fit$parms, transform=transform)
  g
}

Quantile.psm <- if(.newSurvival.) function(fit, ...) {
  dist <- fit$dist
  g <- survreg.auxinfo[[dist]]$Quantile
  formals(g) <- list(q=.5, lp=NULL, parms=logb(fit$scale))
  g
} else function(fit, ...) {
  fam <- fit$family
  dist <- fam["name"]
  transform <- fam[2]
  g <- survreg.auxinfo[[dist]]$quantile
  formals(g) <- list(q=.5, lp=NULL,
                     parms=fit$parms, transform=transform)
  g
}

Mean.psm <- if(.newSurvival.) function(fit, ...) {
 dist <- fit$dist
 g <- survreg.auxinfo[[dist]]$mean
 formals(g) <- list(lp=NULL, parms=logb(fit$scale))
 g
} else function(fit, ...) {
  fam <- fit$family
  dist <- fam["name"]
  transform <- fam[2]
  g <- survreg.auxinfo[[dist]]$mean
  formals(g) <- list(lp=NULL, parms=fit$parms,
                     transform=transform)
  g
}

predict.psm <- 
  function(object, newdata,
           type=c("lp","x","data.frame","terms","adjto","adjto.data.frame",
             "model.frame"),
           se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
           incl.non.slopes, non.slopes, kint=1,
           na.action=na.keep, expand.na=TRUE, center.terms=TRUE, ...)
  predictDesign(object, newdata, type, se.fit, conf.int, conf.type,
                incl.non.slopes, non.slopes, kint,
                na.action, expand.na, center.terms, ...)



residuals.psm <- function(object, type = "censored.normalized", ...)
{
    type <- match.arg(type)
    if(type!='censored.normalized') {
      if(type=='score' && (.newSurvival.))
        stop('score residuals not implemented')
      ## TODO
      return(if(!.R.)residuals.survReg(object, type=type) else
        residuals.survreg(object, type=type))
    }

    y <- object$y
    aty <- attributes(y)
    if(length(y)==0) stop('did not use y=T with fit')
    ncy <- ncol(y)
    if(.newSurvival.) {  ## 17Apr02
      scale <- object$scale
      dist  <- object$dist
    } else {
      scale <- exp(object$parms)
      dist  <- object$family[1]
    }
    r <- (y[,-ncy,drop=FALSE]-object$linear.predictors)/scale
    r <- cbind(r, y[,ncy])
    ## Moved the following line here from bottom
    if(length(object$na.action)) r <- naresid(object$na.action, r)
    attr(r,'dist') <- dist
    attr(r,'type') <- aty$type
    attr(r,'units') <- ' '
    attr(r,'time.label') <- 'Normalized Residual'
    attr(r,'event.label') <- aty$event.label
    oldClass(r) <- c('residuals.psm.censored.normalized','Surv')
    g <- survreg.auxinfo[[dist]]$survival
    formals(g) <- if(.newSurvival.) list(times=NULL, lp=0, parms=0)
    else list(times=NULL, lp=0, parms=0, transform='identify')  ## 17Apr02
    attr(r,'theoretical') <- g
    r
    }

lines.residuals.psm.censored.normalized <- 
  function(x, n=100, lty=1, xlim=range(r[,-ncol(r)],na.rm=TRUE),
           lwd=3, ...) {
    r <- x
	x <- seq(xlim[1], xlim[2], length=n)
    tx <- x
    if(.newSurvival.) {
      dist <- attr(r, 'dist')
      if(dist %in% c('weibull','loglogistic','lognormal')) tx <- exp(x)
      ## $survival functions log x
    }
	lines(x, attr(r,'theoretical')(tx), lwd=lwd, lty=lty)
	invisible()
  }

survplot.residuals.psm.censored.normalized <- 
  function(fit, x, g=4, col, main, ...) {
    r <- fit
    
  if(missing(x)) {
	survplot(survfit(r), conf='none', xlab='Residual', 
			 col=if(missing(col))par('col') else col, ...)
    ## was survfit(r, data=list(r)) 28apr02
	if(!missing(main)) title(main)
  } else {
	if(is.character(x)) x <- as.factor(x)
	if(!is.category(x) && length(unique(x))>5) x <- cut2(x, g=g)
	s <- is.na(r[,1]) | is.na(x)
	if(any(s)) {r <- r[!s,]; x <- x[!s,drop=TRUE]}
	survplot(survfit(r ~ x, data=data.frame(x,r)),  xlab='Residual',
			 conf='none',
			 col=if(missing(col))1:length(levels(x)) else par('col'), ...)
	if(missing(main)) main <- if(length(lab <- attr(x,'label'))) lab 
	  else if(.R.) '' else deparse(substitute(x))
	if(main != '') title(main)
  }
  lines(r, lty=1, lwd=3)
  invisible()
}



#SCCS @(#)residuals.coxph.s	4.5 7/14/92
residuals.cph <-
  function(object, type=c("martingale", "deviance", "score", "schoenfeld",
	"dfbeta","dfbetas", "scaledsch"), collapse=FALSE, weighted=FALSE, ...)
    {
    type <- match.arg(type)
    if(.R.) require('survival')
    otype <- type
    if(type=="dfbeta" | type=="dfbetas") type <- "score"
    if(type=="scaledsch") type <- "schoenfeld"
    n <- length(object$residuals)
    rr <- object$residual
    y <- object$y
    x <- object$x
	weights <- object$weights
    strat <- attr(x,"strata")
    method <- object$method
    if (method=='exact' && (type=='score' || type=='schoenfeld'))
	stop(paste(type, 'residuals are not available for the exact method'))
    if(type!="martingale" & is.null(x))
      stop("you must specify x=T in the fit")
    if(type!="deviance" & type!="martingale" & is.null(y))
      stop("you must specify y=T in the fit")

    if(type!="martingale")				{

	ny <- ncol(y)
	status <- y[,ny,drop=TRUE]

	if (type != 'deviance') {
	    nvar <- ncol(x)
	    if (is.null(strat)) {
		ord <- order(y[,ny-1], -status)
		newstrat <- rep(0,n)
		}
	    else {
		ord <- order(strat, y[,ny-1], -status)
		newstrat <- c(diff(as.numeric(strat[ord]))!=0 ,1)
		}
	    newstrat[n] <- 1

	    # sort the data
	    x <- x[ord,]
	    y <- y[ord,]
	    score <- exp(object$linear.predictor)[ord]
		weights <- if(length(weights)) weights[ord] else rep(1,n)

	    if (ny ==3) subs <- paste("agres", 1:2, sep='')
	    else        subs <- paste("coxres",1:2, sep='')
	    }
	}

    #
    # Now I have gotten the data that I need-- do the work
    #
    if (type=='schoenfeld') {
	if (ny==2)  {
      mintime <- min(y[, 1])
      if(mintime < 0) y <- cbind(2 * mintime - 1, y)
      else y <- cbind(-1, y)
    }
	temp <- if(.R.)
      .C("coxscho",
         n=as.integer(n),
         as.integer(nvar),
         as.double(y),
         resid= x,
         score * weights,
         as.integer(newstrat),
         as.integer(method=='efron'),
         double(3*nvar), PACKAGE="survival") else
      .C(if(.SV4.)'S_coxscho' else "coxscho",  ##14Nov00
         n=as.integer(n),
         as.integer(nvar),
         as.double(y),
         resid= x,
         score * weights,
         as.integer(newstrat),
         as.integer(method=='efron'),
         double(3*nvar))

	deaths <- y[,3]==1

	if (nvar==1) rr <- temp$resid[deaths]
	else rr <- matrix(temp$resid[deaths,], ncol=nvar) #pick rows, and kill attr

	if (length(object$strata)) 
	  attr(rr, "strata")  <- table((strat[ord])[deaths])
	time <- c(y[deaths,2])  # 'c' kills all of the attributes
	if (is.matrix(rr)) dimnames(rr)<- list(time, names(object$coef))
	else               names(rr) <- time

	if (otype=='scaledsch') {
	  ndead <- sum(deaths)
	  coef <- ifelse(is.na(object$coef), 0, object$coef)
	    if (nvar==1) rr <- rr*object$var * ndead + coef
	    else         rr <- rr %*% object$var * ndead + 
		  outer(rep(1,nrow(rr)), coef)
	    }
	return(rr)
	}

    if (type=='score') {
	if (ny==2) {
	    resid <- if(.R.)
          .C("coxscore",
             as.integer(n),
             as.integer(nvar),
             as.double(y),
             x=as.double(x),
             as.integer(newstrat),
             as.double(score),
             as.double(weights),
             as.integer(method=='efron'),
             resid= double(n*nvar),
             double(2*nvar), PACKAGE="survival")$resid else
        .C(if(.SV4.)'S_coxscore' else "coxscore",
           as.integer(n),
           as.integer(nvar),
           as.double(y),
           x=as.double(x),
           as.integer(newstrat),
           as.double(score),
           as.double(weights),
           as.integer(method=='efron'),
           resid= double(n*nvar),
           double(2*nvar))$resid
	    }
	else {
	    resid <- if(.R.)
          .C("agscore",
             as.integer(n),
             as.integer(nvar),
             as.double(y),
             as.double(x),
             as.integer(newstrat),
             as.double(score),
             as.double(weights),
             as.integer(method=='efron'),
             resid=double(n*nvar),
             double(nvar*6), PACKAGE="survival")$resid
        .C(if(.SV4.)'S_agscore' else "agscore",
           as.integer(n),
           as.integer(nvar),
           as.double(y),
           as.double(x),
           as.integer(newstrat),
           as.double(score),
           as.double(weights),
           as.integer(method=='efron'),
           resid=double(n*nvar),
           double(nvar*6))$resid
      }
	if (nvar >1) {
	    rr <- matrix(0, n, nvar)
	    rr[ord,] <- matrix(resid, ncol=nvar)
	    dimnames(rr) <- list(names(object$resid), names(object$coef))
	    }
	else rr[ord] <- resid
	}

    #Expand out the missing values in the result
    if (!is.null(object$na.action)) {
	rr <- naresid(object$na.action, rr)
	if (is.matrix(rr)) n <- nrow(rr)
	else               n <- length(rr)
	if (type=='deviance') status <- naresid(object$na.action, status)
	}

    # Collapse if desired
    if (!missing(collapse)) {
	if (length(collapse) !=n) stop("Wrong length for 'collapse'")
	rr <- rowsum(rr, collapse)
	}

    # Deviance residuals are computed after collapsing occurs
    if (type=='deviance')
	rr <- sign(rr) *sqrt(-2* (rr+
			      ifelse(status==0, 0, status*logb(status-rr))))

    if      (otype=='dfbeta') rr %*% object$var
    else if (otype=='dfbetas') (rr %*% object$var) %*% 
		diag(sqrt(1/diag(object$var)))
    else  rr
    }
residuals.lrm <- function(object, 
  type=c("ordinary","score","score.binary","pearson",
         "deviance","pseudo.dep","partial",
         "dfbeta","dfbetas","dffit","dffits","hat","gof","lp1"),
  pl=FALSE, xlim, ylim, kint=1, label.curves=TRUE, 
  which, ...) {

gotsupsmu <- FALSE
type <- match.arg(type)
dopl <- (is.logical(pl) && pl) || is.character(pl)

k <- object$non.slopes
L <- object$linear.predictors
if(length(L)==0) stop('you did not use linear.predictors=T for the fit')

if(kint<1 | kint>k) stop(paste('kint must be from 1-',k,sep=''))

cof <- object$coef
ordone <- type %in% c('partial','gof','score','score.binary')
 # residuals explicitly handled for ordinal model
if(ordone && !missing(kint)) 
  stop('may not specify kint for partial, score, score.binary, or gof')

if(k>1 && kint!=1 && !ordone) L <- L - cof[1] + cof[kint]
P <- 1/(1+exp(-L))
if(length(Y <- object$y)==0) stop("you did not specify y=T in the fit")
rnam <- names(Y)
cnam <- names(cof)
if(!is.factor(Y)) Y <- factor(Y)  ## 11Apr02
lev  <- levels(Y)
lev2 <- names(cof)[1:k]
Y <- oldUnclass(Y) - 1
if(!ordone && k>1) Y <- Y >= kint
if(k>1 && missing(kint) && !ordone)
  warning(paste('using first intercept and ',
				lev2[kint],
				' to compute residuals or test GOF',
				sep=''))

if(type=="gof") {
  if(length(X <- object$x)==0)
    stop("you did not use x=T in the fit")
  stats <- matrix(NA, nrow=k, ncol=5, dimnames=list(if(k>1)lev2,
                  c("Sum of squared errors","Expected value|H0","SD","Z","P")))
  X <- cbind(1,X)
  for(j in 1:k) {
    y <- Y>=j
    p <- 1/(1+exp(-(L-cof[1]+cof[j])))
    sse <- sum((y-p)^2)
    wt <- p*(1-p)
    d <- 1-2*p
    z <- lm.wfit(X,d,wt,method='qr')
##    res <- summary(lm.wfit(X, d, wt, method="qr"))$residuals  11Apr02
    res <- z$residuals * sqrt(z$weights)
    ## Was without the summary( ).  Thanks to jorge.sirgo@us.pwcglobal.com
    ## Variance was too big    19Mar01
    sd <- sqrt(sum(res^2))
    ev <- sum(wt)
    z <- (sse-ev)/sd
    P <- 2*(1-pnorm(abs(z)))
    stats[j,] <- c(sse,ev,sd,z,P)
  }
  return(drop(stats))
}

naa <- object$na.action

if(type=="ordinary") return(naresid(naa,Y - 1/(1+exp(-L))))

if(type %in% c('score','score.binary','partial')) {
  nc <- length(cof)
  if(missing(which)) which <- if(type=='score')1:nc else 1:(nc-k) else
  if(type=='score') which <- k+which
}

if(type=='score' || type=='score.binary')
  plotit <- function(w, ylim, xlab, ylab, lev=names(w)) {
     statsum <- function(x) {
       n <- length(x)
       xbar <- sum(x)/n
       if(n<2) {low <- hi <- NA} else {
         se   <- 1.96*sqrt(sum((x-xbar)^2)/(n-1)/n)
         low  <- xbar-se; hi <- xbar+se
       }
       c(mean=xbar, lower=low, upper=hi)
     }
     k <- length(w)
     w <- lapply(w, statsum)
     plot(c(1,k), c(0,0), xlab=xlab, ylab=ylab,
          ylim=if(length(ylim)==0) range(unlist(w)) else ylim, type='n', axes=FALSE)
     mgp.axis(2)
     mgp.axis(1, at=1:k, labels=lev)
     abline(h=0, lty=2, lwd=1)
     ii <- 0
     for(ww in w) {
       ii <- ii+1
       points(ii, ww[1])
       errbar(ii, ww[1], ww[3], ww[2], add=TRUE)
     }
    }

if(type=='score.binary') {
  if(k==1) stop('score.binary only applies to ordinal models')
  if(!dopl) stop('score.binary only applies if you are plotting')
  if(!length(X <- oldUnclass(object$x))) stop('you did not specify x=T for the fit')
  xname <- dimnames(X)[[2]]
  yname <- as.character(formula(object))[2]
  for(i in which) {
    xi <- X[,i]
    r <- vector('list',k)
    names(r) <- lev[-1]
    for(j in 1:k) {
      r[[j]] <- xi*((Y>=j)-1/(1+exp(-(L-cof[1]+cof[j]))))
    }
    if(pl!='boxplot') plotit(r, ylim=if(missing(ylim))NULL else ylim,
      xlab=yname, ylab=xname[i]) else
    boxplot(r, varwidth=TRUE, notch=TRUE, err=-1,
            ylim=if(missing(ylim))quantile(unlist(r),c(.1,.9)) else ylim, ...)
    title(xlab=yname, ylab=xname[i])
  }
  invisible()
}

if(type=="score")						{
  if(!length(X <- oldUnclass(object$x)))
    stop("you did not specify x=T for the fit")
  if(k==1) return(naresid(naa,cbind(1,X)*(Y-P)))
  z <- function(i,k,L,coef) 1/(1+exp(-(coef[pmin(pmax(i,1),k)]+L)))
  # Mainly need the pmax - 0 subscript will drop element from vector
#  z$k <- k; z$L <- L-cof[1]; z$coef <- cof
  formals(z) <- list(i=NA,k=k,L=L-cof[1],coef=cof)
  ## set defaults in fctn def'n
  u <- matrix(NA, nrow=length(L), ncol=length(which),
              dimnames=list(names(L),names(cof)[which]))
  pq <- function(x) x*(1-x)
  # Compute probabilities of being in observed cells
  pc <- ifelse(Y==0, 1-z(1), ifelse(Y==k, z(k), z(Y)-z(Y+1)) )
  xname <- dimnames(X)[[2]]
  yname <- as.character(formula(object))[2]
  ii <- 0
  for(i in which) {
	ii <- ii + 1
    di  <- if(i<=k) ifelse(Y==0, if(i==1) 1 else 0, Y==i) else X[,i-k]
    di1 <- if(i<=k) ifelse(Y==0 | Y==k, 0, (Y+1)==i)      else X[,i-k]
    ui  <- ifelse(Y==0, -z(1)*di,
                  ifelse(Y==k, (1-z(k))*di,
                         (pq(z(Y))*di-pq(z(Y+1))*di1)/pc ) )
   u[,ii] <- ui
   if(dopl && i>k) {
    if(pl=='boxplot') {
      boxplot(split(ui, Y), varwidth=TRUE, notch=TRUE, names=lev, err=-1, 
              ylim=if(missing(ylim))quantile(ui,c(.1,.9)) else ylim, ...)
      title(xlab=yname, ylab=paste('Score Residual for',xname[i-k]))
    }
    else plotit(split(ui,Y),ylim=if(missing(ylim))NULL else ylim,lev=lev,
                xlab=yname, ylab=paste('Score Residual for',xname[i-k]))
   }
  }
  return(if(dopl)invisible(naresid(naa, u)) else naresid(naa, u))
								}

if(type=="pearson")						{
  return(naresid(naa,(Y-P)/sqrt(P*(1-P))))
								}

if(type=="deviance")						{
  r <- ifelse(Y==0,-sqrt(2*abs(logb(1-P))),sqrt(2*abs(logb(P))))
  return(naresid(naa,r))
								}

if(type=="pseudo.dep")						{
  r <- L + (Y-P)/P/(1-P)
  return(naresid(naa,r))
								}

if(type=="partial")						{
  if(!length(X <- oldUnclass(object$x)))
    stop("you did not specify x=T in the fit")
  cof.int <- cof[1:k]
  cof     <- cof[-(1:k)]
  if(!missing(which)) {
	X <- X[,which,drop=FALSE]
	cof <- cof[which]
  }
  atx <- attributes(X)
  dx <- atx$dim
  if(k==1) r <- cof.int[1]+X*matrix(cof,nrow=dx[1],ncol=dx[2],byrow=TRUE) +
	 (Y-P)/P/(1-P) else {
		 r <- X*matrix(cof, nrow=dx[1], ncol=dx[2], byrow=TRUE)
		 R <- array(NA, dim=c(dx,k), dimnames=c(atx$dimnames, list(lev2)))
		 for(j in 1:k) {
		   y <- Y>=j
		   p <- 1/(1+exp(-(L-cof.int[1]+cof.int[j])))
		   R[,,j] <- r + (y-p)/p/(1-p)
		 }
	 }
  if(dopl)
  {
	xname <- atx$dimnames[[2]]; X <- oldUnclass(X)
	for(i in 1:dx[2])
	{
	   if(pl=="loess")
	   {
		  if(k>1) stop('pl="loess" not implemented for ordinal response')
          xi <- X[,i]   #17Apr01
          ri <- r[,i]
          ddf <- data.frame(xi,ri)
          if(.R.) {
            plot(xi, ri, xlim=if(missing(xlim)) range(xi) else xlim,
                 ylim=if(missing(ylim)) range(ri) else ylim,
                 xlab=xname[i], ylab='Partial Residual')
            lines(lowess(xi,ri))
          } else {
            g <- loess(ri ~ xi, data=ddf)
            plot(g, coverage=0.95, confidence=7, xlab=xname[i],
                 ylab="Partial Residual", 
                 ylim=if(missing(ylim))range(ri) else ylim,
                 xlim=if(missing(xlim))range(xi) else xlim)
            points(xi, ri)
          }
	   }
	   else if(k==1)
	   {
	      xi <- X[,i]; ri <- r[,i]
	      plot(xi, ri, xlab=xname[i], ylab="Partial Residual",
		   xlim=if(missing(xlim))range(xi) else xlim,
		   ylim=if(missing(ylim))range(ri) else ylim)
          if(.R. && !gotsupsmu && pl!='lowess') {
            require('modreg')
            gotsupsmu <- TRUE
          }
	      if(pl=="lowess") lines(lowess(xi, ri, iter=0, ...))
	      else lines(supsmu(xi, ri, ...))
	   }
	   else {
		 xi <- X[,i]
		 ri <- R[,i,,drop=TRUE]
		 smoothed <- vector('list',k)
		 ymin <- 1e30; ymax <- -1e30

         if(.R. && !gotsupsmu && pl!='lowess') {
           require('modreg')
           gotsupsmu <- TRUE
         }
		 for(j in 1:k) {
		   w <- if(pl!='supsmu') lowess(xi, ri[,j], iter=0, ...) else
		     supsmu(xi, ri[,j], ...)
		   ymin <- min(ymin,w$y)
		   ymax <- max(ymax,w$y)
		   smoothed[[j]] <- w
		 }
		 plot(0, 0, xlab=xname[i], ylab='Partial Residual',
			  xlim=if(missing(xlim))range(xi) else xlim,
			  ylim=if(missing(ylim))range(pretty(c(ymin,ymax))) else ylim,
                          type='n')
		 us <- par('usr')[1:2]
		 for(j in 1:k) {
		   w <- smoothed[[j]]
		   lines(w, lty=j)
		   if(is.character(label.curves)) {
		     xcoord <- us[1]+(us[2]-us[1])*j/(k+1)
		     text(xcoord, approx(w, xout=xcoord, rule=2)$y, lev2[j])
		   }
		 }
		 if(is.list(label.curves) || 
			(is.logical(label.curves) && label.curves))
		   labcurve(smoothed, lev2, opts=label.curves)
	   }
	}
	return(invisible(if(k==1)naresid(naa,r) else R))   
  }
  return(if(k==1) naresid(naa,r) else R)
								}

##if(type=='convexity') {
##  if(missing(p.convexity)) {
##	pq <- quantile(P, c(.01, .99))
##	if(pq[1]==pq[2]) pq <- range(P)
##	p.convexity <- seq(pq[1], pq[2], length=100)
## }
##  lcp <- length(p.convexity)
##  cp <- single(lcp)
##  for(i in 1:lcp) {
##	p <- p.convexity[i]
##	cp[i] <- mean(((p/P)^Y)*(((1-p)/(1-P))^(1-Y)))
##  }
##  if(pl) plot(p.convexity, cp, xlab='p', ylab='C(p)', type='l')
##  return(invisible(cp))
##}
  
if(type=="dfbeta"|type=="dfbetas"|type=="dffit"|type=="dffits"|type=="hat"|
   type=="lp1")
{
  if(length(X <- oldUnclass(object$x))==0)
    stop("you did not specify x=T for the fit")
  v <- P*(1-P)
  if(.R.) {
    g <- lm(L+(Y-P)/v ~ X, weights=v)
    infl <- lm.influence(g)
    dfb <- coef(infl)    ## R already computed differences
  } else {
    xx <- cbind(1,X)
    g <- lm.wfit(xx, L+(Y-P)/v, v, method="qr", qr=TRUE)
    g$x <- xx
    infl <- lm.influence(g)
    dfb <- t(coef(g) - t(coef(infl)))
  }
  dimnames(dfb) <- list(rnam, c(cnam[kint],cnam[-(1:k)]))
  if(type=="dfbeta") return(naresid(naa,dfb))
  if(type=="dfbetas") {
    i <- c(kint, (k+1):length(cof))
    return(naresid(naa,sweep(dfb,2,diag(object$var[i,i])^.5,"/")))
  }
  if(type=="hat") return(naresid(naa,infl$hat))
  if(type=="dffit") return(naresid(naa,
	(infl$hat * g$residuals)/(1 - infl$hat)))
  if(type=="dffits") return(naresid(naa,
	(infl$hat^.5)*g$residuals/(infl$sigma*(1-infl$hat)) ))
  if(type=="lp1") return(naresid(naa,
        L - (infl$hat * g$residuals)/(1 - infl$hat)))
}
}


plot.lrm.partial <- function(..., labels, center=FALSE) {

dotlist <- list(...)
nfit <- length(dotlist)
if(missing(labels)) labels <- (as.character(sys.call())[-1])[1:nfit]

vname <- dimnames(dotlist[[1]]$x)[[2]]
nv <- length(vname)
if(nv==0) stop('you did not specify x=T on the fit')

r <- vector('list', nv)
for(i in 1:nfit) r[[i]] <- resid(dotlist[[i]], 'partial')

for(i in 1:nv) {
  curves <- vector('list',nfit)
  ymin <- 1e10; ymax <- -1e10
  for(j in 1:nfit) {
	xx <- dotlist[[j]]$x[,vname[i]]
	yy <- r[[j]][,vname[i]]
	if(center)yy <- yy - mean(yy)
	curves[[j]] <- lowess(xx, yy, iter=0)
	ymin <- min(ymin, curves[[j]]$y)
	ymax <- max(ymax, curves[[j]]$y)
  }
  for(j in 1:nfit) {
	if(j==1) plot(curves[[1]], xlab=vname[i], ylab='Partial Residual',
		 ylim=c(ymin, ymax), type='l')
	else lines(curves[[j]], lty=j)
  }
  if(nfit>1) labcurve(curves, labels)
}
invisible()
}

residuals.ols <-
  function(object, 
           type=c("ordinary","score","dfbeta","dfbetas","dffit","dffits","hat",
             "hscore"), ...)
{

type <- match.arg(type)
naa <- object$na.action

if(type=="ordinary") return(naresid(naa, object$residuals))

if(!length(object$x))stop("did not specify x=T in fit")

if(type=="score") return(naresid(naa, object$x*object$residuals))

infl <- ols.influence(object)

if(type=="hscore") return(naresid(naa, object$x *
     (object$residuals/(1-infl$hat))))

if(type=="dfbeta"|type=="dfbetas")
{
   r <- t(coef(object) - t(coef(infl)))
   if(type=="dfbetas") r <- sweep(r,2,diag(object$var)^.5,"/")
}
else if(type=="dffit") r <- (infl$hat * object$residuals)/(1 - infl$hat)
else if(type=="dffits") r <- (infl$hat^.5)*object$residuals/
	      (infl$sigma*(1-infl$hat))
else if(type=="hat") r <- infl$hat

naresid(naa, r)
									}

## lm.influence used to work but now it re-computes X for unknown
## reasons  24Nov00
ols.influence <- function(lm, x) {
GET <- function(x, what)
  {
    ## eventually, x[[what, exact=TRUE]]
    if(is.na(n <- match(what, names(x)))) NULL else x[[n]]
  }
wt <- GET(lm, "weights")
## should really test for < 1/BIG if machine pars available
e <- lm$residuals
n <- length(e)
if(length(wt))
  e <- e * sqrt(wt)
beta <- lm$coef
if(is.matrix(beta)) {
  beta <- beta[, 1]
  e <- e[, 1]
  warning("multivariate response, only first y variable used")
}
na <- is.na(beta)
beta <- beta[!na]
p <- GET(lm, "rank")
if(!length(p))
    p <- sum(!na)
R <- if(.R.) lm$qr$qr else lm$R
if(p < max(dim(R)))
  R <- R[1:p, 1:p]
qr <- GET(lm, "qr")
if(!length(qr)) {
  lm.x <- GET(lm, "x")
  if(length(wt))
    lm.x <- lm.x * sqrt(wt)
  if(any(na))
    lm.x <- lm.x[, !na, drop = FALSE]
  Q <- left.solve(R, lm.x)
}
else {
  if(length(wt) && any(zero <- wt == 0)) {
    Q <- matrix(0., n, p)
    dimnames(Q) <- list(names(e), names(beta))
    Q[!zero,  ] <- qr.Q(qr)[, 1:p, drop = FALSE]
  }
  else {
    Q <- qr.Q(qr)
    if(p < ncol(Q))
      Q <- Q[, 1:p, drop = FALSE]
  }
}
h <- as.vector((Q^2 %*% array(1, c(p, 1))))
h.res <- (1 - h)
z <- e/h.res
v1 <- e^2
z <- t(Q * z)
v.res <- sum(v1)
v1 <- (v.res - v1/h.res)/(n - p - 1)
# BKW (2.8)
dbeta <- backsolve(R, z)
list(coefficients = t(beta - dbeta), sigma = sqrt(v1), hat = h)
}
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

if(keep.prop) storage.mode(Propensity) <- 'single'
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



robcov <- function(fit, cluster, method=c('huber','efron')) {

  method <- match.arg(method)

var <- fit$var
vname <- dimnames(var)[[1]]

if(inherits(fit, "ols") ||
   (length(fit$fitFunction) && any(fit$fitFunction=='ols'))) {
   ##14Nov00 22May01
   var <- fit$df.residual * var/sum(fit$residuals^2)  #back to X'X
##   warning("printing the fit object from robcov (with print.ols) will not print the adjusted std. err.\nUse sqrt(diag(fit$var)) to get adjusted std. err.,\nwhere fit is the result of robcov.") 22Dec01
 } else if(method=='efron') stop('method="efron" only works for ols fits')

X <- as.matrix(residuals(fit, type=if(method=='huber')"score" else "hscore"))

n <- nrow(X)
if(missing(cluster)) cluster <- 1:n
else if(any(is.na(cluster))) stop("cluster contains NAs")
if(length(cluster)!=n)
 stop("length of cluster does not match number of observations used in fit")
cluster <- as.factor(cluster)

p <- ncol(var)
j <- is.na(X %*% rep(1, p))
if(any(j))	{
  X <- X[!j,,drop=FALSE]    # 12Apr02
  cluster <- cluster[!j]
  n <- length(cluster)
		}

j <- order(cluster)
X <- X[j,,drop=FALSE]
clus.size <- table(cluster)
clus.start <- c(1,1+cumsum(clus.size))
nc <- length(levels(cluster))
clus.start <- clus.start[-(nc+1)]
storage.mode(clus.start) <- "integer"

#dyn.load("/suserlib/robcovf.o")
W <- matrix(if(.R.)
            .Fortran("robcovf", n, p, nc, clus.start, clus.size, X, 
                     double(p), double(p*p), w=double(p*p),
                     PACKAGE="Design")$w else
            .Fortran("robcovf", n, p, nc, clus.start, clus.size, X, 
                     double(p), double(p*p), w=double(p*p))$w,
            nrow=p)
#The following has a small bug but comes close to reproducing what robcovf does
#W <- tapply(X,list(cluster[row(X)],col(X)),sum)
#W <- t(W) %*% W
#The following logic will also do it, also at great cost in CPU time
#for(j in levels(cluster))		{
#   s <- cluster==j
#   if(sum(s)==1) sx <- X[s,,drop=F]
#   else {sx <- apply(X[s,,drop=F], 2, sum); dim(sx) <- c(1,p)}
#
#   sc <- sc + t(sx) %*% sx
#
#					}

adjvar <- var %*% W %*% var

#var.new <- diag(adjvar)
#deff <- var.new/var.orig; names(deff) <- vname
#eff.n <- n/exp(mean(log(deff)))

#if(pr)			{
#   v <- cbind(var.orig, var.new, deff)
#   dimnames(v) <- list(vname, c("Original Variance","Adjusted Variance",
#			"Design Effect"))
#   .Options$digits <- 4
#   cat("\n\nEffect of Adjustment for Cluster Sampling on Variances of Parameter #Estimates\n\n")
#   print(v)
#   cat("\nEffective sample size:",format(round(eff.n,1)),"\n\n")
#   nn <- n^2/sum(clus.size^2)
#   cat("\nN^2/[sum of Ni^2]    :",format(round(nn,1)),"\n\n")
#			}

fit$orig.var <- fit$var
fit$var <- adjvar
#fit$design.effects <- deff
#fit$effective.n <- eff.n
#oldClass(fit) <- c("robcov",oldClass(fit))  ##14Nov00

fit

								}
sensuc <- function(fit,
				   or.xu=seq(1,6,by=.5), or.u=or.xu, prev.u=.5,
				   constrain.binary.sample=TRUE,
				   or.method=c('x:u y:u','u|x,y'),
				   event=function(y) if(is.matrix(y))y[,ncol(y)] else 1*y) {

  type <- oldClass(fit)[1]
  if(type %nin% c('lrm','cph')) stop('fit must be from lrm or cph')

  or.method <- match.arg(or.method)

  X <- fit$x
  Y <- fit$y
  if(length(X)==0 || length(Y)==0) stop('did not specify x=T, y=T to fit')
  x <- X[,1]
  unq <- sort(unique(x))
  if(length(unq) != 2 || unq[1] != 0 || unq[2] != 1)
	stop('x is not binary')


  event <- event(Y)
  unq <- sort(unique(event))
  if(length(unq) != 2 || unq[1] != 0 || unq[2] != 1)
	stop('Y or event is not binary')

  ##Function to generate Bernoullis with exact proportion p except for roundoff
  bern <- function(n,p,constrain) {
	if(constrain) {
	  sort.random <- function(x) {
		un <- runif(length(x))
		x[order(un)]
	  }
	  ones <- round(n*p)
	  zeros <- n - ones
	  sort.random(c(rep(0,zeros),rep(1,ones)))
	} else sample(0:1, n, replace=TRUE, c(1-p,p))
  }

  a00 <- mean(!event & !x)
  a10 <- mean(event & !x)
  a01 <- mean(!event & x)
  a11 <- mean(event & x)
  p.event <- mean(event)

  b1  <- p.event
  b0  <- 1 - b1
  c1  <- mean(x)
  c0  <- 1 - c1

  n <- length(event)

  n00 <- sum(!event & !x)
  n10 <- sum(event & !x)
  n01 <- sum(!event & x)
  n11 <- sum(event & x)
  m1  <- prev.u * n
  m0  <- n - m1



  m <- length(or.xu) * length(or.u)

  OR.xu <- OR.u <- effect.x <- OOR.xu <- effect.u <- effect.u.adj <- Z <- 
	if(.R.)double(m) else single(m)

  Prev.u <- matrix(NA,nrow=m,ncol=4,
				   dimnames=list(NULL,c('event=0 x=0','event=1 x=0',
					 'event=0 x=1','event=1 x=1')))

  odds <- function(x) {
	p <- mean(x)
	p/(1-p)
  }

  j <- 0
  cat('Current odds ratio for x:u=')
  for(c.or.xu in or.xu) {
	cat(c.or.xu,'')
	for(c.or.u in or.u) {
	  j <- j + 1
	  OR.xu[j] <- c.or.xu
	  OR.u[j]  <- c.or.u

	  if(or.method=='u|x,y') {
		beta  <- logb(c.or.u)
		gamma <- logb(c.or.xu)
		f <- function(alpha,beta,gamma,a00,a10,a01,a11,prev.u) 
		  a00*plogis(alpha)+
			a10*plogis(alpha+beta)+
			  a01*plogis(alpha+gamma)+
				a11*plogis(alpha+beta+gamma) - prev.u

		alpha <- uniroot(f, lower=-10, upper=10,
						 beta=beta, gamma=gamma, 
						 a00=a00, a10=a10, a01=a01, a11=a11, 
						 prev.u=prev.u)$root
		p00 <- plogis(alpha)
		p10 <- plogis(alpha+beta)
		p01 <- plogis(alpha+gamma)
		p11 <- plogis(alpha+beta+gamma)
	  } else {
		## Raking method, thanks to M Conaway
		rake2x2 <- function(prow,pcol,odds) {
		  pstart <- matrix(1, nrow=2, ncol=2)
		  pstart[1,1] <- odds
		  pstart <- pstart/sum(pstart)
		  oldp <- pstart
		  maxdif <- 1
		  while(maxdif > .0001) {
			## Adjust row totals
			obsrow <- oldp[,1]+oldp[,2]
			adjrow <- prow / obsrow
			newp <- oldp * cbind(adjrow,adjrow)
			## Adjust col totals
			obscol <- newp[1,]+newp[2,]
			adjcol <- pcol / obscol
			newp <- newp * rbind(adjcol,adjcol)
			maxdif <- max(abs(newp - oldp))
			oldp <- newp
		  }
		  c(newp[1,],newp[2,])
		}
		
		lambda <- c.or.xu
		theta  <- c.or.u
		prow <- c(1-prev.u, prev.u)
		pcol <- c(n00,n01,n10,n11)/n
		a <- matrix(c(
					  1,0,1,0,0,0,0,0,
					  0,1,0,1,0,0,0,0,
					  0,0,0,0,1,0,1,0,
					  0,0,0,0,0,1,0,1,
					  1,1,0,0,0,0,0,0,
					  0,0,1,1,0,0,0,0,
					  0,0,0,0,1,1,0,0,
					  0,0,0,0,0,0,1,1,
					  1,0,0,0,1,0,0,0,
					  0,1,0,0,0,1,0,0,
					  0,0,1,0,0,0,1,0,
					  0,0,0,1,0,0,0,1),
					nrow=12,byrow=TRUE)
		aindx <- matrix(c(
						  1,3,
						  2,4,
						  5,7,
						  6,8,
						  1,2,
						  3,4,
						  5,6,
						  7,8,
						  1,5,
						  2,6,
						  3,7,
						  4,8),
						ncol=2, byrow=TRUE)
		pcol1 <- c(pcol[1]+pcol[3], pcol[2]+pcol[4])
		u <- rake2x2(prow, pcol1, lambda)

		pcol2 <- c(pcol[1]+pcol[2],pcol[3]+pcol[4])
		w <- rake2x2(prow, pcol2, theta)

		newp8 <- p8start <- rep(1/8, 8)
		targvec <- c(u, w, pcol)
		d <- 1
		while(d > .0001) {
		  for(i in 1:12) {
			adjust <- targvec[i] / sum(a[i,] * newp8)
			newp8[aindx[i,]] <- adjust * newp8[aindx[i,]]
		  }
		  chktarg <- a %*% as.matrix(newp8)
		  d <- max(abs(chktarg - targvec))
		}
		p00 <- newp8[5]/a00
		p01 <- newp8[6]/a01
		p10 <- newp8[7]/a10
		p11 <- newp8[8]/a11

##		  prn(c(newp8[5],newp8[5]*n,newp8[5]/(newp8[1]+newp8[5]),
##				newp8[5]*n/n00,newp8[5]/a00))
##		  w_newp8
##		  A_w[1];B_w[2];C_w[3];D_w[4];E_w[5];FF_w[6];G_w[7];H_w[8]
##		  prn((FF+H)*(A+C)/(B+D)/(E+G))
##		  prn((G+H)*(A+B)/(E+FF)/(C+D))
##		  w1_p01*b0+p11*b1
##		  w2_p00*b0+p10*b1
##		  prn((w1/(1-w1))/(w2/(1-w2)))
##		  z1_p10*c0+p11*c1
##		  z2_p00*c0+p10*c1
##		  prn((z1/(1-z1))/(z2/(1-z2)))
	  }
	  
	  Prev.u[j,] <- c(p00,p10,p01,p11)

	  u <- rep(0, n)
	  i <- !event & !x
	  u[i] <- bern(sum(i), p00, constrain.binary.sample)
	  i <- event & !x
	  u[i] <- bern(sum(i), p10, constrain.binary.sample)
	  i <- !event & x
	  u[i] <- bern(sum(i), p01, constrain.binary.sample)
	  i <- event & x
	  u[i] <- bern(sum(i), p11, constrain.binary.sample)
	  
	  OOR.xu[j] <-  odds(u[x==1])/odds(u[x==0])

	  if(type=='cph') {
		g <- coxphFit(as.matrix(u),Y,rep(1,n),toler.chol=1e-11,
					   iter.max=15,eps=.0001,method='efron')
		effect.u[j] <- exp(g$coefficients)

		g <- coxphFit(cbind(u,X),Y,rep(1,n),toler.chol=1e-11,
					   iter.max=15,eps=.0001,method='efron')
		cof <- g$coefficients
		vr  <- g$var
	  } else {
		effect.u[j] <- odds(event[u==1])/odds(event[u==0])
		g <- lrm.fit(cbind(u,X),event,maxit=20,tol=1E-11)
		ns <- g$non.slopes
		cof <- g$coefficients[-(1:ns)]
		vr  <- g$var[-(1:ns),-(1:ns)]
	  }
	  z   <- cof/sqrt(diag(vr))
	  
	  effect.u.adj[j] <- exp(cof[1])
	  effect.x[j] <- exp(cof[2])
	  Z[j]    <- z[2]
	}
  }
  cat('\n\n')

  structure(list(OR.xu=OR.xu,OOR.xu=OOR.xu,OR.u=OR.u,
				 effect.x=effect.x,effect.u=effect.u,effect.u.adj=effect.u.adj,
				 Z=Z,prev.u=prev.u,cond.prev.u=Prev.u,
				 type=type), class='sensuc')
}

plot.sensuc <- function(x, ylim=c((1+trunc(min(x$effect.u)-.01))/
							   ifelse(type=='numbers',2,1),
							   1+trunc(max(x$effect.u)-.01)),
						xlab='Odds Ratio for X:U',
						ylab=if(x$type=='lrm')'Odds Ratio for Y:U' else
						'Hazard Ratio for Y:U',
						digits=2, cex.effect=.75, cex.z=.6*cex.effect,
						delta=diff(par('usr')[3:4])/40, 
						type=c('symbols','numbers','colors'),
						pch=c(15,18,5,0),col=c(2,3,1,4),alpha=.05,
						impressive.effect=function(x)x > 1,...) {

  type  <- match.arg(type)
  Z     <- abs(x$Z)
  or    <- x$OOR.xu
  eu    <- x$effect.u
  ex    <- x$effect.x
  zcrit <- qnorm(1-alpha/2)

  plot(or, eu, ylim=ylim, xlab=xlab, ylab=ylab, type='n', ...)

  if(type=='numbers') {
	text(or, eu, round(ex,digits), cex=cex.effect)
	text(or, eu - delta, round(Z,2), cex=cex.z)
  } else {
	i <- impressive.effect(ex) & Z >= zcrit
	if(any(i)) if(type=='symbols') 
	  points(or[i], eu[i], pch=pch[1]) else
    text(or[i], eu[i], round(ex[i],digits), cex=cex.effect, col=col[1])
	i <- impressive.effect(ex) & Z < zcrit
	if(any(i)) if(type=='symbols') 
	  points(or[i], eu[i], pch=pch[2]) else
    text(or[i], eu[i], round(ex[i],digits), cex=cex.effect, col=col[2])
	i <- !impressive.effect(ex) & Z < zcrit
	if(any(i)) if(type=='symbols')
	  points(or[i], eu[i], pch=pch[3]) else
    text(or[i], eu[i], round(ex[i],digits), cex=cex.effect, col=col[3])
	i <- !impressive.effect(ex) & Z >= zcrit
	if(any(i)) if(type=='symbols')
	  points(or[i], eu[i], pch=pch[4]) else
	text(or[i], eu[i], round(ex[i],digits), cex=cex.effect, col=col[4])
	
  }
  title(sub=paste('Prevalence of U:',format(x$prev.u)),adj=0)
  invisible()
}
#Print description of specifications.  Can come from individual variable
#created by dx, complete design created by design(), or complete design
#carried forward in fit
#Mod 10Jul91 - print freq table of strata factors
#Mod 28Aug91 - added option long=F to suppress printing limits
#Mod 30Oct91 - changed to specs.Design called from generic specs
#Mod 25Sep92 - print transformations, clean up code
#Mod 22Sep93 - change to new storage format for attributes

if(.R.) specs <- function(fit, ...) UseMethod('specs')

specs.Design<-function(fit, long=FALSE, ...){

Call <- if(length(fit$call))fit$call else
  if(length(attr(fit,'call')))attr(fit,'call') else attr(fit,'formula')

tl <- attr(fit$terms, "term.labels")
ass <- fit$assign
strata <- fit$strata
if(is.null(fit$assume)) {
  #x <- attr(x$terms, "Design")   17Apr01
  d <- fit$Design    ## 30may02
  if(!length(d)) d <- getOldDesign(fit)  ## 30may02
  fit <- d
}
assume <- fit$assume
if(is.null(assume)) stop("fit does not have design information")
parms <- fit$parms
name  <- fit$name
lim   <- fit$limits
ia.order <- fit$ia.order
label <- fit$label
units <- fit$units

if(length(ass)) {
 if(names(ass)[1]=="(Intercept)" | names(ass)[1]=="Intercept") ass[[1]] <- NULL
 names(ass) <- name[assume!="strata"]
}
f <- length(assume)
d<-matrix("",nrow=f,ncol=3)
d[,1]<-assume
iint <- 0
jfact <- 0
trans <- rep("",f)
#Pick off inner transformation of variable. To complete, need to
#evaluate h function
#from <- c("asis","pol","lsp","rcs","catg","scored","strat","matrx","I")
#from <- paste(from,"(\\(.*\\))",sep="")
#tl <- translate(tl, from, "\\1")
#tl <- paste("h(",tl,")",sep="")

from <- c('asis(*)','pol(*)','lsp(*)','rcs(*)','catg(*)','scored(*)',
  'strat(*)','matrx(*)','I(*)')
to   <- rep('*',9)

tl <- paste("h(",sedit(tl, from, to),")",sep="")
#change wrapping function to h()

h <- function(x,...)deparse(substitute(x))
for(i in 1:f)	{
	if(assume[i]=="interaction") iint <- iint+1
	else		{
	  tr <- eval(parse(text=tl[i]))
	  if(tr!=name[i]) trans[i] <- tr
			}
	len <- if(assume[i]=="strata") 0 else length(ass[[name[i]]])
	d[i,3] <- as.character(len)
	parmi <- parms[[name[i]]]
	if(d[i,1]=="transform")d[i,2]<-"function"	else	{
	if(length(parmi))		{
	if(d[i,1]=="interaction")		{
	   i1 <- parmi[1,-1]
	   i2 <- parmi[2,-1]
	   i3 <- parmi[3,-1]
	   if(parmi[3,1]==0)	{   #2nd order interaction
	      iao <- 1*(any(i1) & !any(i2))+
	             2*(!any(i1) & any(i2))+
		     3*(any(i1) & any(i2) & !any(i1&i2))+
		     4*any(i1 & i2)
		d[i,2]<-c("linear x linear - AB",
			"nonlinear x linear - f(A)B",
			"linear x nonlinear - Ag(B)",
			"Af(B) + Bg(A)",
			"f(A,B) - all cross-products")[iao+1]
				}
	   else			#3rd order
	      d[i,2] <- paste(if(any(i1))"nonlinear" else "linear","x",
			if(any(i2))"nonlinear" else "linear","x",
			if(any(i3))"nonlinear" else "linear")
	   if(ncol(parmi)==1)  d[i,2] <- " "
						}
	
	else	{
		lab<-""
		for(z in parmi)
			if(is.character(z))lab<-paste(lab,z) else
				lab<-paste(lab,
					signif(as.single(z),5))
		d[i,2]<-lab	}	}}
									}
collab <- c("Assumption","Parameters","d.f.")
if(any(trans!=""))	{
  collab <- c("Transformation",collab)
  d <- cbind(trans,d)	}

if(any(name!=label))			{
  collab <- c("Label",collab)
  d <- cbind(label,d)		}
if(length(units) && any(units != '')) {  #9Jun99
  collab <- c('Units',collab)
  unitsb <- rep('',length(assume))
  unitsb[assume!='interaction'] <- units
  d <- cbind(unitsb,d)
}
dimnames(d) <- list(name, collab)

structure(list(call=Call,how.modeled=d,limits=if(long)lim,strata=strata),
          class='specs.Design')
}

print.specs.Design <- function(x, ...) {
dput(x$call)
cat('\n')
print(x$how.modeled, quote=FALSE)
if(length(x$limits)) {cat('\n'); print(x$limits)}
if(length(x$strata)) {cat("\n        Strata\n\n");print(strata,quote=FALSE)}
invisible()
}


# Value adjusted to is irrelevant when the factor does not interact with
# other factors.  Form of factors is as follows: factor1=value1,factor2=val2:
# Values:
#	NA	: test factor, use all default settings
#	w	: adjust this factor to w when estimating effects of others
#	c(lo,hi): use range for effect (lo,hi), adjust to default value
#	c(lo,w,hi): use range (lo,hi), adjust to w.  Any of 3 can be NA.
# For categories and strata values can be character
# values that are original values before translating to is.category -
# only enough letters are needed to uniquely identify the category 
# This applies to category and strata vars.  Default adjusted to is
# from second element of limits vector.
# For category factors, all comparisons to reference category are made.
# Reference category is assumed to be adjusted to value.
# est.all is T to estimate effects for all factors, not just those listed
# in ...

summary.Design <- function(object, ..., est.all=TRUE, antilog, conf.int=.95,
			abbrev=FALSE) {	

obj.name <- as.character(sys.call())[2]
#at <- attr(object$terms, "Design")    17Apr01
at <- object$Design
if(!length(at)) at <- getOldDesign(object)

assume <- at$assume.code
if(is.null(assume))stop("fit does not have design information")
if(any(assume==10))
 warning("summary.Design does not currently work with matrix factors in model")
name <- at$name
parms <- at$parms

scale <- object$scale.pred
if(missing(antilog)) antilog <- length(scale)==2
if(antilog & length(scale)<2) scale <- c("","Antilog")

factors <- list(...)
nf <- length(factors)

if(est.all) which <- (1:length(assume))[assume!=9]
if(nf>0)					{
	jw <- charmatch(names(factors),name,0)
	if(any(jw==0))stop(paste("factor name(s) not in the design:",
		paste(names(factors)[jw==0],collapse=" ")))
	if(!est.all) which <- jw
	if(any(assume[which]==9))
	   stop("cannot estimate effects for interaction terms alone")
						}

Limval <- Getlim(at, allow.null=TRUE, need.all=FALSE)
values <- Limval$values
## The next statement (9Jun98) makes limits[1:3,] keep all levels of
## factors.  Problem is that [.data.frame does not pass drop to []
## when first subscripts are specified
oldopt <- options(drop.factor.levels=FALSE)
on.exit(options(oldopt))

lims <- Limval$limits[1:3,,drop=FALSE]

#The following won't work with new data.frame functions - still keeps it
#a factor() object.   6Feb95
#for(i in 1:length(lims))  #the following still preserves class data.frame
#   if(is.factor(lims[[i]]))lims[[i]] <- as.character(lims[[i]])

#li <- vector("list", length(lims))
#for(i in 1:length(lims)) li[[i]] <- if(is.factor(lims[[i]])) 
#	as.character(lims[[i]]) else lims[[i]]
#lims <- structure(li, class="data.frame", row.names=row.names(lims),
#                  names=names(lims))

#Find underlying categorical variables
ucat <- rep(FALSE, length(assume))
for(i in (1:length(assume))[assume!=5 & assume<8])
   ucat[i] <- !is.null(V <- values[[name[i]]]) && is.character(V)

stats <- NULL
lab <- NULL
lc <- length(object$coef)
#Number of non-slopes:
nrp <- num.intercepts(object)
nrp1 <- nrp+1
# Exclude non slopes
beta <- object$coef[nrp1:lc]
var <- Varcov(object, regcoef.only=TRUE)[nrp1:lc,nrp1:lc]

zcrit <- qnorm((1+conf.int)/2)
cll <- paste(signif(as.single(conf.int),3))

jf <- 0
if(nf>0) for(i in jw)						{
	jf <- jf+1
	z <- value.chk(at, i, factors[[jf]], 0, Limval)
	lz <- length(z)
	if(lz==1 && !is.na(z)) lims[2,i] <-  z
	if(lz==2) 				{
		if(!is.na(z[1])) lims[1,i] <- z[1]
		if(!is.na(z[2])) lims[3,i] <- z[2]	}	else
	if(lz==3) lims[!is.na(z),i] <- z[!is.na(z)]
	if(lz<1 | lz>3) stop("must specify 1,2, or 3 values for a factor")
								}
adj <- lims[2,,drop=FALSE]
isna <- sapply(adj, is.na)

if(any(isna)) stop(
   paste("adjustment values not defined here or with datadist for",
	paste(name[assume!=9][isna],collapse=" ")))
k <- which[assume[which]!=8 & assume[which]!=5 & assume[which]!=10 & 
	!ucat[which]]
m <- length(k)
if(m)
{
   isna <- is.na(lims[1,name[k],drop=FALSE]+lims[3,name[k],drop=FALSE]) #added ,drop 7Feb94
   # k was which[k]   4May94 (also 2 lines down)
   # name[k] was k 4Dec00 - don't know why it ever worked
  #note char. excluded from k
   if(any(isna)) stop(paste("ranges not defined here or with datadist for",
	   paste(name[k[isna]], collapse=" ")))
}

xadj <- oldUnclass(Design.levels(adj, at))	# unclass 28jul03

m <- length(k)
if(m)							{
   adj <- xadj
   M <- 2*m
   odd <- seq(1,M,by=2)
   even<- seq(2,M,by=2)
  #Extend data frame
#   adj[2:M,] <- adj  28jul03
   for(i in 1:length(adj)) adj[[i]] <- rep(adj[[i]], M)
   
   i <- 0
   for(l in k) 	{
     i <- i+1
##     adj[[l]][(2*i-1):(2*i)] <-  lims[c(1,3),l]  4Dec00
     adj[[name[l]]][(2*i-1):(2*i)] <- lims[c(1,3),name[l]]
   }
#  adj <- data.frame(adj)
  xx <- predictDesign(object, newdata=adj, type="x", incl.non.slopes=FALSE)
  xd <- matrix(xx[even,]-xx[odd,],nrow=m)
  xb <- xd %*% beta
  se <- drop((((xd %*% var) * xd) %*% rep(1,ncol(xd)))^.5)
  low <- xb - zcrit*se
  up <- xb + zcrit*se
##   lm <- as.matrix(lims[,k,drop=FALSE])  # 4Dec00
   lm <- as.matrix(lims[,name[k],drop=FALSE])
  stats <- cbind(lm[1,],lm[3,],lm[3,]-lm[1,],xb,se,low,up,1)
  lab <- name[k]
  if(antilog)	{
    stats <- rbind(stats,cbind(stats[,1:3,drop=FALSE],exp(xb),NA,exp(low),exp(up),
		2))
    lab <- c(lab,rep(paste("",scale[2]),m))
    w <- integer(M)
    w[odd] <- 1:m
    w[even]<- m+(1:m)
    stats <- stats[w,]
    lab <- lab[w]
		}
							}

for(i in which[assume[which]==5 | ucat[which]])		{
		#All comparisons with reference category
#		xadj[2,] <- xadj[1,]  		#duplicate row 28jul03
  for(j in 1:length(xadj)) xadj[[j]] <- rep(xadj[[j]], 2)
  
		parmi <- if(ucat[i]) values[[name[i]]] else parms[[name[i]]]
		parmi.a <- if(abbrev) abbreviate(parmi) else parmi
		## iref <- as.character(xadj[1,name[i]])   # as.char 2Dec94
  iref <- as.character(xadj[[name[i]]][1])   # 28jul03
		ki <- match(iref, parmi)
		for(j in parmi)					{
			if(j!=iref)			{
				kj <- match(j, parmi)
				adj <- xadj
#				adj[,name[i]] <- c(iref,j)  28jul03
                adj[[name[i]]] <- c(iref,j)
				adj <- as.data.frame(adj)  # 28jul03
				xx <- predictDesign(object,newdata=adj,
				   type="x",incl.non.slopes=FALSE)
				xd <- matrix(xx[2,]-xx[1,],nrow=1)
				xb <- (xd %*% beta)
				se <- sqrt((xd %*% var) %*% t(xd))
				low <- xb - zcrit*se
				up <- xb + zcrit*se
				stats <- rbind(stats,cbind(ki,kj,NA,
					xb,se,low,up,1))
				lab <-c(lab,paste(name[i]," - ",parmi.a[kj],":",
					parmi.a[ki],sep=""))
				if(antilog)				{
					stats <- rbind(stats,cbind(ki,kj,NA,
						exp(xb),NA,exp(low),exp(up),2))
					lab <- c(lab, paste("",scale[2]))}
								}	}
									     }

dimnames(stats) <- list(lab, c("Low","High",
	"Diff.","Effect","S.E.",paste("Lower",cll),paste("Upper",cll),"Type"))

attr(stats,"heading") <- paste("             Effects              Response : ",
 as.character(formula(object))[2], sep='')
## was as.character(fit$formula)[2], sep='') 30may02
## was as.character(attr(fit$terms,"formula")[2]),sep="")  22may02
attr(stats,"class") <- if(.SV4.)'summary.Design' else
  c("summary.Design","matrix")   ##13Nov00
attr(stats,"scale") <- scale
attr(stats,"obj.name") <- obj.name
interact <- at$interactions
adjust <- ""
if(length(interact))
{
   interact <- sort(unique(interact[interact>0]))
   nam <- name[which[match(which, interact, 0)>0]]
   if(length(nam)) for(nm in nam) 
	adjust <- paste(adjust, nm,"=",
#	   if(is.factor(xadj[1,nm])) as.character(xadj[1,nm])  28jul03
#	   else format(xadj[1,nm])," ",sep="")
                    if(is.factor(xadj[[nm]]))
                    as.character(xadj[[nm]])[1] else
                    format(xadj[[nm]][1])," ",sep="")
}
attr(stats,"adjust") <- adjust

stats									}



print.summary.Design <- function(x, ...) {

cstats <- dimnames(x)[[1]]
for(i in 1:3) cstats <- cbind(cstats, format(signif(as.single(x[,i]),5)))
for(i in 4:7) cstats <- cbind(cstats, format(round(x[,i],2)))
dimnames(cstats) <- list(rep("",nrow(cstats)), 
                c("Factor", dimnames(x)[[2]][1:7]))
cat(attr(x,"heading"),"\n\n")
print(cstats,quote=FALSE)
if((A <- attr(x,"adjust"))!="") cat("\nAdjusted to:", A,"\n\n")

invisible()
}


latex.summary.Design <- function(object, 
  title=if(under.unix) paste('summary',attr(object,'obj.name'),sep='.') else
          paste("sum",substring(first.word(attr(object,"obj.name")),
                                1,5),sep=""),
  ...) { 

##expr= in first.word 18Nov00 removed 25May01
#cstats <- dimnames(stats)[[1]]
title <- title   # because of lazy evaluation
caption <- attr(object, "heading")
scale <- attr(object,"scale")
if(.SV4.)object <- matrix(oldUnclass(object), nrow=nrow(object),
                         dimnames=dimnames(object)) ## 14Nov00
object <- object[,-8,drop=FALSE]
rowl <- dimnames(object)[[1]]
rowl <- ifelse(substring(rowl,1,1)==" ",
	paste("~~{\\it ",substring(rowl,2),"}",sep=""), rowl) # preserve leading blank
rowl <- sedit(rowl, "-", "---")   # was translate
cstats <- matrix("", nrow=nrow(object), ncol=ncol(object), 
  dimnames=dimnames(object))
for(i in 1:3) cstats[,i] <- format(signif(as.single(object[,i]),5))
for(i in 4:7) cstats[,i] <- format(round(object[,i],2))
cstats[is.na(object)] <- ""
caption <- sedit(caption, "    Response","~~~~~~Response") #,multichar=TRUE)
cstats <- as.data.frame(cstats)
attr(cstats,"row.names") <- rowl
names(cstats)[3] <- "$\\Delta$"
latex(cstats, caption=caption, title=title, rowlabel="",
      col.just=rep("r",7), ...)
#      n.rgroup=rep(length(scale),nrow(object)/length(scale)),...)
}


plot.summary.Design <- function(x, at, log=FALSE, 
	q=c(0.7, 0.8, 0.9, 0.95, 0.99), xlim, nbar, cex=1, nint=10, cex.c=.5,
	cex.t=1, clip=c(-1e30,1e30), main, ...)
{


scale  <- attr(x, "scale")
adjust <- attr(x, "adjust")
if(.SV4.) x <- matrix(oldUnclass(x), nrow=nrow(x),
                          dimnames=dimnames(x))  ##14Nov00
## so subscripting works

Type   <- x[,"Type"]
x  <- x[Type==1,,drop=FALSE]
lab    <- dimnames(x)[[1]]
effect <- x[,"Effect"]
se     <- x[,"S.E."]
if(!log && any(Type==2))
{
   fun <- exp
   tlab <- scale[2]
}
else
{
   fun <- function(x) x
   if(log)
   {
      if(length(scale)==2) tlab <- scale[2]
      else tlab <- paste("exp(",scale[1],")",sep="")
   }
   else tlab <- scale[1]
}
if(!length(scale)) tlab <- ''  ## 2dec02; mainly for glmD fits
if(!missing(main)) tlab <- main
augment <- if(log | any(Type==2)) c(.1, .5, .75, 1) else 0
n     <- length(effect)
out   <- qnorm((max(q)+1)/2)
if(missing(xlim) && !missing(at)) xlim <- range(if(log)logb(at) else at) else
if(missing(xlim))
{
   xlim <- fun(range(c(effect-out*se,effect+out*se)))
   xlim[1] <- max(xlim[1],clip[1])
   xlim[2] <- min(xlim[2],clip[2])
}
else augment <- c(augment, if(log)exp(xlim) else xlim)   #added 24oct94

fmt <- function(k) {
  m <- length(k)
  f <- character(m)
  for(i in 1:m) f[i] <- format(k[i])
  f
}
lb <- ifelse(is.na(x[,'Diff.']), lab, paste(lab,' - ',
             fmt(x[,'High']),':',fmt(x[,'Low']),sep=''))
## mxlb <-(1+max(nchar(lb)))*cex*par('cin')[1] 30jul02
if(.R.) { plot.new(); par(new=TRUE) }  ## 9apr03
mxlb <- .1+max(strwidth(lb,units='inches',cex=cex))
tmai <- par('mai')
on.exit(par(mai=tmai))
if(.R.) par(mai=c(tmai[1],mxlb,1.5*tmai[3],tmai[4])) else
        par(mai=c(tmai[1],mxlb,tmai[3:4]))

outer.widths <- fun(effect+out*se)-fun(effect-out*se)
if(missing(nbar)) nbar <- n
npage <- ceiling(n/nbar)
is <- 1
for(p in 1:npage) {
  ie <- min(is+nbar-1, n)
  plot(1:nbar, rep(0,nbar), xlim=xlim, ylim=c(1,nbar), type="n", axes=FALSE, 
	xlab="", ylab="")
  if(cex.t>0) title(tlab, cex=cex.t)
  lines(fun(c(0,0)),c(nbar-(ie-is), nbar),lty=2)
  if(log)
  {
     pxlim <- pretty(exp(xlim), n=nint)
     pxlim <- sort(unique(c(pxlim, augment)))
    # For wome weird reason, sometimes duplicates (at xlim[2]) still remain
     pxlim <- pxlim[pxlim>=exp(xlim[1])]    # was > 24oct94
     if(!missing(at)) pxlim <- at
     axis(3, logb(pxlim), lab=format(pxlim))
  }
  else 
  {
     pxlim <- pretty(xlim, n=nint)
     pxlim <- sort(unique(c(pxlim, augment)))
     pxlim <- pxlim[pxlim>=xlim[1]]         # was > 24oct94
     if(!missing(at)) pxlim <- at
     axis(3, pxlim)
  }
  imax <- (is:ie)[outer.widths[is:ie]==max(outer.widths[is:ie])][1]
  for(i in is:ie)
  {
     confbar(nbar-(i-is+1)+1, effect[i], se[i], q=q, type="h", 
             fun=fun, cex=cex.c, labels=i==imax, clip=clip, ...)
     if(.R.) mtext(lb[i], 2, 0, at=nbar-(i-is+1)+1, cex=cex,
                   adj=1, las=1) else
     mtext(lb[i], 2, 0, at=nbar-(i-is+1)+1, srt=0, cex=cex, adj=1)
  }
  if(adjust!="") 
  {
     adjto <- paste("Adjusted to:",adjust,sep="")
     xx <- par('usr')[2]
     if(nbar>ie) text(xx, nbar-(ie-is+1), adjto, adj=1, cex=cex)
     else title(sub=adjto, adj=1, cex=cex)
  }
  is <- ie+1
}
invisible()
}
#SCCS 2/28/95 @(#)summary.survfit.s	1.7
summary.survfit <- function(object, times, censored=FALSE, scale=1, ...) {
  fit <- object  # FEH
  
  if (!inherits(fit, 'survfit')) stop("Invalid data")

    if(.R.) require('survival')
    n <- length(fit$time)
    stime <- fit$time/scale
    if (!length(fit$strata)) {
	stemp <- rep(1,n)
	nstrat <- 1
	}
    else {
	nstrat <- length(fit$strata)
	stemp <- rep(1:nstrat,fit$strata)
	}

    surv <- as.matrix(fit$surv)
    if (!length(fit$std.err)) std.err <- NULL
    else                      std.err <- fit$std.err * surv

    if (length(fit$lower)) {
	lower <- as.matrix(fit$lower)
	upper <- as.matrix(fit$upper)
	}
    if (missing(times)) {
	if (censored) {
	    times <- stime
	    n.risk<- fit$n.risk
	    n.event <- fit$n.event
	    }
	else {
	    who    <- (fit$n.event >0)
	    times  <-  stime[who]
	    n.risk <-  fit$n.risk[who]
	    n.event <- fit$n.event[who]
	    stemp <- stemp[who]
	    surv <- surv[who,,drop=FALSE]
	    if (length(std.err)) std.err <- std.err[who,,drop=FALSE]
	    if (length(fit$lower)) {
		lower <- lower[who,,drop=FALSE]
		upper <- upper[who,,drop=FALSE]
		}
	    }
	}

    else {  #this case is much harder
	if (any(fit$time<0)) stop("Negative times present.\nIf using survplot, don't ask for confidence bars")
        if(max(fit$time) < min(times))
            stop("Requested times are all beyond the end of the survival curve")
	if (length(times) >1 )
	    if (any(diff(times)<0)) stop("Times must be in increasing order")

	temp <- if(.R.)
      .C("survindex2",  ## NEEDS FIXING!!
         as.integer(n),
         as.double(stime),
         as.integer(stemp),
         as.integer(length(times)),
         as.double(times),
         as.integer(nstrat),
         indx = integer(nstrat * length(times)),
         indx2 = integer(nstrat * length(times)),
         PACKAGE="survival") else
      .C(if(.SV4.)'S_survindex2' else "survindex2",  ##14Nov00
               as.integer(n),
				  as.double(stime),
				  as.integer(stemp),
				  as.integer(length(times)),
				  as.double(times),
				  as.integer(nstrat),
				  indx = integer(nstrat*length(times)),
				  indx2= integer(nstrat*length(times)) )
	keep <- temp$indx >=0
	indx <- temp$indx[keep]
	ones <- (temp$indx2==1)[keep]
	ties <- (temp$indx2==2)[keep]  #data set time === requested time
	times <- rep(times, nstrat)[keep]
	n.risk <- fit$n.risk[indx+1 - (ties+ones)]
	surv   <- surv[indx,,drop=FALSE];   surv[ones,] <- 1
	if (length(std.err)) {
	    std.err<- std.err[indx,,drop=FALSE]
	    std.err[ones,] <-0
	    }
	fit$n.event[stime>max(times)] <- 0
	n.event <- (cumsum(c(0,fit$n.event)))[ifelse(ones, indx, indx+1)]
	n.event<-  diff(c(0, n.event))

	if (length(fit$lower)) {
	    lower <- lower[indx,,drop=FALSE];  lower[ones,] <- 1;
	    upper <- upper[indx,,drop=FALSE];  upper[ones,] <- 1;
	    }

	stemp <- stemp[indx]
	}

    ncurve <- ncol(surv)
    temp <- list(surv=surv, time=times, n.risk=n.risk, n.event=n.event,
			conf.int=fit$conf.int)
    if (ncurve==1) {
	temp$surv <- drop(temp$surv)
	if (length(std.err)) temp$std.err <- drop(std.err)
	if (length(fit$lower)) {
	    temp$lower <- drop(lower)
	    temp$upper <- drop(upper)
	    }
	}
    else {
	if (length(std.err)) temp$std.err <- std.err
	if (length(fit$lower)) {
	    temp$lower <- lower
	    temp$upper <- upper
	    }
	}
    if(existsFunction('print.survfit.computations'))
      temp$table <- print.survfit.computations(fit, scale)
    ## Kattan added 8/15/2001
    if (length(fit$strata))
	temp$strata <- factor(stemp,
	    labels = names(fit$strata)[sort(unique(stemp))])
    temp$call <- fit$call
    if (length(fit$na.action)) temp$na.action <- fit$na.action
    oldClass(temp) <- 'summary.survfit'
    temp
    }

##Use x= if input is a design matrix, newdata= if a data frame or data matrix
##or vector.  Can specify (centered) linear predictor values instead (linear.predictors).
##Strata is attached to linear.predictors or x as "strata" attribute.
##data matrix assumes that categorical variables are coded with integer codes
##7Jun99: took away checks on consistency between type and type used
##in fit
##18Sep99: added what='parallel' for val.surv


survest.cph <- function(fit, newdata, linear.predictors, x, times, fun,
                        loglog=FALSE, conf.int=.95,
                        type=NULL, vartype=NULL,
                        conf.type=c("log-log","log","plain","none"),
                        se.fit=TRUE, what=c("survival","parallel"),
                        individual=FALSE, ...) {

  at <- fit$Design
  if(!length(at)) at <- getOldDesign(fit)

  f <- sum(at$assume.code!=8)		#non-strata factors
  nf <- length(at$name)-f
  num.strata <- if(nf==0)1 else length(fit$strata)
  strata.levels <- fit$strata  ## 7may02

#  N <- if(!missing(newdata)) nrow(newdata) else
#  if(!missing(linear.predictors))length(linear.predictors) else
#  if(!missing(x)) nrow(x) else
#  if(ll <- length(fit$linear.predictors)) ll else
#  if(length(fit$x)) nrow(fit$x)

  conf.type <- match.arg(conf.type)
  what <- match.arg(what)
  if(what=='parallel') {conf.int <- FALSE; conf.type <- 'none'}

  if(!se.fit) conf.int <- 0
  ##stype <- attr(fit$surv,"type")
  ##if(length(stype)==0)stype <- "tsiatis"

  if(individual && (length(fit$x)==0 || length(fit$y)==0 || 
                    attr(fit$y,'type')!='counting')) 
    stop('must specify x=T, y=T, and start and stop time to cph when individual=T')

  if(missing(fun)) fun <-
    if(loglog) function(x) logb(-logb(ifelse(x==0|x==1,NA,x)))
    else function(x) x

  naa <- fit$na.action

  ##First see if use method that depends on x and y being stored in fit

  if(!missing(linear.predictors) && length(fit$surv)==0)
    stop('when using linear.predictors= must have specified surv=T to cph')

  if(length(fit$y) && (f==0 || length(fit$x)) &&
     ((conf.int>0 && f>0) | length(fit$surv)==0) &
     (!missing(newdata) | (missing(linear.predictors) && missing(x))))  {

    if(conf.int>0 && conf.type!="log-log" && length(fit$surv))
      warning(paste("Using conf.type=",conf.type,
                    ".\n Survival estimates stored with fit used conf.type=log-log",sep=""))
    if(!missing(linear.predictors) | !missing(x))
      stop(paste("may not specify linear.predictors or x when survival estimation",
                 "is not using underlying survival estimates stored with surv=T"))

    sf <- function(..., type=NULL, vartype=NULL, cphnull=FALSE) {
      g <- list(...)
      if(length(type))    g$type <- type
      if(length(vartype)) g$vartype <- vartype
      do.call(if(cphnull)'survfit.cph.null' else 'survfit.cph', g)
    }
    
    if(f==0) {
      g <- sf(fit, se.fit=se.fit, conf.int=conf.int, conf.type=conf.type,
              type=type, vartype=vartype, cphnull=TRUE)
      if(missing(newdata))sreq <- 1 else sreq <- 
        attr(predict(fit, newdata, type="lp", expand.na=FALSE),"strata")
    }
    else {
      if(missing(newdata)) {
        g <- sf(fit, se.fit=se.fit, conf.int=conf.int, conf.type=conf.type,
                type=type, vartype=vartype)
        sreq <- 1
      }
      else {
        if(nrow(newdata)>1 && !individual && missing(times))
          stop("must specify times= if predicting for >1 observation")
        g <- sf(fit, newdata=newdata, se.fit=se.fit, conf.int=conf.int,
                conf.type=conf.type, individual=individual,
                type=type, vartype=vartype)
        sreq <- g$requested.strata
      }
      naa <- g$na.action
    }
    sreq <- oldUnclass(sreq)

    if(missing(times)) {
      ##delete extra S(t) curves added by survfit for all strata
      ##No newdata -> requested underlying survival for all strata
      if(missing(newdata)) return(g)
      else {
        if(nf==0) j <- TRUE
        else {
          stemp <- rep(1:num.strata, g$strata)
          j <- stemp==sreq
        }

        tim <- c(0,g$time[j]); nr <- c(g$n.risk[j][1],g$n.risk[j])
        ne <- c(0,g$n.event[j]); surv <- c(1, g$surv[j])
        se <- c(NA, g$std.err[j])
        upper <- c(NA, g$upper[j]); lower <- c(NA,g$lower[j])

        yy <- fit$y; ny <- ncol(yy)
        str <- oldUnclass(attr(yy, "strata"))
        if(length(str)) yy <- yy[str==sreq,ny-1] else yy <- yy[,ny-1]
        maxt <- max(yy)
        if(maxt>tim[length(tim)]) {
          tim <- c(tim,maxt); nr <- c(nr, sum(yy>=maxt-1e-6))
          ne <- c(ne, 0); surv <- c(surv, surv[length(surv)])
          se <- c(se, NA); upper <- c(upper, NA); lower <- c(lower, NA)
        }
        se <- -se/logb(ifelse(surv==0|surv==1,NA,surv))

        surv  <- fun(surv);  surv[is.infinite(surv)] <- NA
        lower <- fun(lower); lower[is.infinite(lower)] <- NA
        upper <- fun(upper); upper[is.infinite(upper)] <- NA
        retlist <- list(time=tim,n.risk=nr,
                        n.event=ne,
                        surv=surv, 
                        std.err=se,
                        upper=upper, lower=lower, 
                        conf.type=g$conf.type,
                        conf.int=g$conf.int, call=g$call)
        if(nf>0) retlist$strata <- sreq  #was g$strata[sreq]
        return(retlist)
      }
    }
    else {
      g <- summary.survfit(g, print.it=FALSE, times=times)
      ##g$surv <- exp(-exp((log(-log(g$upper))+log(-log(g$lower)))/2))
      if(nf>0) { #delete extra cells added by survfit for strat
        if(length(g$time) != length(times)*num.strata)
          stop('summary.survfit could not compute estimates for all strata at all times requested.\nYou probably requested times where data are limited.')
          
        d <- dim(g$surv); if(length(d)==0) d <- c(length(g$surv),1)
        strata.col <- matrix(rep(sreq,d[1]),ncol=d[2],byrow=TRUE)
        ## g$strata had dropped unused strata and renumbered codes 7may02
        gs <- factor(g$strata, strata.levels)
        strata.row <- matrix(rep(oldUnclass(gs),d[2]),ncol=d[2]) # was g$strata
        m <- strata.col==strata.row
        g$surv    <- matrix(g$surv[m],   ncol=d[2])[,,drop=TRUE]
        g$lower   <- matrix(g$lower[m],  ncol=d[2])[,,drop=TRUE]
        g$upper   <- matrix(g$upper[m],  ncol=d[2])[,,drop=TRUE]
        g$std.err <- matrix(g$std.err[m],ncol=d[2])[,,drop=TRUE]
      }
    }
    tim  <- g$time
    nr   <- g$n.risk
    ne   <- g$n.event
    surv <- g$surv
    se   <- g$std.err
    low  <- g$lower
    up   <- g$upper
    tim <- unique(tim)
    if(is.matrix(surv)) {
      surv <- t(surv); se <- t(se); low <- t(low); up <- t(up)
      dn <- list(row.names(newdata),format(tim))
      dimnames(surv) <- dn; dimnames(se) <- dn; dimnames(low) <- dn;
      dimnames(up) <- dn
    }

    se   <- -se/logb(ifelse(surv==0|surv==1,NA,surv))
    if(!.R.) {
      storage.mode(surv) <- "single"; storage.mode(se) <- "single"
      storage.mode(low)  <- "single"; storage.mode(up) <- "single"
    }
    surv <- fun(surv); low <- fun(low); up <- fun(up)
    surv[is.infinite(surv)] <- NA
    low[is.infinite(low)] <- NA
    up[is.infinite(up)] <- NA
    retlist <- list(time=tim, surv=naresid(naa,surv), std.err=naresid(naa,se)
                    ,lower=naresid(naa,low), upper=naresid(naa,up))
    if(nf>0) retlist$strata <- naresid(naa,sreq)
    return(retlist)
  }

#  strata.levels <- fit$strata  7may02
  asnum.strata <- function(str, strata.levels) {
    if(!length(str)) return(NULL)  # 4Aug01; thanks Mike Kattan
    if(is.numeric(str) && any(str<1 | str>length(strata.levels)))
      stop('illegal stratum number')

    if(is.category(str) || is.numeric(str)) return(as.integer(str))

    i <- match(str, strata.levels, nomatch=0)
    if(any(i==0)) stop(paste('illegal strata:',
             paste(str[i==0],collapse=' ')))
    i
  }

  ##Instead use the baseline survival computed at fit time with cph(...,surv=T)

  nt <- if(missing(times))0 else length(times)
  if(conf.int>0 && f>0)
    warning(paste("S.E. and confidence intervals are approximate except",
                  "at predictor means.\nUse cph(...,x=T,y=T) (and don't use linear.predictors=) for better estimates."))
#  if(!missing(type) && method!=stype)
#    warning(paste("method=",stype,
#                  " used since\nusing survival estimates stored with fit",sep=""))
  if(conf.type!="log-log" && conf.type!="none")
	stop("only conf.type=log-log can be used since using survival estimates stored with fit")

  if(missing(linear.predictors)) {
    if(missing(x) && missing(newdata)) {
      linear.predictors <- fit$linear.predictors	#assume was centered
      rnam <- names(linear.predictors)
      if(length(linear.predictors)==0) {
        if(length(fit$x)==0)
          stop("newdata, x, linear.predictors not given but x nor linear.predictors stored in fit")
        linear.predictors <- matxv(fit$x, fit$coef)-fit$center
        strata <- fit$strata
        rnam <- dimnames(fit$x)[[1]]
      }
      else strata <- attr(linear.predictors,"strata")
    }
    else {
      if(missing(x)) 
        {x <- predict(fit,newdata,type="x",expand.na=FALSE)
         naa <- attr(x,"na.action")}
      strata <- attr(x,"strata")
      if(f>0) linear.predictors <- matxv(x,fit$coef) - fit$center 
      else linear.predictors <- 0
      rnam <- dimnames(x)[[1]]
    }
  }
  else {
    strata <- asnum.strata(attr(linear.predictors, "strata"),strata.levels)
    rnam <- names(linear.predictors)
  }
  if(length(strata)==0 && nf>0) 
    stop("strata not stored in x or linear.predictors")
  attr(strata, "class") <- NULL

  if(length(fit$surv)==0 && length(fit$x)==0 && length(fit$y)==0)
    stop("you did not specify surv=T or x=T, y=T in cph")

  if(conf.int>0) zcrit <- qnorm((conf.int+1)/2)
  if(length(strata)==0) {
	n <- length(linear.predictors)
	strata <- rep(1,n)
	ns <- 1
  } else {
    ns <- max(strata, na.rm=TRUE)   #na.rm added 8Jun94
    n <- length(strata)
  }

  if(what=='parallel') {
    if(length(times)>1 && length(times) != n)
      stop('length of times must equal 1 or number of subjects being predicted')
    if(!length(fit$surv))stop('must specify surv=T to cph')
    if(diff(range(strata))==0) {
      estsurv <- approx(fit$time, fit$surv, xout=times,
                        method="constant", f=0)$y
      return(estsurv ^ exp(linear.predictors))
    }
    est.surv <- if(.R.)double(n) else single(n)
    for(zs in unique(strata)) {
      this <- strata==zs
      estsurv <- approx(fit$time[[zs]], fit$surv[[zs]],
                        xout=if(length(times)==1)times else times[this],
                        method='constant', f=0)$y
      est.surv[this] <-
        estsurv ^ exp(if(length(linear.predictors)==1)
                      linear.predictors else
                      linear.predictors[this])
    }
    return(est.surv)
  }
  
  if(n>1 && nt==0)
    stop("must specify times if getting predictions for >1 obs.")
  if(nt==0)	{  #Get est for 1 obs
	if(!is.list(fit$time)) {
      times <- fit$time
      surv <- fit$surv^exp(linear.predictors)
      std.err <- fit$std.err
    }	else {
      times <- fit$time[[strata]]
      surv <- fit$surv[[strata]]^exp(linear.predictors)
      std.err <- fit$std.err[[strata]]
    }
	if(conf.int>0) {
      lower <- surv^exp(zcrit*std.err)
      upper <- surv^exp(-zcrit*std.err)
      lower[1] <- 1
      upper[1] <- 1
      attr(lower,"type") <- NULL
      attr(upper,"type") <- NULL
    }
    surv <- fun(surv); surv[is.infinite(surv)] <- NA
    if(conf.int>0) {
      lower <- fun(lower); lower[is.infinite(lower)] <- NA
      upper <- fun(upper); upper[is.infinite(upper)] <- NA
    }
    if(!.R.) {
      storage.mode(times) <- "single"
      storage.mode(surv) <- "single"
    }
	if(conf.int>0 && !.R.) {      
      storage.mode(std.err) <- "single"
      storage.mode(lower) <- "single"
      storage.mode(upper) <- "single"
    }
	if(!.R.) storage.mode(linear.predictors) <- "single"
	if(nf==0) strata <- NULL
	retlist <- list(time=times,surv=surv,linear.predictors=linear.predictors)
	if(conf.int>0) retlist <- 
      c(retlist,list(lower=lower,upper=upper,std.err=std.err))
	if(nf>0) {
      retlist$strata <- strata
      retlist$requested.strata <- oldUnclass(strata)
    }
	return(retlist)
  }

  ##Selected times for >=1 obs
  ##if(conf.int>0)
  ## warning("for >1 observation being predicted, does not store confidence intervals")
  ##First get survival at times "times" for each stratum

  surv <- matrix(if(.R.)double(1) else single(1),nrow=ns,ncol=nt)
  serr <- matrix(if(.R.)double(1) else single(1),nrow=ns,ncol=nt)
  for(i in 1:ns) {
	if(!is.list(fit$time)) {
      tim <- fit$time
      se  <- fit$std.err
      srv <- fit$surv
    }	else {
      tim <- fit$time[[i]]
      se  <- fit$std.err[[i]]
      srv <- fit$surv[[i]]
    }
	m <- length(tim)
	j <- 0
	for(u in times)	{
      j <- j + 1
      ##		tm <- max((1:length(tim))[max(tim[tim<=u])==tim])
      tm <- max((1:length(tim))[tim<=u+1e-6])
      s <- srv[tm]; Se <- se[tm]
      ##		if(s != min(srv[tim<=u])) warning()
      if(u > tim[m] && srv[m]>0) {s <- NA; Se <- NA}
      surv[i,j] <- s; serr[i,j] <- Se
    }
  }
  srv <- surv[strata,]^exp(linear.predictors)		#was surv[strata,1:nt]
  ft <- format(times)
  if(is.matrix(srv)) {
    dn <- list(rnam, ft)
    dimnames(srv) <- dn	} else names(srv) <- if(n==1) ft else rnam
  if(conf.int>0) {
    serr <- serr[strata,]
    lower <- srv^exp(zcrit*serr); upper <- srv^exp(-zcrit*serr)
    if(is.matrix(lower)) {
      dimnames(serr) <- dn; dimnames(lower) <- dn; dimnames(upper) <- dn
    }
    else {
      names(serr) <- names(lower) <- names(upper) <- if(n==1) ft else rnam
    }
    lower <- fun(lower); lower[is.infinite(lower)] <- NA
    upper <- fun(upper); upper[is.infinite(upper)] <- NA
    if(!.R.) {
      storage.mode(serr)  <- "single"
      storage.mode(lower) <- "single"
      storage.mode(upper) <- "single"
    }
  }

  srv <- fun(srv); srv[is.infinite(srv)] <- NA
  if(!.R.) storage.mode(srv) <- "single"
  if(conf.int==0)
    return(list(time=times, surv=naresid(naa,srv))) #was return(srv)
  retlist <- list(time=times, surv=naresid(naa,srv), lower=naresid(naa,lower),
                  upper=naresid(naa,upper), std.err=naresid(naa,serr))
  if(nf>0) retlist$requested.strata <- naresid(naa,oldUnclass(strata))
  retlist
}
#Use x= if input is a design matrix, newdata= if a data frame or data matrix
#or vector.  Can specify (centered) linear predictor values instead (linear.predictors).
#Strata is attached to linear.predictors or x as "strata" attribute.
#data matrix assumes that categorical variables are coded with integer codes

survest.psm <- function(fit, newdata, linear.predictors, x, times, fun, 
                        loglog=FALSE, conf.int=.95,
                        what=c("survival","hazard","parallel"), 
                        ...) {   # ... so survplot will work

  what <- match.arg(what)
  if(what=='parallel') conf.int <- FALSE

  trans <- switch(what, survival=Survival(fit), hazard=Hazard(fit),
                  parallel=Survival(fit))

  if(missing(fun)) fun <- if(loglog) function(x) logb(ifelse(x==0|x==1,NA,x))
  else function(x) x

  if(what=="hazard" & conf.int>0) {
    warning('conf.int ignored for what="hazard"')
    conf.int <- FALSE
  }

  if(conf.int>0)					{
    cov <- fit$var
    i <- 1:length(fit$coef); cov <- cov[i,i,drop=FALSE] #ignore scale for now
    if(!missing(linear.predictors))	{
      warning("conf.int set to 0 since linear.predictors specified")
      conf.int <- 0
    }
  }

  if(any(attr(fit,'class')=='pphsm'))
    stop("fit should not have been passed thru pphsm")

  nvar <- length(fit$coef)-num.intercepts(fit)
  
  if(missing(linear.predictors))					{
    if(nvar>0 & missing(x) & missing(newdata))		{
      linear.predictors <- fit$linear.predictors
      if(conf.int>0) stop("may not specify conf.int unless x or newdata given")
      rnam <- names(linear.predictors)
    } else				{
      if(nvar==0)	{
        x <- as.matrix(1)	#no predictors
        linear.predictors <- fit$coef[1]  #changed 15May97
      }
      else	{
        if(missing(x)) x <- predict(fit,newdata,type="x")
        linear.predictors <- matxv(x,fit$coef)
      }
      if(conf.int>0)				{
        g1 <- drop(((x %*% cov) * x) %*% rep(1, ncol(x)))
        last <- if(.newSurvival.) {
          nscale <- length(fit$icoef) - 1
          ncol(fit$var)-(1:nscale)+1
        } else ncol(fit$var)
        g2 <- drop(x %*% fit$var[-last,last,drop=FALSE])
      }
      rnam <- dimnames(x)[[1]]
    }
  }
  else  rnam <- names(linear.predictors)

  if(what=='parallel') {   ## 18Sep99
    if(length(times)>1 && (length(times) != length(linear.predictors)))
      stop('length of times must = 1 or number of subjects when what="parallel"')
    return(trans(times,linear.predictors))
  }

  if(missing(times)) times <- seq(0,fit$maxtime,length=200)
  nt <- length(times)
  n <- length(linear.predictors)
  
  if(n>1 & missing(times))
    warning("should specify times if getting predictions for >1 obs.")
  
  if(conf.int>0) zcrit <- qnorm((conf.int+1)/2)
  
  comp <- function(a, b, Trans) Trans(b, a)
  surv <- drop(outer(linear.predictors, times, FUN=comp, Trans=trans)) #1Apr95

  if(conf.int>0 && (nt==1 | n==1))	{
    if(.newSurvival.) {
      dist <- fit$dist
      link <- if(!.R.) survReg.distributions[[dist]]$trans else
                       survreg.distributions[[dist]]$trans
      z <- if(length(link)) link(times) else times
      sc <- fit$scale   ## TODO: generalize for vector of scale params
      logtxb <- outer(linear.predictors,z,function(a,b)b-a)
      se <- sqrt(g1 + logtxb * (2*g2 + logtxb * fit$var[last,last]))/sc
      prm <- 0
      tm <- if(length(link)) 1 else 0
    } else {
      z <- glm.links["link", fit$family[2]]$link(times)
      sc <- exp(fit$parms[1])
      logtxb <- outer(linear.predictors,z,function(a,b)b-a)
      se <- sqrt(g1 + logtxb * (2*g2 + logtxb * fit$var[last,last]))/sc
      prm <- fit$parms   ## fool survreg.auxinfo$survival into
      ## using -linear.predictor as whole inner argument
      prm[1] <- 0        # anti-log of 0 = 1
      tm <- glm.links["inverse", fit$family[2]]$inverse(0)
      ## e.g. tm=1 if using log(t), 0 if using t
    }
    lower <- trans(tm,-drop(logtxb/sc + zcrit*se),parms=prm)
    upper <- trans(tm,-drop(logtxb/sc - zcrit*se),parms=prm)
    if(what=='survival') {
      lower[times==0] <- 1
      upper[times==0] <- 1
    }
    std.err <- drop(se)
    if(!.R.) {
      storage.mode(lower) <- "single"
      storage.mode(upper) <- "single"
      storage.mode(std.err) <- "single"
    }
  }

  if(!.R.) {
    storage.mode(times)  <- "single"
    storage.mode(surv)   <- "single"
    storage.mode(linear.predictors) <- "single"
  }

  if(nt==1 | n==1){
    surv <- fun(surv); surv[is.infinite(surv)] <- NA
    if(conf.int>0) {
      lower <- fun(lower); lower[is.infinite(lower)] <- NA
      upper <- fun(upper); upper[is.infinite(upper)] <- NA
      retlist <- list(time=times,surv=surv,
                      lower=lower,upper=upper,
                      std.err=std.err,
                      linear.predictors=linear.predictors)
    }
    else retlist <- list(time=times,surv=surv,
                         linear.predictors=linear.predictors)
    retlist <- structure(c(retlist,
                           list(conf.int=conf.int,units=fit$units,
                                n.risk=fit$stats["Obs"],
                                n.event=fit$stats["Events"], what=what)),
                         class='survest.psm')
    return(retlist)
  }
  
  if(n==1) names(surv) <- format(times) else		{
    if(is.matrix(surv))
      dimnames(surv) <- list(rnam, format(times))
    else names(surv) <- rnam
  }
  
  if(!.R.) storage.mode(surv) <- "single"
  surv
}
		

print.survest.psm <- function(x, ...) {
  cat('\nN:',x$n.risk,'\tEvents:',x$n.event)
  z <- if(length(unique(x$time)) > 1) data.frame(Time=x$time) else {
    cat('\tTime:',x$time[1],' ',x$units,'s',sep='')
    data.frame(LinearPredictor=x$linear.predictors)
  }
  cat('\n\n')
  z$whatever  <- x$surv
  names(z)[2] <- x$what
  if(x$conf.int) {
    z$Lower <- x$lower
    z$Upper <- x$upper
    z$SE    <- x$std.err
  }
  print.data.frame(z)
  invisible()
}


	


survest <- function(fit, ...) UseMethod("survest")
#SCCS  @(#)survfit.coxph.null.s	4.5 % G%
#Set score to 0 if linear predictors not there - FEH 28 Sep 93
survfit.cph.null <-
  function(object, newdata, se.fit=TRUE, conf.int=.95, individual=FALSE,
	    type, vartype,
	    conf.type=c('log-log', 'log', 'plain', 'none'), ...) {

    if(.R.) require('survival')
    ## Sense whether survival5 is in effect and if so use this later version
    s5 <- exists('coxpenal.fit')

    # May have strata and/or offset terms, linear predictor = offset
    #  newdata doesn't make any sense
    #  This is survfit.coxph with lots of lines removed
    y <- object$y
    if(is.null(y))stop("must use y=T with fit")
    n <- nrow(y)
    Strata <- attr(y,"strata")
    if(length(object$linear.predictor)==0) score <- rep(1, n)
    else score <- exp(object$linear.predictor)
    
    temp <- c("aalen", "kalbfleisch-prentice", "efron", "tsiatis", "breslow", 
              "kaplan-meier", "fleming-harringon", "greenwood", "exact")
    temp2 <- c(2, 1, 3, 2, 2, 1, 3, 1, 1)
    if(missing(type)) type <- object$method
    if(missing(vartype)) vartype <- type
    method <- temp2[match(match.arg(type, temp), temp)]
    if(is.na(method)) stop("Invalid survival curve type")
    vartype <- temp2[match(match.arg(vartype, temp), temp)]
    if(is.na(vartype)) stop("Invalid variance type specified")
    km <- method==1
    
    if (!se.fit) conf.type <- 'none'
    else conf.type <- match.arg(conf.type)

    ny <- ncol(y)
    type <- attr(y, 'type')
    if (type=='counting') {
      ord <- order(Strata, y[,2], -y[,3])
      ##if (method=='kaplan-meier')   bug correction FEH 6Jun99
      if(km) stop ("KM method not valid for counting type data")
	}
    else if (type=='right') {
      ord <- order(Strata, y[,1], -y[,2])
      miny <- min(y[, 1])
      if(miny < 0) y <- cbind(2 * miny - 1, y)
      else y <- cbind(-1, y)
	}
    else stop("Cannot handle \"", type, "\" type survival data")

    if (!is.null(Strata)) {
	newstrat <- (as.numeric(Strata))[ord]
	newstrat <- as.integer(c(1*(diff(newstrat)!=0), 1))
	}
    else newstrat <- as.integer(c(rep(0,n-1),1))

    if ( !missing(newdata))
	stop("A newdata argument does not make sense for a null model")
    dimnames(y) <- NULL   #I only use part of Y, so names become invalid
    storage.mode(y) <- 'double'
    surv <- if(.R.)
      .C('agsurv2',
         as.integer(n),
         as.integer(0),
         y = y[ord,],
         as.double(score[ord]),
         strata = as.integer(newstrat),
         surv = double(n),
         varhaz = double(n),
         double(1),
         as.double(0),
         nsurv = as.integer(c(method,vartype)),
         double(2),
         as.integer(1),
         double(1),
         newrisk= as.double(1), PACKAGE="survival") else
      .C(if(.SV4.)'S_agsurv2' else 'agsurv2', ##14Nov00
         as.integer(n),
         as.integer(0),
         y = y[ord,],
         as.double(score[ord]),
         strata = as.integer(newstrat),
         surv = double(n),
         varhaz = double(n),
         double(1),
         as.double(0),  ##14Nov00
         nsurv = as.integer(if(s5) c(method,vartype) else km),
         double(2),
         as.integer(1),
         double(1),
         newrisk= as.double(1))

    nsurv <- surv$nsurv[1]
    ntime <- 1:nsurv
    tsurv <- surv$surv[ntime]
    tvar  <- surv$varhaz[ntime]
    if (surv$strata[1] <=1)
	temp <- list(time=surv$y[ntime,1],
		 n.risk=surv$y[ntime,2],
		 n.event=surv$y[ntime,3],
		 surv=tsurv )
    else {
	temp <- surv$strata[1:(1+surv$strata[1])]
	tstrat <- diff(c(0, temp[-1])) #n in each strata
	names(tstrat) <- levels(Strata)
	temp <- list(time=surv$y[ntime,1],
		 n.risk=surv$y[ntime,2],
		 n.event=surv$y[ntime,3],
		 surv=tsurv,
		 strata= tstrat)
	}
    if (se.fit) temp$std.err <- sqrt(tvar)

    zval <- qnorm(1- (1-conf.int)/2, 0,1)
    if (conf.type=='plain') {
	temp1 <- temp$surv + zval* temp$std * temp$surv
	temp2 <- temp$surv - zval* temp$std * temp$surv
	temp <- c(temp, list(upper=pmin(temp1,1), lower=pmax(temp2,0),
			conf.type='plain', conf.int=conf.int))
	}
    if (conf.type=='log') {
	xx <- ifelse(temp$surv==0,1,temp$surv)  #avoid some "log(0)" messages
	temp1 <- ifelse(temp$surv==0, 0*temp$std, exp(logb(xx) + zval* temp$std))
	temp2 <- ifelse(temp$surv==0, 0*temp$std, exp(logb(xx) - zval* temp$std))
	temp <- c(temp, list(upper=pmin(temp1,1), lower=temp2,
			conf.type='log', conf.int=conf.int))
	}
    if (conf.type=='log-log') {
	who <- (temp$surv==0 | temp$surv==1) #special cases
	xx <- ifelse(who, .1,temp$surv)  #avoid some "log(0)" messages
	temp1 <- exp(-exp(logb(-logb(xx)) + zval*temp$std/logb(xx)))
	temp1 <- ifelse(who, temp$surv + 0*temp$std, temp1)
	temp2 <- exp(-exp(logb(-logb(xx)) - zval*temp$std/logb(xx)))
	temp2 <- ifelse(who, temp$surv + 0*temp$std, temp2)
	temp <- c(temp, list(upper=temp1, lower=temp2,
			conf.type='log-log', conf.int=conf.int))
	}

    temp$call <- call
#    if(!is.null(strata)) attr(temp, "strata") <- Strata
    oldClass(temp) <- if(.SV4.)'survfit.cox' else
     c("survfit.cph", "survfit.cox", "survfit")
    temp
    }
#SCCS @(#)survfit.coxph.s	4.5 7/14/92
survfit.cph <-
  function(object, newdata, se.fit=TRUE, conf.int=.95, 
           individual=FALSE,
           type, vartype,
           conf.type=c('log-log', 'log', 'plain', 'none'))
  {
    if(.R.) require('survival')
    ## Sense whether survival5 is in effect and if so use this later version
    s5 <- exists('coxpenal.fit')

	if(length(object$weights) || length(object$call$weights)) 
      stop('survfit cannot yet handle a weighted model')

    nvar <- length(object$coef)
    score <- exp(object$linear.predictors)
    temp <- c("aalen", "kalbfleisch-prentice", "efron", "tsiatis", "breslow", 
              "kaplan-meier", "fleming-harringon", "greenwood", "exact")
    temp2 <- c(2, 1, 3, 2, 2, 1, 3, 1, 1)
    if(missing(type)) type <- object$method
    if(missing(vartype)) vartype <- type
    method <- temp2[match(match.arg(type, temp), temp)]
    if(is.na(method)) stop("Invalid survival curve type")
    vartype <- temp2[match(match.arg(vartype, temp), temp)]
    if(is.na(vartype)) stop("Invalid variance type specified")

    km <- method==1
    coxmethod <- object$method
    if (!se.fit) conf.type <- 'none'
    else conf.type <- match.arg(conf.type)

    x <- object$x
    y <- object$y
    if(is.null(x) | is.null(y)) stop("you must specify x=T and y=T in the fit")
    n <- nrow(y)
    stratum <- attr(x, "strata")
    offset <- attr(x, "offset")
    if(length(offset)==0) offset <- rep(0,n)
    type <- attr(y, 'type')
    if (type=='counting') {
      if(km) stop ("KM method not valid for counting type data")
      ## if (method=='kaplan-meier')  bug correction FEH 6Jun99
      ord <- if(length(stratum)) order(stratum, y[,2], -y[,3]) else
             order(y[,2], -y[,3])
	}
    else if (type=='right') {
      ord <- if(length(stratum)) order(stratum, y[, 1],  - y[, 2])
             else order(y[, 1],  - y[, 2])
      ## length(stratum) was length(strat) 14Nov00
    miny <- min(y[, 1])
    if(miny < 0)
      y <- cbind(2 * miny - 1, y)
    else y <- cbind(-1, y)
	}
    else stop("Cannot handle \"", type, "\" type survival data")

    x <- x[ord,]
    weights <- if(length(object$weights)) object$weights[ord] else rep(1,n)

    if (length(stratum)) {
      newstrat <- (as.numeric(stratum))[ord]
      newstrat <- as.integer(c(1*(diff(newstrat)!=0), 1))
	}
    else newstrat <- as.integer(c(rep(0,n-1),1))

    if (individual && !missing(newdata)) stype <- 1
    else stype <- 2
    if(stype == 1 && method != vartype)
      stop("The type and vartype args must agree for individual=T")
    if(stype == 1 && method == 1)
      stop("Only Aalen and F-H estimates available for individual=T")
    
    if (!missing(newdata)) {
      x2 <- predictDesign(object,newdata,type="x",expand.na=FALSE)
      n2 <- nrow(x2)
      offset2 <- attr(x2, "offset")
      if(length(offset2)==0) offset2 <- rep(0, n2)
      strata2 <- attr(x2,"strata")
      if(length(strata2)==0) strata2 <- rep(1, n2)
      naa <- attr(x2, "na.action")
      if(stype==1) {
		Terms <- terms(attr(object$terms,'formula'))
		m2 <- model.frame(Terms,newdata)    #,Des=F) 17Jul01
		y2 <- model.extract(m2, 'response')
		if(length(y2)==0) stop('no Surv object in newdata')
		if(attr(y2,'type')!=type)
		  stop('Survival type of newdata does not match fitted model')
		if(nrow(y2) != n2) stop('wrong # rows in Y')
      }
    }

    else 	{
      x2 <- matrix(object$means, nrow=1)
      offset2 <- 0; strata2 <- 0; n2 <- 1
    }
    coef <- ifelse(is.na(object$coef), 0, object$coef)
    newrisk <- exp(c(x2 %*% coef) + offset2 - object$center)

    dimnames(y) <- NULL   #I only use part of Y, so names become invalid
    storage.mode(x) <- 'double'
    storage.mode(y) <- 'double'
    ndead <- sum(y[,3])
    
    if (stype==1) {
      surv <- if(.R.)
        .C("agsurv1",
           as.integer(n),
           as.integer(nvar),
           y[ord,],
           as.double(score[ord]),
           strata=as.integer(newstrat),
           surv=double(ndead*n2),
           varhaz=double(ndead*n2),
           nsurv=if(s5)as.integer(method==3) else
           as.integer(2+1*(coxmethod=="efron")),
           x,
           double(3*nvar),
           as.double(object$var),
           y = double(3*n*n2),
           as.integer(n2),
           as.double(y2),
           as.double(x2),
           as.double(newrisk),
           as.integer(strata2), PACKAGE="survival") else
        .C(if(.SV4.)'S_agsurv1' else "agsurv1", ##14Nov00
           as.integer(n),
           as.integer(nvar),
           y[ord,],
           as.double(score[ord]),   ## was a bug in survival4 as.double14Nov00
           strata=as.integer(newstrat),
           surv=double(ndead*n2),   # was bug in surv4
           varhaz=double(ndead*n2), # was bug in surv4
           nsurv=if(s5)as.integer(method==3) else
           as.integer(2+1*(coxmethod=="efron")),
           x,
           double(3*nvar),
           as.double(object$var),
           y = double(3*n*n2),
           as.integer(n2),
           as.double(y2),
           as.double(x2),
           as.double(newrisk),
           as.integer(strata2) )
      ntime <- 1:surv$nsurv
      temp <- (matrix(surv$y, ncol=3))[ntime,]
      temp <- list(time = temp[,1],
                   n.risk= temp[,2],
                   n.event=temp[,3],
                   surv = surv$surv[ntime])
      if (se.fit) temp$std.err <- sqrt(surv$varhaz[ntime])
    }
    else {
      if(!s5) temp <- ifelse(km, 1, 2+as.integer(coxmethod=="efron"))

      surv <- .C(if(.SV4.)'S_agsurv2' else 'agsurv2',  ##14Nov00
                 as.integer(n),
                 as.integer(nvar* se.fit),
                 y = y[ord,],
                 as.double(score[ord]),
                 strata = as.integer(newstrat),
                 surv   = double(ndead*n2),   # was bug in surv4
                 varhaz = double(ndead*n2),   # was bug in surv4
                 x,
                 as.double(object$var),
                 nsurv = if(s5) as.integer(c(method,vartype))
                         else as.integer(temp),
                 double(3*nvar),
                 as.integer(n2),
                 as.double(x2),
                 as.double(newrisk))

      nsurv <- surv$nsurv[1]
      ntime <- 1:nsurv
      if (n2>1) {
	    tsurv <- matrix(surv$surv[1:(nsurv*n2)], ncol=n2)
	    tvar  <- matrix(surv$varhaz[1:(nsurv*n2)], ncol=n2)
	    dimnames(tsurv) <- list(NULL, dimnames(x2)[[1]])
      }
      else {
	    tsurv <- surv$surv[ntime]
	    tvar  <- surv$varhaz[ntime]
      }
      if (surv$strata[1] <=1)
	    temp <- list(time=surv$y[ntime,1],
                    n.risk=surv$y[ntime,2],
                    n.event=surv$y[ntime,3],
                    surv=tsurv )
      else {
	    temp <- surv$strata[1:(1+surv$strata[1])]
	    tstrat <- diff(c(0, temp[-1])) #n in each strata
	    names(tstrat) <- levels(stratum)
	    temp <- list(time=surv$y[ntime,1],
                    n.risk=surv$y[ntime,2],
                    n.event=surv$y[ntime,3],
                    surv=tsurv,
                    strata= tstrat)
      }
      if (se.fit) temp$std.err <- sqrt(tvar)
	}
    zval <- qnorm(1- (1-conf.int)/2, 0,1)
    if (conf.type=='plain') {
      temp1 <- temp$surv + zval* temp$std * temp$surv
      temp2 <- temp$surv - zval* temp$std * temp$surv
      temp <- c(temp, list(upper=pmin(temp1,1), lower=pmax(temp2,0),
                           conf.type='plain', conf.int=conf.int))
	}
    if (conf.type=='log') {
      xx <- ifelse(temp$surv==0,1,temp$surv)  #avoid some "log(0)" messages
      temp1 <- ifelse(temp$surv==0, 0*temp$std, exp(logb(xx) + zval* temp$std))
      temp2 <- ifelse(temp$surv==0, 0*temp$std, exp(logb(xx) - zval* temp$std))
      temp <- c(temp, list(upper=pmin(temp1,1), lower=temp2,
                           conf.type='log', conf.int=conf.int))
	}
    if (conf.type=='log-log') {
      who <- (temp$surv==0 | temp$surv==1) #special cases
      xx <- ifelse(who, .1,temp$surv)  #avoid some "log(0)" messages
      temp1 <- exp(-exp(logb(-logb(xx)) + zval*temp$std/logb(xx)))
      temp1 <- ifelse(who, temp$surv + 0*temp$std, temp1)
      temp2 <- exp(-exp(logb(-logb(xx)) - zval*temp$std/logb(xx)))
      temp2 <- ifelse(who, temp$surv + 0*temp$std, temp2)
      temp <- c(temp, list(upper=temp1, lower=temp2,
                           conf.type='log-log', conf.int=conf.int))
	}
    
    temp$call <- call
    if(!missing(newdata)) 
      {
        temp$requested.strata <- strata2
        attr(temp, "na.action") <- naa
      }
    oldClass(temp) <- if(.SV4.) 'survfit.cox' else
     c("survfit.cph", "survfit.cox", "survfit")  ##14Nov00
    temp
  }
##modified version of Therneau's survfit - keeps attributes of Surv
##object and uses interaction() to form strata labels

survfit <- function (formula, data, weights, subset, na.action=na.delete, 
                     conf.type=c("log-log","log","plain","none"),...) {
  call <- match.call()
  conf.type <- match.arg(conf.type)
  if(.R.) require('survival')
  ## Real tricky -- find out if the first arg is "Surv(...)" without
  ##  evaluating it.  If this is so, or it is a survival object, turn it
  ##  into a formula
  if ((mode(call[[2]]) == 'call' &&  call[[2]][[1]] == as.name('Surv'))
      || inherits(formula, 'Surv'))  {
	## The dummy function stops an annoying warning message "Looking for
	##  'formula' of mode function, ignored one of mode ..."
	xx <- function(x) formula(x)
	formula <- xx(paste(deparse(call[[2]]), 1, sep="~"))
  }

  if (!inherits(formula, 'formula')) temp <- UseMethod("survfit")
  else {
	m <- match.call(expand=FALSE)
	m$conf.type <- m$... <- NULL

	Terms <- terms(formula, 'strata')
	ord <- attr(Terms, 'order')
	if (length(ord) & any(ord !=1))
      stop("Interaction terms are not valid for this function")
	m$formula <- Terms
#    m$Des <- F   #turn off Design()   3Jun99  17Jul01
	m[[1]] <- as.name("model.frame")
	m <- eval(m, sys.parent())

	n <- nrow(m)
	Y <- model.extract(m, response)
	casewt <- model.extract(m, "weights")
	## The second line below works around a bug in Splus 3.0.1, which later
	##    went away, i.e., casewt is returned as an unevaluated arg.
	if (is.null(casewt)) casewt <- rep(1,n)
	else if (mode(casewt)=='argument') casewt <- eval(casewt[[1]])

	if(length(attr(Terms,'offset'))) warning('Offset term ignored')

	ll <- attr(Terms, 'term.labels')
	if (length(ll) == 0) X <- factor(rep(1,n))
	else {
      temp <-  rep(1, length(ll))
      strat <- untangle.specials(Terms, 'strata',1)$terms
      if (length(strat)) temp[strat] <- 0
      X <- m
      X[[1]] <- NULL
      X <- interaction(X,drop=TRUE,sep=" ")
    }

	temp <- survfit.km(X, Y, casewt, conf.type=conf.type, ...)
	attr(temp, "class") <- "survfit"
	if (!is.null(attr(m, 'na.action'))) temp$na.action <- attr(m, 'na.action')
	ny <- ncol(Y)
	temp$maxtime <- max(Y[,ny-1])
	temp$units <- attr(Y,"units")
	temp$time.label <- attr(Y,"time.label")
  }
  temp$call <- call
  temp
}
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

survplot <- function(fit, ...) UseMethod("survplot")
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


# SCCS @(#)survreg.distributions.s	4.3 11/19/92
#
# Create the survreg.distributions object
#
# Infinite mean in log logistic courtesy of Victor Moreno
# SERC, Institut Catala d'Oncologia  (V.Moreno@ico.scs.es)  9Feb98

# survival package defines basic quantile function ignoring link
# Actual quantile function called Quantile here, for SV4 or R

survreg.auxinfo <- if(.SV4. || .R.) list(
exponential = list(
    survival = function(times, lp, parms) exp(-times/exp(lp)),
    hazard = function(times, lp, parms) exp(-lp),
    Quantile = function(q=.5, lp, parms) {
		names(parms) <- NULL
		f <- function(lp, q, parms) -logb(1-q)*exp(lp)
		names(q) <- format(q)
		drop(outer(lp, q, FUN=f, parms=parms))
		},
    mean = function(lp, parms) exp(lp),
    latex = function(...) '\\exp(-t/\\exp(X\\beta))'
  ),
  
extreme = list(
    survival = function(times, lp, parms) { 
		exp(-exp((times-lp)/exp(parms)))
		},
    hazard = function(times, lp, parms) {
		scale <- exp(parms[1])   #14Jun97
		exp((times-lp)/scale)/scale
		},
    Quantile = function(q=.5, lp, parms) {
		names(parms) <- NULL
		f <- function(lp, q, parms) lp + exp(parms)*logb(-logb(1-q))
		names(q) <- format(q)
		drop(outer(lp, q, FUN=f, parms=parms))
		},
    mean = function(lp, parms) {
		names(parms) <- NULL
		lp-.57722*exp(parms)
		},
    latex = function(scale) {
		yvar <- "t"
		z <- if(scale==1) paste(yvar,"-X\\beta") else paste(
			"\\frac{", yvar, "-X\\beta}{",format(scale),"}",sep="")
		z <- paste("\\exp[-\\exp(",z,")]")
		z
		}
    ),

weibull = list(
    survival = function(times, lp, parms) { 
		t.trans <- logb(times)
		names(t.trans) <- format(times)
		exp(-exp((t.trans-lp)/exp(parms)))
		},
    hazard = function(times, lp, parms) {
		t.trans <- logb(times)
		t.deriv <- 1/times
		names(t.trans) <- format(times)
		scale <- exp(parms[1])   #14Jun97
		ifelse(times==0,exp(-lp/scale)/scale,
                        exp((t.trans-lp)/scale)*t.deriv/scale)
		},
    Quantile = function(q=.5, lp, parms) {
		names(parms) <- NULL
		f <- function(lp, q, parms) lp + exp(parms)*logb(-logb(1-q))
		names(q) <- format(q)
		drop(exp(outer(lp, q, FUN=f, parms=parms)))
		},
    mean = function(lp, parms, transform) {
		names(parms) <- NULL
		exp(lp)*gamma(exp(parms)+1)
		},
    latex = function(scale) {
		yvar <- "\\log(t)"
		z <- if(scale==1) paste(yvar,"-X\\beta") else paste(
			"\\frac{", yvar, "-X\\beta}{",format(scale),"}",sep="")
		z <- paste("\\exp[-\\exp(",z,")]")
		z
		}
    ),
                    
logistic = list(
    survival = function(times, lp, parms) { 
		1/(1+exp((times-lp)/exp(parms)))
		},
    hazard = function(times, lp, parms) {
		scale <- exp(parms)
		1/scale/(1+exp(-(times-lp)/scale))
		},
    Quantile = function(q=.5, lp, parms) {
		names(parms) <- NULL
		f <- function(lp, q, parms) lp + exp(parms)*logb(q/(1-q))
		names(q) <- format(q)
		drop(outer(lp, q, FUN=f, parms=parms))
		},
    mean = function(lp, parms) lp,
    latex = function(scale){
		yvar <- "t"
		z <- if(scale==1) paste(yvar,"-X\\beta") else paste(
			"\\frac{", yvar, "-X\\beta}{",format(scale),"}",sep="")
		z <- paste("[1+\\exp(",z,")]^{-1}")
		z
		}
    ),

loglogistic = list(
    survival = function(times, lp, parms) { 
		1/(1+exp((logb(times)-lp)/exp(parms)))
		},
    hazard = function(times, lp, parms) {
		t.trans <- logb(times)
		t.deriv <- 1/times
		scale <- exp(parms)
		names(t.trans) <- format(times)
		t.deriv/scale/(1+exp(-(t.trans-lp)/scale))
		},
    Quantile = function(q=.5, lp, parms) {
		names(parms) <- NULL
		f <- function(lp, q, parms) lp + exp(parms)*logb(q/(1-q))
		names(q) <- format(q)
		drop(exp(outer(lp, q, FUN=f, parms=parms)))
		},
    mean = function(lp, parms) {
		names(parms) <- NULL
		if(exp(parms)>1) rep(Inf,length(lp)) else
			   exp(lp)*pi*exp(parms)/sin(pi*exp(parms))
		},
    latex = function(scale) {
		yvar <- "\\log(t)"
		z <- if(scale==1) paste(yvar,"-X\\beta") else paste(
			"\\frac{", yvar, "-X\\beta}{",format(scale),"}",sep="")
		z <- paste("[1+\\exp(",z,")]^{-1}")
		z
		}),
    
gaussian = list(
    survival = function(times, lp, parms) 1-pnorm((times-lp)/exp(parms)),
    hazard = function(times, lp, parms) {
		scale <- exp(parms)
		z <- (times-lp)/scale
		dnorm(z)/scale/(1-pnorm(z))
		},
    Quantile = function(q=.5, lp, parms) {
		names(parms) <- NULL
		f <- function(lp, q, parms) lp + exp(parms)*qnorm(q)
		names(q) <- format(q)
		drop(outer(lp, q, FUN=f, parms=parms))
		},
    mean = function(lp, parms) lp,
    latex = function(scale) {
		yvar <- "t"
		z <- if(scale==1) paste(yvar,"-X\\beta") else paste(
			"\\frac{", yvar, "-X\\beta}{",format(scale),"}",sep="")
		z <- paste("1-\\Phi(",z,")")
		z
		}
    ),

lognormal = list(
    survival = function(times, lp, parms) { 
		t.trans <- logb(times)
		names(t.trans) <- format(times)
		1-pnorm((t.trans-lp)/exp(parms))
		},
    hazard = function(times, lp, parms) {
		t.trans <- logb(times)
		t.deriv <- 1/times
		scale <- exp(parms)
		names(t.trans) <- format(times)
		z <- (t.trans-lp)/scale
		t.deriv*dnorm(z)/scale/(1-pnorm(z))
		},
    Quantile = function(q=.5, lp, parms) {
		names(parms) <- NULL
		f <- function(lp, q, parms) lp + exp(parms)*qnorm(q)
		names(q) <- format(q)
		drop(exp(outer(lp, q, FUN=f, parms=parms)))
		},
    mean = function(lp, parms) {
		names(parms) <- NULL
		exp(lp+exp(2*parms)/2)
		},
    latex = function(scale) {
		yvar <- "\\log(t)"
		z <- if(scale==1) paste(yvar,"-X\\beta") else paste(
			"\\frac{", yvar, "-X\\beta}{",format(scale),"}",sep="")
		z <- paste("1-\\Phi(",z,")")
		z
		}
    ),
  
t = list(
    survival = function(times, lp, parms) {
		scale <- exp(parms[1])
		df <- parms[2]
		1-pt((times-lp)/scale,df)
		},
    hazard = function(times, lp, parms) {
		scale <- exp(parms[1])
		df <- parms[2]
		z <- (times-lp)/scale
		dt(z,df)/scale/(1-pt(z,df))
		},
    Quantile = function(q=.5, lp, parms) {
		names(parms) <- NULL
		f <- function(lp, q, parms) lp + exp(parms[1])*qt(q, parms[2])
		names(q) <- format(q)
		drop(outer(lp, q, FUN=f, parms=parms))
		},
    mean = function(lp, parms) lp,
    latex = function(scale,df) {
		yvar <- "t"
		z <- if(scale==1) paste(yvar,"-X\\beta") else paste(
			"\\frac{", yvar, "-X\\beta}{",format(scale),"}",sep="")
		z <- paste("1-T_{",df,"}(",z,")", sep="")
		z
      }
  )
 ) else list(
'extreme' = list(
    survival = function(times, lp, parms, transform) { 
		t.trans <- glm.links["link",transform]$link(times)
		names(t.trans) <- format(times)
		exp(-exp((t.trans-lp)/exp(parms)))
		},
    survival.inverse = function(q, parms=0) logb(-logb(q))*exp(parms),
    hazard = function(times, lp, parms, transform) {
		t.trans <- glm.links["link",transform]$link(times)
		t.deriv <- glm.links["deriv", transform]$deriv(times)
		names(t.trans) <- format(times)
		scale <- exp(parms[1])   #14Jun97
		exp((t.trans-lp)/scale)*t.deriv/scale
		},
    quantile = function(q=.5, lp, parms, transform) {
		names(parms) <- NULL
		inv <- glm.links["inverse", transform]$inverse
		f <- function(lp, q, parms) lp + exp(parms)*logb(-logb(1-q))
		names(q) <- format(q)
		drop(inv(outer(lp, q, FUN=f, parms=parms)))
		},
    mean = function(lp, parms, transform) {
		names(parms) <- NULL
		switch(transform, identity=lp-.57722*exp(parms), 
		  log=exp(lp)*gamma(exp(parms)+1), 
		  stop(paste(transform,"not implemented")))
		},
    latex = function(parms, transform) {
		yvar <- switch(transform, identity="t", log="\\log(t)",
			paste(transform,"(t)",sep=""))
		scale <- exp(parms)
		z <- if(scale==1) paste(yvar,"-X\\beta") else paste(
			"\\frac{", yvar, "-X\\beta}{",format(scale),"}",sep="")
		z <- paste("\\exp[-\\exp(",z,")]")
		z
		}
    ),

logistic = list(
    survival = function(times, lp, parms, transform) { 
		t.trans <- glm.links["link",transform]$link(times)
		names(t.trans) <- format(times)
		1/(1+exp((t.trans-lp)/exp(parms)))
		},
    survival.inverse = function(q, parms=0) logb((1-q)/q) * exp(parms),
    hazard = function(times, lp, parms, transform) {
		t.trans <- glm.links["link",transform]$link(times)
		t.deriv <- glm.links["deriv", transform]$deriv(times)
		scale <- exp(parms)
		names(t.trans) <- format(times)
		t.deriv/scale/(1+exp(-(t.trans-lp)/scale))
		},
    quantile = function(q=.5, lp, parms, transform) {
		names(parms) <- NULL
		inv <- glm.links["inverse", transform]$inverse
		f <- function(lp, q, parms) lp + exp(parms)*logb(q/(1-q))
		names(q) <- format(q)
		drop(inv(outer(lp, q, FUN=f, parms=parms)))
		},
    mean = function(lp, parms, transform) {
		names(parms) <- NULL
		switch(transform, identity=lp, 
		  log=if(exp(parms)>1) rep(Inf,length(lp)) else
			   exp(lp)*pi*exp(parms)/sin(pi*exp(parms)),
		  stop(paste(transform,"not implemented")))
		},
    latex = function(parms, transform) {
		yvar <- switch(transform, identity="t", log="\\log(t)",
			paste(transform,"(t)",sep=""))
		scale <- exp(parms)
		z <- if(scale==1) paste(yvar,"-X\\beta") else paste(
			"\\frac{", yvar, "-X\\beta}{",format(scale),"}",sep="")
		z <- paste("[1+\\exp(",z,")]^{-1}")
		z
		}
    ),

gaussian = list(
    survival = function(times, lp, parms, transform) { 
		t.trans <- glm.links["link",transform]$link(times)
		names(t.trans) <- format(times)
		1-pnorm((t.trans-lp)/exp(parms))
		},
    survival.inverse = function(q, parms=0) qnorm(q)*exp(parms),
    hazard = function(times, lp, parms, transform) {
		t.trans <- glm.links["link",transform]$link(times)
		t.deriv <- glm.links["deriv", transform]$deriv(times)
		scale <- exp(parms)
		names(t.trans) <- format(times)
		z <- (t.trans-lp)/scale
		t.deriv*dnorm(z)/scale/(1-pnorm(z))
		},
    quantile = function(q=.5, lp, parms, transform) {
		names(parms) <- NULL
		inv <- glm.links["inverse", transform]$inverse
		f <- function(lp, q, parms) lp + exp(parms)*qnorm(q)
		names(q) <- format(q)
		drop(inv(outer(lp, q, FUN=f, parms=parms)))
		},
    mean = function(lp, parms, transform) {
		names(parms) <- NULL
		switch(transform, identity=lp, log=exp(lp+exp(2*parms)/2), 
		  stop(paste(transform,"not implemented")))
		},
    latex = function(parms, transform) {
		yvar <- switch(transform, identity="t", log="\\log(t)",
			paste(transform,"(t)",sep=""))
		scale <- exp(parms)
		z <- if(scale==1) paste(yvar,"-X\\beta") else paste(
			"\\frac{", yvar, "-X\\beta}{",format(scale),"}",sep="")
		z <- paste("1-\\Phi(",z,")")
		z
		}
    ),

t = list(
    survival = function(times, lp, parms, transform) { 
		t.trans <- glm.links["link",transform]$link(times)
		scale <- exp(parms[1])
		df <- parms[2]
		names(t.trans) <- format(times)
		1-pt((t.trans-lp)/scale,df)
		},
    survival.inverse = function(q, parms=c("log(scale)"=0, df=2))
		qt(q, parms[2])*exp(parms[1]),
    hazard = function(times, lp, parms, transform) {
		t.trans <- glm.links["link",transform]$link(times)
		t.deriv <- glm.links["deriv", transform]$deriv(times)
		scale <- exp(parms[1])
		df <- parms[2]
		names(t.trans) <- format(times)
		z <- (t.trans-lp)/scale
		t.deriv*dt(z,df)/scale/(1-pt(z,df))
		},
    quantile = function(q=.5, lp, parms, transform) {
		names(parms) <- NULL
		inv <- glm.links["inverse", transform]$inverse
		f <- function(lp, q, parms) lp + exp(parms[1])*qt(q, parms[2])
		names(q) <- format(q)
		drop(inv(outer(lp, q, FUN=f, parms=parms)))
		},
    mean = function(lp, parms, transform) {
		names(parms) <- NULL
		switch(transform, identity=lp,  
		  log=stop("mean of log-t distribution does not exist"),
		  stop(paste(transform,"not implemented")))
		},
    latex = function(parms, transform) {
		yvar <- switch(transform, identity="t", log="\\log(t)",
			paste(transform,"(t)",sep=""))
		scale <- exp(parms[1])
		df <- format(parms[2])
		z <- if(scale==1) paste(yvar,"-X\\beta") else paste(
			"\\frac{", yvar, "-X\\beta}{",format(scale),"}",sep="")
		z <- paste("1-\\Phi(",z,")")
		z <- paste("1-T_{",df,"}(",z,")", sep="")
		z
		}
    )
 )

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
                          xlab, ylab, xlim, ylim, datadensity=TRUE, ...) {
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
    f <- survfit.km(factor(rep(1,length(est.surv))),
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
    f <- survfit.km(factor(rep(1,sum(s))),
                    Surv(1-est.surv[s],S[s,2]),
                    se.fit = FALSE, conf.type = "none")
    curves[[i]] <- list(x=c(0,f$time), y=1-c(1,f$surv))
}  
  labcurve(curves, pl=TRUE, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
  abline(a=0,b=1,lty=2)
  invisible()
}



validate.cph <- function(fit,method="boot",
                         B=40,bw=FALSE,rule="aic",type="residual",
                         sls=.05,aics=0,pr=FALSE,
                         dxy=FALSE,u,tol=1e-9, ...) {

  atr <- fit$Design
  if(!length(atr)) atr <- getOldDesign(fit)

  need.surv <- dxy & any(atr$assume.code==8)
  if(need.surv & missing(u))
    stop("Presence of strata -> survival estimates needed for dxy; u omitted")

  modtype <- fit$method

  discrim <- function(x,y,fit,iter,evalfit=FALSE,dxy=FALSE,need.surv=FALSE,u,
                      modtype,pr=FALSE,...) {
    stra <- y[,3]
    y <- y[,1:2]
    n <- length(stra)
    ## length(unique()) 25apr03
    if(!length(x) || length(unique(x))==1) {
      Dxy <- 0
      slope <- 1
      D <- 0
      U <- 0
      R2 <- 0
    }	else {
      x <- as.matrix(x)
      dimnames(x) <- list(as.character(1:nrow(x)),as.character(1:ncol(x)))
      if(evalfit) {	#Fit was for training sample
        lr <- -2 * (fit$loglik[1]-fit$loglik[2])
        ll0 <- -2 * fit$loglik[1]
        slope <- 1
        D <- (lr - 1)/ll0
        U <- -2/ll0
        R2.max <- 1 - exp(-ll0/n)
        R2 <- (1 - exp(-lr/n))/R2.max
      }	else	{
        storage.mode(x) <- "double"
        f <- coxphFit(x,y,stra,iter.max=10,eps=.0001,method=modtype)

#        if(is.character(f) || any(is.na(f$coef)))
#          stop(paste("fit failure in discrim:",f))  1may02
        if(f$fail) stop('fit failure in discrim,coxphFit')
        ##x is x*beta from training sample
        lr <- -2 * (f$loglik[1]-f$loglik[2])
        ll0 <- -2 * f$loglik[1]
        slope <- f$coef[1]
        D <- (lr - 1)/ll0
        R2.max <- 1 - exp(-ll0/n)
        R2 <- (1 - exp(-lr/n))/R2.max
        f.frozen <- coxphFit(x,y,stra,iter.max=0,eps=.0001,
                             method=modtype, init=1)
#        if(is.character(f.frozen) || any(is.na(f.frozen$coef)))
#          stop(paste("fit failure in discrim:",
#                     f.frozen))
        if(f.frozen$fail) stop('fit failure in discrim for f.frozen')
        U <- -2 * (f.frozen$loglik[2] - f$loglik[2]) / ll0
      }	}
    
    Q <- D - U
    z <- c(R2, slope, D, U, Q)
    nam <- c("R2","Slope", "D", "U", "Q")
    if(dxy)	{
      if(need.surv) {
        attr(x, "strata") <- stra
        x <- survest(fit, linear.predictors=x, times=u, conf.int=FALSE)$surv }
      Dxy <- rcorr.cens(x, y)["Dxy"]
      z <- c(Dxy, z)
      nam <- c("Dxy", nam)
	}
    names(z) <- nam
    z
  }

  cox.fit <- function(x, y, u, need.surv=FALSE, modtype, tol=1e-9,
                      ...) {
    if(!length(x)) return(list(fail=FALSE,coefficients=numeric(0))) ##25apr03
	if(!need.surv) u <- 0
	stra <- y[,3]
	y <- y[,1:2]
    ##	coxph(x,y,e,pr=F,surv=need.surv)
	if(!need.surv) {
      storage.mode(x) <- "double"
      x <- as.matrix(x)
      dimnames(x) <- list(as.character(1:nrow(x)),as.character(1:ncol(x)))
      f <- coxphFit(x,y,stra,iter.max=10,eps=.0001,
                    method=modtype, toler.chol=tol)
#      if(is.character(f)) {  1may02
#		cat("Fit error in coxph.fit:",f,sep="\n")
#		return(list(fail=T))	}
      if(f$fail) return(f)
      if(any(is.na(f$coef))) {
		cat('Singularity in coxph.fit. Coefficients:\n'); print(f$coef)
        return(list(fail=TRUE))
      }
#      f$fail <- FALSE
      return(f)
    }

	x <- x	#Don't want lazy evaluation of complex expression
	attributes(y) <- c(attributes(y), list(class="Surv", type="right"))
	f <- cph(y ~ x + strat(stra), surv=TRUE, method=modtype)
	f$non.slopes <- f$assume.code <- f$assign <- f$name <- f$assume <- NULL
    ##Don't fool fastbw called from predab.resample
	f
  }	

  predab.resample(fit,method=method,fit=cox.fit,measure=discrim,pr=pr,
                  B=B,bw=bw,rule=rule,type=type,sls=sls,aics=aics,dxy=dxy,
                  u=u,need.surv=need.surv,strata=TRUE,modtype=modtype,tol=tol,
                  ...)
}


#Resampling optimism of discrimination and reliability of a logistic 
#regression model
#B: # reps
#bw=T to incorporate backward stepdown (using fastbw) with params rule,type,sls
#pr=T to print results of each bootstrap rep

validate.lrm <- function(fit,method="boot",
	B=40, bw=FALSE, rule="aic", type="residual",
	sls=.05, aics=0, pr=FALSE,
    kint,
	Dxy.method=if(k==1)"somers2" else "lrm",
	emax.lim=c(0,1), ...)						{

k <- fit$non.slopes
y <- fit$y
if(length(y)==0) stop("fit did not use x=T,y=T")
if(!is.factor(y)) y <- factor(y)   ## was category 11Apr02
fit$y <- oldUnclass(y)-1  #mainly for Brier score (B)

if(missing(kint)) kint <- floor((k+1)/2)

penalty.matrix <- fit$penalty.matrix

discrim <- function(x,y,fit,iter,evalfit=FALSE,pr=FALSE,Dxy.method="somers2",
                    penalty.matrix, kint, ...) {
if(evalfit)		{	#Fit was for bootstrap sample
	stats <- fit$stats
	lr <- stats["Model L.R."]
	if(Dxy.method=="lrm") Dxy <- stats["Dxy"] else
	Dxy <- somers2(x,y)["Dxy"]
	intercept <- 0
	shrink <- 1
	n <- stats["Obs"]
	D <- (lr-1)/n
	U <- -2/n
	Q <- D - U
	R2 <- stats["R2"]
			}
else			{	
	k <- fit$non.slopes
	null.model <- length(fit$coef)==k
	refit <- if(null.model) lrm.fit(y=y) else lrm.fit(x,y,tol=1e-13)
	kr <- refit$non.slopes   # 9Nov98
	#Model with no variables = null model
	stats <- refit$stats
	lr <- stats["Model L.R."]
	if(Dxy.method=="lrm")Dxy <- stats["Dxy"] else
	Dxy <- somers2(x,y)["Dxy"]
	intercept <- refit$coef[kint]  # was [floor((k+1)/2)] was [1] 1Apr02
	if(null.model) shrink <- 1   else   shrink <- refit$coef[kr+1] # 9Nov98
	n <- stats["Obs"]
	D <- (lr-1)/n
	L01 <- -2 * sum( (y>=kint)*x - logb(1 + exp(x)), na.rm=TRUE)  ## was y 1Apr02
	U <- (L01 - refit$deviance[2] -2)/n
	Q <- D - U
	R2 <- stats["R2"]
			}
	P <- plogis(x)  # 1/(1+exp(-x)) 14Sep00
	B <- sum(((y>=kint)-P)^2)/n  ## was y 1Apr02
z <- c(Dxy, R2, intercept, shrink, D, U, Q, B)
names(z) <- c("Dxy", "R2", "Intercept", "Slope", "D", "U", "Q", "B")
z									}

lrmfit <- function(x, y, maxit=12, tol=1e-7, penalty.matrix=NULL, 
                   xcol=NULL, ...) {
  if(length(xcol) && length(penalty.matrix)>0)
    penalty.matrix <- penalty.matrix[xcol,xcol,drop=FALSE]
  lrm.fit(x, y, maxit=maxit, penalty.matrix=penalty.matrix, tol=tol)
}

z <- predab.resample(fit,method=method,fit=lrmfit,measure=discrim,pr=pr,
	B=B,bw=bw,rule=rule,type=type,sls=sls,aics=aics,Dxy.method=Dxy.method,
	non.slopes.in.x=FALSE, penalty.matrix=penalty.matrix, kint=kint, ...)

calib <- z[3:4,5]
p <- seq(emax.lim[1],emax.lim[2],.0005)
L <- logb(p/(1-p))
P <- plogis(calib[1]+calib[2]*L)  # 1/(1+exp(-calib[1]-calib[2]*L)) 14Sep00
emax <- max(abs(p-P), na.rm=TRUE)  # 14Sep00
z <- rbind(z[1:4,],c(0,0,emax,emax,emax,z[1,6]),z[5:8,])
dimnames(z) <- list(c("Dxy", "R2","Intercept", "Slope", "Emax", "D", "U", "Q",
			"B"),
		    c("index.orig","training","test","optimism",
		      "index.corrected","n"))
if(k>1) z <- z[-(7:9),,drop=FALSE]
z									}

	
#Resampling optimism of discrimination and reliability of an ols regression
#B: # reps
#bw=T to incorporate backward stepdown (using fastbw) with params rule,type,sls
#pr=T to print results of each bootstrap rep
#Requires: predab.resample, fastbw, ols
#Frank Harrell 11 June 91

validate.ols <- function(fit, method="boot",
	B=40,bw=FALSE,rule="aic",type="residual",
	sls=.05,aics=0,pr=FALSE, u=NULL, rel=">", tolerance=1e-7, ...) {

  fit.orig <- fit
  
penalty.matrix <- fit.orig$penalty.matrix

discrim <- function(x,y,fit,iter,evalfit=FALSE,u=NULL,rel=NULL,pr=FALSE,...)	{

#if(evalfit) 			{
#	resid <- fit$residuals
#	df <- length(fit$coef)	}  else	{
#	resid <- y - x
#	df <- 2		}

resid <- if(evalfit) fit$residuals else y - x

n <- length(resid)
sst <- sum(y^2) - (sum(y)^2)/n
mse <- sum(resid^2)
rsquare <- 1 - mse/sst
mse <- mse/n

#if(df>=n) 
#{
#   warning("too few observations to evaluate MSE")
#   mse <- NA
#}
#else mse <- mse/(n-df)

if(evalfit)		{	#Fit being examined on sample used to fit
	intercept <- 0
	slope <- 1	}	else	{
        if(length(fit$coef)==1) {intercept <- mean(y)-mean(x); slope <- 1}
	else {coef <- lsfit(x, y)$coef   #Note x is really x*beta from other fit
	intercept <- coef[1]
	slope <- coef[2]}		}


z <- c(rsquare, mse, intercept, slope)
nam <- c("R-square", "MSE", "Intercept", "Slope")
if(length(u))						{
   if(rel==">") yy <- ifelse(y>u,1,0)
   else if(rel==">=") yy <- ifelse(y>=u,1,0)
   else if(rel=="<") yy <- ifelse(y<u,1,0)
   else yy <- ifelse(y<=u,1,0)
   z <- c(z, somers2(x,yy)["Dxy"])
   nam <- c(nam, paste("Dxy Y",rel,format(u),sep=""))
   if(rel==">"|rel==">=") P <- 1-pnorm((u-x)/sqrt(mse))
   else P <- pnorm((u-x)/sqrt(mse))
   P0 <- sum(yy)/n
   L <- -2*sum(yy*logb(P)+(1-yy)*logb(1-P))
   L0<- -2*sum(yy*logb(P0)+(1-yy)*logb(1-P0))
   R2 <- (1-exp(-(L0-L)/n))/(1-exp(-L0/n))
   z <- c(z, R2)
   nam <- c(nam, paste("R2 Y",rel,format(u),sep=""))
 }
names(z) <- nam
z
}

ols.fit <- function(x,y,tolerance=1e-7,backward, 
                    penalty.matrix=NULL, xcol=NULL, ...) {
  if(!length(x)) {
    ybar <- mean(y)
    n <- length(y)
    residuals <- y-ybar
    v <- sum(residuals^2)/(n-1)
    return(list(coef=ybar,var=v/n,residuals=residuals,fail=FALSE))
							}
  if(length(penalty.matrix)>0) {
    if(length(xcol)) {
      xcol <- xcol[-1]-1   # remove position for intercept
      penalty.matrix <- penalty.matrix[xcol,xcol,drop=FALSE]
    }
    fit <- lm.pfit(x, y, penalty.matrix=penalty.matrix,
                   tol=tolerance)
  } else {
    fit <- lm.fit.qr.bare(x,as.vector(y),tolerance=tolerance,
                          intercept=FALSE, xpxi=TRUE)
    if(backward) 
      fit$var <- sum(fit$residuals^2)*fit$xpxi/
        (length(y) - length(fit$coefficients))
  }
  c(fit,fail=FALSE)
}

predab.resample(fit.orig,method=method,fit=ols.fit,measure=discrim,pr=pr,
	B=B,bw=bw,rule=rule,type=type,sls=sls,aics=aics,tolerance=tolerance,
	backward=bw,u=u, penalty.matrix=penalty.matrix,
        rel=rel, ...)
}

	
if(.newSurvival.) {
validate.psm <- function(fit,method="boot",B=40,
		   bw=FALSE,rule="aic",type="residual",sls=.05,aics=0,pr=FALSE,
		   dxy=FALSE,tol=1e-12, rel.tolerance=1e-5, maxiter=15, ...) {

xb <- fit$linear.predictors
ny <- dim(fit$y)
nevents <- sum(fit$y[,ny[2]])

#Note: fit$y already has been transformed by the link function by psm

dist  <- fit$dist
scale <- fit$scale
parms <- fit$parms

distance <- function(x,y,fit,iter,evalfit=FALSE,fit.orig,dxy=FALSE,dist,parms,
	tol=1e-12, maxiter=15, rel.tolerance=1e-5, ...){
	#Assumes y is matrix with 1st col=time, 2nd=event indicator
if(evalfit)		{	#Fit was for training sample
  lr <- 2*diff(fit$loglik)
  ll0 <- -2*fit$loglik[1]
  R2.max <- 1 - exp(-ll0/length(x))
  R2 <- (1 - exp(-lr/length(x)))/R2.max
  intercept <- 0
  slope <- 1
  D <- (lr-1)/ll0
  U <- -2/ll0		}
else			{
  f <- survreg.fit2(cbind(1,x),y,iter=iter,dist=dist,parms=parms,tol=tol,
	maxiter=maxiter,rel.tolerance=rel.tolerance)
  if(f$fail) stop("survreg.fit2 failed in distance")
  lr <- 2*diff(f$loglik)
  ll0 <- -2*f$loglik[1]
  R2.max <- 1 - exp(-ll0/length(x))
  R2 <- (1 - exp(-lr/length(x)))/R2.max
  intercept <- f$coefficients[1]
  slope <- f$coefficients[2]
  D <- (lr-1)/ll0
  init <- c(0, 1, if(length(f$scale)) log(f$scale) else NULL)
  f.frozen <- survreg.fit2(cbind(1,x),y,dist=dist,parms=parms,tol=tol,
        maxiter=0, init=init)
  if(f.frozen$fail) stop('survreg.fit2 failed for frozen coefficient re-fit')
  ll0 <- -2*f.frozen$loglik[1]
  frozen.lr <- 2*diff(f.frozen$loglik)
  U <- (frozen.lr - lr)/ll0
}

Q <- D - U
z <- c(R2, intercept, slope, D, U, Q)
nam <- c("R2", "Intercept", "Slope", "D", "U", "Q")
if(dxy)		{
  Dxy <- rcorr.cens(x,y)["Dxy"]; z <- c(Dxy, z); nam <- c("Dxy", nam)   }
names(z) <- nam
z
}

predab.resample(fit, method=method,
		fit=survreg.fit2, measure=distance,
		pr=pr, B=B, bw=bw, rule=rule, type=type,  
		dxy=dxy,
		dist=dist, inverse=inverse, parms=parms,
		sls=sls, aics=aics, strata=FALSE, tol=tol, maxiter=maxiter, 
		rel.tolerance=rel.tolerance, ...)
}


survreg.fit2 <- function(x,y,iter=0,dist,parms=NULL,tol,maxiter=15, 
	init=NULL, rel.tolerance=1e-5, fixed=NULL, ...) {
	e <- y[,2]
	if(sum(e)<5)return(list(fail=TRUE))
	x <- x	#Get around lazy evaluation creating complex expression
    if(!.R.) survreg.distributions <- survReg.distributions
    dlist <- survreg.distributions[[dist]]
    logcorrect <- 0
    if (length(dlist$trans)) {
      exactsurv <- y[,ncol(y)] ==1
      if(any(exactsurv)) {
        ytrans <- if(length(dlist$itrans)) dlist$itrans(y[exactsurv,1]) else
                 y[exactsurv,1]
        logcorrect <- sum(logb(dlist$dtrans(ytrans)))
      }
    }
    if (length(dlist$dist)) dlist <- survreg.distributions[[dlist$dist]]

	f <- if(!.R.)
      survReg.fit(as.matrix(x),y,dist=dlist,parms=parms,
                  controlvals=survReg.control(maxiter=maxiter,
                    rel.tolerance=rel.tolerance,
                    failure=2),
                  offset=rep(0,length(e)),init=init) else
    survreg.fit(as.matrix(x),y,dist=dlist,parms=parms,
                controlvals=survreg.control(maxiter=maxiter,
                  rel.tolerance=rel.tolerance,
                  failure=2),
                offset=rep(0,length(e)),init=init)
	if(is.character(f)) { warning(f); return(list(fail=TRUE)) }
	f$fail <- FALSE
    
    # TODO: fetch scale properly if fixed
    nstrata <- length(f$icoef) - 1
    if (nstrata > 0) {
      nvar <- length(f$coef) - nstrata
      f$scale <- exp(f$coef[-(1:nvar)])
      names(f$scale) <- NULL  # get rid of log( ) in names
      f$coefficients  <- f$coefficients[1:nvar]
    }
	else f$scale <- scale
    f$loglik <- f$loglik + logcorrect
    
#	f$var <- solvet(f$imat, tol=tol)
#	sd <- survreg.distributions[[dist]]
#	f$deviance <- sum(sd$deviance(y,f$parms, f$deriv[,1]))
#	f$null.deviance <- f$deviance + 2*(f$loglik[2] - f$ndev[2])
#	f$loglik <- c(f$ndev[2], f$loglik[2])
    f
  }

NULL
} else {

  validate.psm <- function(fit,method="boot",dxy=FALSE,B=40,
		   bw=FALSE,rule="aic",type="residual",sls=.05,aics=0,pr=FALSE,
		   tol=1e-12, rel.tolerance=1e-5, maxiter=15, ...) {

## f$coef -> f$coefficients 14Aug01
#call <- match.call()
xb <- fit$linear.predictors
ny <- dim(fit$y)
nevents <- sum(fit$y[,ny[2]])

#Note: fit$y already has been transformed by the link function by psm

link <- fit$family["link"]
inverse <- glm.links["inverse",link][[1]]
family <- fit$family
fixed <- fit$fixed
fixed <- if(length(fixed)==1 && is.logical(fixed) && !fixed) list() else
   list(scale=TRUE)
dist <- fit$family["name"]


distance <- function(x,y,fit,iter,evalfit=FALSE,fit.orig,dxy=FALSE,dist,fixed,
	family, tol=1e-12, maxiter=15, rel.tolerance=1e-5, ...){
	#Assumes y is matrix with 1st col=time, 2nd=event indicator
if(evalfit)		{	#Fit was for training sample
  ll0 <- fit$null.deviance
  lr <- ll0-fit$deviance
  R2.max <- 1 - exp(-ll0/length(x))
  R2 <- (1 - exp(-lr/length(x)))/R2.max
  intercept <- 0
  slope <- 1
  D <- (lr-1)/ll0
  U <- -2/ll0		}
else			{
  f <- survreg.fit2(cbind(1,x),y,iter=iter,dist=dist,fixed=fixed,tol=tol,
	family=family,maxiter=maxiter,rel.tolerance=rel.tolerance)
  if(f$fail) stop("survreg.fit2 failed in distance")
  ll0 <- f$null.deviance
  lr <- ll0-f$deviance
  R2.max <- 1 - exp(-ll0/length(x))
  R2 <- (1 - exp(-lr/length(x)))/R2.max
  intercept <- f$coefficients[1]
  slope <- f$coefficients[2]
  D <- (lr-1)/ll0
  init <- if(length(fixed)==0 || (is.logical(fixed) && !fixed)) 
	c(0,1,f$parms) else c(0, 1)

  f.frozen <- survreg.fit2(cbind(1,x),y,dist=dist,fixed=fixed,tol=tol,
        maxiter=0, init=init,family=family)
  if(f.frozen$fail) stop('survreg.fit2 failed for frozen coefficient re-fit')
  ll0 <- f.frozen$null.deviance
  U <- (f.frozen$deviance-f$deviance)/ll0
			}

Q <- D - U
z <- c(R2, intercept, slope, D, U, Q)
nam <- c("R2", "Intercept", "Slope", "D", "U", "Q")
if(dxy)		{
  Dxy <- rcorr.cens(x,y)["Dxy"]; z <- c(Dxy, z); nam <- c("Dxy", nam)   }
names(z) <- nam
z								}

predab.resample(fit, method=method,
		fit=survreg.fit2, measure=distance,
		pr=pr, B=B, bw=bw, rule=rule, type=type,  
		dxy=dxy,
		dist=dist, inverse=inverse, fixed=fixed, family=family,
		sls=sls, aics=aics, strata=FALSE, tol=tol, maxiter=maxiter, 
		rel.tolerance=rel.tolerance, ...)
								}


survreg.fit2 <- function(x,y,iter=0,dist,fixed=NULL,family,tol,maxiter=15, 
	init=NULL, rel.tolerance=1e-5, parms=NULL, ...) {
	e <- y[,2]
	if(sum(e)<5)return(list(fail=TRUE))
	x <- x	#Get around lazy evaluation creating complex expression
	f <- survreg.fit(as.matrix(x),y,dist=dist,
		fixed=fixed,
		controlvals=survreg.control(maxiter=maxiter,
		  rel.tolerance=rel.tolerance,
		  failure=2),
		offset=rep(0,length(e)),init=init)
	if(is.character(f)) { warning(f); return(list(fail=TRUE)) }
    if(!length(f$coefficients)) {
      f$coefficients <- f$coef # 14Aug01
      f$coef <- NULL
    }
	if(length(fixed)==0 || !fixed$scale)
      f$coefficients <- f$coefficients[-length(f$coefficients)]
	f$family <- family
	f$fail <- FALSE
	f$var <- solvet(f$imat, tol=tol)
	sd <- survreg.distributions[[dist]]
	f$deviance <- sum(sd$deviance(y,f$parms, f$deriv[,1]))
	f$null.deviance <- f$deviance + 2*(f$loglik[2] - f$ndev[2])
	f$loglik <- c(f$ndev[2], f$loglik[2])
	f
  }
NULL
}

validate <-
  function(fit,  method="boot", B=40,
           bw=FALSE, rule="aic", type="residual", sls=0.05, aics=0, 
           pr=FALSE,...) UseMethod("validate")

validate.tree <- function(fit, method, B, bw, rule, type, sls, aics,
                          pr=TRUE, k, rand, xval = 10, 
                          FUN, ...) {

  wm <- if(inherits(fit,'tree'))'tree' else
    if(inherits(fit,'rpart'))'rpart' else 'none'
  if(wm=='none') stop("Not legitimate tree")
#  if(.R. && wm=='rpart') require('maptree')  # for prune.Rpart
  if(missing(FUN))
    FUN <- switch(wm,
                  tree=prune.tree,
                  rpart=function(...,k)prune.rpart(...,cp=k))
  act <- (fit$call)$na.action
  if(!length(act))
    act <- function(x) x
  m <- model.frame(fit, na.action = act)
  if(!is.data.frame(m)) stop('you must specify model=T in the fit')
  y <- model.extract(m, response)
  binary <- is.logical(y) ||
   ((length(un <- sort(unique(y[!is.na(y)]))) == 
                               2) && un[1] == 0 && un[2] == 1)
  if(binary && is.factor(y)) y <- as.numeric(y) - 1  ## 12dec02
  call <- match.call()
  method <- call$method
  size <- NULL
  if(missing(k)) {
    if(wm=='tree') {
      ff <- FUN(fit, ...)
      k <- ff$k
      size <- ff$size[is.finite(k)]
      k <- k[is.finite(k)]  # tree makes first k -Inf   7dec02
    } else {
      k    <- fit$cptable[,'CP']
      size <- fit$cptable[,'nsplit']
    }
  }
  if(missing(rand))
    rand <- sample(xval, length(m[[1]]), replace = TRUE)
  which <- unique(rand)
  pdyx.app <- pdyx.val <- pb.app <- pb.val <- double(length(k))
  l <- 0
  for(kk in k) {
    l <- l + 1
    dyx.val <- dyx.app <- b.val <- b.app <- double(length(which))
    j <- 0
    for(i in which) {
      j <- j + 1
      s <- rand != i
      tlearn <- switch(wm,tree=tree(model=m[s,]),rpart=rpart(model=m[s,]))
      papp <- if(kk == 0) tlearn else FUN(tlearn, k = kk, ...)
      if(nrow(papp$frame) == 1) {
        dyx.app[j] <- dyx.val[j] <- 0	#no splits
        b.app[j] <- b.val[j] <- mean((y - mean(y))^2, na.rm = TRUE)
      }
      else {
        yhat <- predict(papp, newdata = m[s,  ])
        if(is.matrix(yhat) && ncol(yhat) > 1)
          yhat <- yhat[,ncol(yhat),drop=TRUE]
        ## tree with factor binary y  7dec02
        b.app[j] <- mean((yhat - y[s])^2)
        dyx.app[j] <- if(binary) somers2(yhat, y[s])["Dxy"] else
          rcorr.cens(yhat, y[s])["Dxy"]
        s <- rand == i
        yhat <- predict(papp, newdata = m[s,  ])
        b.val[j] <- mean((yhat - y[s])^2)
        dyx.val[j] <- if(binary) somers2(yhat, y[s])["Dxy"] else
          rcorr.cens(yhat, y[s])["Dxy"]
      }
    }
    pdyx.app[l] <- mean(dyx.app)
    pdyx.val[l] <- mean(dyx.val)
    pb.app[l] <- mean(b.app)
    pb.val[l] <- mean(b.val)
    if(pr) {
      dyx.app <- c(dyx.app, pdyx.app[l])
      dyx.val <- c(dyx.val, pdyx.val[l])
      b.app <- c(b.app, pb.app[l])
      b.val <- c(b.val, pb.val[l])
      cat("\n\nk=", format(kk), ":\n\n")
      dyx <- cbind(dyx.app, dyx.val, b.app, b.val)
      dimnames(dyx) <- list(c(as.character(1:j), "Mean"),
                            c("Dxy Training", "Dxy Test", "MSE Training", 
                              "MSE Test"))
      print(dyx)
    }
  }
  structure(list(k = k, size = size, dxy.app = pdyx.app, dxy.val = 
                 pdyx.val, mse.app = pb.app, mse.val = pb.val,
                 binary = binary, xval = xval),
            class = "validate.tree")
}

print.validate.tree <- function(x, ...) {
  cat(x$xval, "-fold cross-validation\n\n", sep = "")
  w <- cbind(k = x$k, size = x$size, Dxy.apparent = x$dxy.app, 
             Dxy.val = x$dxy.val, MSE.apparent = x$mse.app, MSE.val = 
             x$mse.val)
  if(x$binary)
    dimnames(w) <- list(NULL, c("k", if(length(x$size)) "size", 
                                "Dxy.apparent", "Dxy.val", "Brier.apparent", 
                                "Brier.val"))
  invisible(print(w))
}

plot.validate.tree <- function(x, what = c("mse", "dxy"),
                               legendloc = locator, ...) {
  obj <- x
  if(length(obj$size)) {
    x <- obj$size
    xlab <- "Number of Nodes"
  }
  else {
    x <- obj$k
    xlab <- "Cost/Complexity Parameter"
  }
  if("mse" %in% what) {
    blab <- if(obj$binary) "Brier Score" else "Mean Squared Error"
    ylim <- range(c(obj$mse.app, obj$mse.val))
    plot(x, obj$mse.app, xlab = xlab, ylab = blab, ylim = ylim, 
         type = "n")
    lines(x, obj$mse.app, lty = 3)
    lines(x, obj$mse.val, lty = 1)
    title(sub = paste(obj$xval, "-fold cross-validation", sep = ""),
          adj = 0)
    if(is.function(legendloc))
      legend(legendloc(1), c("Apparent", "Cross-validated"),
             lty = c(3, 1), bty = "n") else {
               par(usr=c(0,1,0,1))
               legend(legendloc[1],legendloc[2],
                      c("Apparent", "Cross-validated"),
                      lty = c(3, 1), bty = "n")
             }
  }
  if("dxy" %in% what) {
    ylim <- range(c(obj$dxy.app, obj$dxy.val))
    plot(x, obj$dxy.app, xlab = xlab, ylab = "Somers' Dxy", ylim = 
         ylim, type = "n")
    lines(x, obj$dxy.app, lty = 3)
    lines(x, obj$dxy.val, lty = 1)
    title(sub = paste(obj$xval, "-fold cross-validation", sep = ""),
          adj = 0)
    if(is.function(legendloc))
      legend(legendloc(1), c("Apparent", "Cross-validated"),
             lty = c(3, 1), bty = "n") else {
               par(usr=c(0,1,0,1))
               legend(legendloc[1],legendloc[2],
                      c("Apparent", "Cross-validated"),
                      lty = c(3, 1), bty = "n")
             }
  }
  invisible()
}

validate.rpart <- function(fit, ...) validate.tree(fit, ...)
vif <- function(fit) {

v <- Varcov(fit, regcoef.only=TRUE)
nam <- dimnames(v)[[1]]
ns <- num.intercepts(fit)
v <- solve(v)
if(ns>0) {v <- v[-(1:ns),-(1:ns),drop=FALSE]; nam <- nam[-(1:ns)]}
d <- diag(v)^.5
v <- diag(solve(v/(d %o% d)))
names(v) <- nam
v

}
which.influence <- function(fit, cutoff=.2)				{

  cox <- inherits(fit,"cph") || (length(fit$fitFunction) &&
                                 any(fit$fitFunction=='cph'))
                                 ##14Nov00 22May01

  stats <- resid(fit, "dfbetas")
  stats <- stats[!is.na(stats[,1]), ]   ##delete rows added back due to NAs
  rnam <- dimnames(stats)[[1]]
  if(!length(rnam)) rnam <- 1:nrow(stats)

  at <- fit$Design
  if(!length(at)) at <- getOldDesign(fit)

  
  w <- list()
  namw <- NULL
  k <- 0
#  .Options$warn <- -1   14Sep00
  oldopt <- options(warn=-1)
  on.exit(options(oldopt))

  if(!cox)			{
	ww <- rnam[abs(stats[,1])>=cutoff]
	if(length(ww))	{
	  k <- k+1
	  w[[k]] <- ww
	  namw <- "Intercept"
	}
  }

  Assign <- fit$assign
  nm <- names(Assign)[1]
  if(nm=="Intercept" | nm=="(Intercept)") Assign[[1]] <- NULL
  ##remove and re-number

  j <- 0
  for(i in (1:length(at$name))[at$assume.code!=8])	{
	j <- j+1
	as <- Assign[[j]]
	if(length(as)==1) ww <- rnam[abs(stats[,as])>=cutoff]
	  else	{
		z <- rep(FALSE,length(rnam))
		for(r in as)
		  z <- z | abs(stats[,r])>=cutoff
		ww <- rnam[z]
	  }
	if(length(ww))	{
	  k <- k+1
	  w[[k]] <- ww
	  namw <- c(namw, at$name[i])
	}
	TRUE
  }
  if(length(w))names(w) <- namw

  w								}


##show.influence was written by:
##Jens Oehlschlaegel-Akiyoshi
##oehl@psyres-stuttgart.de
##Center for Psychotherapy Research
##Christian-Belser-Strasse 79a
##D-70597 Stuttgart Germany

show.influence <- function(object, dframe, report=NULL, sig=NULL, id=NULL) {
  who <- unlist(object)
  nam <- names(object)  # was names(w) 24Nov00
  ## In future parse out interaction components in case main effects
  ## not already selected 24Nov00
  ia <- grep('\\*',nam)              # remove interactions   28may02
  if(length(ia)) nam <- nam[-ia]
  nam <- nam[nam %nin% 'Intercept']  # remove Intercept
  rnam <- dimnames(dframe)[[1]]
  if(!length(rnam)) rnam <- 1:nrow(dframe)
  if (length(report)) col <- c(nam,
	dimnames(dframe[,report,drop=FALSE])[[2]] )
	else col <- nam
  row <- rnam %in% who
  if(any(col %nin% names(dframe)))
    stop(paste('needed variables not in dframe:',
               paste(col[col %nin% names(dframe)],collapse=' ')))
  dframe <- dframe[row,col,drop=FALSE]
  rnam <- rnam[row]
  Count <- table(who)
  Count <- as.vector(Count[match(rnam,names(Count))])
  for (i in 1:length(nam)){
    ni <- nam[i]        # 24Nov00
	val <- dframe[,ni]  #i]
	if (length(sig) && is.numeric(val)) val <- signif(val, sig) else
    val <- format(val)
	dframe[,ni] <- paste(ifelse(rnam %in% object[[ni]],"*",""), val, sep  = "")
    ## In future change i to also find any object containing the
    ## variable (e.g., interaction)   was object[[i]] dframe[,i] 24Nov00
  }
  if (length(sig) && length(report))
	for (i in (length(nam)+1):dim(dframe)[2])
			if(is.numeric(dframe[,i]))
				dframe[,i] <- signif(dframe[,i],sig)
  dframe <- data.frame(Count,dframe)
  if(length(id)) row.names(dframe) <- id[as.numeric(row.names(dframe))]
  ## 24Nov00
  print(dframe, quote=FALSE)
  invisible(dframe)
}







