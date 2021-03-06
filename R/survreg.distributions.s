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

