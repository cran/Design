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

