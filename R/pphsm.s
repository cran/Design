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
