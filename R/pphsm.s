pphsm <- function(fit)
{

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

print.pphsm <- function(x, digits = max(options()$digits - 4, 3),
                        correlation = TRUE, ...)
{
  if (length(f <- x$fail) && f) {
    stop(" Survreg failed.  No summary provided")
    return(invisible(x))
  }

  cat("Parametric Survival Model Converted to PH Form\n\n")

  stats <- x$stats
  stats[3] <- round(stats[3],2)
  stats[5] <- round(stats[5],4)
  stats[6] <- round(stats[6],2)
  if(.R.) print(format.sep(stats),quote=FALSE) else print(stats)
  cat("\n")

  print(c(x$coef, x$icoef[2]), digits=digits)

  correl <- x$correl
  if (correlation && !is.null(x$correl)) {  ## FEH
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

  invisible()
}

Varcov.pphsm <- function(object, ...)
{
  .NotYetImplemented()
}
