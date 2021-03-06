\name{validate}
\alias{validate}
\title{
  Resampling Validation of a Fitted Model's Indexes of Fit
}
\description{
  The \code{validate} function when used on an object created by one of the
  \code{Design} series does resampling validation of a 
  regression model, with or without backward step-down variable deletion.
}
\usage{
# fit <- fitting.function(formula=response ~ terms, x=TRUE, y=TRUE)
validate(fit, method="boot", B=40,
         bw=FALSE, rule="aic", type="residual", sls=0.05, aics=0, 
         pr=FALSE, \dots)
}
\arguments{
  \item{fit}{
    a fit derived by e.g. \code{lrm}, \code{cph}, \code{psm}, \code{ols}. The options \code{x=TRUE} and \code{y=TRUE}
    must have been specified.
  }
  \item{method}{
    may be \code{"crossvalidation"}, \code{"boot"} (the default), \code{".632"}, or
    \code{"randomization"}.
    See \code{predab.resample} for details.  Can abbreviate, e.g.
    \code{"cross", "b", ".6"}.
  }
  \item{B}{
    number of repetitions.  For \code{method="crossvalidation"}, is the
    number of groups of omitted observations.
  }
  \item{bw}{
    \code{TRUE} to do fast step-down using the \code{fastbw} function,
    for both the overall model and for each repetition. \code{fastbw}
    keeps parameters together that represent the same factor.
  }
  \item{rule}{
    Applies if \code{bw=TRUE}.  \code{"aic"} to use Akaike's information criterion as a
    stopping rule (i.e., a factor is deleted if the \eqn{\chi^2}{chi-square} falls below
    twice its degrees of freedom), or \code{"p"} to use \eqn{P}-values.
  }
  \item{type}{
    \code{"residual"} or \code{"individual"} - stopping rule is for individual factors or
    for the residual \eqn{\chi^2}{chi-square} for all variables deleted
  }
  \item{sls}{
    significance level for a factor to be kept in a model, or for judging the
    residual \eqn{\chi^2}{chi-square}.
  }
  \item{aics}{
    cutoff on AIC when \code{rule="aic"}.
  }
  \item{pr}{
    \code{TRUE} to print results of each repetition
  }
  \item{\dots}{
    parameters for each specific validate function, and parameters to
    pass to \code{predab.resample} (note especially the \code{group},
    \code{cluster}, amd \code{subset} parameters).
    
    For \code{psm}, you can pass the \code{maxiter} parameter here (passed to 
    \code{survreg.control}, default is 15 iterations) as well as a \code{tol} parameter 
    for judging matrix singularity in \code{solvet} (default is 1e-12)
    and a \code{rel.tolerance} parameter that is passed to
    \code{survreg.control} (default is 1e-5).
  }
}
\details{
  It provides bias-corrected indexes that are specific to each type
  of model. For \code{validate.cph} and \code{validate.psm}, see \code{validate.lrm},
  which is similar. \cr
  For \code{validate.cph} and \code{validate.psm}, there is
  an extra argument \code{dxy}, which if \code{TRUE} causes the \code{rcorr.cens}
  function to be invoked to compute the Somers' \eqn{D_{xy}}{Dxy} rank correlation
  to be computed at each resample (this takes a bit longer than
  the likelihood based statistics). The values corresponting to the row
  \eqn{D_{xy}}{Dxy} are equal to \eqn{2 * (C - 0.5)} where C is the
  C-index or concordance probability. \cr
  
  For \code{validate.cph} with \code{dxy=TRUE},
  you must specify an argument \code{u} if the model is stratified, since
  survival curves can then cross and \eqn{X\beta}{X beta} is not 1-1 with
  predicted survival. \cr
  There is also \code{validate} method for
  \code{tree}, which only does cross-validation and which has a different
  list of arguments.   
}
\value{
a matrix with rows corresponding to the statistical indexes and
columns for columns for the original index, resample estimates, 
indexes applied to
the whole or omitted sample using the model derived from the resample,
average optimism, corrected index, and number of successful re-samples.
}
\section{Side Effects}{
prints a summary, and optionally statistics for each re-fit
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link{validate.ols}}, \code{\link{validate.cph}}, \code{\link{validate.lrm}}, \code{\link{validate.tree}},
\code{\link{predab.resample}}, \code{\link{fastbw}}, \code{\link{Design}}, \code{\link{Design.trans}}, \code{\link{calibrate}}
}
\examples{
# See examples for validate.cph, validate.lrm, validate.ols
# Example of validating a parametric survival model:


n <- 1000
set.seed(731)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
sex <- factor(sample(c('Male','Female'), n, TRUE))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
dt <- -log(runif(n))/h
e <- ifelse(dt <= cens,1,0)
dt <- pmin(dt, cens)
units(dt) <- "Year"
S <- Surv(dt,e)


f <- psm(S ~ age*sex, x=TRUE, y=TRUE)  # Weibull model
# Validate full model fit
validate(f, B=10)                # usually B=150


# Validate stepwise model with typical (not so good) stopping rule
# bw=TRUE does not preserve hierarchy of terms at present
validate(f, B=10, bw=TRUE, rule="p", sls=.1, type="individual")
}
\keyword{models}
\keyword{regression}
\keyword{methods}
\keyword{survival}
\concept{model validation}
\concept{predictive accuracy}
\concept{bootstrap}
