\name{validate.cph}
\alias{validate.cph}
\alias{validate.psm}
\title{
Validation of a Fitted Cox or Parametric Survival Model's Indexes of Fit
}
\description{
This is the version of the \code{validate} function specific to models
fitted with \code{cph} or \code{psm}. Statistics validated include the
Nagelkerke \eqn{R^2}, 
\eqn{D_{xy}}, slope shrinkage,  the
discrimination index\eqn{D} [(model L.R. \eqn{\chi^2}{chi-square} - 1)/L], the unreliability index
\eqn{U} = (difference in -2 log likelihood between uncalibrated
\eqn{X\beta}{X beta} and  
\eqn{X\beta}{X beta} with overall slope calibrated to test sample) / L,
and the overall quality index \eqn{Q = D - U}. 
L is -2 log likelihood with beta=0.  The "corrected" slope
can be thought of as shrinkage factor that takes into account overfitting.
See \code{predab.resample} for the list of resampling methods.
}
\usage{
# fit <- cph(formula=Surv(ftime,event) ~ terms, x=TRUE, y=TRUE, \dots)
\method{validate}{cph}(fit,method="boot",
                         B=40,bw=FALSE,rule="aic",type="residual",
                         sls=.05,aics=0,pr=FALSE,dxy=FALSE,u,tol=1e-9, \dots)

\method{validate}{psm}(fit, method="boot",B=40,
                   bw=FALSE,rule="aic",type="residual",sls=.05,aics=0,pr=FALSE,
                   dxy=FALSE,tol=1e-12, rel.tolerance=1e-5, maxiter=15, \dots)
}
\arguments{
\item{fit}{
a fit derived \code{cph}. The options \code{x=TRUE} and \code{y=TRUE}
must have been specified. If the model contains any stratification factors
and dxy=TRUE,
the options \code{surv=TRUE} and \code{time.inc=u} must also have been given,
where \code{u} is the same value of \code{u} given to \code{validate}.
}
\item{method}{see \code{\link{validate}}}
\item{B}{}
\item{rel.tolerance}{}
\item{maxiter}{}
\item{bw}{}
\item{rule}{}
\item{type}{}
\item{sls}{}
\item{aics}{}
\item{pr}{}
\item{tol}{}
\item{...}{see \code{\link{validate}} or \code{\link{predab.resample}}}
\item{dxy}{
set to \code{TRUE} to validate Somers' \eqn{D_{xy}}  using
\code{rcorr.cens}, which takes longer.
}
\item{u}{
must be specified if the model has any stratification factors and \code{dxy=TRUE}.
In that case, strata are not included in \eqn{X\beta}{X beta} and the
survival curves may cross.  Predictions at time \code{t=u} are
correlated with observed survival times.  Does not apply to
\code{validate.psm}.
}
}
\value{
matrix with rows corresponding to \eqn{D_{xy}}, Slope, \eqn{D},
\eqn{U}, and \eqn{Q}, and columns for the original index, resample estimates, 
indexes applied to whole or omitted sample using model derived from resample, average optimism, corrected index, and number
of successful resamples.
}
\section{Side Effects}{
prints a summary, and optionally statistics for each re-fit (if \code{pr=TRUE})
}
\author{
Frank Harrell\cr
Division of Biostatistics and Epidemiology, University of Virginia\cr
fharrell@virginia.edu
}
\seealso{
\code{\link{validate}}, \code{\link{predab.resample}}, \code{\link{fastbw}}, \code{\link{Design}}, \code{\link{Design.trans}}, \code{\link{calibrate}},
\code{\link[Hmisc]{rcorr.cens}}, \code{\link{cph}}, \code{\link[survival]{coxph.fit}}
}
\examples{
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


f <- cph(S ~ age*sex, x=TRUE, y=TRUE)
# Validate full model fit
validate(f, B=10)               # normally B=150


# Validate a model with stratification.  Dxy is the only
# discrimination measure for such models, by Dxy requires
# one to choose a single time at which to predict S(t|X)
f <- cph(S ~ rcs(age)*strat(sex), 
         x=TRUE, y=TRUE, surv=TRUE, time.inc=2)
validate(f, dxy=TRUE, u=2, B=10)   # normally B=150
# Note u=time.inc
}
\keyword{validate}
\keyword{accuracy}
\keyword{bootstrap}
\keyword{models}
\keyword{regression}
\keyword{survival}
