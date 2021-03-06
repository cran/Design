\name{calibrate}
\alias{calibrate}
\alias{calibrate.default}
\alias{calibrate.cph}
\alias{calibrate.psm}
\alias{print.calibrate}
\alias{print.calibrate.default}
\alias{plot.calibrate}
\alias{plot.calibrate.default}
\title{
Resampling Model Calibration
}
\description{
Uses bootstrapping or cross-validation to get bias-corrected (overfitting-
corrected) estimates of predicted vs. observed values based on
subsetting predictions into intervals (for survival models) or on
nonparametric smoothers (for other models). There are calibration
functions for Cox (\code{cph}), parametric survival models (\code{psm}),
binary and ordinal logistic models (\code{lrm}) and ordinary least squares (\code{ols}).
For survival models,
"predicted" means predicted survival probability at a single
time point, and "observed" refers to the corresponding Kaplan-Meier 
survival estimate, stratifying on intervals of predicted survival.
For logistic and linear models, a nonparametric calibration curve is
estimated over a sequence of predicted values.
The fit must have specified \code{x=TRUE, y=TRUE}.  The \code{print} and \code{plot} methods
for \code{lrm} and \code{ols} models (which use \code{calibrate.default}) print the mean
absolute error in predictions, the mean squared error, and the 0.9 quantile
of the absolute error.  Here, error refers to the difference between the
predicted values and the corresponding bias-corrected calibrated values.

Below, the second, third, and fourth invocations of \code{calibrate}
are, respectively, for \code{ols} and \code{lrm}, \code{cph}, and
\code{psm}.  The first and second \code{plot} invocation are
respectively for \code{lrm} and \code{ols} fits or all other fits.
}
\usage{
calibrate(fit, \dots)
\method{calibrate}{default}(fit, predy, 
  method=c("boot","crossvalidation",".632","randomization"),
  B=40, bw=FALSE, rule=c("aic","p"),
  type=c("residual","individual"),
  sls=.05, pr=FALSE, kint, smoother="lowess", \dots)
\method{calibrate}{cph}(fit, method="boot", u, m=150, cuts, B=40, 
  bw=FALSE, rule="aic", type="residual", sls=0.05, aics=0, 
  pr=FALSE, what="observed-predicted", tol=1e-12, \dots)
\method{calibrate}{psm}(fit, method="boot", u, m=150, cuts, B=40,
  bw=FALSE,rule="aic",
  type="residual",sls=.05,aics=0,
  pr=FALSE,what="observed-predicted",tol=1e-12, maxiter=15, 
  rel.tolerance=1e-5, \dots)

\method{print}{calibrate}(x, \dots)
\method{print}{calibrate.default}(x, \dots)

\method{plot}{calibrate}(x, xlab, ylab, subtitles=TRUE, conf.int=TRUE,
 cex.subtitles=.75, \dots)

\method{plot}{calibrate.default}(x, xlab, ylab, xlim, ylim,
  legend=TRUE, subtitles=TRUE, scat1d.opts=NULL, \dots)
}
\arguments{
\item{fit}{
a fit from \code{ols}, \code{lrm}, \code{cph} or \code{psm}
}
\item{x}{an object created by \code{calibrate}}
\item{method}{}
\item{B}{}
\item{bw}{}
\item{rule}{}
\item{type}{}
\item{sls}{}
\item{aics}{see \code{\link{validate}}}
\item{u}{
the time point for which to validate predictions for survival models. For \code{cph} fits,
you must have specified \code{surv=TRUE, time.inc=u}, where \code{u} is
the constant specifying the time to predict.
}
\item{m}{
group predicted \code{u}-time units survival into intervals containing
\code{m} subjects on the average (for survival models only)
}
\item{cuts}{
actual cut points for predicted survival probabilities. You may
specify only one of \code{m} and \code{cuts} (for survival models only)
}
\item{pr}{
set to \code{TRUE} to print intermediate results for each re-sample
}
\item{what}{
The default is \code{"observed-predicted"}, meaning to estimate optimism
in this difference. This is preferred as it accounts for skewed
distributions of predicted probabilities in outer intervals. You can
also specify \code{"observed"}.  This argument applies to survival models only.
}
\item{tol}{criterion for matrix singularity (default is \code{1e-12})}
\item{maxiter}{for \code{psm}, this is passed to
  \code{\link[survival]{survreg.control}} (default is 15 iterations)
}
\item{rel.tolerance}{parameter passed to
  \code{\link[survival]{survreg.control}} for \code{psm} (default is 1e-5).
  }
\item{predy}{
a scalar or vector of predicted values to calibrate (for \code{lrm},
\code{ols}).  Default is 50 equally spaced points between the 5th
smallest and the 5th largest  predicted values.  For \code{lrm} the
predicted values are probabilities (see \code{kint}).
}
\item{kint}{
For an ordinal logistic model the default predicted
probability that \eqn{Y\geq} the middle level.  Specify \code{kint} to specify the
intercept to use, e.g., \code{kint=2} means to calibrate \eqn{Prob(Y\geq
  b)}, where \eqn{b} is the second level of \eqn{Y}.
}
\item{smoother}{
a function in two variables which produces \eqn{x}- and
\eqn{y}-coordinates by smoothing the input \code{y}.  The default is to
use \code{lowess(x, y, iter=0)}. 
}
\item{\dots}{
other arguments to pass to \code{predab.resample}, such as \code{group},
\code{cluster}, and \code{subset}.
Also, other arguments for \code{plot}.
}
\item{xlab}{
defaults to "Predicted x-units Survival" or to a suitable label for
other models
}
\item{ylab}{
defaults to "Fraction Surviving x-units" or to a suitable label for
other models
}
\item{xlim}{}
\item{ylim}{2-vectors specifying x- and y-axis limits, if not using defaults}
\item{subtitles}{
set to \code{FALSE} to suppress subtitles in plot describing method and for \code{lrm}
and \code{ols} the mean absolute error and original sample size
}
\item{conf.int}{
set to \code{FALSE} to suppress plotting 0.95 confidence intervals for
Kaplan-Meier estimates
}
\item{cex.subtitles}{character size for plotting subtitles}
\item{legend}{
set to \code{FALSE} to suppress legends (for \code{lrm}, \code{ols}
only) on the calibration plot, or specify a list with elements \code{x}
and \code{y} containing the coordinates of the upper left corner of the
legend.  By default, a legend will be drawn in the lower right 1/16th of
the plot.
}
\item{scat1d.opts}{a list containing additional arguments to
  \code{\link[Hmisc]{scat1d}}}
}
\value{
matrix specifying mean predicted survival in each interval, the
corresponding estimated bias-corrected Kaplan-Meier estimates,
number of subjects, and other statistics.  For linear and logistic models,
the matrix instead has rows corresponding to the prediction points, and
the vector of predicted values being validated is returned as an attribute.
The returned object has class \code{"calibrate"} or \code{"calibrate.default"}.
}
\section{Side Effects}{
prints, and stores an object \code{pred.obs} or \code{.orig.cal}
}
\details{
If the fit was created using penalized maximum likelihood estimation,
the same \code{penalty} and \code{penalty.scale} parameters are used during
validation.
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link{validate}}, \code{\link{predab.resample}},
  \code{\link{groupkm}}, \code{\link[Hmisc]{errbar}},
  \code{\link[Hmisc]{scat1d}}, \code{\link{cph}}, \code{\link{psm}},
\code{\link{lowess}}
}
\examples{
set.seed(1)
d.time <- rexp(200)
x1 <- runif(200)
x2 <- factor(sample(c('a','b','c'),200,TRUE))
f <- cph(Surv(d.time) ~ pol(x1,2)*x2, x=TRUE, y=TRUE, surv=TRUE, time.inc=2)
#or f <- psm(S ~ \dots)
cal <- calibrate(f, u=2, m=50, B=20)  # usually B=200 or 300
plot(cal)


y <- sample(0:2, 200, TRUE)
x1 <- runif(200)
x2 <- runif(200)
x3 <- runif(200)
x4 <- runif(200)
f <- lrm(y ~ x1+x2+x3*x4, x=TRUE, y=TRUE)
cal <- calibrate(f, kint=2, predy=seq(.2,.8,length=60), 
                 group=y)
# group= does k-sample validation: make resamples have same 
# numbers of subjects in each level of y as original sample


plot(cal)
#See the example for the validate function for a method of validating
#continuation ratio ordinal logistic models.  You can do the same
#thing for calibrate
}
\keyword{methods}
\keyword{models}
\keyword{regression}
\keyword{survival}
\keyword{hplot}
\concept{bootstrap}
\concept{model validation}
\concept{calibration}
\concept{model reliability}
\concept{predictive accuracy}
