\name{val.prob}
\alias{val.prob}
\alias{val.surv}
\alias{print.val.prob}
\alias{plot.val.prob}
\alias{plot.val.surv}
\title{
Validate Predicted Probabilities
}
\description{
The \code{val.prob} and \code{val.surv} functions are useful for validating
predicted probabilities against binary events and predicted survival
probabilities against right-censored failure times, respectively.
First \code{val.prob} is described.

Given a set of predicted probabilities \code{p} or predicted log odds
\code{logit}, and a vector of binary outcomes \code{y} that were not
used in developing the predictions \code{p} or \code{logit},
\code{val.prob} computes the following indexes and statistics: Somers'
\eqn{D_{xy}} rank correlation between \code{p} and \code{y}
[\eqn{2(C-.5)}, \eqn{C}=ROC area], Nagelkerke-Cox-Snell-Maddala-Magee
R-squared index, Discrimination index \code{D} [ (Logistic model
L.R. \eqn{\chi^2}{chi-square} - 1)/n], L.R. \eqn{\chi^2}{chi-square},
its \eqn{P}-value, Unreliability index \eqn{U}, \eqn{\chi^2}{chi-square}
with 2 d.f.  for testing unreliability (H0: intercept=0, slope=1), its
\eqn{P}-value, the quality index \eqn{Q}, \code{Brier} score (average
squared difference in \code{p} and \code{y}), \code{Intercept}, and
\code{Slope}, and \eqn{E_{max}}=maximum absolute difference in predicted
and calibrated probabilities.  If \code{pl=TRUE}, plots fitted logistic
calibration curve and optionally a smooth nonparametric fit using
\code{lowess(p,y,iter=0)} and grouped proportions vs.  mean predicted
probability in group.  If the predicted probabilities or logits are
constant, the statistics are returned and no plot is made.

When \code{group} is present, different statistics are computed,
different graphs are made, and the object returned by \code{val.prob} is
different.  \code{group} specifies a stratification variable.
Validations are done separately by levels of group and overall.  A
\code{print} method prints summary statistics and several quantiles of
predicted probabilities, and a \code{plot} method plots calibration
curves with summary statistics superimposed, along with selected
quantiles of the predicted probabilities (shown as tick marks on
calibration curves).  Only the \code{lowess} calibration curve is
estimated.  The statistics computed are the average predicted
probability, the observed proportion of events, a 1 d.f. chi-square
statistic for testing for overall mis-calibration (i.e., a test of the
observed vs. the overall average predicted probability of the event)
(\code{ChiSq2}), and a 2 d.f. chi-square statistic for testing
simultaneously that the intercept of a linear logistic calibration curve
is zero and the slope is one (\code{B ChiSq}), average absolute
calibration error (average absolute difference between the
\code{lowess}-estimated calibration curve and the line of identity,
labeled \code{Eavg}), \code{Eavg} divided by the difference between the
0.95 and 0.05 quantiles of predictive probabilities (\code{Eavg/P90}), a
"median odds ratio", i.e., the anti-log of the median absolute
difference between predicted and calibrated predicted log odds of the
event (\code{Med OR}), the C-index (ROC area), the Brier quadratic error
score (\code{B}), a chi-square test of goodness of fit based on the
Brier score, and the Brier score computed on calibrated rather than raw
predicted probabilities (\code{B cal}).  The first chi-square test is a
test of overall calibration accuracy ("calibration in the large"), and
the second will also detect errors such as slope shrinkage caused by
overfitting or regression to the mean.  See Cox (1970) for both of these
score tests.  The goodness of fit test based on the (uncalibrated) Brier
score is due to Hilden, Habbema, and Bjerregaard (1978) and is discussed
in Spiegelhalter (1986).  When \code{group} is present you can also
specify sampling \code{weights} (usually frequencies), to obtained
weighted calibration curves.

To get the behavior that results from a grouping variable being present
without having a grouping variable, use \code{group=TRUE}.  In the
\code{plot} method, calibration curves are drawn and labeled by default
where they are maximally separated using the \code{labcurve} function.
The following parameters do not apply when \code{group} is present:
\code{pl}, \code{smooth}, \code{logistic.cal}, \code{m}, \code{g},
\code{cuts}, \code{emax.lim}, \code{legendloc}, \code{riskdist},
\code{mkh}, \code{connect.group}, \code{connect.smooth}.  The following
parameters apply to the \code{plot} method but not to \code{val.prob}:
\code{xlab}, \code{ylab}, \code{lim}, \code{statloc}, \code{cex}.

\code{val.surv} uses Cox-Snell (1968) residuals on the cumulative
probability scale to check on the calibration of a survival model
against right-censored failure time data.  If the predicted survival
probability at time \eqn{t} for a subject having predictors \eqn{X} is
\eqn{S(t|X)}, this method is based on the fact that the predicted
probability of failure before time \eqn{t}, \eqn{1 - S(t|X)}, when
evaluated at the subject's actual survival time \eqn{T}, has a uniform
(0,1) distribution.  The quantity \eqn{1 - S(T|X)} is right-censored
when \eqn{T} is.  By getting one minus the Kaplan-Meier estimate of the
distribution of \eqn{1 - S(T|X)} and plotting against the 45 degree line
we can check for calibration accuracy.  A more stringent assessment can
be obtained by stratifying this analysis by an important predictor
variable.  The theoretical uniform distribution is only an approximation
when the survival probabilities are estimates and not population values.

When \code{censor} is specified to \code{val.surv}, a different
validation is done that is more stringent but that only uses the
uncensored failure times.  This method is used for type I censoring when
the theoretical censoring times are known for subjects having uncensored
failure times.  Let \eqn{T}, \eqn{C}, and \eqn{F} denote respectively
the failure time, censoring time, and cumulative failure time
distribution (\eqn{1 - S}).  The expected value of \eqn{F(T | X)} is 0.5
when \eqn{T} represents the subject's actual failure time.  The expected
value for an uncensored time is the expected value of \eqn{F(T | T \leq
C, X) = 0.5 F(C | X)}.  A smooth plot of \eqn{F(T|X) - 0.5 F(C|X)} for
uncensored \eqn{T} should be a flat line through \eqn{y=0} if the model
is well calibrated.  A smooth plot of \eqn{2F(T|X)/F(C|X)} for
uncensored \eqn{T} should be a flat line through \eqn{y=1.0}. The smooth
plot is obtained by smoothing the (linear predictor, difference or
ratio) pairs.
}
\usage{
val.prob(p, y, logit, group, weights=NULL, normwt=FALSE, 
         pl=TRUE, smooth=TRUE, logistic.cal=TRUE,
         xlab="Predicted Probability", ylab="Actual Probability",
         lim=c(0, 1), m, g, cuts, emax.lim=c(0,1), legendloc, 
         statloc=lim, riskdist="calibrated", cex=.75, mkh=.02,
         connect.group=FALSE, connect.smooth=TRUE, g.group=4, 
         evaluate=100, nmin=0)

\method{print}{val.prob}(x, \dots)

\method{plot}{val.prob}(x, xlab="Predicted Probability", 
     ylab="Actual Probability",
     lim=c(0,1), statloc=lim, stats=1:12, cex=.5, 
     lwd.overall=4, quantiles=c(.05,.95), flag, \dots)

val.surv(fit, newdata, S, est.surv, censor)

\method{plot}{val.surv}(x, group, g.group=4, what=c('difference','ratio'),
     type=c('l','b','p'),
     xlab, ylab, xlim, ylim, datadensity=TRUE, \dots)
}
\arguments{
\item{p}{
predicted probability
}
\item{y}{
vector of binary outcomes
}
\item{logit}{
predicted log odds of outcome.  Specify either \code{p} or \code{logit}.
}
\item{fit}{
a fit object created by \code{cph} or \code{psm}
}
\item{group}{
a grouping variable.  If numeric this variable is grouped into
\code{g.group} quantile groups (default is quartiles).  Set \code{group=TRUE} to
use the \code{group} algorithm but with a single stratum for \code{val.prob}.
}
\item{weights}{
an optional numeric vector of per-observation weights (usually frequencies),
used only if \code{group} is given.
}
\item{normwt}{
set to \code{TRUE} to make \code{weights} sum to the number of non-missing observations.
}
\item{pl}{
TRUE to plot calibration curves and optionally statistics
}
\item{smooth}{
plot smooth fit to \code{(p,y)} using \code{lowess(p,y,iter=0)}
}
\item{logistic.cal}{
plot linear logistic calibration fit to \code{(p,y)}
}
\item{xlab}{
x-axis label, default is \code{"Predicted Probability"} for \code{val.prob}.
Other defaults are used by \code{val.surv}.
}
\item{ylab}{
y-axis label, default is \code{"Actual Probability"} for \code{val.prob}.  Other
defaults are used by \code{val.surv}.
}
\item{lim}{
limits for both x and y axes  
}
\item{m}{
If grouped proportions are desired, average no. observations per group
}
\item{g}{
If grouped proportions are desired, number of quantile groups
}
\item{cuts}{
If grouped proportions are desired, actual cut points for constructing
intervals, e.g. \code{c(0,.1,.8,.9,1)} or \code{seq(0,1,by=.2)}
}
\item{emax.lim}{
Vector containing lowest and highest predicted probability over which to
compute \code{Emax}.
}
\item{legendloc}{
If \code{pl=TRUE}, list with components \code{x,y} or vector \code{c(x,y)} for upper left corner
of legend for curves and points.  Default is \code{c(.55, .27)} scaled to
\code{lim}.  Use \code{locator(1)} to use the mouse, \code{FALSE} to suppress legend.
}
\item{statloc}{
\eqn{D_{xy}}, \eqn{C}, \eqn{R^2}, \eqn{D}, \eqn{U}, \eqn{Q}, \code{Brier} score, \code{Intercept}, \code{Slope}, and \eqn{E_{max}}
will be added to plot, using
\code{statloc} as the upper left corner of a box (default is \code{c(0,.9)}).
You can specify a list or a vector.  Use \code{locator(1)}
for the mouse, \code{FALSE} to suppress statistics.  
This is plotted after the curve legends.
}
\item{riskdist}{
Defaults to \code{"calibrated"} to plot the relative frequency distribution of
calibrated robabilities after dividing into 101 bins from \code{lim[1]} to
\code{lim[2]}.
Set to \code{"predicted"} to use raw assigned risk, \code{FALSE} to omit risk distribution.
Values are scaled so that highest bar is \code{0.15*(lim[2]-lim[1])}.
}
\item{cex}{
Character size for legend or for table of statistics when \code{group} is given
}
\item{mkh}{
Size of symbols for legend.   Default is 0.02 (see \code{par()}).
}
\item{connect.group}{
Defaults to \code{FALSE} to only represent group fractions as triangles.
Set to \code{TRUE} to also connect with a solid line.
}
\item{connect.smooth}{
Defaults to \code{TRUE} to draw smoothed estimates using a dashed line.
Set to \code{FALSE} to instead use dots at individual estimates.
}
\item{g.group}{
number of quantile groups to use when \code{group} is given and variable is
numeric.
}
\item{evaluate}{
number of points at which to store the \code{lowess}-calibration curve.
Default is 100.  If there are more than \code{evaluate} unique predicted
probabilities, \code{evaluate} equally-spaced quantiles of the unique
predicted probabilities, with linearly interpolated calibrated values,
are retained for plotting (and stored in the object returned by
\code{val.prob}.
}
\item{nmin}{
applies when \code{group} is given.  When \code{nmin} \eqn{> 0}, \code{val.prob} will not
store coordinates of smoothed calibration curves in the outer tails,
where there are fewer than \code{nmin} raw observations represented in
those tails.  If for example \code{nmin}=50, the \code{plot} function will only
plot the estimated calibration curve from \eqn{a} to \eqn{b}, where there are
50 subjects with predicted probabilities \eqn{< a} and \eqn{> b}.
\code{nmin} is ignored when computing accuracy statistics.
}
\item{x}{result of \code{val.prob} (with \code{group} in effect) or
  \code{val.surv}}
\item{\dots}{
optional arguments for \code{labcurve} (through \code{plot}).  Commonly used
options are \code{col} (vector of colors for the strata plus overall) and
\code{lty}.  Ignored for \code{print}.
}
\item{stats}{
vector of column numbers of statistical indexes to write on plot
}
\item{lwd.overall}{
line width for plotting the overall calibration curve
}
\item{quantiles}{
a vector listing which quantiles should be indicated on each
calibration curve using tick marks.  The values in \code{quantiles} can be
any number of values from the following: .01, .025, .05, .1, .25, .5, .75, .9, .95, .975, .99.
By default the 0.05 and 0.95 quantiles are indicated.
}
\item{flag}{
a function of the matrix of statistics (rows representing groups)
returning a vector of character strings (one value for each group, including
"Overall").  \code{plot.val.prob}
will print this vector of character values to the left of the
statistics.  The \code{flag} function 
can refer to columns of the matrix used as input to the function by
their names given in the description above.  The default function
returns \code{"*"} if either \code{ChiSq2} or \code{B ChiSq} is significant at the
0.01 level and \code{" "} otherwise.
}
\item{newdata}{
a data frame for which \code{val.surv} should obtain predicted survival
probabilities.  If omitted, survival estimates are made for all of the
subjects used in \code{fit}.
}
\item{S}{an \code{\link[survival]{Surv}} object}
\item{est.surv}{
a vector of estimated survival probabilities corresponding to times in
the first column of \code{S}.
}
\item{censor}{
a vector of censoring times.  Only the censoring times for uncensored
observations are used.
}
\item{what}{
the quantity to plot when \code{censor} was in effect.  The default is to
show the difference between cumulative probabilities and their
expectation given the censoring time.  Set \code{what="ratio"} to show the
ratio instead.
}
\item{type}{
Set to the default (\code{"l"}) to plot the trend line only, \code{"b"} to plot
both individual subjects ratios and trend lines, or \code{"p"} to plot only points.
}
\item{xlim}{
}
\item{ylim}{
axis limits for \code{plot.val.surv} when the \code{censor} variable was used.
}
\item{datadensity}{
By default, \code{plot.val.surv} will show the data density on each curve
that is created as a result of \code{censor} being present.  Set
\code{datadensity=FALSE} to suppress these tick marks drawn by \code{scat1d}.
}}
\value{
\code{val.prob} without \code{group} returns a vector with the following named
elements: \code{Dxy}, \code{R2}, \code{D}, \code{D:Chi-sq}, \code{D:p},
\code{U}, \code{U:Chi-sq}, \code{U:p}, \code{Q}, \code{Brier}, \code{Intercept}, \code{Slope}, \code{Emax}.
When \code{group} is present \code{val.prob} returns an object of class
\code{val.prob} containing a list with summary statistics and calibration
curves for all the strata plus \code{"Overall"}.  
}
\details{
The 2 d.f. \eqn{\chi^2}{chi-square} test and \code{Med OR} exclude predicted or
calibrated predicted probabilities \eqn{\leq 0} to zero or \eqn{\geq 1},
adjusting the sample size as needed.
}
\author{
Frank Harrell\cr
Division of Biostatistics and Epidemiology, University of Virginia\cr
fharrell@virginia.edu
}
\section{Side Effects}{
plots calibration curve
}
\references{
Harrell FE, Lee KL, Mark DB (1996): Multivariable prognostic models:
Issues in developing models, evaluating assumptions and adequacy, and
measuring and reducing errors.  Stat in Med 15:361--387.


Harrell FE, Lee KL (1987):  Using logistic calibration to assess the
accuracy of probability predictions (Technical Report).


Miller ME, Hui SL, Tierney WM (1991): Validation techniques for
logistic regression models.  Stat in Med 10:1213--1226.


Harrell FE, Lee KL (1985):  A comparison of the
\emph{discrimination}
of discriminant analysis and logistic regression under multivariate
normality.  In Biostatistics: Statistics in Biomedical, Public Health,
and Environmental Sciences.  The Bernard G. Greenberg Volume, ed. PK
Sen. New York: North-Holland, p. 333--343.


Cox DR (1970): The Analysis of Binary Data, 1st edition, section 4.4.
London: Methuen.


Spiegelhalter DJ (1986):Probabilistic prediction in patient management.
Stat in Med 5:421--433.


Cox DR, Snell EJ (1968):A general definition of residuals (with
discussion).  JRSSB 30:248--275.
}
\seealso{
\code{\link{validate.lrm}}, \code{\link{lrm.fit}}, \code{\link{lrm}}, \code{\link[Hmisc]{labcurve}}, \code{\link[Hmisc]{wtd.rank}},
\code{\link[Hmisc]{wtd.loess.noiter}}, \code{\link[Hmisc]{scat1d}}, \code{\link{cph}}, \code{\link{psm}}
}
\examples{
# Fit logistic model on 100 observations simulated from the actual 
# model given by Prob(Y=1 given X1, X2, X3) = 1/(1+exp[-(-1 + 2X1)]),
# where X1 is a random uniform [0,1] variable.  Hence X2 and X3 are 
# irrelevant.  After fitting a linear additive model in X1, X2,
# and X3, the coefficients are used to predict Prob(Y=1) on a
# separate sample of 100 observations.


set.seed(1)
n <- 200
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
logit <- 2*(x1-.5)
P <- 1/(1+exp(-logit))
y <- ifelse(runif(n)<=P, 1, 0)
d <- data.frame(x1,x2,x3,y)
f <- lrm(y ~ x1 + x2 + x3, subset=1:100)
pred.logit <- predict(f, d[101:200,])
phat <- 1/(1+exp(-pred.logit))
val.prob(phat, y[101:200], m=20, cex=.5)  # subgroups of 20 obs.


# Validate predictions more stringently by stratifying on whether
# x1 is above or below the median


v <- val.prob(phat, y[101:200], group=x1[101:200], g.group=2)
v
plot(v)
plot(v, flag=function(stats) ifelse(
  stats[,'ChiSq2'] > qchisq(.95,2) |
  stats[,'B ChiSq'] > qchisq(.95,1), '*', ' ') )
# Stars rows of statistics in plot corresponding to significant
# mis-calibration at the 0.05 level instead of the default, 0.01


plot(val.prob(phat, y[101:200], group=x1[101:200], g.group=2), 
              col=1:3) # 3 colors (1 for overall)


# Weighted calibration curves
# plot(val.prob(pred, y, group=age, weights=freqs))


# Survival analysis examples
# Generate failure times from an exponential distribution
set.seed(123)              # so can reproduce results
n <- 2000
age <- 50 + 12*rnorm(n)
sex <- factor(sample(c('Male','Female'), n, rep=TRUE, prob=c(.6, .4)))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
t <- -log(runif(n))/h
label(t) <- 'Time to Event'
ev <- ifelse(t <= cens, 1, 0)
t <- pmin(t, cens)
S <- Surv(t, ev)


# First validate true model used to generate data
w <- val.surv(est.surv=exp(-h*t), S=S)
plot(w)
plot(w, group=sex)  # stratify by sex


# Now fit an exponential model and validate
# Note this is not really a validation as we're using the
# training data here
f <- psm(S ~ age + sex, dist='exponential', y=TRUE)
w <- val.surv(f)
plot(w, group=sex)


# We know the censoring time on every subject, so we can
# compare the predicted Pr[T <= observed T | T>c, X] to
# its expectation 0.5 Pr[T <= C | X] where C = censoring time
# We plot a ratio that should equal one
w <- val.surv(f, censor=cens)
plot(w)


plot(w, group=age, g=3)   # stratify by tertile of age
}
\keyword{validate}
\keyword{models}
\keyword{regression}
\keyword{htest}
\keyword{accuracy}
\keyword{lrm}
\keyword{smooth}
\keyword{sampling}
\keyword{survival}

