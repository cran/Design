\name{rm.impute}
\alias{rm.impute}
\alias{pbind}
\title{
Imputation of Repeated Measures
}
\description{
NOTE: This function is under development and is not correct at present.
Uses the method of Lavori, Dawson, and Shera (1995) to analyze
uniformly (over subjects) collected repeated measurement data subject
to non-random dropout.  Separately for each imputation and for each
time period, a binary logistic model is developed (using the Design
\code{lrm} function) to predict the probability that each subject remains
in the study at that time period.  Predictors for the first time
period are those listed in the \code{pformula} formula.  These are assumed
to be baseline variables that are never missing.  For later time
periods, predictors include the baseline predictors plus the matrix of
response (\code{y}) values for all earlier periods.  These "previous
responses" will have missing values imputed from the earlier steps.


Missing responses for time period i are imputed, for one of the
\code{n.impute} multiple imputations, as follows.  The period i fitted
propensity model described above is evaluated to obtain the predicted
probability that each subject remained in the study until at least
period i.  The estimated propensity is divided into \code{g} quantile
groups.  If for period i within a propensity quantile group there are
\code{a} subjects still in the study and \code{b} subjects who have dropped out,
Rubin's approximate Bayesian bootstrap is used to estimate the
predictive distribution of the response values for the \code{b} dropouts,
given that the propensity for remaining in the study is approximately
constant for all subjects (dropouts and non-dropouts) in the group.  A
sample of size \code{a} is selected with replacement from the \code{a} subjects
still in the study from the propensity group.  Then a sample of size
\code{b} with replacement is selected from this sample of size \code{a}.  These
\code{b} responses are used to fill-in the responses for the \code{b} dropouts
in the quantile group for the current imputation and current time
period.


If the right-hand-side of a formula is specified for a univariate
response summary (which may be the last response, mean, or area under
the time-response curve), \code{rm.impute} goes on to fit \code{rformula} to
this response summary for each of the multiple imputations using a
fitting function \code{fitter}.  After all \code{n.impute} imputations have been
done, the average "apparent" covariance matrix and the
between-imputation covariance matrix are computed to derive Rubin's
multiple-imputation-corrected covariance matrix for the average of
\code{n.impute} sets of regression coefficients.  See \code{fit.mult.impute} for
more details.


The response variable \code{y} may be an array to handle multiple responses
at each time period.  This array has number of rows equal to the
number of subjects, number of columns equal to the number of periods,
and number of "pages" equal to the number of different response
measurements.  A utility function \code{pbind} is supplied for creating
such arrays from a series of matrices.  When multiple responses are
present, all responses are used in the current propensity model, and
the \code{which}, \code{nk}, \code{rinteraction}, and \code{rint.with} arguments will
apply equally to all responses.
}
\usage{
rm.impute(pformula, y, last, 
          rformula, fitter=ols, which=c("last", "mean", "auc"), 
          data=sys.parent(1), n.impute=10, g=5, 
          nk=0, rinteraction, rint.with=c('all','recent'),
          pr=FALSE, pra=FALSE, npr, 
          keep.prop=FALSE, keep.pfits=FALSE)


pbind(\dots)
}
\arguments{
\item{pformula}{
right-hand-side propensity formula, e.g., ~ treatment + x1 + x2 +
x3*x4.  This formula (as well as \code{rformula} if \code{fitter} is one of the Design
library fitting functions) can contain any of the Design library's
transformation functions such as \code{rcs}, \code{pol}, etc.
}
\item{y}{
matrix of responses over time periods.  To use \code{which="auc"}, column
names of \code{y} must contain numeric measurement times.
}
\item{last}{
an integer vector whose value for the jth subject is the last period
before the subject dropped out.  A subject who never had a follow-up
response measured will have \code{last=0}.
}
\item{...}{
a series of matrices to "page bind" into an array.  New names may be
supplied using the notation \code{pbind(newname1=y1,newname2=y2)}.  The
\code{dimnames} of the first argument (which will be converted to a matrix
if it is a vector, for the unusual one-period case) will be used as
the first two \code{dimnames} of the resulting array, and the names of the
matrices will form the third vector of \code{dimnames}.
}
\item{rformula}{
right-hand-side response formula, e.g., ~ x1 + pol(x2) + treatment.
If omitted, \code{rm.impute} will return only the multiple response imputations.
}
\item{fitter}{
any S-Plus or Design library fitting function for a univariate
response summary.  The default is \code{ols}.  If there are multiple
response variables at each time period and you want to use a different
fitter for different response variables, specify a list of \code{nr}
fitting functions as this argument, where \code{nr} is the number of
response variables.
}
\item{which}{
which response summary is used if \code{rformula} is given.  The default is
the last column of the response matrix.
}
\item{data}{
usually a data frame, if the variables in \code{pformula} and \code{rformula}
are not already available via \code{attach()}
}
\item{n.impute}{
number of imputations.  The more missing data, the higher \code{n.impute}
should be.
}
\item{g}{
number of propensity quantile groups
}
\item{nk}{
number of knots to use in expanding each previous response into a
restricted cubic spline in the propensity model.  Default is zero
(assume linearity).
}
\item{rinteraction}{
a character vector specifying the names of baseline variables that
should be interacted with each response in the propensity model.
Default is no such interactions.
}
\item{rint.with}{
set to \code{"recent"} to allow the variables in \code{rinteraction} to only
interact with the response for the most recent time period, and not with
the most recent and all previous responses (the default)
}
\item{pr}{
set to \code{TRUE} to print each logistic propensity model fit.
}
\item{pra}{
if \code{pr=TRUE}, you can also set \code{pra=TRUE} to print the Design \code{anova()}
results for each propensity model fit.
}
\item{npr}{
if \code{pr=TRUE}, printing will be done for the first \code{npr} imputations
}
\item{keep.prop}{
set to \code{TRUE} to store the array \code{propensity} in the returned list.  The
dimensions for \code{propensity} are the same as \code{Y}.
}
\item{keep.pfits}{
set to \code{TRUE} to store all propensity model fits from \code{lrm} in the result
returned by \code{rm.impute}
}}
\value{
a list with elements \code{Y} and optionally \code{fit} (if \code{rformula} is given)
and \code{propensity} (if \code{keep.prop=TRUE}).  \code{Y} and \code{propensity} are arrays
whose last dimension
corresponds to the multiple imputations and whose first two dimensions
correspond to \code{y}.  \code{Y} is the multiply-imputed response
array and \code{fit} is the imputation-corrected fit object.  Note: Aside
from the regression coefficient vector and covariance matrix, this fit
object will have parameters from the fit of the response summary for
the last imputation.  If \code{keep.pfits=TRUE}, the returned list will also
have an array of propensity fit objects (\code{lrm} objects) for all
response periods and imputations.  If there is more than one response
variable at each time period, \code{fit} will be a list of \code{nr} fit objects
for \code{nr} response variables.
}
\section{Side Effects}{
prints, and creates variables such as y.1, y.2, \dots and in.period.i in
the session database (frame 0)
}
\details{
The algorithm used here will not correct for non-random dropout due to
variables that are not included in the propensity model.  A worst-case
would be having dropouts at period i due to unmeasured responses at
period i.


Ironically, there must be a sufficient number of dropouts for the
propensity score method to work, as the propensity models must have
adequate numbers of dropouts and non-dropouts at each time period.
}
\author{
Frank Harrell
\cr
Department of Biostatistics
\cr
Vanderbilt University School of Medicine
\cr
f.harrell@vanderbilt.edu
\cr
Much valuable input was received from Chris Barker (Roche
Pharmaceuticals) and Phil Lavori (Stanford University).
}
\references{
Lavori PW, Dawson R, Shera, D: A multiple imputation strategy for
clinical trials with truncation of patient data.  Stat in Med
14:1913--1925, 1995. 


Rubin D, Shenker N: Multiple imputation in health-care data bases: An
overview and some applications.  Stat in Med 10:585--598, 1991.


Engels JM, Diehr P: Imputation of missing longitudinal data: a
comparison of methods.  J Clin Epi 56:968--976, 2003.
}
\seealso{
\code{\link[Hmisc]{transcan}}, \code{\link[Hmisc]{fit.mult.impute}}, \code{\link{lrm}}, \code{\link[Hmisc]{rm.boot}}, \code{\link[Hmisc]{reShape}}
}
\examples{
\dontrun{
# Generate multiple imputes of the response matrix for later use
Y <- rm.impute(~treatment + pol(age,2)*sex, responses, 
               last=lastvisit, data=mydata)$Y
# Do some analysis for each imputation
fits <- vector('list',10)
for(i in 1:10) {
  y <- Y[,,i]
  fits[[i]] <- my.analysis(X,y)
}


# Function to generate a 4-variate equal correlation pattern response
# with missing-at-random responses; missingness is a function of x and
# previous responses.
#
# pna is a function that computes the probability that a subject
# drops out at the current visit.  For visit 1 pna is a function
# of treatment and baseline covariable x.  For visits > 1 pna is
# a function of the matrix of responses for all previous visits.
#
# If second=TRUE we also generate a second response variable having
# NAs in the same positions as this first one.  y2 is generated
# so that its NAs are completely unrelated to any y2 values if
# y2B.effect=0, as the pna function is only given the first 
# response variable.
# y2 is N(0,1) for treat='A' and N(y2.treat.effect,1) for treat='B'.


testdf <- function(n=1500, seed=7, pna, second=FALSE, y2.treat.effect=0) {


set.seed(seed)
treat <- sample(c('A','B'),n,TRUE)
x <- runif(n)
nt <- 4


mvrnorm <- function(n, p = 1, u = rep(0, p), S = diag(p)) {
  Z <- matrix(rnorm(n * p), p, n)
  t(u + t(chol(S)) \%*\% Z)
}


# Generate multivariate normal errors for n subjects at nt times
# Assume equal correlations of rho=.5, independent subjects


rho   <- .5
y <- mvrnorm(n, p=nt, S=diag(rep(1-rho,nt))+rho)


y[treat=='B',] <- y[treat=='B',] + 1


cat('\n\nTreatment-specific means for last period in response variable 1 before generating NAs:\n')
print(tapply(y[,4], treat, mean, na.rm=TRUE))


y[runif(n) < pna(treat, x), 1] <- NA
y[is.na(y[,1]) | runif(n) < pna(treat, x, y[,1]),   2] <- NA
y[is.na(y[,2]) | runif(n) < pna(treat, x, y[,1:2]), 3] <- NA
y[is.na(y[,3]) | runif(n) < pna(treat, x, y[,1:3]), 4] <- NA


last <- rep(4, n)
last[is.na(y[,4])] <- 3
last[is.na(y[,3])] <- 2
last[is.na(y[,2])] <- 1
last[is.na(y[,1])] <- 0


cat('\nNumber of NAs for each time period:\n')
print(apply(y, 2, function(x)sum(is.na(x))))
cat('\n\nTreatment-specific means for last period in response variable 1 after excluding NAs:\n')
print(tapply(y[,4], treat, mean, na.rm=TRUE))
cat('\n\nNaive complete-case analysis:\n\n')
prn(ols(y[,4] ~ pol(x,2) + treat))


if(second) {
  y2 <- matrix(rnorm(n*4),ncol=4)
  y2[treat=='B',] <- y2[treat=='B',] + y2.treat.effect
  cat('\n\nTreatment-specific means for last period in response variable 2 before generating NAs:\n')
  print(tapply(y2[,4], treat, mean, na.rm=TRUE))


  y2[is.na(y[,1]),1] <- NA
  y2[is.na(y[,2]),2] <- NA
  y2[is.na(y[,3]),3] <- NA
  y2[is.na(y[,4]),4] <- NA
  cat('\n\nTreatment-specific means for last period in response variable 2 after excluding NAs:\n')
  print(tapply(y2[,4], treat, mean, na.rm=TRUE))


  y <- pbind(y1=y, y2=y2)
}


list(x=x, treat=treat, y=y, last=last)
}


pna <- function(treat, x, yprev) {
# In this model for the probability of dropout just before the
# current visit, the probability does not depend on the baseline
# covariable x.  For treat='B' the probability of dropout is a
# constant 0.1.  For treat='A' it is a curtailed quadratic
# function of the previous visit's response.
#
# If no previous responses available, we are at first follow-up visit


if(missing(yprev)) 0 else {
  if(is.matrix(yprev)) yprev <- yprev[,ncol(yprev)]
  ifelse(treat=='B', .1,
         pmax(0, pmin(1, .124 +.0835*yprev + .020868*yprev^2)))
  }
}


df <- testdf(pna = pna, second=TRUE)


g <- rm.impute(~ pol(x,2) + treat, df$y, last=df$last, 
               rformula=~ pol(x,2) + treat,
               n.impute=10, g=4, nk=3, 
               rinteraction='treat', rint.with='all',
               pr=TRUE, pra=TRUE, data=df, keep.prop=TRUE, keep.pfits=TRUE)
# Base propensity model is in.study ~ pol(x,2) + treat
# for visits 2,3,4, filled-in y's from previous visits will also be
# used as predictors, and these interact with treat.  
# Restricted cubic spline with 3 knots is assumed for the propensity models
# To fit the multiply-imputed last (4th) response an additive model
# in quadratic x and treat is used


g$fit[[1]]       # shows response fit for first response variable
                 # (y1), with variances adj. for imputation
page(g$Y)        # show all 10 imputations for both responses x 4 periods


# Check for the first imputation how well propensity matching achieved 
# balance in baseline and period 3 filled-in responses for
# dropouts and non-dropouts.  For continuous variables show ECDFs
# using the Hmisc ecdf function, for first 4 imputations.  Do this
# with and without stratifying on quintiles of propensity, and also
# show the estimated 3rd period response  vs. propensity stratified 
# by dropout status.  Use only first response (y1) for all of this.


for(imp in 1:4) {
  y3     <- g$Y[,3,1,imp]
  prop3  <- g$propensity[,3,imp]
  prop3g <- cut2(prop3,g=5)
  ti <- paste('Imputation',imp)
  print(ecdf(~ y3, groups=df$last >= 3, subset=unclass(prop3g)<5))
  title(ti)
  print(ecdf(~ y3 | prop3g, groups=df$last >= 3, 
             subset=unclass(prop3g)<5))
  # Not enough dropouts in highest quintile of propensity completing
  # visit 3
  title(ti)
  plsmo(prop3, y3, group=df$last >= 3, datadensity=TRUE, col=1:2)
  title(ti)
}


# Examine propensity fit for sixth imputation, 4th response
f <- g$pfits[4,6][[1]]
dfr <- as.data.frame(df)
# Edit names of dfr so that responses called y.1, y.2, etc.
# For this example, these are already OK
dd <- datadist(dfr)
options(datadist='dd')   
# datadist makes plot below work without specifying variable settings
plot(f, y.3=NA, treat=NA, conf.int=FALSE)


# Analyze multiple response variables.  Both systolic.bp and
# diastolic.bp are matrices (columns = time periods)


f <- rm.impute(~treatment + pol(age,2)*sex,
               pbind(systolic.bp, diastolic.bp),
               last=lastvisit, data=mydata)


# To deal with a continuous and a binary endpoint you can specify
# pbind(sysbolic.bp, stroke), fitter=list(ols, lrm)
}
}
\keyword{regression}
\keyword{htest}
\keyword{multivariate}
\keyword{array}
\concept{bootstrap}
\concept{repeated measures}
\concept{longitudinal data}
