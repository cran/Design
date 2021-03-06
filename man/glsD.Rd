\name{glsD}
\alias{glsD}
\alias{print.glsD}
\title{Fit Linear Model Using Generalized Least Squares}
\usage{
glsD(model, data, correlation, weights, subset, method, na.action,
    control, verbose, B=0, dupCluster=FALSE, pr=FALSE,
    opmeth=c('optimize','optim'))

\method{print}{glsD}(x, digits=4, \dots)
}
\arguments{
 \item{model}{a two-sided linear formula object describing the
    model, with the response on the left of a \code{~} operator and the
    terms, separated by \code{+} operators, on the right.}
 \item{data}{an optional data frame containing the variables named in
   \code{model}, \code{correlation}, \code{weights}, and
   \code{subset}. By default the variables are taken from the
   environment from which \code{gls} is called.}
 \item{correlation}{an optional \code{corStruct} object describing the
   within-group correlation structure. See the documentation of
   \code{corClasses} for a description of the available \code{corStruct}
   classes. If a grouping variable is to be used, it must be specified in
   the \code{form} argument to the \code{corStruct}
   constructor. Defaults to \code{NULL}, corresponding to uncorrelated 
   errors.}  
 \item{weights}{an optional \code{varFunc} object or one-sided formula
   describing the within-group heteroscedasticity structure. If given as
   a formula, it is used as the argument to \code{varFixed},
   corresponding to fixed variance weights. See the documentation on
   \code{varClasses} for a description of the available \code{varFunc}
   classes. Defaults to \code{NULL}, corresponding to homoscesdatic
   errors.} 
 \item{subset}{an optional expression indicating which subset of the rows of
   \code{data} should  be  used in the fit. This can be a logical
   vector, or a numeric vector indicating which observation numbers are
   to be included, or a  character  vector of the row names to be
   included.  All observations are included by default.}
 \item{method}{a character string.  If \code{"REML"} the model is fit by
   maximizing the restricted log-likelihood.  If \code{"ML"} the
   log-likelihood is maximized.  Defaults to \code{"REML"}.}
 \item{na.action}{a function that indicates what should happen when the
   data contain \code{NA}s.  The default action (\code{na.fail}) causes
   \code{gls} to print an error message and terminate if there are any
   incomplete observations.}
 \item{control}{a list of control values for the estimation algorithm to
   replace the default values returned by the function \code{glsControl}.
   Defaults to an empty list.}
 \item{verbose}{an optional logical value. If \code{TRUE} information on
   the evolution of the iterative algorithm is printed. Default is
   \code{FALSE}.}
 \item{B}{number of bootstrap resamples to fit and store, default is
   none}
 \item{dupCluster}{set to \code{TRUE} to have \code{glsD} when
   bootstrapping to consider multiply-sampled clusters as if they were
   one large cluster when fitting using the \code{gls} algorithm}
 \item{pr}{set to \code{TRUE} to show progress of bootstrap resampling}
 \item{opmeth}{specifies whether the \code{optimize} or the \code{optim}
   function is to be used for optimization}
 \item{x}{the result of \code{glsD}}
 \item{digits}{number of significant digits to print}
 \item{\dots}{ignored}
}
\description{
  This function fits a linear model using generalized least
  squares. The errors are allowed to be correlated and/or have unequal
  variances.  \code{glsD} is a slightly enhanced version of the
  Pinheiro and Bates \code{glsD} function in the \code{nlme} package to
  make it easy to use with the Design library and to implement cluster
  bootstrapping (primarily for nonparametric estimates of the
  variance-covariance matrix of the parameter estimates and for
  nonparametric confidence limits of correlation parameters).
}
\value{
  an object of classes \code{glsD}, \code{Design}, and \code{gls}
  representing the linear model
  fit. Generic functions such as \code{print}, \code{plot}, and 
  \code{summary} have methods to show the results of the fit. See
  \code{glsObject} for the components of the fit. The functions
  \code{resid}, \code{coef}, and \code{fitted} can be used to extract
  some of its components.  \code{glsD} returns the following components
  not returned by \code{gls}: \code{Design}, \code{assign},
  \code{formula}, \code{opmeth} (see arguments), \code{B} (see
  arguments), \code{bootCoef} (matrix of \code{B} bootstrapped
  coefficients), \code{boot.Corr} (vector of bootstrapped correlation
  parameters), \code{Nboot} (vector of total sample size used in each
  bootstrap (may vary if have unbalanced clusters), and \code{var}
  (sample variance-covariance matrix of bootstrapped coefficients).
}
\references{
  Pinheiro J, Bates D (2000): Mixed effects models in S and S-Plus.  New
  York: Springer-Verlag.
}
\author{Jose Pinheiro \email{jcp@research.bell-labs.com},
  Douglas Bates \email{bates@stat.wisc.edu},
  Frank Harrell \email{f.harrell@vanderbilt.edu},
  Patrick Aboyoun \email{aboyoun@insightful.com}
} 
\seealso{
  \code{\link{gls}}
  \code{\link{glsControl}}, \code{\link{glsObject}},
  \code{\link{varFunc}}, \code{\link{corClasses}}, \code{\link{varClasses}}
}
\examples{
\dontrun{
ns  <- 20  # no. subjects
nt  <- 10  # no. time points/subject
B   <- 10  # no. bootstrap resamples
           # usually do 100 for variances, 1000 for nonparametric CLs
rho <- .5  # AR(1) correlation parameter
V <- matrix(0, nrow=nt, ncol=nt)
V <- rho^abs(row(V)-col(V))   # per-subject correlation/covariance matrix

d <- expand.grid(tim=1:nt, id=1:ns)
d$trt <- factor(ifelse(d$id <= ns/2, 'a', 'b'))
true.beta <- c(Intercept=0,tim=.1,'tim^2'=0,'trt=b'=1)
d$ey  <- true.beta['Intercept'] + true.beta['tim']*d$tim +
  true.beta['tim^2']*(d$tim^2) +  true.beta['trt=b']*(d$trt=='b')
set.seed(13)
library(MASS)   # needed for mvrnorm
d$y <- d$ey + as.vector(t(mvrnorm(n=ns, mu=rep(0,nt), Sigma=V)))

dd <- datadist(d); options(datadist='dd')
# library(nlme)  # S-Plus: library(nlme3) or later
f <- glsD(y ~ pol(tim,2) + trt, correlation=corCAR1(form= ~tim | id),
          data=d, B=B)
f
f$var      # bootstrap variances
f$varBeta  # original variances
summary(f)
anova(f)
plot(f, tim=NA, trt=NA)
# v <- Variogram(f, form=~tim|id, data=d)
}
}
\keyword{models}
