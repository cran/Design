\name{lrm.fit}
\alias{lrm.fit}
\title{Logistic Model Fitter}
\description{
Fits a binary or ordinal logistic model for a given design matrix and response
vector with no missing values in either.  Ordinary or penalized maximum
likelihood estimation is used.
}
\usage{
lrm.fit(x, y, offset, initial, est, maxit=12, eps=.025,
        tol=1E-7, trace=FALSE, penalty.matrix, weights, normwt)
}
\arguments{
\item{x}{
design matrix with no column for an intercept
}
\item{y}{
response vector, numeric, categorical, or character
}
\item{offset}{optional numeric vector containing an offset on the logit scale}
\item{initial}{
vector of initial parameter estimates, beginning with the
intercept
}
\item{est}{
indexes of \code{x} to fit in the model (default is all columns of \code{x}).
Specifying \code{est=c(1,2,5)} causes columns 1,2, and 5 to have
parameters estimated. The score vector \code{u} and covariance matrix \code{var}
can be used to obtain score statistics for other columns
}
\item{maxit}{
maximum no. iterations (default=\code{12}). Specifying \code{maxit=1}
causes logist to compute statistics at initial estimates.
}
\item{eps}{
difference in \preformatted{-2  log} likelihood for declaring convergence.
Default is \code{.025}.
}
\item{tol}{
Singularity criterion. Default is 1E-7
}
\item{trace}{
set to \code{TRUE} to print -2 log likelihood, step-halving
fraction, and rank of variance matrix at each iteration
}
\item{penalty.matrix}{
  a self-contained ready-to-use penalty matrix - see \code{lrm}
}
\item{weights}{
  a vector (same length as \code{y}) of possibly fractional case weights
}
\item{normwt}{
 set to \code{TRUE} to scale \code{weights} so they sum to the length of
 \code{y}; useful for sample surveys as opposed to the default of
 frequency weighting 
 }
}

\value{
a list with the following components:

\item{call}{
calling expression
}
\item{freq}{
table of frequencies for \code{y} in order of increasing \code{y}
}
\item{stats}{
vector with the following elements: number of observations used in the
fit, maximum absolute value of first
derivative of log likelihood, model likelihood ratio chi-square, d.f.,
P-value,
\eqn{c} index (area under ROC curve), Somers' \eqn{D_{xy}},
Goodman-Kruskal \eqn{\gamma}{gamma}, and Kendall's \eqn{\tau_a}{tau-a}
rank correlations 
between predicted probabilities and observed response, the
Nagelkerke \eqn{R^2} index, and the Brier probability score with
respect to computing the probability that \eqn{y >} lowest level. 
Probabilities are rounded to the nearest 0.002
in the computations or rank correlation indexes.
When \code{penalty.matrix} is present, the \eqn{\chi^2}{chi-square},
d.f., and P-value are not corrected for the effective d.f.
}
\item{fail}{
set to \code{TRUE} if convergence failed (and \code{maxiter>1})
}
\item{coefficients}{
estimated parameters
}
\item{var}{
estimated variance-covariance matrix (inverse of information matrix).
Note that in the case of penalized estimation, \code{var} is not the
improved sandwich-type estimator (which \code{lrm} does compute).
}
\item{u}{
vector of first derivatives of log-likelihood
}
\item{deviance}{
-2 log likelihoods. 
When an offset variable is present, three
deviances are computed: for intercept(s) only, for
intercepts+offset, and for intercepts+offset+predictors.
When there is no offset variable, the vector contains deviances for
the intercept(s)-only model and the model with intercept(s) and predictors.
}
\item{est}{
vector of column numbers of \code{X} fitted (intercepts are not counted)
}
\item{non.slopes}{
number of intercepts in model
}
\item{penalty.matrix}{
see above
}}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link{lrm}}, \code{\link{glm}}, \code{\link{matinv}}, \code{\link[Hmisc]{solvet}}, \code{\link{cr.setup}}
}
\examples{
#Fit an additive logistic model containing numeric predictors age, 
#blood.pressure, and sex, assumed to be already properly coded and 
#transformed
#
# fit <- lrm.fit(cbind(age,blood.pressure,sex), death)
}
\keyword{models}
\keyword{regression}
\concept{logistic regression model}
