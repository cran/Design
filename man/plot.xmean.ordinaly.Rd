\name{plot.xmean.ordinaly}
\alias{plot.xmean.ordinaly}
\title{
Plot Mean X vs. Ordinal Y
}
\description{
Separately for each predictor variable \eqn{X} in a formula, plots the mean of
\eqn{X} vs. levels of \eqn{Y}.  Then under the proportional odds assumption,
the expected value of the predictor for each \eqn{Y} value is also plotted (as
a dotted line).  This plot is useful for assessing the ordinality assumption 
for \eqn{Y} separately for each \eqn{X}, and for assessing the proportional odds
assumption in a simple univariable way.  If several predictors do not
distinguish adjacent categories of \eqn{Y}, those levels may need to be 
pooled.  This display assumes
that each predictor is linearly related to the log odds of each event in
the proportional odds model.  There is also an option to plot the
expected means assuming a forward continuation ratio model.
}
\usage{
plot.xmean.ordinaly(x, data, subset, na.action, subn=TRUE, cr=FALSE, \dots)
}
\arguments{
\item{x}{
an S formula.  Response variable is treated as ordinal.  Predictor
variables must be binary or continuous.  Interactions or non-linear
effects are not allowed.
}
\item{data}{
a data frame or frame number
}
\item{subset}{
vector of subscripts or logical vector describing subset of data to
analyze
}
\item{na.action}{
defaults to \code{na.keep} so all NAs are initially retained.  Then NAs
are deleted only for each predictor currently being plotted.
Specify \code{na.action=na.delete} to remove observations that are missing
on any of the predictors (or the response).
}
\item{subn}{
set to \code{FALSE} to suppress a left bottom subtitle specifying the sample size
used in constructing each plot
}
\item{cr}{
set to \code{TRUE} to plot expected values by levels of the response, assuming a
forward continuation ratio model holds.  The function is fairly slow
when this option is specified.
}
\item{...}{
other arguments passed to \code{plot} and \code{lines}
}}
\section{Side Effects}{
plots
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\references{
Harrell FE et al. (1998): Development of a clinical prediction model for
an ordinal outcome. Stat in Med 17:909--44.
}
\seealso{
\code{\link{lrm}}, \code{\link{residuals.lrm}}, \code{\link{cr.setup}}, \code{\link[Hmisc]{cumcategory}}
}
\examples{
# Simulate data from a population proportional odds model
set.seed(1)
n <- 400
age <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
L <- .2*(age-50) + .1*(blood.pressure-120)
p12 <- plogis(L)    # Pr(Y>=1)
p2  <- plogis(L-1)  # Pr(Y=2)
p   <- cbind(1-p12, p12-p2, p2)   # individual class probabilites
# Cumulative probabilities:
cp  <- matrix(cumsum(t(p)) - rep(0:(n-1), rep(3,n)), byrow=TRUE, ncol=3)
y   <- (cp < runif(n)) \%*\% rep(1,3)
# Thanks to Dave Krantz <dhk@paradox.psych.columbia.edu> for this trick


par(mfrow=c(1,2))
plot.xmean.ordinaly(y ~ age + blood.pressure, cr=TRUE)
par(mfrow=c(1,1))
}
\keyword{category}
\keyword{models}
\keyword{regression}
\keyword{hplot}
\concept{model validation}
\concept{logistic regression model}