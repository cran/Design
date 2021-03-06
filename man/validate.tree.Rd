\name{validate.tree}
\alias{validate.tree}
\alias{validate.rpart}
\alias{print.validate.tree}
\alias{plot.validate.tree}
\title{
Dxy and Mean Squared Error by Cross-validating a Tree Sequence
}
\description{
Uses \code{xval}-fold cross-validation of a sequence of trees to derive
estimates of the mean squared error and Somers' \code{Dxy} rank correlation
between predicted and observed responses.  In the case of a binary response
variable, the mean squared error is the Brier accuracy score.
This function is a modification of \code{cv.tree} which should be
consulted for details.  There are \code{print} and \code{plot} methods for
objects created by \code{validate.tree}.
}
\usage{
# f <- tree(formula=y ~ x1 + x2 + \dots) # or rpart
\method{validate}{tree}(fit, method, B, bw, rule, type, sls, aics, pr=TRUE,
    k, rand, xval=10, FUN, \dots)
\method{validate}{rpart}(fit, \dots)
\method{print}{validate.tree}(x, \dots)
\method{plot}{validate.tree}(x, what=c("mse","dxy"), legendloc=locator, \dots)
}
\arguments{
\item{fit}{
an object created by \code{tree} or \code{rpart} or having the same
attributes as one created by \code{tree}.  If it was created by
\code{rpart} you must have specified the \code{model=TRUE} argument to
\code{rpart}. 
}
\item{method,B,bw,rule,type,sls,aics}{are there only for consistency
  with the generic \code{validate} function; these are ignored}
\item{x}{the result of \code{validate.tree}}
\item{k}{
a sequence of cost/complexity values.  By default these are obtained
from calling \code{FUN} with no optional arguments (if \code{tree}) or
from the \code{rpart} \code{cptable} object in the original fit object.
You may also specify a scalar or vector.
}
\item{rand}{
see \code{cv.tree}
}
\item{xval}{
number of splits
}
\item{FUN}{
the name of a function which produces a sequence of trees, such
as \code{prune.tree} or \code{shrink.tree} or \code{prune.rpart}.  Default is
\code{prune.tree} for fits from \code{tree} and \code{prune.rpart} for fits from \code{rpart}.
}
\item{\dots}{
additional arguments to \code{FUN} (ignored by \code{print,plot}).  For
\code{validate.rpart}, \dots can be the same arguments used in
\code{validate.tree}.
}
\item{pr}{
set to \code{FALSE} to prevent intermediate results for each \code{k} to be printed
}
\item{what}{
a vector of things to plot.  By default, 2 plots will be done, one for
\code{mse} and one for \code{Dxy}.
}
\item{legendloc}{
a function that is evaluated with a single argument equal to \code{1} to
generate a list with components \code{x, y} specifying coordinates of the
upper left corner of a legend, or a 2-vector.  For the latter,
\code{legendloc} specifies the relative fraction of the plot at which to
center the legend.
}}
\value{
a list of class \code{"validate.tree"} with components named \code{k, size, dxy.app},
\code{dxy.val, mse.app, mse.val, binary, xval}.  \code{size} is the number of nodes,
\code{dxy} refers to Somers' \code{D}, \code{mse} refers to mean squared error of prediction,
\code{app} means apparent accuracy on training samples, \code{val} means validated
accuracy on test samples, \code{binary} is a logical variable indicating whether
or not the response variable was binary (a logical or 0/1 variable is
binary).  \code{size} will not be present if the user specifies \code{k}.
}
\section{Side Effects}{
prints if \code{pr=TRUE}
}
\author{
Frank Harrell
\cr
Department of Biostatistics
\cr
Vanderbilt University
\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link[rpart]{rpart}}, \code{\link[Hmisc]{somers2}},
\code{\link[Hmisc]{rcorr.cens}}, \code{\link[tree]{cv.tree}},
\code{\link{locator}}, \code{\link{legend}}
}
\examples{
\dontrun{
n <- 100
set.seed(1)
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
y  <- 1*(x1+x2+rnorm(n) > 1)
table(y)
library(rpart)
f <- rpart(y ~ x1 + x2 + x3, model=TRUE)
v <- validate(f)
v    # note the poor validation
par(mfrow=c(1,2))
plot(v, legendloc=c(.2,.5))
par(mfrow=c(1,1))
}
}
\keyword{models}
\keyword{tree}
\keyword{category}
\concept{model validation}
\concept{predictive accuracy}
