\name{pphsm}
\alias{pphsm}
\alias{print.pphsm}
\title{
Parametric Proportional Hazards form of AFT Models
}
\description{
Translates an accelerated failure time (AFT) model fitted by
\code{psm} to proportional hazards form, if the fitted model was
a Weibull or exponential model (extreme value distribution with
"log" link).
}
\usage{
pphsm(fit)
\method{print}{pphsm}(x, correlation=TRUE, \dots)
}
\arguments{
  \item{fit}{fit object created by \code{psm}}
  \item{x}{result of \code{psm}}
  \item{correlation}{set to \code{FALSE} to suppress printing of
	correlation matrix of parameter estimates}
\item{\dots}{ignored}
}
\value{
a new fit object with transformed parameter estimates
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link{psm}}, \code{\link{summary.Design}}, \code{\link{print.pphsm}}
}
\examples{
set.seed(1)
S <- Surv(runif(100))
x <- runif(100)
dd <- datadist(x); options(datadist='dd')
f <- psm(S ~ x, dist="exponential")
summary(f)        # effects on log(T) scale
f.ph <- pphsm(f)
summary(f.ph)     # effects on hazard ratio scale
options(datadist=NULL)
}
\keyword{models}
\keyword{survival}
\keyword{regression}
% Converted by Sd2Rd version 1.21.
