\name{survfit}
\alias{survfit}
\title{Modified Version of survival Package survfit Function}
\description{
This modification of \code{survfit} keeps attributes of the \code{Surv}
object, users \code{Hmisc}'s \code{interaction()} to form strata labels,
and uses a default confidence interval basis of log-log.
}
\usage{
survfit(formula, data, weights, subset, na.action = na.delete, conf.type = c("log-log", "log", "plain", "none"), ...)
}
\arguments{
  \item{formula}{}
  \item{data}{}
  \item{weights}{}
  \item{subset}{}
  \item{na.action}{}
  \item{conf.type}{}
  \item{\dots}{see \code{\link[survival]{survfit}}}
}
\value{see \code{\link[survival]{survfit}}}
\seealso{\code{\link[survival]{survfit}}}
\examples{
\dontrun{
#fit a Kaplan-Meier and print the results
data(aml)
survfit(Surv(time, status) ~ x, data=aml)
}
}
\keyword{survival}

