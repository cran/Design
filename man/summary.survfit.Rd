\name{summary.survfit}
\alias{summary.survfit}
\title{Design version of survival Package summary.survfit}}
\description{
This is a modified version of the survival packages's
\code{summary.survfit} function.  It returns a list containing the
survival curve, confidence limits, and other information.
}
\usage{
\method{summary}{survfit}(object, times, censored = FALSE, scale = 1, \dots)
}
\arguments{
  \item{object}{result of \code{survfit}}
  \item{times}{}
  \item{censored}{}
  \item{scale}{}
  \item{\dots}{see \code{\link[survival]{summary.survfit}}}
}
\value{a list with components \code{time,surv,n.risk,n.event,std.err,conf.int,lower,upper,strata,call,na.action}}
}
\seealso{\code{\link[survival]{summary.survfit}}}
\examples{
}
\keyword{survival}
