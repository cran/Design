\name{survfit.cph}
\alias{survfit.cph}
\title{
Cox Predicted Survival
}
\description{
This is a slightly modified version of Therneau's \code{survfit.coxph}
function. The difference is that \code{survfit.cph} assumes that \code{x=TRUE,y=TRUE}
were specified to the fit. This assures that the environment in effect
at the time of the fit (e.g., automatic knot estimation for spline functions)
is the same one used for basing predictions. Unlike \code{survfit.coxph},
the default basis for confidence intervals is \code{"log-log"}.
}
\usage{
survfit.cph(object, newdata, se.fit=TRUE, conf.int=0.95, 
            individual=FALSE, type, vartype,
            conf.type=c("log-log", "log", "plain", "none"))
}
\arguments{
\item{object}{
a fit object from \code{cph} or \code{coxph}
see \code{\link[survival]{survfit.coxph}}
}
\item{newdata}{}
\item{se.fit}{}
\item{conf.int}{}
\item{individual}{}
\item{type}{}
\item{vartype}{}
\item{conf.type}{see \code{\link[survival]{survfit}}}
}
\value{
see \code{survfit.coxph}
}
\seealso{
\code{\link{survest.cph}}
}
\keyword{survival}
