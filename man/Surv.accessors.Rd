\name{Surv.accessors}
\alias{Surv.event}
\alias{Surv.strata}
\alias{Surv.time}

\title{ Accessors for survivals Surv object. }
\description{
  These function allow code to access the components of a survival Surv
  object in a implementation independant way.
}
\usage{
Surv.event(surv)
Surv.strata(surv)
Surv.time(surv)
}

\arguments{
  \item{surv}{ a survival object }
}
\value{
  returns the part of the object requested
}
\author{ Charles Dupont }
\seealso{ \code{\link{Surv}} }
\keyword{ manip }
