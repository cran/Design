\name{glmD}
\alias{glmD}
\alias{print.glmD}
\title{Design Version of glm}
\description{
This function saves \code{Design} attributes with the fit object so that
\code{anova.Design}, \code{plot.Design}, etc. can be used just as with
\code{ols} and other fits.  No \code{validate} or \code{calibrate}
methods exist for \code{glmD} though.
}
\usage{
glmD(formula, family = gaussian, data = list(), weights = NULL, subset =
NULL, na.action = na.fail, start = NULL, offset = NULL, control =
glm.control(...), model = TRUE, method = "glm.fit", x = FALSE, y = TRUE,
contrasts = NULL, \dots)

\method{print}{glmD}(x, digits=4, \dots)
}
\arguments{
  \item{formula}{}
  \item{family}{}
  \item{data}{}
  \item{weights}{}
  \item{subset}{}
  \item{na.action}{}
  \item{start}{}
  \item{offset}{}
  \item{control}{}
  \item{model}{}
  \item{method}{}
  \item{x}{}
  \item{y}{}
  \item{contrasts}{see \code{\link{glm}}; for \code{print}, \code{x} is
	the result of \code{glmD}}
  \item{\dots}{ignored for \code{print}}
  \item{digits}{number of significant digits to print}
}
\value{a fit object like that produced by \code{\link{glm}} but with
  \code{Design} attributes and a \code{class} of \code{"Design"},
  \code{"glmD"}, and \code{"glm"} or \code{"glm.null"}.
}
\seealso{\code{\link{glm}},\code{\link{Design}}}
\examples{
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
f <- glm(counts ~ outcome + treatment, family=poisson())
f
anova(f)
summary(f)
f <- glmD(counts ~ outcome + treatment, family=poisson())
# could have had rcs( ) etc. if there were continuous predictors
f
anova(f)
summary(f, outcome=c('1','2','3'), treatment=c('1','2','3'))
}
\keyword{models}
\keyword{regression}
