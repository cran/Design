\name{predict.lrm}
\alias{predict.lrm}
\title{
Predicted Values for Binary and Ordinal Logistic Models
}
\description{
Computes a variety of types of predicted values for fits from
\code{lrm}, either from the original dataset or for new observations.
}
\usage{
\method{predict}{lrm}(object, \dots, type=c("lp", "fitted", "fitted.ind", "mean", "x", 
            "data.frame", "terms", "adjto","adjto.data.frame", 
            "model.frame"), se.fit=FALSE, codes=FALSE)
}
\arguments{
\item{object}{
a object created by \code{lrm}
}
\item{...}{
arguments passed to \code{predict.Design}, such as \code{kint} and \code{newdata}
(which is used if you are predicting \code{out of data}).  See
\code{predict.Design} to see how NAs are handled.
}
\item{type}{
See \code{predict.Design} for \code{"x", "data.frame", "terms", "adjto", 
"adjto.data.frame"} and \code{"model.frame"}. \code{type="lp"} is used to get
linear predictors (always using the first intercept). \code{type="fitted"}
is used to get all the probabilities \eqn{Y\geq
  j}. \code{type="fitted.ind"} gets all the individual probabilities
\eqn{Y=j}. For an ordinal response variable, \code{type="mean"} computes
the estimated mean \eqn{Y} by summing values of \eqn{Y} 
multiplied by the estimated \eqn{Prob(Y=j)}. If \eqn{Y} was a character or
\code{factor} object, the levels are the character values or factor levels,
so these must be translatable to numeric, unless \code{codes=TRUE}.
See the Hannah and Quigley reference below for the method of estimating
(and presenting) the mean score.  If you specify
\code{type="fitted","fitted.ind","mean"} you may not specify \code{kint}.
}
\item{se.fit}{
applies only to \code{type="lp"}, to get standard errors.
}
\item{codes}{
if \code{TRUE}, \code{type="mean"} uses the integer codes
\eqn{1,2,\ldots,k} for the \eqn{k}-level response in computing the
predicted mean response.
}
}
\value{
a vector (\code{type="lp"} with \code{se.fit=FALSE}, or \code{type="mean"} or only one
observation being predicted), a list (with elements \code{linear.predictors}
and \code{se.fit} if \code{se.fit=TRUE}), a matrix (\code{type="fitted"} or \code{type="fitted.ind"}),
a data frame, or a design matrix.  
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\references{
Hannah M, Quigley P: Presentation of ordinal regression analysis on the original
scale.  Biometrics 52:771--5; 1996.
}
\seealso{
\code{\link{lrm}}, \code{\link{predict.Design}}, \code{\link{naresid}}, \code{\link{contrast.Design}}
}
\examples{
# See help for predict.Design for several binary logistic
# regression examples


# Examples of predictions from ordinal models
set.seed(1)
y <- factor(sample(1:3, 400, TRUE), 1:3, c('good','better','best'))
x1 <- runif(400)
x2 <- runif(400)
f <- lrm(y ~ rcs(x1,4)*x2)
predict(f, type="fitted.ind")[1:10,]   #gets Prob(better) and all others
d <- data.frame(x1=.5,x2=.5)
predict(f, d, type="fitted")        # Prob(Y>=j) for new observation
predict(f, d, type="fitted.ind")    # Prob(Y=j)
predict(f, d, type='mean', codes=TRUE) # predicts mean(y) using codes 1,2,3
}
\keyword{models}
\keyword{regression}
\concept{logistic regression model}
