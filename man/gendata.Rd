\name{gendata}
\alias{gendata}
\alias{gendata.Design}
\alias{gendata.default}
\title{
Generate Data Frame with Predictor Combinations
}
\description{
If \code{nobs} is not specified, allows user to specify predictor settings
by e.g. \code{age=50, sex="male"}, and any omitted predictors are set to
reference values (default=median for continuous variables, first level
for categorical ones - see \code{datadist}).  If any predictor has more than one
value given, \code{expand.grid} is called to generate all possible combinations
of values.  If \code{nobs} is given, a data frame is first generated which has
\code{nobs} of adjust-to values duplicated.  Then an editor window is opened
which allows the user to subset the variable names down to ones which she
intends to vary (this streamlines the \code{data.ed} step).  Then, if any
predictors kept are discrete and \code{viewvals=TRUE}, a window (using \code{page})
is opened defining the possible values of this subset, to facilitate
data editing.  Then the \code{data.ed} function is invoked to allow interactive
overriding of predictor settings in the \code{nobs} rows.  The subset of
variables are combined with the other predictors which were not
displayed with \code{data.ed}, and a final full data frame is returned.
\code{gendata} is most useful for creating a \code{newdata} data frame to pass
to \code{predict}.
}
\usage{
gendata(fit, \dots)
\method{gendata}{Design}(fit, nobs, viewvals=FALSE,
  editor=.Options$editor, \dots, factors)
\method{gendata}{default}(fit, \dots)
}
\arguments{
\item{fit}{
a fit object created with \code{Design} in effect
}
\item{nobs}{
number of observations to create if doing it interactively using X-windows
}
\item{viewvals}{
if \code{nobs} is given, set \code{viewvals=TRUE} to open a window displaying the
possible value of categorical predictors
}
\item{editor}{
editor to use to edit the list of variable names to consider.
Default is \code{options(editor=)} value (\code{"xedit"} is this is not specified
by \code{using.X()==TRUE}.
}
\item{...}{
predictor settings, if \code{nobs} is not given. 
}
\item{factors}{
a list containing predictor settings with their names.  This is an
alternative to specifying the variables separatey in \dots .
}}
\value{
a data frame with all predictors, and an attribute \code{names.subset} if
\code{nobs} is specified.  This attribute contains the vector of variable
names for predictors which were passed to \code{data.ed} and hence were
allowed to vary.  If neither \code{nobs} nor any predictor settings were
given, returns a data frame with adjust-to values.
}
\section{Side Effects}{
optionally writes to the terminal, opens X-windows, and generates a
temporary file using \code{sink}.
}
\details{
if you have a variable in \code{\dots} that is named \code{n, no, nob, nob}, add
\code{nobs=FALSE} to the invocation to prevent that variable from being misrecognized
as \code{nobs}
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link{predict.Design}}, \code{\link{survest.cph}}, \code{\link{survest.psm}}, \code{\link{Design.Misc}}, \code{\link{expand.grid}}, \code{\link{data.entry}}, \code{\link{page}},
\code{\link{print.datadist}}, \code{\link{edit.data.frame}},
\code{\link{plot.Design}}
}
\examples{
set.seed(1)
age <- rnorm(200, 50, 10)
sex <- factor(sample(c('female','male'),200,TRUE))
race <- factor(sample(c('a','b','c','d'),200,TRUE))
y <- sample(0:1, 200, TRUE)
dd <- datadist(age,sex,race)
options(datadist="dd")
f <- lrm(y ~ age*sex + race)
gendata(f)
gendata(f, age=50)
d <- gendata(f, age=50, sex="female")  # leave race=reference category
d <- gendata(f, age=c(50,60), race=c("b","a"))  # 4 obs.
d$Predicted <- predict(f, d, type="fitted")
d      # Predicted column prints at the far right
options(datadist=NULL)
\dontrun{
d <- gendata(f, nobs=5, view=TRUE)        # 5 interactively defined obs.
d[,attr(d,"names.subset")]             # print variables which varied
predict(f, d)
}
}
\keyword{methods}
\keyword{models}
\keyword{regression}
\keyword{manip}
