\name{Design}
\alias{Design}
\title{
Design Methods and Generic Functions
}
\description{
This is a series of special transformation functions (\code{asis}, \code{pol},
\code{lsp}, \code{rcs}, \code{catg}, \code{scored}, \code{strat}, \code{matrx}), fitting functions
(e.g.,  \code{lrm},\code{cph}, \code{psm}, or \code{ols}), and generic analysis functions
(\code{anova.Design}, \code{summary.Design}, \code{predict.Design}, 
\code{plot.Design}, \code{survplot},
\code{fastbw}, \code{validate}, \code{calibrate}, \code{specs.Design},
\code{which.influence}, \code{latex.Design}, \code{nomogram.Design}, \code{datadist}, \code{gendata}) 
that help automate many
analysis steps, e.g. fitting restricted interactions and multiple
stratification variables, analysis of variance (with tests of linearity
of each factor and pooled tests), plotting effects of variables in the
model, estimating and graphing effects of variables that appear non-linearly in the
model using e.g. inter-quartile-range hazard ratios, bootstrapping
model fits, and constructing nomograms for obtaining predictions manually. 
Behind the scene is the \code{Design} function, 
called by a modified version of \code{model.frame.default} to
store extra attributes. \code{Design()} is not intended to be
called by users.  
\code{Design} causes detailed design attributes
and descriptions of the distribution of predictors to be stored 
in an attribute of the \code{terms} component called \code{Design}.
In addition to \code{model.frame.default} being replaced by a modified
version, \code{[.} and \code{[.factor} are replaced by versions which carry
along the \code{label} attribute of a variable.  In this way, when an
\code{na.action} function is called to subset out NAs, labels are still
defined for variables in the model.
}
\synopsis{
Design(mf, allow.offset=TRUE, intercept=1)
}
\usage{
Design(mf, allow.offset=TRUE, intercept=1)
# not to be called by the user; called by fitting routines
# dist <- datadist(x1,x2,sex,age,race,bp)   
# or dist <- datadist(my.data.frame)
# Can omit call to datadist if not using summary.Design, plot.Design, 
# survplot.Design, or if all variable settings are given to them
# options(datadist="dist")
# f <- fitting.function(formula = y ~ rcs(x1,4) + rcs(x2,5) + x1\%ia\%x2 +
#                       rcs(x1,4)\%ia\%rcs(x2,5) +
#                       strat(sex)*age + strat(race)*bp)
# See Design.trans for rcs, strat, etc.
# \%ia\% is restricted interaction - not doubly nonlinear
# for x1 by x2 this uses the simple product only, but pools x1*x2
# effect with nonlinear function for overall tests
# specs(f)
# anova(f)
# summary(f)
# fastbw(f)
# pred <- predict(f, newdata=expand.grid(x1=1:10,x2=3,sex="male",
#                 age=50,race="black"))
# pred <- predict(f, newdata=gendata(f, x1=1:10, x2=3, sex="male"))
# This leaves unspecified variables set to reference values from datadist
# pred.combos <- gendata(f, nobs=10)   # Use X-windows to edit predictor settings
# predict(f, newdata=pred.combos)
# plot(f, x1=NA)
# latex(f)
# nomogram(f)
}
\arguments{
  \item{mf}{a model frame}
  \item{allow.offset}{set to \code{TRUE} if model fitter allows an
	offset term}
  \item{intercept}{1 if an ordinary intercept is present, 0 otherwise}
}
\value{
  a data frame augmented with additional information about the
  predictors and model formulation
  }
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link{Design.trans}}, \code{\link{Design.Misc}}, \code{\link{cph}}, \code{\link{lrm}}, \code{\link{ols}}, \code{\link{specs.Design}}, \code{\link{anova.Design}},
\code{\link{summary.Design}}, \code{\link{predict.Design}}, \code{\link{gendata}}, \code{\link{plot.Design}}, \code{\link{fastbw}}, 
\code{\link{validate}}, \code{\link{calibrate}}, \code{\link{which.influence}},
\code{\link[Hmisc]{latex}}, \code{\link{latex.Design}}, \code{\link{model.frame.default}}, \code{\link{datadist}}, \code{\link[Hmisc]{describe}},
\code{\link{nomogram}}, \code{\link{vif}}, \code{\link[Hmisc]{dataRep}}
}
\examples{
\dontrun{
library(Design, first=TRUE)  # omit first for R
dist <- datadist(data=2)     # can omit if not using summary, plot, survplot,
                             # or if specify all variable values to them. Can
                             # also  defer.  data=2: get distribution summaries
                             # for all variables in search position 2
                             # run datadist once, for all candidate variables
dist <- datadist(age,race,bp,sex,height)   # alternative
options(datadist="dist")
f <- cph(Surv(d.time, death) ~ rcs(age,4)*strat(race) +
         bp*strat(sex)+lsp(height,60),x=TRUE,y=TRUE)
anova(f)
anova(f,age,height)          # Joint test of 2 vars
fastbw(f)
summary(f, sex="female")     # Adjust sex to "female" when testing
                             # interacting factor bp
plot(f, age=NA, height=NA)   # 3-D plot
plot(f, age=10:70, height=60)
latex(f)                     # LaTeX representation of fit


f <- lm(y ~ x)               # Can use with any fitting function that
                             # calls model.frame.default, e.g. lm, glm
specs.Design(f)              # Use .Design since class(f)="lm"
anova(f)                     # Works since Varcov(f) (=Varcov.lm(f)) works
fastbw(f)
options(datadist=NULL)
f <- ols(y ~ x1*x2)          # Saves enough information to do fastbw, anova
anova(f)                     # Will not do plot.Design since distributions
fastbw(f)                    # of predictors not saved
plot(f, x1=seq(100,300,by=.5), x2=.5) 
                             # all values defined - don't need datadist
dist <- datadist(x1,x2)      # Equivalent to datadist(f)
options(datadist="dist")
plot(f, x1=NA, x2=.5)        # Now you can do plot, summary
nomogram(f, interact=list(x2=c(.2,.7)))
}
}
\keyword{models}
\keyword{regression}
\keyword{survival}
\keyword{math}
\keyword{manip}
\keyword{methods}
\concept{logistic regression model}
