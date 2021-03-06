\name{predict.Design}
\alias{predict.Design}
\alias{predictDesign}
\alias{predict.bj}
\alias{predict.cph}
\alias{predict.glmD}
\alias{predict.glsD}
\alias{predict.ols}
\alias{predict.psm}

\title{Predicted Values from Model Fit}
\description{
The \code{predict} function is used to obtain a variety of values or
predicted values from either the data used to fit the model (if
\code{type="adjto"} or \code{"adjto.data.frame"} or if \code{x=TRUE} or
\code{linear.predictors=TRUE} were specified to the modeling function), or from
a new dataset. Parameters such as knots and factor levels used in creating 
the design
matrix in the original fit are "remembered".
See the \code{Function} function for another method for computing the
linear predictors.
}
\usage{
\method{predict}{bj}(object, newdata,
        type=c("lp", "x", "data.frame",
                 "terms", "adjto", "adjto.data.frame", "model.frame"),
        se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
        incl.non.slopes,
        non.slopes, kint=1, na.action=na.keep, expand.na=TRUE,
        center.terms=TRUE, \dots) # for bj

\method{predict}{cph}(object, newdata,
        type=c("lp", "x",
                 "data.frame", "terms", "adjto", "adjto.data.frame",
                 "model.frame"),
        se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
        incl.non.slopes=NULL,
        non.slopes=NULL, kint=1, na.action=na.keep, expand.na=TRUE,
        center.terms=TRUE, \dots) # cph

\method{predict}{glmD}(object, newdata,
        type= c("lp", "x", "data.frame",
                 "terms", "adjto", "adjto.data.frame", "model.frame"),
        se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
        incl.non.slopes,
        non.slopes, kint=1, na.action=na.keep, expand.na=TRUE,
        center.terms=TRUE, \dots) # glmD

\method{predict}{glsD}(object, newdata,
        type=c("lp", "x", "data.frame",
                 "terms", "adjto", "adjto.data.frame", "model.frame"),
        se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
        incl.non.slopes,
        non.slopes, kint=1, na.action=na.keep, expand.na=TRUE,
        center.terms=TRUE, \dots) # glsD

\method{predict}{ols}(object, newdata,
        type=c("lp", "x", "data.frame",
                 "terms", "adjto", "adjto.data.frame", "model.frame"),
        se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
        incl.non.slopes,
        non.slopes, kint=1, na.action=na.keep, expand.na=TRUE,
        center.terms=TRUE, \dots) # ols

\method{predict}{psm}(object, newdata,
        type=c("lp", "x", "data.frame",
                 "terms", "adjto", "adjto.data.frame", "model.frame"),
        se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
        incl.non.slopes,
        non.slopes, kint=1, na.action=na.keep, expand.na=TRUE,
        center.terms=TRUE, \dots) # psm
}
\arguments{
\item{object}{a fit object with a \code{Design} fitting function}
\item{newdata}{
An S data frame, list or a matrix specifying new data for which predictions
are desired.  If \code{newdata} is a list, it is converted to a matrix first.
A matrix is converted to a data frame.  For the matrix form, categorical
variables (\code{catg} or \code{strat}) must be coded as integer category
numbers corresponding to the order in which value labels were stored.
For list or matrix forms, \code{matrx} factors must be given a single
value.  If this single value is the S missing value \code{NA}, the adjustment
values of matrx (the column medians) will later replace this value.
If the single value is not \code{NA}, it is propagated throughout the columns
of the \code{matrx} factor.  For \code{factor} variables having numeric levels,
you can specify the numeric values in \code{newdata} without first converting
the variables to factors.  These numeric values are checked to make sure
they match a level, then the variable is converted internally to a \code{factor}.
It is most typical to use a data frame
for newdata, and the S function \code{expand.grid} is very handy here.
For example, one may specify 
\cr
\code{newdata=expand.grid(age=c(10,20,30),}
\cr
   \code{race=c("black","white","other"),}
\cr
   \code{chol=seq(100,300,by=25))}.
}
\item{type}{
Type of output desired.  The default is \code{"lp"} to get the linear predictors -
predicted \eqn{X\beta}{X beta}.  For Cox models, these predictions are centered.
You may specify \code{"x"} to get an expanded design matrix
at the desired combinations of values, \code{"data.frame"} to get an
S data frame of the combinations, \code{"model.frame"} to get a data frame
of the transformed predictors, \code{"terms"} to get a matrix with
each column being the linear combination of variables making up
a factor, \code{"adjto"} to return a vector of \code{limits[2]} (see \code{datadist}) in coded
form, and \code{"adjto.data.frame"} to return a data frame version of these
central adjustment values.  If \code{newdata} is not given, \code{predict}
will attempt to return information stored with the fit object if the
appropriate options were used with the modeling function (e.g., \code{x, y, linear.predictors, se.fit}).
}
\item{se.fit}{
Defaults to \code{FALSE}.  If \code{type="linear.predictors"}, set \code{se.fit=TRUE} to return
a list with components \code{linear.predictors} and \code{se.fit} instead of just
a vector of fitted values.
}
\item{conf.int}{
Specify \code{conf.int} as a positive fraction to obtain upper and lower
confidence intervals (e.g., \code{conf.int=0.95}).  The \eqn{t}-distribution is
used in the calculation for \code{ols} fits.  Otherwise, the normal
critical value is used.
}
\item{conf.type}{
specifies the type of confidence interval.  Default is for the mean.
For \code{ols} fits there is the option of obtaining confidence limits for
individual predicted values by specifying \code{conf.type="individual"}.
}
\item{incl.non.slopes}{
Default is \code{TRUE} if \code{non.slopes} or \code{kint} is specified, the model has a scale
parameter (e.g., a parametric survival model), or \code{type!="x"}.
Otherwise the default is \code{FALSE}.
Set to \code{TRUE} to use an intercept in the prediction if the model has
any intercepts (except for \code{type="terms"} which doesn't need
intercepts).  Set to \code{FALSE} to get predicted \eqn{X\beta}{X beta} ignoring intercepts.
}
\item{non.slopes}{
For models such as the ordinal logistic models containing more than
one intercept, this specifies dummy variable values to pick off intercept(s)
to use in computing predictions.  For example, if there are 3 intercepts,
use \code{non.slopes=c(0,1,0)} to use the second.  Default is
\code{c(1,0,\dots,0)}.  You may alternatively specify \code{kint}.
}
\item{kint}{
a single integer specifying the number of the intercept to use in
multiple-intercept models
}
\item{na.action}{
Function to handle missing values in \code{newdata}.  For predictions
"in data", the same \code{na.action} that was used during model fitting is
used to define an \code{naresid} function to possibly restore rows of the data matrix
that were deleted due to NAs.  For predictions "out of data", the default
\code{na.action} is \code{na.keep}, resulting in NA predictions when a row of
\code{newdata} has an NA.  Whatever \code{na.action} is in effect at the time
for "out of data" predictions, the corresponding \code{naresid} is used also.
}
\item{expand.na}{
set to \code{FALSE} to keep the \code{naresid} from having any effect, i.e., to keep
from adding back observations removed because of NAs in the returned object.
If \code{expand.na=FALSE}, the \code{na.action} attribute will be added to the returned
object.
}
\item{center.terms}{
set to \code{FALSE} to suppress subtracting the mean from columns of the design
matrix before computing terms with \code{type="terms"}.
}
\item{\dots}{ignored}
}
\details{
\code{datadist} and \code{options(datadist=)} should be run before \code{predict.Design}
if using \code{type="adjto"}, \code{type="adjto.data.frame"}, or \code{type="terms"},
or if the fit is a Cox model fit and you are requesting \code{se.fit=TRUE}.
For these cases, the adjustment values are needed (either for the
returned result or for the correct covariance matrix computation).
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link{plot.Design}}, \code{\link{summary.Design}}, \code{\link{Design}}, \code{\link{Design.trans}}, \code{\link{predict.lrm}},
\code{\link{residuals.cph}}, \code{\link{naresid}}, \code{\link{datadist}}, \code{\link{gendata}}, \code{\link{Function.Design}}, \code{\link[Hmisc]{reShape}},
\code{\link[Hmisc]{xYplot}}, \code{\link{contrast.Design}}
}
\examples{
n <- 1000    # define sample size
set.seed(17) # so can reproduce the results
age            <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
cholesterol    <- rnorm(n, 200, 25)
sex            <- factor(sample(c('female','male'), n,TRUE))
treat          <- factor(sample(c('a','b','c'), n,TRUE))


# Specify population model for log odds that Y=1
L <- .4*(sex=='male') + .045*(age-50) +
  (log(cholesterol - 10)-5.2)*(-2*(sex=='female') + 2*(sex=='male')) +
  .3*sqrt(blood.pressure-60)-2.3 + 1*(treat=='b')
# Simulate binary y to have Prob(y=1) = 1/[1+exp(-L)]
y <- ifelse(runif(n) < plogis(L), 1, 0)


ddist <- datadist(age, blood.pressure, cholesterol, sex, treat)
options(datadist='ddist')


fit <- lrm(y ~ rcs(blood.pressure,4) + 
           sex * (age + rcs(cholesterol,4)) + sex*treat*age)


# Use xYplot to display predictions in 9 panels, with error bars,
# with superposition of two treatments


dat <- expand.grid(treat=levels(treat),sex=levels(sex),
                   age=c(20,40,60),blood.pressure=120,
                   cholesterol=seq(100,300,length=10))
# Add variables linear.predictors and se.fit to dat
dat <- cbind(dat, predict(fit, dat, se.fit=TRUE))
# xYplot in Hmisc extends xyplot to allow error bars
xYplot(Cbind(linear.predictors,linear.predictors-1.96*se.fit,
             linear.predictors+1.96*se.fit) ~ cholesterol | sex*age,
       groups=treat, data=dat, type='b')




# Since blood.pressure doesn't interact with anything, we can quickly and
# interactively try various transformations of blood.pressure, taking
# the fitted spline function as the gold standard. We are seeking a
# linearizing transformation even though this may lead to falsely
# narrow confidence intervals if we use this data-dredging-based transformation


bp <- 70:160
logit <- predict(fit, expand.grid(treat="a", sex='male', age=median(age),
                 cholesterol=median(cholesterol),
                 blood.pressure=bp), type="terms")[,"blood.pressure"]
#Note: if age interacted with anything, this would be the age
#      "main effect" ignoring interaction terms
#Could also use
#   logit <- plot(f, age=ag, \dots)$x.xbeta[,2]
#which allows evaluation of the shape for any level of interacting
#factors.  When age does not interact with anything, the result from
#predict(f, \dots, type="terms") would equal the result from
#plot if all other terms were ignored


plot(bp^.5, logit)               # try square root vs. spline transform.
plot(bp^1.5, logit)              # try 1.5 power
plot(sqrt(bp-60), logit)


#Some approaches to making a plot showing how predicted values
#vary with a continuous predictor on the x-axis, with two other
#predictors varying


combos <- gendata(fit, age=seq(10,100,by=10), cholesterol=c(170,200,230),
                  blood.pressure=c(80,120,160))
#treat, sex not specified -> set to mode
#can also used expand.grid


combos$pred <- predict(fit, combos)
xyplot(pred ~ age | cholesterol*blood.pressure, data=combos, type='l')
xYplot(pred ~ age | cholesterol, groups=blood.pressure, data=combos, type='l')
Key()   # Key created by xYplot
xYplot(pred ~ age, groups=interaction(cholesterol,blood.pressure),
       data=combos, type='l', lty=1:9)
Key()


#Add upper and lower 0.95 confidence limits for individuals
combos <- cbind(combos, predict(fit, combos, conf.int=.95))
xYplot(Cbind(linear.predictors, lower, upper) ~ age | cholesterol,
       groups=blood.pressure, data=combos, type='b')
Key()


# Plot effects of treatments (all pairwise comparisons) vs.
# levels of interacting factors (age, sex)


d <- gendata(fit, treat=levels(treat), sex=levels(sex), age=seq(30,80,by=10))
x <- predict(fit, d, type="x")
betas <- fit$coef
cov   <- fit$var


i <- d$treat=="a"; xa <- x[i,]; Sex <- d$sex[i]; Age <- d$age[i]
i <- d$treat=="b"; xb <- x[i,]
i <- d$treat=="c"; xc <- x[i,]


doit <- function(xd, lab) {
  xb <- xd\%*\%betas
  se <- apply((xd \%*\% cov) * xd, 1, sum)^.5
  q <- qnorm(1-.01/2)   # 0.99 confidence limits
  lower <- xb - q * se; upper <- xb + q * se
  #Get odds ratios instead of linear effects
  xb <- exp(xb); lower <- exp(lower); upper <- exp(upper)
  #First elements of these agree with 
  #summary(fit, age=30, sex='female',conf.int=.99))
  for(sx in levels(Sex)) {
    j <- Sex==sx
    errbar(Age[j], xb[j], upper[j], lower[j], xlab="Age", 
           ylab=paste(lab,"Odds Ratio"), ylim=c(.1,20), log='y')
    title(paste("Sex:",sx))
    abline(h=1, lty=2)
  }
}


par(mfrow=c(3,2), oma=c(3,0,3,0))
doit(xb - xa, "b:a")
doit(xc - xa, "c:a")
doit(xb - xa, "c:b")


# NOTE: This is much easier to do using contrast.Design


\dontrun{
#A variable state.code has levels "1", "5","13"
#Get predictions with or without converting variable in newdata to factor
predict(fit, data.frame(state.code=c(5,13)))
predict(fit, data.frame(state.code=factor(c(5,13))))


#Use gendata function (gendata.Design) for interactive specification of
#predictor variable settings (for 10 observations)
df <- gendata(fit, nobs=10, viewvals=TRUE)
df$predicted <- predict(fit, df)  # add variable to data frame
df


df <- gendata(fit, age=c(10,20,30))  # leave other variables at ref. vals.
predict(fit, df, type="fitted")


# See reShape (in Hmisc) for an example where predictions corresponding to 
# values of one of the varying predictors are reformatted into multiple
# columns of a matrix
}
options(datadist=NULL)
}
\keyword{models}
\keyword{regression}


