\name{contrast.Design}
\alias{contrast}
\alias{contrast.Design}
\alias{print.contrast.Design}
\title{
General Contrasts of Regression Coefficients
}
\description{
This function computes one or more contrasts of the estimated
regression coefficients in a fit from one of the functions in Design,
along with standard errors, confidence limits, t or Z statistics, P-values.
General contrasts are handled by obtaining the design matrix for two
sets of predictor settings (\code{a}, \code{b}) and subtracting the
corresponding rows of the two design matrics to obtain a new contrast
design matrix for testing the \code{a} - \code{b} differences.  This allows for
quite general contrasts (e.g., estimated differences in means between
a 30 year old female and a 40 year old male).
This can also be used
to obtain a series of contrasts in the presence of interactions (e.g.,
female:male log odds ratios for several ages when the model contains
age by sex interaction).  Another use of \code{contrast} is to obtain
center-weighted (Type III test) and subject-weighted (Type II test)
estimates in a model containing treatment by center interactions.  For
the latter case, you can specify \code{type="average"} and an optional
\code{weights} vector to average the within-center treatment contrasts.
The design contrast matrix computed by \code{contrast.Design} can be used
by the \code{bootplot} and \code{confplot} functions to obtain bootstrap
nonparametric confidence intervals for contrasts.

By omitting the \code{b} argument, \code{contrast} can be used to obtain
an average or weighted average of a series of predicted values, along
with a confidence interval for this average.  This can be useful for
"unconditioning" on one of the predictors (see the next to last
example).

When more than one contrast is computed, the list created by
\code{contrast.Design} is suitable for plotting (with error bars or bands)
with \code{xYplot} or \code{Dotplot} (see the last example).
}
\usage{
contrast(fit, \dots)
\method{contrast}{Design}(fit, a, b, cnames=NULL, 
         type=c("individual", "average"), 
         weights="equal", conf.int=0.95, \dots)

\method{print}{contrast.Design}(x, X=FALSE, fun=function(u)u, ...)
}
\arguments{
\item{fit}{
a fit of class \code{"Design"}
}
\item{a}{
a list containing settings for all predictors that you do not wish to
set to default (adjust-to) values.  Usually you will specify two
variables in this list, one set to a constant and one to a sequence of
values, to obtain contrasts for the sequence of values of an
interacting factor.  The \code{gendata} function will generate the
necessary combinations and default values for unspecified predictors.
}
\item{b}{
another list that generates the same number of observations as \code{a},
unless one of the two lists generates only one observation.  In that
case, the design matrix generated from the shorter list will have its
rows replicated so that the contrasts assess several differences
against the one set of predictor values.  This is useful for comparing
multiple treatments with control, for example.  If \code{b} is missing, the
design matrix generated from \code{a} is analyzed alone.
}
\item{cnames}{
vector of character strings naming the contrasts when
\code{type="individual"}.  Usually \code{cnames} is not necessary as
\code{contrast.Design} tries to name the contrasts by examining which
predictors are varying consistently in the two lists.  \code{cnames} will
be needed when you contrast "non-comparable" settings, e.g., you
compare \code{list(treat="drug", age=c(20,30))} with
\code{list(treat="placebo"), age=c(40,50))}
}
\item{type}{
set \code{type="average"} to average the individual contrasts (e.g., to
obtain a Type II or III contrast)
}
\item{weights}{
a numeric vector, used when \code{type="average"}, to obtain weighted contrasts
}
\item{conf.int}{
confidence level for confidence intervals for the contrasts
}
\item{\dots}{unused}
\item{x}{result of \code{contrast}}
\item{X}{
set \code{X=TRUE} to  print design matrix used in computing the contrasts (or
the average contrast)
}
\item{fun}{
a function to transform the contrast, SE, and lower and upper
confidence limits before printing.  For example, specify \code{fun=exp} to
anti-log them for logistic models.
}}
\value{
a list of class \code{"contrast.Design"} containing the elements
\code{Contrast}, \code{SE}, \code{Z}, \code{var}, \code{df.residual}
\code{Lower}, \code{Upper}, \code{Pvalue}, \code{X}, \code{cnames}, which denote the contrast
estimates, standard errors, Z or t-statistics, variance matrix,
residual degrees of freedom (this is \code{NULL} if the model was not
\code{ols}), lower and upper confidence limits, 2-sided P-value, design
matrix, and contrast names (or \code{NULL}).
}
\author{
Frank Harrell\cr
Division of Biostatistics and Epidemiology\cr
Department of Health Evaluation Sciences\cr
University of Virginia Health Sciences Center\cr
fharrell@virginia.edu
}
\seealso{
\code{\link{predict.Design}}, \code{\link{gendata}}, \code{\link{bootcov}}, \code{\link{summary.Design}},
\code{\link{anova.Design}}, \code{\link{plot.Design}}
}
\examples{
set.seed(1)
age <- rnorm(200,40,12)
sex <- factor(sample(c('female','male'),200,TRUE))
logit <- (sex=='male') + (age-40)/5
y <- ifelse(runif(200) <= plogis(logit), 1, 0)
f <- lrm(y ~ pol(age,2)*sex)
# Compare a 30 year old female to a 40 year old male
# (with or without age x sex interaction in the model)
contrast(f, list(sex='female', age=30), list(sex='male', age=40))


# For a model containing two treatments, centers, and treatment
# x center interaction, get 0.95 confidence intervals separately
# by cente
center <- factor(sample(letters[1:8],500,TRUE))
treat  <- factor(sample(c('a','b'),  500,TRUE))
y      <- 8*(treat=='b') + rnorm(500,100,20)
f <- ols(y ~ treat*center)


lc <- levels(center)
contrast(f, list(treat='b', center=lc),
            list(treat='a', center=lc))


# Get 'Type III' contrast: average b - a treatment effect over
# centers, weighting centers equally (which is almost always
# an unreasonable thing to do)
contrast(f, list(treat='b', center=lc),
            list(treat='a', center=lc),
         type='average')


# Get 'Type II' contrast, weighting centers by the number of
# subjects per center.  Print the design contrast matrix used.
k <- contrast(f, list(treat='b', center=lc),
                 list(treat='a', center=lc),
              type='average', weights=table(center))
print(k, X=TRUE)
# Note: If other variables had interacted with either treat 
# or center, we may want to list settings for these variables
# inside the list()'s, so as to not use default settings


# For a 4-treatment study, get all comparisons with treatment 'a'
treat  <- factor(sample(c('a','b','c','d'),  500,TRUE))
y      <- 8*(treat=='b') + rnorm(500,100,20)
dd     <- datadist(treat,center); options(datadist='dd')
f <- ols(y ~ treat*center)
lt <- levels(treat)
contrast(f, list(treat=lt[-1]),
            list(treat=lt[ 1]),
         cnames=paste(lt[-1],lt[1],sep=':'), conf.int=1-.05/3)


# Compare each treatment with average of all others
for(i in 1:length(lt)) {
  cat('Comparing with',lt[i],'\n\n')
  print(contrast(f, list(treat=lt[-i]),
                    list(treat=lt[ i]), type='average'))
}
options(datadist=NULL)

# Six ways to get the same thing, for a variable that
# appears linearly in a model and does not interact with
# any other variables.  We estimate the change in y per
# unit change in a predictor x1.  Methods 4, 5 also
# provide confidence limits.  Method 6 computes nonparametric
# bootstrap confidence limits.  Methods 2-6 can work
# for models that are nonlinear or non-additive in x1.
# For that case more care is needed in choice of settings
# for x1 and the variables that interact with x1.


\dontrun{
coef(fit)['x1']                            # method 1
diff(predict(fit, gendata(x1=c(0,1))))     # method 2
g <- Function(fit)                         # method 3
g(x1=1) - g(x1=0)
summary(fit, x1=c(0,1))                    # method 4
k <- contrast(fit, list(x1=1), list(x1=0)) # method 5
print(k, X=TRUE)
fit <- update(fit, x=TRUE, y=TRUE)               # method 6
b <- bootcov(fit, B=500, coef.reps=TRUE)
bootplot(b, X=k$X)    # bootstrap distribution and CL


# In a model containing age, race, and sex,
# compute an estimate of the mean response for a
# 50 year old male, averaged over the races using
# observed frequencies for the races as weights


f <- ols(y ~ age + race + sex)
contrast(f, list(age=50, sex='male', race=levels(race)),
         type='average', weights=table(race))
}


# Plot the treatment effect (drug - placebo) as a function of age
# and sex in a model in which age nonlinearly interacts with treatment
# for females only
set.seed(1)
n <- 800
treat <- factor(sample(c('drug','placebo'), n,TRUE))
sex   <- factor(sample(c('female','male'),  n,TRUE))
age   <- rnorm(n, 50, 10)
y     <- .05*age + (sex=='female')*(treat=='drug')*.05*abs(age-50) + rnorm(n)
f     <- ols(y ~ rcs(age,4)*treat*sex)
d     <- datadist(age, treat, sex); options(datadist='d')
# show separate estimates by treatment and sex
plot(f, age=NA, treat=NA, sex='female')
plot(f, age=NA, treat=NA, sex='male')
ages  <- seq(35,65,by=5); sexes <- c('female','male')
w     <- contrast(f, list(treat='drug',    age=ages, sex=sexes),
                     list(treat='placebo', age=ages, sex=sexes))
xYplot(Cbind(Contrast, Lower, Upper) ~ age | sex, data=w,
       ylab='Drug - Placebo')
xYplot(Cbind(Contrast, Lower, Upper) ~ age, groups=sex, data=w,
       ylab='Drug - Placebo', method='alt bars')
options(datadist=NULL)
}
\keyword{htest}
\keyword{models}
\keyword{regression}
% Converted by Sd2Rd version 1.21.
