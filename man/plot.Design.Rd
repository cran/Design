\name{plot.Design}
\alias{plot.Design}
\alias{print.plot.Design}
\alias{perimeter}
\alias{lines.perimeter}
\alias{Legend}
\alias{Legend.default}
\alias{Legend.plot.Design}
\alias{datadensity}
\alias{datadensity.plot.Design}
\title{
Plot Effects of Variables
}
\description{
Plots the effect of one or two predictors on the linear
predictor or X beta scale, or on some transformation of that scale.
The predictor is always plotted in its original coding on the \eqn{x} or \eqn{y}-axis.
\code{perimeter} is a function used to generate the boundary of data to plot
when a 3-d plot is made.  It finds the area where there are sufficient
data to generate believable interaction fits.

\code{Legend} is a generic function for adding legends to an existing graph
according to the specific plot made by \code{plot.Design}.  The specific
\code{Legend} method for \code{plot.Design} is \code{Legend.plot.Design}.  It handles
legends for \code{image} plots.  For other plots with one or more curves,
make legends using the \code{label.curves} parameter.

\code{datadensity} is a function for showing the data density (raw data) on
each curve generated for curve-type plots.  This is a rug plot showing
the location/density of data values for the \eqn{x}-axis variable.  If
there was a second variable specified to \code{plot} that generated separate
curves, the data density specific to each class of points is shown.
This assumes that the second variable was categorical.  The rug plots
are drawn by \code{scat1d}.

To plot effects instead of estimates (e.g., treatment differences as a
function of interacting factors) see \code{contrast.Design} and \code{summary.Design}.
}
\usage{
perimeter(x, y, xinc=diff(range(x))/10, 
          n=10, lowess.=FALSE)

\method{plot}{Design}(x, \dots, xlim, ylim, fun, xlab, ylab, 
     conf.int=.95, conf.type=c('mean','individual'),
     add=FALSE, label.curves=TRUE,
     eye=c(-6,-6,9), lty, col=1, lwd, lwd.conf=1, pch=1, 
     adj.zero=FALSE, ref.zero=FALSE, adj.subtitle, cex.adj,
     non.slopes, time=NULL, loglog=FALSE, val.lev=FALSE,
     digits=4, log="", perim,
     method=c("persp","contour","image","dotchart","default"),
     sortdot=c('neither','ascending','descending'),
     nlevels=10, name, zlim=range(zmat,na.rm=TRUE),
     vnames=c('labels','names'), abbrev=FALSE)
# Just say plot(fit, ...)

\method{print}{plot.Design}(x, \dots)

\method{lines}{perimeter}(x, \dots)

Legend(object, \dots)
\method{Legend}{plot.Design}(object, x, y, size=c(1,1), horizontal=TRUE,
       nint=50, fun, at, zlab, \dots)

\method{datadensity}{plot.Design}(object, x1, x2, \dots)
}
\arguments{
\item{fit}{
a fit object created with \code{Design()} in effect.  \code{options(datadist="d")}
must have been specified (where \code{d} was created by \code{datadist}), or
it must have been in effect with \code{fit} was created.
}
\item{\dots}{
The first variable in this list is displayed on the \eqn{x}-axis. Specify
\code{x=NA} to use the 
default display range, or any range you choose (e.g. 
\code{seq(0,100,by=2),c(2,3,7,14)}). 
The default list of values for which predictions are made
is taken as the list of unique values of the variable if they number fewer
than 11. For variables with \eqn{>10} unique values, 100 equally spaced
values in the range are used for plotting if the range is not specified.
If there is a second variable listed, and its range is \code{NA} or a single value,
that variable is displayed on the y-axis. If the second variable's range
has fewer than 40 levels, separate curves are generated for each value
of the variable. Otherwise, a three dimensional perspective plot is drawn
using 40 equally-spaced values of \code{y}. Names may be abbreviated.
\code{plot} senses that a variable is not to be displayed by checking if the list of values for
the variable is a scalar instead of a vector. Variables not specified are set to the default
adjustment value \code{limits[2]}, i.e. the median for continuous variables and a reference category for
non-continuous ones.  Due to a bug in S, the first variable mentioned
may not be named \code{x}.  This would cause the general scatterplot function
\code{plot} to be invoked by mistake.
Variables after the first or second specified to \code{plot} define adjustment settings.
For categorical variables, specify the class labels in quotes when specifying variable values.  If the levels of a categorical variable are numeric,
you may omit the quotes.  For variables not described using \code{datadist},
you must specify explicit ranges and adjustment settings for predictors
that were in the model.  Note that you can omit \code{variables} entirely.  In
that case, all non-interaction effects will be plotted automatically as
if you said \code{plot(fit, age=NA); plot(fit, sex=NA); \dots}.  In this case
you have no control over the settings of the variables for the x-axis,
i.e., \code{NA} is always assumed.

For a plot made up of multiple curves, these are extra graphical arguments
will be passed to \code{key} from \code{Legend}.  For \code{image} plots, these
arguments are passed to \code{par} and have temporary effect.
For \code{datadensity} these extra arguments are passed along to \code{scat1d}.
}
\item{x}{
first variable of a pair of predictors forming a 3-d plot, to specify
to \code{perim}.  For \code{Legend}, is either a vector of 1 or 2
\eqn{x}-coordinates or a list with elements \code{x} and \code{y} each
with 1 or 2 coordinates. For \code{method="image"} plots, 1 or 2
coordinates may be given, and for other plot types, 1 coordinate is
given.  A single coordinate represents the upper left corner of the
legend box.  For \code{Legend}, \code{x} and \code{y} are optional.  If
omitted, \code{locator} is used to position legends with the mouse.
For \code{lines.perimeter}, \code{x} is the result of \code{perimeter}.
For \code{print.plot.Design}, \code{x} is the result of
\code{plot.Design}.
}
\item{y}{
second variable of the pair for \code{perim}, or \eqn{y}-coordinates for
\code{Legend}.  If omitted, \code{x} is assumed to be a list with both
\code{x} and \code{y} components.
}
\item{xinc}{
increment in \code{x} over which to examine the density of \code{y} in \code{perimeter}
}
\item{n}{
within intervals of \code{x} for \code{perimeter}, takes the informative range of \code{y} to be
the \eqn{n}th smallest to the \eqn{n}th largest values of \code{y}.  If there aren't
at least 2\eqn{n} \code{y} values in the \code{x} interval, no \code{y} ranges are used
for that interval.
}
\item{lowess.}{
set to \code{TRUE} to have \code{lowess} smooth the data perimeters
}
\item{xlim}{
This parameter is seldom used, as limits are usually controlled with the
\code{variables} specifications.  One reason to use \code{xlim} is to plot a
\code{factor} variable on the x-axis that was created with the \code{cut2} function
with the \code{levels.mean} option, with \code{val.lev=TRUE} specified to \code{plot.Design}. 
In this case you may want the axis to
have the range of the original variable values given to \code{cut2} rather
than the range of the means within quantile groups.
}
\item{ylim}{
Range for plotting on response variable axis. Computed by default.
}
\item{fun}{
Function used to transform \eqn{X\beta}{X beta} and its confidence interval before plotting.
For example, to transform from a logit to a probability scale, use
\code{fun=function(x)1/(1+exp(-x))} or \code{fun=plogis}, 
and to take the anti-log, specify \code{fun=exp}.
for \code{Legend}, \code{fun} is
a function for transforming tick mark labels for color or gray scale
legends for \code{method="image"}.   For example, if \code{plot.Design} is used
to make an image plot of log odds ratios, specifying \code{fun=plogis} will
cause the color legend to be labeled with probability values rather
than log odds.
}
\item{xlab}{
Label for \code{x}-axis. Default is one given to \code{asis, rcs}, etc., which may have been
the \code{"label"} attribute of the variable.
}
\item{ylab}{
Label for \code{y}-axis (\code{z}-axis if perspective plot). If \code{fun} is not given,
default is \code{"log Odds"} for
\code{lrm}, \code{"log Relative Hazard"} for \code{cph}, name of the response
variable for \code{ols}, \code{TRUE} or \code{log(TRUE)} for \code{psm}, or \code{"X * Beta"} otherwise.
If \code{fun} is given, the default is \code{""}.
If \code{time} is given, the default is
\code{"(time) (units) Survival Probability"} or
\code{"log[-log S(time)]"} depending on the \code{loglog} parameter.
}
\item{conf.int}{
Default is \code{.95}.  Specify \code{FALSE} to suppress confidence bands.
}
\item{conf.type}{
specifies the type of confidence interval.  Default is for the mean.
For \code{ols} fits there is the option of obtaining confidence limits for
individual predicted values by specifying \code{conf.type="individual"}.
}
\item{add}{
Set to \code{TRUE} to add to an existing plot without drawing new axes. Default is \code{FALSE}.
See the warning note under \code{sortdot}.
}
\item{label.curves}{
Set to \code{FALSE} to suppress labeling of separate curves.
Default is \code{TRUE}, which
causes \code{labcurve} to be invoked to place labels at positions where the
curves are most separated, labeling each curve with the full curve label.
Set \code{label.curves} to a \code{list} to specify options to
\code{labcurve}, e.g., \code{label.curves=} \code{list(method="arrow", cex=.8)}.
These option names may be abbreviated in the usual way arguments
are abbreviated.  Use for example \code{label.curves=list(keys=letters[1:5])}
to draw single lower case letters on 5 curves where they are most
separated, and automatically position a legend
in the most empty part of the plot.  The \code{col}, \code{lty}, and \code{lwd} parameters
are passed automatically to \code{labcurve} although they may be overridden
here.
}
\item{eye}{
Argument to S \code{persp} function for defining perspective in 3-d plots.
Default is \code{(-6, -6, 9)}.
}
\item{lty}{
Vector of line types to use in plotting separate curves. Default is 1,2, \dots
}
\item{lwd}{
Vector of line widths corresponding to separate curves, 
default is \code{par("lwd")}.
}
\item{lwd.conf}{
scalar width of lines for confidence bands.  Default is 1.
}
\item{pch}{
symbol to use when plotting unconnected points when a categorical
variable is on the x-axis or when \code{method="dotchart"}.  Default is 1
(open circle).  See \code{points} for other values, or use the \code{show.pch} function in Hmisc.
}
\item{col}{
S color number for displaying curves.  Default is \code{1} (black).  Specify
a vector of integers to assign different colors to different curves.
}
\item{adj.subtitle}{
Set to \code{FALSE} to suppress subtitling the graph with the list of settings of non-graphed adjustment values. Default is \code{TRUE} if \code{ <= 6} non-plotted factors.
}
\item{cex.adj}{
\code{cex} parameter for size of adjustment settings in subtitles.  Default is
0.75 times \code{par("cex")}.
}
\item{adj.zero}{
Set to \code{TRUE} to adjust all non-plotted variables to 0 (or reference cell for
categorical variables) and to omit intercept(s) from consideration. Default
is \code{FALSE}.
}
\item{ref.zero}{
Subtract a constant from \eqn{X\beta}{X beta} before plotting so that
the reference value of the \code{x}-variable yields \code{y=0}.  This is
done before applying function \code{fun}.
}
\item{non.slopes}{
This is only useful in a multiple intercept model such as the ordinal
logistic model. There to use to second of three intercepts, for example,
specify \code{non.slopes=c(0,1,0)}. The default is \code{non.slopes=rep(0,k)}
if \code{adj.zero=TRUE}, where \code{k} is the number of intercepts in the model.
If \code{adj.zero=FALSE}, the default is \code{(1,0,0,\dots,0)}.
}
\item{time}{
Specify a single time \code{u} to cause function \code{survest} to be invoked
to plot the probability of surviving until time \code{u} when the fit
is from \code{cph} or \code{psm}.
}
\item{loglog}{
Specify \code{loglog=TRUE} to plot \code{log[-log(survival)]} instead of survival,
when \code{time} is given.
}
\item{val.lev}{
When plotting a categorical or strata factor with category labels that
are strings of legal numeric values, set to \code{TRUE} to use these values in
plotting.  An ordinary axis with uniform spacing will be used rather than
spacing dictated by the value labels.  When \code{val.lev=FALSE}, category
labels dictate how axis tick marks are made.  \code{val.lev} is used
typically when the variable being plotted is a categorical variable
that was collapsed into intervals, with the value label for a category
representing interval means or midpoints.  Such variables are created
for example by the \code{cut2} function, specifying \code{levels.mean=TRUE}.  For
plotting a discrete numeric variable you can specify \code{val.lev=TRUE} to
force plotting of the variable as if it were continuous.
}
\item{digits}{
Controls how ``adjust-to'' values are plotted.  The default is 4 significant
digits.
}
\item{log}{
Set \code{log="x", "y"} or \code{"xy"} to plot log scales on one or both axes.
}
\item{perim}{
names a matrix created by \code{perimeter} when used for 3-d plots of
two continuous predictors.  When the combination of variables is outside
the range in \code{perim}, that section of the plot is suppressed.  If \code{perim}
is omitted, 3-d plotting will use the marginal distributions of the
two predictors to determine the plotting region, when the grid is
not specified explicitly in \code{variables}.  When instead a series of
curves is being plotted, \code{perim} specifies a function having two
arguments.  The first is the vector of values of the first variable that
is about to be plotted on the x-axis.  The second argument is the single
value of the variable representing different curves, for the current
curve being plotted.  The function's returned value must be a logical
vector whose length is the same as that of the first argument, with
values \code{TRUE} if the corresponding point should be plotted for the
current curve, \code{FALSE} otherwise.  See one of the latter examples.
}
\item{method}{
For 3-d plots, use \code{method="persp"} for perspective plots (\code{persp()}, 
the default), \code{method="contour"} to use \code{contour()}, or \code{method="image"}
to use \code{image()}. Specify \code{method="dotchart"} to make a horizontal dot
chart to represent predicted values associated with categorical predictors.
The \code{log} argument does not apply to these plot types.
You can specify \code{method="default"} to get the default behaviour.
For \code{"dotchart"}, the \code{dotchart2} function in the Hmisc library is used.
}
\item{sortdot}{
applies when \code{method="dotchart"}.  The default is to plot the points in
the order requested for predictions.  Specify \code{method="ascending"} or
an abbreviation such as \code{method="a"} to sort in ascending order before
plotting the dot chart.  You may also specify \code{method="descending"}.
Unless \code{method="neither"}, specifying \code{add=TRUE} may not work properly.
}
\item{nlevels}{
number of contour levels if \code{method="contour"}
}
\item{name}{
Instead of specifying the variable to plot on the x-axis in the
\code{variables} list, you can specify one or more variables to plot by
specifying a vector of character string variable names in the
\code{name} argument.  Using this mode you cannot specify a list of
variable values to use; plotting is done as if you had said e.g.
\code{age=NA}.  Also, interacting factors can only be set to their reference values
using this notation.
}
\item{zlim}{
If 'type="persp"' controls the range for plottin in the
z-axis. Computed by default.
}
\item{vnames}{
applies when no x-variable is specified (i.e., when all predictors are
being plotted).  To override the default x-axis label in that case
(variable \code{"label"} attributes) to instead use variable names, specify
\code{vnames="names"}.
}
\item{object}{an object created by \code{plot.Design}}
\item{abbrev}{
Set to \code{TRUE} to use the \code{abbreviate} function to abbreviate levels of
categorical factors for labeling tick marks on the x-axis.
}
\item{size}{
}
\item{horizontal}{
}
\item{nint}{
see \code{image.legend}
}
\item{at}{
If \code{fun} is specified to \code{Legend}, \code{at} may be given.  \code{at} is a vector
of values at which to evaluate \code{fun} for drawing tick marks in the
color legend.  For example, if you want to show the median survival time
for a log-normal survival model whereas the linear predictor (log median)
was used in constructing the image plot, and if you want to place tick
marks at nice median values, specify \code{fun=exp, at=log(c(1,10,100,1000))}.
}
\item{zlab}{
label for \code{image} color axis legend.  Default is from the model
(e.g., \code{"Log Odds"}), but \code{zlab} will often be specified if
\code{fun} was specified to \code{plot.Design} or \code{Legend}.
}
\item{x1}{
data vector for first variable in \code{plot} (\code{x}-axis variable)
}
\item{x2}{
data vector for second variable in \code{plot} if it was not constant
(curve-generating variable)
}
}
\value{
\code{perimeter} returns a matrix of class \code{perimeter}.  This outline can be
conveniently plotted by \code{lines.perimeter}.
\code{Legend.plot.Design} invisibly returns the position of the legend.
\code{plot.Design} invisibly returns an invisible object of class \code{"plot.Design"}
with the following components (use \code{print.plot.Design}
to get a nice printout of the object):

\item{x.xbeta}{
data frame of values plotted. First column is sequence of \code{x}-values. If a
second variable was plotted, second column is sequence of \code{y}-values.
Next column is estimated \eqn{X\beta}{X beta}, followed by a column of lower confidence
limits and upper confidence limits. If fun was specified, these last three
columns are transformed by the specified function.
}
\item{adjust}{
character string of the form \code{"sex=male age=50"} containing settings of
non-plotted factors.
}
\item{curve.labels}{
character vector containing values of \code{y}-variable if it determined different
curves on the plot. E.g., \code{c("female","male")} or \code{c("10","20","30")}.
This vector is useful as an argument to the S \code{key} function if
\code{label.curves=FALSE}.
}
\item{plot.type}{
\code{"curves"} or \code{"3d"}
}
\item{method}{
from \code{plot.Design} call
}
\item{lty}{
vector of line types used for curves (if the plot used a few curves to
represent a second variable plotted)
}
\item{lwd}{
vector of line widths used
}
\item{col}{
vector of color codes
}}
\details{
When there are no intercepts in the fitted model, plot subtracts adjustment values from
each factor while computing variances for confidence limits.


\code{perimeter} is a kind of generalization of \code{datadist} for 2 continuous
variables.  First, the \code{n} smallest and largest \code{x} values are determined.
These form the lowest and highest possible \code{x}s to display.  Then \code{x}
is grouped into intervals bounded by these two numbers, with the interval
widths defined by \code{xinc}.  Within each interval, \code{y} is sorted and the
\eqn{n}th smallest and largest \code{y} are taken as the interval containing
sufficient data density to plot interaction surfaces.  The interval
is ignored when there are insufficient \code{y} values.  When \code{plot.Design}
readies the data for \code{persp}, it uses the \code{approx} function to do
linear interpolation of the \code{y}-boundaries as a function of the
\code{x} values actually used in forming the grid (the values of the
first variable specified to \code{plot}).  To make the perimeter smooth,
specify \code{lowess.=TRUE} to \code{perimeter}.


Specifying \code{time} will not work for Cox models with time-dependent
covariables.  Use \code{survest} or \code{survfit} for that purpose.


Use \code{ps.slide}, \code{win.slide}, \code{gs.slide} to set up nice defaults for
plotting.  These also set a system option \code{mgp.axis.labels} to allow x
and y-axes to have differing \code{mgp} graphical parameters (see \code{par}).
This is important when labels for y-axis tick marks are to be written
horizontally (\code{par(las=1)}), as a larger gap between the labels and
the tick marks are needed.  You can set the axis-specific 2nd
components of \code{mgp} using \code{mgp.axis.labels(c(xvalue,yvalue))}.

Note that because the generic \code{plot} method has the variable
\code{x} as its first argument, you cannot explicitly specify that you
want to plot the effect of a predictor named \code{x}.
}
\author{
Frank Harrell\cr
Division of Biostatistics and Epidemiology, University of Virginia\cr
fharrell@virginia.edu
}
\seealso{
\code{\link{datadist}}, \code{\link{predict.Design}}, \code{\link{contrast.Design}}, \code{\link{summary.Design}}, 
\code{\link{persp}}, \code{\link{Design}}, 
\code{\link{Design.trans}}, \code{\link{survest}}, \code{\link{survplot}}, \code{\link{Design.Misc}}, 
\code{\link{contour}}, \code{\link{image}}, \code{\link[Hmisc]{labcurve}}, \code{\link[Hmisc]{scat1d}}, \code{\link[Hmisc]{dotchart2}}, 
\code{\link[Hmisc]{mgp.axis.labels}} \code{\link[Hmisc]{Overview}}, \code{\link{par}},
\code{\link[Hmisc]{ps.slide}}, \code{\link[Hmisc]{xYplot}}, \code{\link[Hmisc]{smearingEst}}
}
\examples{
n <- 1000    # define sample size
set.seed(17) # so can reproduce the results
age            <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
cholesterol    <- rnorm(n, 200, 25)
sex            <- factor(sample(c('female','male'), n,TRUE))
label(age)            <- 'Age'      # label is in Hmisc
label(cholesterol)    <- 'Total Cholesterol'
label(blood.pressure) <- 'Systolic Blood Pressure'
label(sex)            <- 'Sex'
units(cholesterol)    <- 'mg/dl'   # uses units.default in Hmisc
units(blood.pressure) <- 'mmHg'


# Specify population model for log odds that Y=1
L <- .4*(sex=='male') + .045*(age-50) +
  (log(cholesterol - 10)-5.2)*(-2*(sex=='female') + 2*(sex=='male'))
# Simulate binary y to have Prob(y=1) = 1/[1+exp(-L)]
y <- ifelse(runif(n) < plogis(L), 1, 0)


ddist <- datadist(age, blood.pressure, cholesterol, sex)
options(datadist='ddist')


fit <- lrm(y ~ blood.pressure + sex * (age + rcs(cholesterol,4)),
               x=TRUE, y=TRUE)


par(mfrow=c(2,2))
plot(fit)                # Plot effects of all 4 predictors
par(mfrow=c(1,2))
plot(fit, name=c('age','cholesterol'))   # Make 2 plots
par(mfrow=c(1,1))
plot(fit, age=seq(20,80,length=100), sex=NA, conf.int=FALSE)
                         # Plot relationship between age and log
                         # odds, separate curve for each sex,
                         # no C.I.
z <- plot(fit, age=NA, sex=NA, label.curves=FALSE)
                         # use label.curves=list(keys=c('a','b'))'
                         # to use 1-letter abbreviations
datadensity(z, age, sex) # rug plots (1-dimensional scatterplots)
                         # on each treatment curve, with treatment-
                         # specific density of age
plot(fit, age=seq(20,80,length=100), sex='male')  # works if datadist not used
plot(fit, age=NA, cholesterol=NA)# 3-dimensional perspective plot for age,
                         # cholesterol, and log odds using default
                         # ranges for both variables
boundaries <- perimeter(age, cholesterol, lowess=TRUE)
plot(age, cholesterol)   # show bivariate data density
lines(boundaries)        # and perimeter that will be used for 3-D plot
z <- plot(fit, age=NA, cholesterol=NA, perim=boundaries, method='image')
                         # draws image() plot
                         # don't show estimates where data are sparse
                         # doesn't make sense here since vars don't interact
if(!.R.)Legend(z, fun=plogis, at=qlogis(c(.01,.05,.1,.2,.3,.4,.5)),
               zlab='Probability')   # gray scale or color legend for prob.
plot(fit, age=NA, fun=function(x) 1/(1+exp(-x)) , # or fun=plogis
     ylab="Prob", conf.int=.9)    # Plot estimated probabilities instead of
                                  # log odds


# Plot the age effect as an odds ratio
# comparing the age shown on the x-axis to age=30 years


ddist$limits$age[2] <- 30    # make 30 the reference value for age
# Could also do: ddist$limits["Adjust to","age"] <- 30
fit <- update(fit)   # make new reference value take effect
plot(fit, age=NA, ref.zero=TRUE, fun=exp, ylab='Age=x:Age=30 Odds Ratio')
abline(h=1, lty=2, col=2); abline(v=30, lty=2, col=2)


# Make two curves, and plot the predicted curves as two trellis panels
w <- plot(fit, age=NA, sex=NA)   # Would be nice if a pl=FALSE option was avail.
z <- data.frame(w$x.xbeta)     # Makes variable names legal
if(.R.) library(lattice)
xyplot(log.odds ~ age | sex, data=z, type='l')
# To add confidence bands we need to use the Hmisc xYplot function in
# place of xyplot
xYplot(Cbind(log.odds,lower,upper) ~ age | sex, data=z, 
       method='bands', type='l')
# If non-displayed variables were in the model, add a subtitle to show
# their settings using title(sub=paste('Adjusted to',w$adjust),adj=0)
# See predict.Design for an example using predict and xYplot without plot()




# Plots for a parametric survival model
n <- 1000
set.seed(731)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
sex <- factor(sample(c('Male','Female'), n, 
              rep=TRUE, prob=c(.6, .4)))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
t <- -log(runif(n))/h
label(t) <- 'Follow-up Time'
e <- ifelse(t<=cens,1,0)
t <- pmin(t, cens)
units(t) <- "Year"
ddist <- datadist(age, sex)
Srv <- Surv(t,e)


# Fit log-normal survival model and plot median survival time vs. age
f <- psm(Surv(t, e) ~ rcs(age), dist=if(.R.)'lognormal' else 'gaussian')
med <- Quantile(f)       # Creates function to compute quantiles
                         # (median by default)
plot(f, age=NA, fun=function(x)med(lp=x), ylab="Median Survival Time")
# Note: This works because med() expects the linear predictor (X*beta)
#       as an argument.  Would not work if use 
#       plot(\dots, ref.zero=TRUE or adj.zero=TRUE)
# Also, confidence intervals from this method are approximate since
# they don't take into account estimation of scale parameter


# Fit an ols model to log(y) and plot the relationship between x1
# and the predicted mean(y) on the original scale without assuming
# normality of residuals; use the smearing estimator
set.seed(1)
x1 <- runif(300)
x2 <- runif(300)
ddist <- datadist(x1,x2)
y  <- exp(x1+x2-1+rnorm(300))
f <- ols(log(y) ~ pol(x1,2)+x2)
r <- resid(f)
smean <- function(yhat)smearingEst(yhat, exp, res, statistic='mean')
formals(smean) <- list(yhat=numeric(0), res=r[!is.na(r)])
#smean$res <- r[!is.na(r)]   # define default res argument to function
plot(f, x1=NA, fun=smean, ylab='Predicted Mean on y-scale')


options(datadist=NULL)


\dontrun{
# Example in which separate curves are shown for 4 income values
# For each curve the estimated percentage of voters voting for
# the democratic party is plotted against the percent of voters
# who graduated from college.  scat1d is used to indicate
# the income-interval-specific data density for college.  For
# this purpose show the distribution of percent in college for
# those having an income level within +/- the half-width of
# the income interval.  scat1d shows the rug plot superimposed
# on the estimated curve.  Data are county-level percents.
# This can't be done automatically using datadensity on the object
# returned by plot.Design, as the variable representing different
# curves (income) is a continuous variable.


incomes <- seq(22900, 32800, length=4)  
# equally spaced to outer quintiles
pl <- plot(f, college=NA, income=incomes, 
           conf.int=FALSE, xlim=c(0,35), ylim=c(30,55),
           lty=1, lwd=c(.25,1.5,3.5,6), col=c(1,1,2,2))
graph.points <- pl$x.xbeta
for(i in 1:4) {
  college.in.income.group <- college[abs(income-incomes[i]) < 1650]
  this.income <- graph.points[,'income']==incomes[i]
  scat1d(college.in.income.group,
         curve=list(x=graph.points[this.income,'college'],
           y=graph.points[this.income,'democrat']))
}


# Instead of showing a rug plot on each curve, erase end portions
# of each curve where there are fewer than 10 counties having
# \% college graduates to the left of the x-coordinate being plotted,
# for the subset of counties having median family income with 1650
# of the target income for the curve


show.pts <- function(college.pts, income.pt) {
  s <- abs(income - income.pt) < 1650  #assumes income known to top frame
  x <- college[s]
  x <- sort(x[!is.na(x)])
  n <- length(x)
  low <- x[10]; high <- x[n-9]
  college.pts >= low & college.pts <= high
}


plot(f, college=NA, income=incomes,
     conf.int=FALSE, xlim=c(0,35), ylim=c(30,55),
     lty=1, lwd=c(.25,1.5,3.5,6), col=c(1,1,2,2),
     perim=show.pts)
}
}
\keyword{models}
\keyword{hplot}
\keyword{htest}
% Converted by Sd2Rd version 1.21.


