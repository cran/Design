\name{anova.Design}
\alias{anova.Design}
\alias{print.anova.Design}
\alias{text.anova.Design}
\alias{plot.anova.Design}
\alias{latex.anova.Design}
\title{
  Analysis of Variance (Wald and F Statistics)
}
\description{
  The \code{anova} function automatically tests most meaningful hypotheses in
  a design. For example, suppose that age and cholesterol are
  predictors, and that a general interaction is modeled using a
  restricted spline surface. \code{anova} prints Wald statistics (\eqn{F}
  statistics for an \code{ols} fit) for testing
  linearity of age, linearity of cholesterol, age effect (age + age by
  cholesterol interaction), cholesterol effect (cholesterol + age by
  cholesterol interaction), linearity of the age by cholesterol
  interaction (i.e., adequacy of the simple age * cholesterol 1
  d.f. product), linearity of the interaction in age alone, and
  linearity of the interaction in cholesterol alone. Joint tests of all
  interaction terms in the model and all nonlinear terms in the model
  are also performed.  For any multiple d.f. effects for continuous
  variables that were not modeled through \code{rcs}, \code{pol},
  \code{lsp}, etc., 
  tests of linearity will be omitted.  This applies to matrix predictors
  produced by e.g.  \code{poly} or \code{ns}.  \code{print.anova.Design} is the
  printing method.  \code{text.anova.Design} is the \code{text} method for
  inserting anova tables on graphs.  \code{plot.anova.Design} draws dot
  charts depicting the importance of variables in the model, as measured
  by Wald \eqn{\chi^2}{chi-square}, \eqn{\chi^2}{chi-square} minus d.f., AIC, \eqn{P}-values, partial \eqn{R^2},
  \eqn{R^2} for the whole model after deleting the effects in question, or
  proportion of overall model \eqn{R^2} that is due to each predictor. 
  \code{latex.anova.Design} is the \code{latex} method.  It substitutes
  Greek/math symbols in column headings, uses boldface for \samp{TOTAL}
  lines, and constructs a caption.  Then it passes the result to
  \code{latex.default} for conversion to LaTeX.
}
\usage{
\method{anova}{Design}(object, \ldots, main.effect=FALSE, tol=1e-9, 
      test=c('F','Chisq'), ss=TRUE)

\method{print}{anova.Design}(x, which=c('none','subscripts','names','dots'), \dots)

\method{plot}{anova.Design}(x, 
     what=c("chisqminusdf","chisq","aic","P","partial R2","remaining R2",
            "proportion R2"), 
     xlab=NULL, pch=16, 
     rm.totals=TRUE, rm.ia=FALSE, rm.other=NULL, newnames,
     sort=c("descending","ascending","none"), pl=TRUE, \dots)

\method{text}{anova.Design}(x, at, cex=.5, font=2, \dots)

\method{latex}{anova.Design}(object, title, psmall=TRUE, 
      dec.chisq=2, dec.F=2, dec.ss=NA, dec.ms=NA, dec.P=4, \dots)
}
\arguments{
  \item{object}{
    a \code{Design} fit object.  \code{object} must
    allow \code{Varcov} to return the variance-covariance matrix.  For
    \code{latex}, is the result of \code{anova}.
  }
  \item{\dots}{
    If omitted, all variables are tested, yielding tests for individual factors
    and for pooled effects. Specify a subset of the variables to obtain tests
    for only those factors, with a pooled Wald tests for the combined effects
    of all factors listed. Names may be abbreviated.  For example, specify
    \code{anova(fit,age,cholesterol)} to get a Wald statistic for testing the joint
    importance of age, cholesterol, and any factor interacting with them.

    Can be optional graphical parameters to send to \code{text} or
    \code{dotchart2}, or other parameters to send to \code{latex.default}.
    Ignored for \code{print}.
  }
  \item{main.effect}{
    Set to \code{TRUE} to print the (usually meaningless) main effect tests even when
    the factor is involved in an interaction. The default is \code{FALSE}, to print only
    the effect of the main effect combined with all interactions involving that
    factor.
  }
  \item{tol}{
    singularity criterion for use in matrix inversion
  }
  \item{test}{
    For an \code{ols} fit, set \code{test="Chisq"} to use Wald \eqn{\chi^2}{chi-square} tests rather than F-tests.
  }
  \item{ss}{
    For an \code{ols} fit, set \code{ss=FALSE} to suppress printing partial sums of squares, mean
    squares, and the Error SS and MS.
  }
  \item{x}{for \code{print},\code{plot},\code{text} is the result of \code{anova}.
  }
  \item{which}{
    If \code{which} is not \code{"none"} (the default), \code{print.anova.Design} will
    add to the rightmost column of the output the list of parameters being
    tested by the hypothesis being tested in the current row.  Specifying
    \code{which="subscripts"} causes the subscripts of the regression
    coefficients being tested to be printed (with a subscript of one for
    the first non-intercept term).  \code{which="names"} prints the names of
    the terms being tested, and \code{which="dots"} prints dots for terms being
    tested and blanks for those just being adjusted for.
  }
  \item{at}{
    for \code{text} is a list containing the x- and y-coordinates for the
    upper left corner of the anova table to be drawn on an existing plot,
    e.g. \code{at=locator(1)}
  }
  \item{cex}{
    character expansion size for \code{text.anova.Design}
  }
  \item{font}{
    font for \code{text.anova.Design}.  Default is \code{2} (usually Courier).
  }
  \item{what}{
    what type of statistic to plot.  The default is the Wald
    \eqn{\chi^2}{chi-square} 
    statistic for each factor (adding in the effect of higher-ordered
    factors containing that factor) minus its degrees of freedom.  The
    last three choice for \code{what} only apply to \code{ols} models.
  }
  \item{xlab}{
    x-axis label, default is constructed according to \code{what}.
    \code{plotmath} symbols are used for \R, by default.
  }
  \item{pch}{
    character for plotting dots in dot charts.  Default is \code{16} (solid dot).
  }
  \item{rm.totals}{
    set to \code{FALSE} to keep total \eqn{\chi^2}{chi-square}s (overall, nonlinear, interaction totals)
    in the chart.
  }
  \item{rm.ia}{
    set to \code{TRUE} to omit any effect that has \code{"*"} in its name
  }
  \item{rm.other}{
    a list of other predictor names to omit from the chart
  }
  \item{newnames}{
    a list of substitute predictor names to use, after omitting any.
  }
  \item{sort}{
    default is to sort bars in descending order of the summary statistic
  }
  \item{pl}{
    set to \code{FALSE} to suppress plotting.  This is useful when you only wish to
    analyze the vector of statistics returned.
  }
  \item{title}{
    title to pass to \code{latex}, default is name of fit object passed to \code{anova}
    prefixed with \code{"anova."}.  For Windows, the default is \code{"ano"} followed
    by the first 5 letters of the name of the fit object.
  }
  \item{psmall}{
    The default is \code{psmall=TRUE}, which causes \code{P<0.00005} to print as \code{<0.0001}.
    Set to \code{FALSE} to print as \code{0.0000}.
  }
  \item{dec.chisq}{
    number of places to the right of the decimal place for typesetting
    \eqn{\chi^2}{chi-square} values (default is \code{2}).  Use zero for integer, \code{NA} for
    floating point.
  }
  \item{dec.F}{
    digits to the right for \eqn{F} statistics (default is \code{2})
  }
  \item{dec.ss}{
    digits to the right for sums of squares (default is \code{NA}, indicating
    floating point)
  }
  \item{dec.ms}{
    digits to the right for mean squares (default is \code{NA})
  }
  \item{dec.P}{
    digits to the right for \eqn{P}-values
  }
}
\value{
  \code{anova.Design} returns a matrix of class \code{anova.Design} containing factors 
  as rows and \eqn{\chi^2}{chi-square}, d.f., and \eqn{P}-values as
  columns (or d.f., partial \eqn{SS, MS, F, P}).
  \code{plot.anova.Design} invisibly returns the vector of quantities
  plotted.  This vector has a names attribute describing the terms for
  which the statistics in the vector are calculated.
}
\details{
  If the statistics being plotted with \code{plot.anova.Design} are few in
  number and one of them is negative or zero, \code{plot.anova.Design}
  will quit because of an error in \code{dotchart2}.
}
\author{
  Frank Harrell\cr
  Department of Biostatistics, Vanderbilt University\cr
  \email{f.harrell@vanderbilt.edu}
}
\section{Side Effects}{
  \code{print} prints, \code{text} uses \code{tempfile} to get a temporary Unix file name,
  \code{sink}, and \code{unix} (to remove the temporary file).  \code{latex} creates a
  file with a name of the form \file{title.tex} (see the \code{title} argument above).
}
\seealso{
  \code{\link[Design]{Design}}, \code{\link[Design]{Design.Misc}}, \code{\link[Design:Design.Misc]{lrtest}}, \code{\link[Design]{Design.trans}}, \code{\link[Design]{summary.Design}}, \code{\link[Hmisc]{solvet}}, 
  \code{\link[graphics]{text}}, \code{\link[graphics]{locator}}, \code{\link[Hmisc]{dotchart2}}, \code{\link[Hmisc]{latex}}, 
  \code{\link[Hmisc:xYplot]{Dotplot}}, \code{\link[stats]{anova.lm}}, \code{\link[Design:contrast]{contrast.Design}}
}
\examples{
n <- 1000    # define sample size
set.seed(17) # so can reproduce the results
treat <- factor(sample(c('a','b','c'), n,TRUE))
num.diseases <- sample(0:4, n,TRUE)
age <- rnorm(n, 50, 10)
cholesterol <- rnorm(n, 200, 25)
weight <- rnorm(n, 150, 20)
sex <- factor(sample(c('female','male'), n,TRUE))
label(age) <- 'Age'      # label is in Hmisc
label(num.diseases) <- 'Number of Comorbid Diseases'
label(cholesterol) <- 'Total Cholesterol'
label(weight) <- 'Weight, lbs.'
label(sex) <- 'Sex'
units(cholesterol) <- 'mg/dl'   # uses units.default in Hmisc


# Specify population model for log odds that Y=1
L <- .1*(num.diseases-2) + .045*(age-50) +
     (log(cholesterol - 10)-5.2)*(-2*(treat=='a') +
     3.5*(treat=='b')+2*(treat=='c'))
# Simulate binary y to have Prob(y=1) = 1/[1+exp(-L)]
y <- ifelse(runif(n) < plogis(L), 1, 0)


fit <- lrm(y ~ treat + scored(num.diseases) + rcs(age) +
               log(cholesterol+10) + treat:log(cholesterol+10))
anova(fit)                       # Test all factors
anova(fit, treat, cholesterol)   # Test these 2 by themselves
                                 # to get their pooled effects
g <- lrm(y ~ treat*rcs(age))
dd <- datadist(treat, num.diseases, age, cholesterol)
options(datadist='dd')
plot(g, age=NA, treat="b")
s <- anova(g)
print(s)
#p <- locator(1)                  # click mouse at upper left corner of table
p <- list(x=32,y=2.1)
text(s, at=p)                    # add anova table to regression plot
plot(s)                          # new plot - dot chart of chisq-d.f.
# latex(s)                       # nice printout - creates anova.g.tex
options(datdist=NULL)



# Simulate data with from a given model, and display exactly which
# hypotheses are being tested


set.seed(123)
age <- rnorm(500, 50, 15)
treat <- factor(sample(c('a','b','c'), 500,TRUE))
bp  <- rnorm(500, 120, 10)
y   <- ifelse(treat=='a', (age-50)*.05, abs(age-50)*.08) + 3*(treat=='c') +
       pmax(bp, 100)*.09 + rnorm(500)
f   <- ols(y ~ treat*lsp(age,50) + rcs(bp,4))
print(names(coef(f)), quote=FALSE)
specs(f)
anova(f)
an <- anova(f)
options(digits=3)
print(an, 'subscripts')
print(an, 'dots')


an <- anova(f, test='Chisq', ss=FALSE)
plot(0:1)                        # make some plot
text(an, at=list(x=1.5,y=.6))    # add anova table to plot
plot(an)                         # new plot - dot chart of chisq-d.f.
# latex(an)                      # nice printout - creates anova.f.tex


# Suppose that a researcher wants to make a big deal about a variable 
# because it has the highest adjusted chi-square.  We use the
# bootstrap to derive 0.95 confidence intervals for the ranks of all
# the effects in the model.  We use the plot method for anova, with
# pl=FALSE to suppress actual plotting of chi-square - d.f. for each
# bootstrap repetition.  We rank the negative of the adjusted
# chi-squares so that a rank of 1 is assigned to the highest.
# It is important to tell plot.anova.Design not to sort the results,
# or every bootstrap replication would have ranks of 1,2,3 for the stats.


mydata <- data.frame(x1=runif(200), x2=runif(200),
                     sex=factor(sample(c('female','male'),200,TRUE)))
set.seed(9)  # so can reproduce example
mydata$y <- ifelse(runif(200)<=plogis(mydata$x1-.5 + .5*(mydata$x2-.5) + 
                   .5*(mydata$sex=='male')),1,0)


if(.R.) {
library(boot)
b <- boot(mydata, function(data, i, ...) rank(-plot(anova(
                lrm(y ~ rcs(x1,4)+pol(x2,2)+sex,data,subset=i)), 
                sort='none', pl=FALSE)),
                R=25)  # should really do R=500 but will take a while
Rank <- b$t0
lim <- t(apply(b$t, 2, quantile, probs=c(.025,.975)))
} else {
b <- bootstrap(mydata, rank(-plot(anova(
                lrm(y ~ rcs(x1,4)+pol(x2,2)+sex,mydata)), sort='none', pl=FALSE)),
               B=25)  # should really do B=500 but will take a while
Rank <- b$observed
lim <- limits.emp(b)[,c(1,4)]  # get 0.025 and 0.975 quantiles
}


# Use the Hmisc Dotplot function to display ranks and their confidence
# intervals.  Sort the categories by descending adj. chi-square, for ranks
original.chisq <- plot(anova(lrm(y ~ rcs(x1,4)+pol(x2,2)+sex,data=mydata)),
                       sort='none', pl=FALSE)
predictor <- as.factor(names(original.chisq))
predictor <- reorder.factor(predictor, -original.chisq)

Dotplot(predictor ~ Cbind(Rank, lim), pch=3, xlab='Rank', 
		main=if(.R.) expression(paste(
'Ranks and 0.95 Confidence Limits for ',chi^2,' - d.f.')) else
'Ranks and 0.95 Confidence Limits for Chi-square - d.f.')
}
\keyword{models}
\keyword{regression}
\keyword{htest}
\keyword{aplot}
\concept{bootstrap}


