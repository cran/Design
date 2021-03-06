# $Date: 2006-08-28 18:12:46 $ $Id: survfit.formula.S 11207 2009-02-09 04:25:50Z therneau $
survfit.formula <- function(formula, data, weights, subset, 
			    na.action, etype, id,
                            conf.type=c("log-log", "log", "plain", "none"), ...)
{
  Call <- match.call()
  Call[[1]] <- as.name('survfit')       #make nicer printout for the user
                                        # create a copy of the call that has only the arguments we want,
                                        #  and use it to call model.frame()
  mfnames <- c('formula', 'data', 'weights', 'subset','na.action',
               'etype', 'id')  #legal args for model.frame
  temp <- Call[c(1, match(mfnames, names(Call), nomatch=0))]
  temp[[1]] <- as.name("model.frame")
  if (is.R()) m <- eval.parent(temp)
  else        m <- eval(temp, sys.parent())
  
  Terms <- terms(formula, 'strata')
  ord <- attr(Terms, 'order')
  if (length(ord) & any(ord !=1))
    stop("Interaction terms are not valid for this function")

  n <- nrow(m)
  Y <- model.extract(m, 'response')
  if (!is.Surv(Y)) stop("Response must be a survival object")

  casewt <- model.extract(m, "weights")
  if (is.null(casewt)) casewt <- rep(1,n)

  if (!is.null(attr(Terms, 'offset'))) warning("Offset term ignored")

  ll <- attr(Terms, 'term.labels')
  if (length(ll) == 0) X <- factor(rep(1,n))  # ~1 on the right
  else X <- interaction(m[ll], drop=TRUE, sep=' ')
  
  etype <- model.extract(m, 'etype')
  id    <- model.extract(m, 'id')

  if (length(etype) >0)
    temp <- survival:::survfitCI(X, Y, weights=casewt, etype=etype, id=id,  ...)
  else if (attr(Y, 'type') != 'right' && attr(Y, 'type') != 'counting')
    temp <- survival:::survfitTurnbull(X, Y, casewt, ...)
  else
    temp <- survival:::survfitKM(X, Y, casewt, ...)

  if (is.R()) class(temp) <- 'survfit'
  else        oldClass(temp) <- "survfit"
  if (!is.null(attr(m, 'na.action')))
    temp$na.action <- attr(m, 'na.action')

  temp$maxtime <- max(Y[, ncol(Y)-1])
  temp$units <- valueUnit(Y)
  temp$time.label <- attr(Y, "time.label")

  temp$call <- Call
  temp
}
