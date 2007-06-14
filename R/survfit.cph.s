#SCCS @(#)survfit.coxph.s	4.5 7/14/92
survfit.cph <- function(object, newdata, se.fit=TRUE, conf.int=.95, 
                        individual=FALSE,
                        type, vartype,
                        conf.type=c('log-log', 'log', 'plain', 'none')) {
  if(.R.) {
    require('survival')
    coxpenal.fit <- getFromNamespace("coxpenal.fit", "survival")
  }
  
  ## Sense whether survival5 is in effect and if so use this later version
  s5 <- exists('coxpenal.fit')

  if(.R. || version$major < 6 || (version$major == 6 && version$minor < 2))
    if(length(object$weights) || length(object$call$weights)) 
      stop('survfit cannot yet handle a weighted model')

  nvar <- length(object$coef)
  score <- exp(object$linear.predictors)
  temp <- c("aalen", "kalbfleisch-prentice", "efron", "tsiatis", "breslow", 
            "kaplan-meier", "fleming-harringon", "greenwood", "exact")
  temp2 <- c(2, 1, 3, 2, 2, 1, 3, 1, 1)

  if(missing(type))
    type <- object$method

  if(missing(vartype))
    vartype <- type

  method <- temp2[match(match.arg(type, temp), temp)]

  if(is.na(method))
    stop("Invalid survival curve type")

  vartype <- temp2[match(match.arg(vartype, temp), temp)]

  if(is.na(vartype))
    stop("Invalid variance type specified")

  km <- method==1
  coxmethod <- object$method

  if (!se.fit)
    conf.type <- 'none'
  else
    conf.type <- match.arg(conf.type)

  x <- object$x
  y <- object$y

  if(is.null(x) | is.null(y))
    stop("you must specify x=T and y=T in the fit")

  n <- nrow(y)
  stratum <- attr(x, "strata")
  offset <- attr(x, "offset")

  if(length(offset)==0)
    offset <- rep(0,n)

  type <- attr(y, 'type')

  if (type=='counting') {
    if(km)
      stop ("KM method not valid for counting type data")

    ## if (method=='kaplan-meier')  bug correction FEH 6Jun99
    ord <- if(length(stratum))
      order(stratum, y[,2], -y[,3])
    else
      order(y[,2], -y[,3])
  }
  else if (type=='right') {
    ord <- if(length(stratum))
      order(stratum, y[, 1],  - y[, 2])
    else
      order(y[, 1],  - y[, 2])

    ## length(stratum) was length(strat) 14Nov00
    miny <- min(y[, 1])
    if(miny < 0)
      y <- cbind(2 * miny - 1, y)
    else
      y <- cbind(-1, y)
  }
  else
    stop("Cannot handle \"", type, "\" type survival data")

  x <- x[ord,]
  weights <- if(length(object$weights))
    object$weights[ord]
  else
    rep(1,n)

  if (length(stratum)) {
    newstrat <- (as.numeric(stratum))[ord]
    newstrat <- as.integer(c(1*(diff(newstrat)!=0), 1))
  }
  else
    newstrat <- as.integer(c(rep(0,n-1),1))

  if (individual && !missing(newdata))
    stype <- 1
  else
    stype <- 2
  
  if(stype == 1 && method != vartype)
    stop("The type and vartype args must agree for individual=T")

  if(stype == 1 && method == 1)
    stop("Only Aalen and F-H estimates available for individual=T")
  
  if (!missing(newdata)) {
    x2 <- predictDesign(object,newdata,type="x",expand.na=FALSE)
    n2 <- nrow(x2)
    offset2 <- attr(x2, "offset")
    if(length(offset2)==0)
      offset2 <- rep(0, n2)
    
    strata2 <- attr(x2,"strata")
    if(length(strata2)==0)
      strata2 <- rep(1, n2)
    
    naa <- attr(x2, "na.action")
    if(stype==1) {
      Terms <- terms(object$terms)
      m2 <- model.frame(Terms,newdata)
      y2 <- model.extract(m2, 'response')
      if(length(y2)==0)
        stop('no Surv object in newdata')

      if(attr(y2,'type')!=type)
        stop('Survival type of newdata does not match fitted model')

      if(nrow(y2) != n2)
        stop('wrong # rows in Y')
    }
  }
  else {
    x2 <- matrix(object$means, nrow=1)
    offset2 <- 0
    strata2 <- 0
    n2 <- 1
  }

  coef <- ifelse(is.na(object$coef), 0, object$coef)
  newrisk <- exp(c(x2 %*% coef) + offset2 - object$center)

  dimnames(y) <- NULL   #I only use part of Y, so names become invalid
  storage.mode(x) <- 'double'
  storage.mode(y) <- 'double'
  ndead <- sum(y[,3])
  
  if (stype==1) {
    surv <- if(.R.)
      .C("agsurv1",
         as.integer(n),
         as.integer(nvar),
         y[ord,],
         as.double(score[ord]),
         strata=as.integer(newstrat),
         surv=double(ndead*n2),
         varhaz=double(ndead*n2),
         nsurv=if(s5) as.integer(method==3)
               else as.integer(2+1*(coxmethod=="efron")),
         x,
         double(3*nvar),
         as.double(object$var),
         y = double(3*n*n2),
         as.integer(n2),
         as.double(y2),
         as.double(x2),
         as.double(newrisk),
         as.integer(strata2), PACKAGE="survival")
    else
      .C(if(.SV4.) 'S_agsurv1'
         else "agsurv1",
         as.integer(n),
         as.integer(nvar),
         y[ord,],
         as.double(score[ord]),   ## was a bug in survival4 as.double14Nov00
         strata=as.integer(newstrat),
         surv=double(ndead*n2),   # was bug in surv4
         varhaz=double(ndead*n2), # was bug in surv4
         nsurv=if(s5)as.integer(method==3)
               else as.integer(2+1*(coxmethod=="efron")),
         x,
         double(3*nvar),
         as.double(object$var),
         y = double(3*n*n2),
         as.integer(n2),
         as.double(y2),
         as.double(x2),
         as.double(newrisk),
         as.integer(strata2))

    ntime <- 1:surv$nsurv
    temp <- (matrix(surv$y, ncol=3))[ntime,]
    temp <- list(time = temp[,1],
                 n.risk= temp[,2],
                 n.event=temp[,3],
                 surv = surv$surv[ntime])

    if (se.fit)
      temp$std.err <- sqrt(surv$varhaz[ntime])
  }
  else {
    if(!s5)
      temp <- ifelse(km, 1, 2+as.integer(coxmethod=="efron"))

    ## 8feb04
    if(.SV4. && (version$major > 6 || (version$major == 6 &&
                                       version$minor >= 2)))
      surv <- .C('S_agsurv2',
                 as.integer(n),
                 as.integer(nvar* se.fit),
                 y = y[ord,],
                 as.double(score[ord]),
                 strata = as.integer(newstrat),
                 wt = as.double(weights),
                 surv   = double(ndead*n2),
                 varhaz = double(ndead*n2),
                 x,
                 as.double(object$var),
                 nsurv = if(s5|| .R.) as.integer(c(method,vartype))
                         else as.integer(temp),
                 double(3*nvar),
                 as.integer(n2),
                 as.double(x2),
                 as.double(newrisk))
    else
      surv <- .C(if(.SV4.) 'S_agsurv2'
                 else 'agsurv2',  ##14Nov00
                 as.integer(n),
                 as.integer(nvar* se.fit),
                 y = y[ord,],
                 as.double(score[ord]),
                 strata = as.integer(newstrat),
                 surv   = double(ndead*n2),   # was bug in surv4
                 varhaz = double(ndead*n2),   # was bug in surv4
                 x,
                 as.double(object$var),
                 nsurv = if(s5) as.integer(c(method,vartype))
                         else as.integer(temp),
                 double(3*nvar),
                 as.integer(n2),
                 as.double(x2),
                 as.double(newrisk))

    nsurv <- surv$nsurv[1]
    ntime <- 1:nsurv
    if (n2>1) {
      tsurv <- matrix(surv$surv[1:(nsurv*n2)], ncol=n2)
      tvar  <- matrix(surv$varhaz[1:(nsurv*n2)], ncol=n2)
      dimnames(tsurv) <- list(NULL, dimnames(x2)[[1]])
    }
    else {
      tsurv <- surv$surv[ntime]
      tvar  <- surv$varhaz[ntime]
    }
    if (surv$strata[1] <=1)
      temp <- list(time=surv$y[ntime,1],
                   n.risk=surv$y[ntime,2],
                   n.event=surv$y[ntime,3],
                   surv=tsurv)
    else {
      temp <- surv$strata[1:(1+surv$strata[1])]
      tstrat <- diff(c(0, temp[-1])) #n in each strata
      names(tstrat) <- levels(stratum)
      temp <- list(time=surv$y[ntime,1],
                   n.risk=surv$y[ntime,2],
                   n.event=surv$y[ntime,3],
                   surv=tsurv,
                   strata= tstrat)
    }

    if (se.fit)
      temp$std.err <- sqrt(tvar)
  }

  zval <- qnorm(1- (1-conf.int)/2, 0,1)
  if (conf.type=='plain') {
    temp1 <- temp$surv + zval* temp$std * temp$surv
    temp2 <- temp$surv - zval* temp$std * temp$surv
    temp <- c(temp, list(upper=pmin(temp1,1), lower=pmax(temp2,0),
                         conf.type='plain', conf.int=conf.int))
  }
  if (conf.type=='log') {
    xx <- ifelse(temp$surv==0,1,temp$surv)  #avoid some "log(0)" messages
    temp1 <- ifelse(temp$surv==0, 0*temp$std, exp(logb(xx) + zval* temp$std))
    temp2 <- ifelse(temp$surv==0, 0*temp$std, exp(logb(xx) - zval* temp$std))
    temp <- c(temp, list(upper=pmin(temp1,1), lower=temp2,
                         conf.type='log', conf.int=conf.int))
  }
  if (conf.type=='log-log') {
    who <- (temp$surv==0 | temp$surv==1) #special cases
    xx <- ifelse(who, .1,temp$surv)  #avoid some "log(0)" messages
    temp1 <- exp(-exp(logb(-logb(xx)) + zval*temp$std/logb(xx)))
    temp1 <- ifelse(who, temp$surv + 0*temp$std, temp1)
    temp2 <- exp(-exp(logb(-logb(xx)) - zval*temp$std/logb(xx)))
    temp2 <- ifelse(who, temp$surv + 0*temp$std, temp2)
    temp <- c(temp, list(upper=temp1, lower=temp2,
                         conf.type='log-log', conf.int=conf.int))
  }
  
  temp$call <- call
  if(!missing(newdata)) {
    temp$requested.strata <- strata2
    attr(temp, "na.action") <- naa
  }

  oldClass(temp) <- if(.SV4.)
    'survfit.cox'
  else
    c("survfit.cph", "survfit.cox", "survfit")  ##14Nov00

  temp
}
