#SCCS  @(#)survfit.coxph.null.s	4.5 % G%
#Set score to 0 if linear predictors not there - FEH 28 Sep 93
survfit.cph.null <-
  function(object, newdata, se.fit=TRUE, conf.int=.95, individual=FALSE,
	    type, vartype,
	    conf.type=c('log-log', 'log', 'plain', 'none'), ...) {

    if(.R.) require('survival')
    ## Sense whether survival5 is in effect and if so use this later version
    s5 <- exists('coxpenal.fit')

    # May have strata and/or offset terms, linear predictor = offset
    #  newdata doesn't make any sense
    #  This is survfit.coxph with lots of lines removed
    y <- object$y
    if(is.null(y))stop("must use y=T with fit")
    n <- nrow(y)
    Strata <- Surv.strata(y)
    if(length(object$linear.predictor)==0) score <- rep(1, n)
    else score <- exp(object$linear.predictor)
    
    temp <- c("aalen", "kalbfleisch-prentice", "efron", "tsiatis", "breslow", 
              "kaplan-meier", "fleming-harringon", "greenwood", "exact")
    temp2 <- c(2, 1, 3, 2, 2, 1, 3, 1, 1)
    if(missing(type)) type <- object$method
    if(missing(vartype)) vartype <- type
    method <- temp2[match(match.arg(type, temp), temp)]
    if(is.na(method)) stop("Invalid survival curve type")
    vartype <- temp2[match(match.arg(vartype, temp), temp)]
    if(is.na(vartype)) stop("Invalid variance type specified")
    km <- method==1
    
    if (!se.fit) conf.type <- 'none'
    else conf.type <- match.arg(conf.type)

    ny <- ncol(y)
    type <- attr(y, 'type')
    if (type=='counting') {
      ord <- order(Strata, y[,2], -y[,3])
      ##if (method=='kaplan-meier')   bug correction FEH 6Jun99
      if(km) stop ("KM method not valid for counting type data")
	}
    else if (type=='right') {
      ord <- order(Strata, y[,1], -y[,2])
      miny <- min(y[, 1])
      if(miny < 0) y <- cbind(2 * miny - 1, y)
      else y <- cbind(-1, y)
	}
    else stop("Cannot handle \"", type, "\" type survival data")

    if (!is.null(Strata)) {
	newstrat <- (as.numeric(Strata))[ord]
	newstrat <- as.integer(c(1*(diff(newstrat)!=0), 1))
	}
    else newstrat <- as.integer(c(rep(0,n-1),1))

    if ( !missing(newdata))
	stop("A newdata argument does not make sense for a null model")
    dimnames(y) <- NULL   #I only use part of Y, so names become invalid
    storage.mode(y) <- 'double'
    surv <- if(.R.)
      .C('agsurv2',
         as.integer(n),
         as.integer(0),
         y = y[ord,],
         as.double(score[ord]),
         strata = as.integer(newstrat),
         surv = double(n),
         varhaz = double(n),
         double(1),
         as.double(0),
         nsurv = as.integer(c(method,vartype)),
         double(2),
         as.integer(1),
         double(1),
         newrisk= as.double(1), PACKAGE="survival") else
    if(.SV4. && (version$major > 6 || (version$major == 6 &&
                                       version$minor >= 2))) {
      weights <- if(length(object$weights)) object$weights[ord] else rep(1,n)
           .C('S_agsurv2',
              as.integer(n),
              as.integer(0),
              y = y[ord,],
              as.double(score[ord]),
              strata = as.integer(newstrat),
              wt = as.double(weights),
              surv = double(n),
              varhaz = double(n),
              double(1),
              as.double(0),
              nsurv = as.integer(if(s5) c(method,vartype) else km),
              double(2),
              as.integer(1),
              double(1),
              newrisk= as.double(1))
         } else
    .C(if(.SV4.)'S_agsurv2' else 'agsurv2', ##14Nov00
       as.integer(n),
       as.integer(0),
       y = y[ord,],
       as.double(score[ord]),
       strata = as.integer(newstrat),
       surv = double(n),
       varhaz = double(n),
       double(1),
       as.double(0),  ##14Nov00
       nsurv = as.integer(if(s5) c(method,vartype) else km),
       double(2),
       as.integer(1),
       double(1),
       newrisk= as.double(1))

    nsurv <- surv$nsurv[1]
    ntime <- 1:nsurv
    tsurv <- surv$surv[ntime]
    tvar  <- surv$varhaz[ntime]
    if (surv$strata[1] <=1)
	temp <- list(time=surv$y[ntime,1],
		 n.risk=surv$y[ntime,2],
		 n.event=surv$y[ntime,3],
		 surv=tsurv )
    else {
	temp <- surv$strata[1:(1+surv$strata[1])]
	tstrat <- diff(c(0, temp[-1])) #n in each strata
	names(tstrat) <- levels(Strata)
	temp <- list(time=surv$y[ntime,1],
		 n.risk=surv$y[ntime,2],
		 n.event=surv$y[ntime,3],
		 surv=tsurv,
		 strata= tstrat)
	}
    if (se.fit) temp$std.err <- sqrt(tvar)

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
#    if(!is.null(strata)) attr(temp, "strata") <- Strata
    oldClass(temp) <- if(.SV4.)'survfit.cox' else
     c("survfit.cph", "survfit.cox", "survfit")
    temp
    }
