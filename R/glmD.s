glmD <- if(.R.)
  function(formula, family = gaussian, data = list(), weights = NULL,
    subset = NULL, na.action = na.fail, start = NULL, offset = NULL,
    control = glm.control(...), model = TRUE, method = "glm.fit",
    x = FALSE, y = TRUE, contrasts = NULL, ...)
{
    call <- match.call()
    if (is.character(family))
        family <- get(family)
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("`family' not recognized")
    }
    mt <- terms(formula, data = data)
    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    mf$family <- mf$start <- mf$control <- mf$maxit <- NULL
    mf$model <- mf$method <- mf$x <- mf$y <- mf$contrasts <- NULL
    mf$... <- NULL
    mf$drop.unused.levels <- TRUE   # FEH 31jul02
    mf[[1]] <- as.name("model.frame")

    dul <- .Options$drop.unused.levels  # FEH 31jul02
    if(!length(dul) || dul) {
      on.exit(options(drop.unused.levels=dul))
      options(drop.unused.levels=FALSE)
    }

    mf <- Design(eval(mf, parent.frame()))   # FEH 13Apr01
    desatr <- attr(mf,'Design')
    attr(mf,'Design') <- NULL
    
    switch(method, model.frame = return(mf), glm.fit = 1, glm.fit.null = 1,
        stop(paste("invalid `method':", method)))
    xvars <- as.character(attr(mt, "variables"))[-1]
    if ((yvar <- attr(mt, "response")) > 0)
        xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(mf[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts)
    Y <- model.response(mf, "numeric")
    weights <- model.weights(mf)
    offset <- model.offset(mf)
    if (!is.null(weights) && any(weights < 0))
        stop("Negative wts not allowed")
    if (!is.null(offset) && length(offset) != NROW(Y))
        stop(paste("Number of offsets is", length(offset), ", should equal",
            NROW(Y), "(number of observations)"))
    fit <- (if (is.empty.model(mt))
        glm.fit.null
    else glm.fit)(x = X, y = Y, weights = weights, start = start,
        offset = offset, family = family, control = control,
        intercept = attr(mt, "intercept") > 0)
    if (any(offset) && attr(mt, "intercept") > 0) {
        fit$null.deviance <- if (is.empty.model(mt))
            fit$deviance
        else glm.fit(x = X[, "(Intercept)", drop = FALSE], y = Y,
            weights = weights, start = start, offset = offset,
            family = family, control = control, intercept = TRUE)$deviance
    }
    if (model)
        fit$model <- mf
    if (x)
        fit$x <- X
    if (!y)
        fit$y <- NULL
    fit <- c(fit, list(call = call, formula = formula, terms = mt,
        data = data, offset = offset, control = control, method = method,
        contrasts = attr(X, "contrasts"), xlevels = xlev,
                       Design=desatr,
                       assign=DesignAssign(desatr,1,mt)))
    ##FEH 13Apr01 24nov02 above
    class(fit) <- c("Design", 'glmD',
                    if (is.empty.model(mt)) "glm.null", "glm",
        "lm") # FEH 13Apr01  glmD 26nov02
    fit
} else function(formula = formula(data), family = gaussian,
                data = sys.parent(),
        weights, subset, na.action, start = eta, control = glm.control(...),
        method = "glm.fit", model = FALSE, x = FALSE, y = TRUE, contrasts = NULL, ...)
{
        call <- match.call()
        m <- match.call(expand = FALSE)
        m$family <- m$method <- m$model <- m$x <- m$y <- m$control <- m$
                contrasts <- m$... <- NULL
        m$drop.unused.levels <- TRUE
        m[[1]] <- as.name("model.frame")

        m <- Design(eval(m, sys.parent()))  # FEH 13Apr01
        desatr <- attr(m,'Design')
        attr(m,'Design') <- NULL
        
        Terms <- attr(m, "terms")
        if(method == "model.frame")
                return(m)
        xvars <- as.character(attr(Terms, "variables"))
        if(length(xvars) > 0) {
                xlevels <- lapply(m[xvars], levels)
                xlevels <- xlevels[!sapply(xlevels, is.null)]
                if(length(xlevels) == 0)
                        xlevels <- NULL
        }
        else xlevels <- NULL
        a <- attributes(m)
        Y <- model.extract(m, response)
        X <- model.matrix(Terms, m, contrasts)
        w <- model.extract(m, weights)
        if(!length(w))
                w <- rep(1, nrow(m))
        else if(any(w < 0))
                stop("negative weights not allowed")
        start <- model.extract(m, start)
        offset <- model.extract(m, offset)
        family <- as.family(family)
        if(missing(method))
                method <- attr(family, "method")
        if(!is.null(method)) {
                if(!existsFunction(method))
                        stop(paste("unimplemented method:", method))
        }
        else method <- "glm.fit"
        glm.fitter <- get(method)
        fit <- glm.fitter(x = X, y = Y, w = w, start = start, offset = offset,
                family = family, maxit = control$maxit, epsilon = control$
                epsilon, trace = control$trace, null.dev = TRUE, ...)
        #
        # If an offset and intercept is present, iterations are needed to
        # compute the Null deviance; these are done here, unless the model
        # is NULL, in which case the computations have been done already
        #
        if(any(offset) && attr(Terms, "intercept")) {
                null.deviance <- if(length(Terms)) glm.fitter(X[, "(Intercept)",                                drop = FALSE], Y, w, offset = offset, family =
                                family, maxit = control$maxit, epsilon =
                                control$epsilon, null.dev = NULL)$deviance
                         else fit$deviance
                fit$null.deviance <- null.deviance
        }
        oldClass(fit) <- if(.SV4.) 'Design' else
            c("Design","glmD","glm","lm")  # FEH 13Apr01 16aug02
        ## glmD 2dec02 8p
        if(!is.null(xlevels))
                attr(fit, "xlevels") <- xlevels
        fit$terms <- Terms
        fit$formula <- as.vector(attr(Terms, "formula"))
        fit$call <- call
        fit$Design <- desatr   # FEH 13Apr01
        fit$assign <- DesignAssign(desatr,1,Terms) ## 24nov02
        if(model)
                fit$model <- m
        if(x)
                fit$x <- X
        if(!y)
                fit$y <- NULL
        fit$control <- control
        if(!is.null(attr(m, "na.action")))
                fit$na.action <- attr(m, "na.action")
        fit$fitFunction <- c('glmD','glm','lm')  ## glmD 26nov02
        fit
}

## 26nov02
print.glmD <- function(x, digits=4, ...) {
  cat('General Linear Model\n\n')
  dput(x$call); cat('\n\n')
  cof <- coef(x)
  lr <- x$null.deviance - x$deviance
  names(cof) <- ifelse(names(cof)=='(Intercept)','Intercept',names(cof))
  dof <- x$rank - (names(cof)[1]=='Intercept')
  pval <- 1 - pchisq(lr, dof)
  print(c('Model L.R.'=format(lr,digits=2), 'd.f.'=format(dof),
          'P'=format(pval,digits=4)), quote=FALSE)
  cat('\n')
  se <- sqrt(diag(Varcov(x)))
  z <- cof/se
  p <- 1 - pchisq(z^2, 1)
  w <- cbind(format(cof, digits=digits),
             format(se,  digits=digits),
             format(z,   digits=2),
             format(p,   digits=4))
  dimnames(w) <- list(names(cof), c('Coef','S.E.','Wald Z','P'))
  print(w, quote=FALSE)
  invisible()
}

## 26nov02
summary.glmD <- function(...) summary.Design(...)
## 2dec02
Varcov.glmD <- function(object, regcoef.only=FALSE, ...)
  Varcov.glm(object, regcoef.only, ...)
## 6dec02
predict.glmD <- 
  function(object, newdata,
           type=c("lp","x","data.frame","terms","adjto","adjto.data.frame",
             "model.frame"),
           se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
           incl.non.slopes, non.slopes, kint=1,
           na.action=na.keep, expand.na=TRUE, center.terms=TRUE, ...)
  predictDesign(object, newdata, type, se.fit, conf.int, conf.type,
                incl.non.slopes, non.slopes, kint,
                na.action, expand.na, center.terms, ...)


latex.glmD <- function(...) latexDesign(...)
