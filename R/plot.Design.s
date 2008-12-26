## $Id: plot.Design.s 125 2008-12-26 20:47:05Z dupontct $
## First variable in ... cannot be named x - S methods try to call plot.default

plot.Design <-
  function(x, ..., xlim, ylim, fun, xlab, ylab,
           conf.int=.95, conf.type=c('mean','individual'), add=FALSE,
           label.curves=TRUE, eye, theta=0, phi=15, perspArgs=NULL, lty, col=1,
           lwd=par("lwd"), lwd.conf=1, pch=1, adj.zero=FALSE, ref.zero=FALSE,
           adj.subtitle, cex.adj, non.slopes, time=NULL, loglog=FALSE,
           val.lev=FALSE, digits=4, log="", perim,
           method=c("persp","contour","image", "dotchart","default"),
           sortdot=c('neither','ascending','descending'), nlevels=10, name,
           zlim=range(zmat,na.rm=TRUE),
           vnames=c('labels','names'), abbrev=FALSE, forceLines=FALSE)
{

  fit <- x
  conf.type <- match.arg(conf.type)
  vnames <- match.arg(vnames)
  
  dotlist <- list(...)
  fname  <- if(missing(name)) '' else name
  at     <- fit$Design
  if(!length(at)) at <- getOldDesign(fit)
  assume <- at$assume.code


  name <- at$name	##interactions are placed at end by design

  if(any(name == 'time'))
    {
      dotlist$time <- time
      time <- NULL
    }

  if(length(fname)>1 || (length(dotlist)==0 && fname==''))
    {
      m <- match.call(expand=FALSE)
      m[[1]] <- as.name('plot.Design')
      for(nam in if(length(fname)>1)fname else name[assume!=9])
        {
          m$name <- nam
          if(vnames=='names') m$xlab <- nam
          lv <- eval(m)
        }
      return(invisible(lv))
    }

  oldopt <- options(digits=digits)
  on.exit(options(oldopt))

  method <- match.arg(method)
  sortdot <- match.arg(sortdot)

  cex <- par('cex')
  if(missing(cex.adj)) cex.adj <- .75*cex

  Pretty <- function(x)
    { ##handles chron objects
      if(inherits(x,"dates") | inherits(x,"times"))
        structure(pretty(oldUnclass(x)), class=oldClass(x))
      else pretty(x)
    }

  f <- sum(assume!=9)	##limit list to main effects factors
  parms  <- at$parms
  label  <- at$label
  values <- at$values
  yunits <- fit$units
  units  <- at$units
  scale  <- fit$scale.pred
  if(!length(scale)) scale <- "X * Beta"

  Center <- fit$center
  if(length(Center)==0) Center <- 0

  if(missing(ylab))
    {
      if(!missing(fun))
        ylab <- ""
      else if(length(time))
        {
          if(loglog)
            ylab <- paste("log[-log S(",format(time),")]",sep="")
          else
            ylab <- paste(format(time),yunits,"Survival Probability")
        }
      else
        ylab <- scale[1]
    }

  trlab <- if(.R. && ylab=='X * Beta') function(k) expression(X*hat(beta))
  else function(k) k

  if(ref.zero & length(time))
    stop("ref.zero=T does not make sense with time given")

  labelc <- is.list(label.curves) || label.curves

  if(fname=='') factors <- dotlist else
  {
    factors <- list(NA)
    names(factors) <- fname
  }
  
  nf <- length(factors)
  
  if(nf<1) stop("must specify 1 or 2 factors to plot")

  which <- charmatch(names(factors), name, 0)
  if(any(which==0))
    stop(paste("factor(s) not in design:",
               paste(names(factors)[which==0],collapse=" ")))

  if(any(assume[which]==9))
    stop("cannot plot interaction terms")

  ##Number of non-slopes
  nrp <- num.intercepts(fit)
  if(missing(non.slopes))
    {
      non.slopes <- rep(0, nrp)
      if(!adj.zero)
        non.slopes[1] <- 1
    }
  
  if(nrp>0 && length(non.slopes)!=nrp)
    stop("wrong # values in non.slopes")

  if(is.logical(conf.int))
    {
      if(conf.int) conf.int <- .95
      else conf.int <- 0
    }

  if(conf.int)
    {
      vconstant <- 0
      if(conf.type=='individual')
        {
          vconstant <- fit$stats['Sigma']^2
          if(is.na(vconstant))
            stop('conf.type="individual" requires that fit be from ols')
        }
      zcrit <- if(length(idf <- fit$df.residual))
        qt((1 + conf.int) / 2, idf)
      else
        qnorm((1 + conf.int) / 2)
    }

  ix <- which[1]
  ixt <- assume[ix]
  Limval <- Getlim(at, allow.null=TRUE, need.all=FALSE)
  parmx <- parms[[name[ix]]]
  if(ixt!=7 && is.numeric(parmx)) parmx <- NULL

  if(ixt != 5 && ixt < 7)
    {
      Lv <- Limval$values
      if(any(names(Lv)==name[ix]))
        parmx <- Lv[[name[ix]]]
    }

  if(ixt>8)
    stop("matrix or interaction factor may not be displayed")

  xadj <- Limval$limits[2,,drop=FALSE]

  if(adj.zero) for(i in 1:length(xadj))
    {
      if(assume[i]==5 | assume[i]==8)
        xadj[[i]] <- parms[[name[i]]][1]
      else if(length(V <- Limval$values[[name[i]]]) & is.character(V))
        xadj[[i]] <- V[1]
      else xadj[[i]] <- 0
    }

  ##Use default adjusted to, replace some later.  Will be NA if 
  ## datadist doesn't have the variable

  nv <- 1
  xseq <- factors[[1]]

  if(nf>1)
    {
      y <- factors[[2]]
      ny <- length(y)
      if(ny>1 || (ny==1 && is.na(y))) nv <- 2	
    }

  if(missing(adj.subtitle))
    {
      if(add)adj.subtitle <- FALSE
      else adj.subtitle <- f-nv <= 6	
    }
  
  jf <- nv
  if(nf>nv) for(i in which[(nv+1):nf])
    {
      jf <- jf+1
      z <- factors[[jf]]
      if(!is.na(z))
        z <- value.chk(at, i, z, 0, Limval)
      
      if(length(z)!=1)
        stop("must specify single value for adjustment factors")

      if(!is.na(z))
        xadj[[name[i]]] <- z
    }

  beta <- fit$coef
  if(length(beta) & conf.int>0) cov <- Varcov(fit, regcoef.only=TRUE)

  plot.type <- "curves"
  curve.labels <- NULL
  xval <- parmx
  if(nv>1)
    {
      iy <- which[2]
      isna <- sapply(xadj,is.na); isna[c(ix,iy)] <- FALSE
      if(any(isna))
        stop(paste("variables not set to values or defined with datadist:",
                   paste(name[isna],collapse=" ")))

      iyt <- assume[iy]
      parmy <- parms[[name[iy]]]
      if(iyt!=5 && iyt<8 && length(V <- Limval$values[[name[iy]]]))
        parmy <- V

      if(iyt>8)
        stop("matrix or interaction factor may not be displayed")

      if(((length(xval)) | ixt==5 | ixt==7 | ixt==8) &
         (iyt!=5 & iyt!=7 & iyt!=8))
        warning("plot more effective if continuous factor on x-axis")

      y <- value.chk(at, iy, y, 40, Limval)
      ny <- length(y)
      if(!(ixt==5 | ixt==7 | ixt==8 | iyt==5 | iyt==7 | iyt==8 | 
           ny<=10)) plot.type <- "3d"
    }

  ## Use 40 x 40 grid if two continuous factors (for perspective or contour plot)

  if(plot.type=="3d") conf.int <- 0

  xseq <- value.chk(at, ix, xseq, 
                    if(plot.type=="curves")100
                    else 40,
                    Limval)

  if(is.factor(xseq)) xseq <- as.character(xseq)

  if(missing(xlab))
    xlab <- labelPlotmath(label[ix],units[name[ix]],
                          plotmath=!(plot.type=='3d' && method=='persp'))

  ##No intercept in model -> expand factors at adjustment values for later
  ##subtraction in variance.  Also needed for ref.zero option.

  if(ref.zero | (nrp==0 & conf.int>0))
    {
      if(any(sapply(xadj,is.na))) 
        stop("with conf.int>0 for ref.zero=T or Cox model, all variables must have been defined with datadist")

      xadjdf <- structure(xadj, class="data.frame", row.names="1")
      xsub <- predictDesign(fit,xadjdf,type="x",non.slopes=non.slopes)
    }

  if(ref.zero) ycenter <- matxv(xsub, beta) - Center
  else ycenter <- 0

  xseqn <- if(length(parmx) && is.character(parmx)) 1:length(xseq)
  else xseq

  if(val.lev) xseqn <- as.numeric(xseq)

  if(nv==1)
    {
      nna <- sapply(xadj, is.na); nna[ix] <- FALSE  ##ignore this one
      if(any(nna))
        stop(paste("variables not set to values here or defined with datadist:"
                   , paste(name[nna],collapse=" ")))

      if(missing(lty)) lty <- 1

      ##Expand xadj to ##rows=length(xseq), replace col. ix with xseq
      adj <- oldUnclass(xadj)  ##so can expand one column
      adj[[name[ix]]] <- xseq
      adj <- expand.grid(adj)
      if(!length(time))
        {
          xx <- predictDesign(fit, adj, type="x",non.slopes=non.slopes)
          if(length(attr(xx,"strata")) && any(is.na(attr(xx,"strata"))))
            warning("Computed stratum NA.  Requested stratum may not\nexist or reference values may be illegal strata combination\n")

          if(length(xx)==0)
            stop("model has no covariables and survival not plotted")

          xb <- matxv(xx, beta) - Center

          if(conf.int>0)
            {
              ##Subtract from rows if need to center for variance
              if(nrp==0) xxx <- sweep(xx,2,xsub)
              else xxx <- xx

              var <- ((xxx %*% cov) * xxx) %*% rep(1,ncol(xxx))
              lower <- xb - zcrit*sqrt(vconstant + var)
              upper <- xb + zcrit*sqrt(vconstant + var)
            }
          if(ref.zero)
            {
              xb <- xb-ycenter
              if(conf.int>0)
                {
                  lower <- lower-ycenter;
                  upper <- upper-ycenter
                }
            }
        }
      else
        {
          xb <- survest(fit, adj, times=time,loglog=loglog, 
                        conf.int=conf.int)
          lower <- xb$lower; upper <- xb$upper; xb <- xb$surv
        }
      
      if(!missing(fun))
        {
          xb <- fun(xb)
          if(conf.int>0)
            {
              lower <- fun(lower)
              upper <- fun(upper)
            }	
        }
      if(missing(ylim))
        {
          if(conf.int>0) ylim <- range(c(lower,upper),na.rm=TRUE)
          else ylim <- range(xb,na.rm=TRUE)
        }

      if((ixt==5 | ixt==8 | (!forceLines & length(xval))) & !val.lev)
        {
          if(method=='dotchart')
            {
              w <- xb
              names(w) <- if(is.numeric(xseq)) format(xseq)
              else if(abbrev)abbreviate(xseq) else xseq

              isd <- switch(sortdot,
                            ascending=order(w),
                            descending=order(-w),
                            neither=1:length(w))
            }
          
          if(add)
            {
              if(method=='dotchart') 
                dotchart2(w[isd], add=TRUE, pch=pch, reset.par=FALSE)
              else points(xseqn,xb,col=col[1],pch=pch)
            }
          else
            {
              labs <- format(xseq)
              if(is.character(xseq) ||
                 ((ixt==5 | ixt==8) || 
                  (length(xseq)==length(xval) &&
                   all(abs(xseq-as.numeric(xval))<.00001)))) {
                att <- 1
                at <- xseqn
              }
              else
                {
                  att <- 2
                  at <- Pretty(xseqn)
                }
              
              xlm <- if(missing(xlim)) range(at) else xlim
              
              if(method=='dotchart')
                dotchart2(w[isd], pch=pch, reset.par=FALSE,
                          xlim=ylim, xlab=trlab(ylab), ylab='')
              else
                {
                  plot(xseqn,xb,xlab=xlab, xlim=xlm, ##Handle NAs in Y
                       axes=FALSE, ylim=ylim, ylab=trlab(ylab), log=log, type='n')
                  points(xseqn,xb,col=col[1],pch=pch)
                  if(att==1)
                    mgp.axis(1,at=at,
                             labels=if(abbrev && ixt %in% c(5,8))
                             abbreviate(labs)
                             else labs)
                  else mgp.axis(1,at=at)
                  mgp.axis(2,at=pretty(ylim))
                }
            }
          
          if(conf.int>0)
            {
              if(method=='dotchart')
                {
                  dotchart2(lower[isd], add=TRUE, pch='[', reset.par=FALSE)
                  dotchart2(upper[isd], add=TRUE, pch=']', reset.par=FALSE)
                }
              else
                {
                  points(xseqn,lower,pch="-",col=col[1])
                  points(xseqn,upper,pch="-",col=col[1])
                }
            }
        }
      else
        {
          if(!add)
            {
              xlm <- if(missing(xlim)) range(Pretty(xseqn))
              else xlim
              
              plot(xseqn,xb,xlab=xlab, xlim=xlm, ylim=ylim,
                   ylab=trlab(ylab), log=log,
                   type='n', axes=FALSE)
              mgp.axis(1, at=pretty(xlm))
              mgp.axis(2, at=pretty(ylim))
            }
          lines(xseqn,xb,lty=lty[1],lwd=lwd[1],col=col[1])
          
          if(conf.int>0)
            {
              lines(xseqn,lower,lty=2,col=col[1],lwd=lwd.conf)
              lines(xseqn,upper,lty=2,col=col[1],lwd=lwd.conf)
            }
        }
      
      xx <- list(); xx[[name[ix]]] <- xseq
      xx <- structure(xx,class="data.frame",
                      row.names=as.character(1:length(xseq)))	
    }  ## end nv=1
  else
    {
      ##Expand xadj into length(xseq)*ny rows, replace columns
      ##ix and iy with xseq and y
      adj <- oldUnclass(xadj)   ##so can expand a column
      adj[[name[ix]]] <- NULL	##Guarantee y moves fastest (expand.grid moves
      ##first factor first
      adj[[name[iy]]] <- y
      adj[[name[ix]]] <- xseq
      adj <- expand.grid(adj)
      xx <- predictDesign(fit,adj,type="x",non.slopes=non.slopes)
      if(length(attr(xx,"strata")) && any(is.na(attr(xx,"strata"))))
        warning("Computed stratum NA.  Requested stratum may not\nexist or reference values may be illegal strata combination\n")
      
      if(length(xx)==0)
        {
          xb <- 0
          if(!length(time))
            stop("model has no covariables and survival not plotted") 
        }
      else
        xb <- matxv(xx, beta) - Center

      if(length(time))
        {
          xb <- survest(fit, adj, times=time, loglog=loglog, 
                        conf.int=conf.int)
          lower <- xb$lower; upper <- xb$upper; xb <- xb$surv
        }
      else {
        if(conf.int>0)
          {
            if(nrp==0)
              xxx <- sweep(xx,2,xsub)
            else xxx <- xx
            
            var <- ((xxx %*% cov) * xxx) %*% rep(1,ncol(xxx))
            lower <- xb - zcrit*sqrt(vconstant + var)
            upper <- xb + zcrit*sqrt(vconstant + var)			
          }
        
        if(ref.zero)
          {
            xb <- xb-ycenter
            if(conf.int>0)
              {
                lower <- lower-ycenter;
                upper <- upper-ycenter
              }	
          }
      }
      
      if(!missing(fun))
        {
          xb <- fun(xb)
          if(conf.int>0)
            {
              lower <- fun(lower)
              upper <- fun(upper)
            }
        }

      if(!missing(perim) && plot.type=="3d")
        {
          if(.SV4.)
            perim <- matrix(oldUnclass(perim), nrow=nrow(perim),
                            dimnames=dimnames(perim))

          Ylo <- approx(perim[,1],perim[,2],adj[[name[ix]]])$y
          Yhi <- approx(perim[,1],perim[,3],adj[[name[ix]]])$y
          Ylo[is.na(Ylo)] <-  1e30
          Yhi[is.na(Yhi)] <- -1e30
          xb[adj[[name[iy]]] < Ylo] <- NA
          xb[adj[[name[iy]]] > Yhi] <- NA
        }
      

      if(missing(ylim))
        {
          if(conf.int>0) ylim <- range(c(lower,upper),na.rm=TRUE)
          else ylim <- range(xb,na.rm=TRUE)
        }
      
      xx <- adj[,c(name[ix],name[iy])]
      
      if(plot.type=="3d")
        {
          zmat <- matrix(pmin(ylim[2],pmax(ylim[1],xb)),
                         nrow=length(xseqn),
                         ncol=ny, byrow=TRUE)
          
          laby <- labelPlotmath(label[iy],units[name[iy]],
                                plotmath=method!='persp')
          if(method=="contour")
            {
              contour(xseqn, y, zmat,
                      nlevels=nlevels, xlab=xlab, ylab=laby, col=col)
            }
          else if(method=="persp")
            {
              a <- c(list(xseqn, y, zmat, zlim=zlim, xlab=xlab, ylab=laby,
                          zlab=trlab(ylab), box=TRUE, col=col), perspArgs)
              a <- c(a,
                     if(.R.) list(theta=theta, phi=phi)
                     else if(!missing(eye)) list(eye=eye))
              
              do.call('persp', a)
            }
          else
            {
              if(length(col) > 1)
                image(xseqn, y, zmat, xlab=xlab, ylab=laby, col=col)
              else image(xseqn, y, zmat, xlab=xlab, ylab=laby)
            }
        }
      else
        {
          ## One curve for each value of y, excl style used for C.L.
          
          lty <- if(missing(lty)) seq(ny+1)[-2]
          else rep(lty, length=ny)
          
          i <- 0
          if(labelc)
            curves <- vector('list',ny)

          col <- rep(col, length=ny)
          lwd <- rep(lwd, length=ny)
          for(ay in y)
            {
              i <- i+1
              index.y <- (1:nrow(xx))[xx[,2]==ay]
              xseq <- xx[index.y, name[ix]]
              if(is.factor(xseq))
                xseq <- as.character(xseq)
              
              xl <- if(val.lev) as.single(xseq)
              else if(is.character(xseq)) match(xseq, parmx)
              else xseq
              
              curve.labels <-c(curve.labels, format(ay))
              if(!missing(perim))
                {
                  if(!is.function(perim))
                    stop('perim must be a function')
                  
                  show.pt <- perim(xl, ay)
                  xb[index.y][!show.pt] <- NA
                  if(conf.int)
                    {
                      lower[index.y][!show.pt] <- NA
                      upper[index.y][!show.pt] <- NA
                    }
                }
              if(labelc) curves[[i]] <- list(xl, xb[index.y])
              
              if(i>1 | add)
                lines(xl,xb[index.y],
                      lty=lty[i],col=col[i],lwd=lwd[i])
              else
                {
                  if((ixt==5 | ixt==8 | (length(xval))) &
                     !val.lev)
                    {
                      labs <- format(xseq)
                      if(is.character(xseq) || 
                         ((ixt==5 | ixt==8) || (length(xseq)== length(xval) &&
                                      all(abs(xseq-as.numeric(xval))<.00001))))
                        {
                          att <- 1
                          at <- xl
                        }
                      else
                        {
                          att <- 2
                          at <- Pretty(xl)
                        }
                      
                      xlm <- if(missing(xlim)) range(at) else xlim
                      
                      plot(xl,xb[index.y],log=log,
                           xlab=xlab,ylab=trlab(ylab),
                           xlim=xlm, ylim=ylim,
                           type='n', axes=FALSE)
                      lines(xl,xb[index.y],lty=lty[i],lwd=lwd[i],
                            col=col[i])
                      if(att==1) mgp.axis(1,at=at,
                           label=if(abbrev && ixt %in% c(5,8))
                           abbreviate(labs)
                           else labs)
                      else mgp.axis(1,at=at)
                      
                      mgp.axis(2,at=pretty(ylim))
                    }
                  else
                    {
                      xlm <- if(missing(xlim)) range(Pretty(xl)) else xlim
                      
                      plot(xl,xb[index.y],xlab=xlab,
                           ylab=trlab(ylab),
                           xlim=xlm,ylim=ylim,log=log,type='n',
                           axes=FALSE)
                      mgp.axis(1,at=pretty(xlm))
                      mgp.axis(2,at=pretty(ylim))
                      lines(xl,xb[index.y],col=col[i],
                            lwd=lwd[i],lty=lty[i])
                    }
                }
              if(conf.int>0)
                {
                  lines(xl,lower[index.y],
                        lty=2,col=col[i],lwd=lwd.conf)
                  lines(xl,upper[index.y],
                        lty=2,col=col[i],lwd=lwd.conf) 
                }
            }
          
          if(labelc)
            labcurve(curves, curve.labels, 
                     opts=label.curves, lty=lty, lwd=lwd, col=col)
        }
    }
  
  xx[[if(ylab=="")"Z" else ylab]] <- xb
  
  if(conf.int>0)
    {
      xx$lower <- lower;
      xx$upper <- upper
    }

  adjust <- ""
  if(f>nv) for(i in 1:f)
    {
      if(!any(i==which[1:nv]))
        adjust <- paste(adjust, name[i], "=", 
                        if(is.factor(xadj[[i]])) as.character(xadj[[i]])
                        else format(xadj[[i]])," ",sep="")
    }

  if(adjust!="" & adj.subtitle)
    title(sub=paste("Adjusted to:",adjust),adj=0,cex=cex.adj)

  R <- list(x.xbeta=xx, adjust=adjust, curve.labels=curve.labels, 
            plot.type=plot.type, method=method, col=col, 
            lty=if(plot.type=='curves') lty, lwd=lwd)
  oldClass(R) <- "plot.Design"
  invisible(R)
}

print.plot.Design <- function(x, ...)
{
  print(x$x.xbeta)
  if(x$adjust!="")
    cat("Adjust to:",x$adjust,"\n")
  
  if(length(cl <- x$curve.labels))
    cat("Curves:",cl,"\n")

  invisible()
}

Legend <- function(object, ...) UseMethod("Legend")

Legend.plot.Design <- function(object, x, y, size = c(1, 1), 
                               horizontal = TRUE, nint = 50, fun, at, 
                               zlab,  ...)
{
  
  if(missing(x)) 
    if(object$method=="image" && missing(size))
      {
        cat("Using function \"locator(2)\" to place opposite corners of image.legend\n")
        x <- locator(2)
      }
    else
      {
        cat("Using function \"locator(1)\" to place upper left corner of legend\n")
        x <- locator(1)
      }
  
  if(!missing(y)) x <- list(x=x,y=y)
  xb <- object$x.xbeta
  
  if(object$method!="image")
    stop('expecting to use results from method="image"')
  
  z <- xb[,3]
  irgz <- seq(min(z,na.rm=TRUE), max(z,na.rm=TRUE), length = nint)
  lirgz <- length(irgz)
  if(horizontal)
    {
      f <- expression({
        if(length(list(...)))
          par(...) ##axis() does not respect mgp
        
        image(x=irgz, y=1:lirgz, z=matrix(irgz, lirgz, lirgz), yax="n",
              xaxt=if(missing(fun))"s"
              else "n")
        
        if(!missing(fun))
          mgp.axis(1,
                   if(missing(at)) pretty(irgz)
                   else at,
                   labels=format(fun(if(missing(at)) pretty(irgz)
                   else at)))
        
        title(xlab=if(missing(zlab)) names(xb)[3]
        else zlab)
      })
      subplot(x = x, y = y, size = size, fun = f, hadj=0, vadj=1)
    }
  else
    {
      f <- expression({
        if(length(list(...)))
          par(...)
        
        image(x = 1:lirgz, y = irgz, z = matrix(irgz, lirgz, lirgz,byrow=TRUE), 
              xaxt = "n", yaxt=if(missing(fun))"s" else "n")
        
        if(!missing(fun))
          mgp.axis(2, if(missing(at)) pretty(irgz) else at,
                   labels=format(fun(if(missing(at)) pretty(irgz) else at)))
        
        title(ylab=if(missing(zlab)) names(xb)[3] else zlab)
      })
      
      subplot(x = x, y = y, size = size, fun = f, hadj=0, vadj=1)
    }
  
  invisible(if(missing(y))x else list(x=x,y=y))
}

perimeter <- function(x, y, xinc=diff(range(x))/10, n=10,
                      lowess.=TRUE)
{

  s <- !is.na(x+y)
  x <- x[s]
  y <- y[s]
  m <- length(x)
  if(m<n)
    stop("number of non-NA x must be >= n")

  i <- order(x)
  x <- x[i]
  y <- y[i]
  s <- n:(m-n+1)
  x <- x[s]
  y <- y[s]

  x <- round(x/xinc)*xinc

  g <- function(y, n)
    {
      y <- sort(y)
      m <- length(y)
      if(n>(m-n+1)) c(NA,NA)
      else c(y[n], y[m-n+1])
    }

  r <- unlist(tapply(y, x, g, n=n))
  i <- seq(1, length(r), by=2)
  rlo <- r[i]
  rhi <- r[-i]
  s <- !is.na(rlo+rhi)
  if(!any(s))
    stop("no intervals had sufficient y observations")

  x <- sort(unique(x))[s]
  rlo <- rlo[s]
  rhi <- rhi[s]
  if(lowess.)
    {
      rlo <- lowess(x, rlo)$y
      rhi <- lowess(x, rhi)$y
    }

  structure(cbind(x, rlo, rhi), dimnames=list(NULL,
                                  c("x","ymin","ymax")), class='perimeter')
}

lines.perimeter <- function(x, ...)
{
  lines(x[,'x'], x[,'ymin'],...)
  lines(x[,'x'], x[,'ymax'],...)
  invisible()
}

datadensity.plot.Design <- function(object, x1, x2, ...)
{
  if(missing(x1))
    stop('must specify x1')
  
  r <- object$x.xbeta
  nam <- names(r)
  x1.name <- deparse(substitute(x1))
  if(x1.name!=nam[1])
    warning(paste(x1.name,
                  ' is not first variable mentioned in plot() (',nam[1],')',sep=''))
  
  if(missing(x2))
    {
      x <- r[[1]]
      y <- r[[2]]
      x1 <- x1[!is.na(x1)]
      y.x1 <- approx(x, y, xout=x1)$y
      scat1d(x1, y=y.x1, ...)
      return(invisible())
    }
  
  x2.name <- deparse(substitute(x2))
  if(x2.name!=nam[2])
    warning(paste(x2.name,
                  ' is not second variable mentioned in plot() (',nam[2],')',sep=''))

  x     <- r[[1]]
  curve <- r[[2]]
  y     <- r[[3]]
  for(s in if(is.factor(curve))levels(curve) else unique(curve))
    {
      i <- curve==s
      xs <- x[i]
      ys <- y[i]
      x1s <- x1[x2==s]
      x1s <- x1s[!is.na(x1s)]
      y1s <- approx(xs, ys, xout=x1s)$y
      scat1d(x1s, y=y1s, ...)
    }

  invisible()
}
