nomogram <- function(fit, ...) UseMethod("nomogram")

nomogram.Design <- function(fit, ..., adj.to, 
			lp=TRUE, lp.at, lplabel="Linear Predictor",
			fun, fun.at, fun.lp.at, funlabel="Predicted Value",
			fun.side,
			interact=NULL, intercept=1,
			conf.int=FALSE, 
			col.conf=c(1,if(under.unix).3 else 12),
			conf.space=c(.08,.2), 
			conf.lp=c("representative", "all", "none"),
			est.all=TRUE, abbrev=FALSE, minlength=4, maxscale=100, nint=10, 
			label.every=1, force.label=FALSE, 
			xfrac=.35, cex.axis=.85, cex.var=1,
			col.grid=FALSE, vnames=c("labels","names"),
			varname.label=TRUE, varname.label.sep="=", ia.space=.7, 
			tck=-.009, lmgp=.4, omit=NULL, naxes,
			points.label='Points', total.points.label='Total Points',
			total.sep.page=FALSE, total.fun, verbose=FALSE) {	

conf.lp <- match.arg(conf.lp)
vnames  <- match.arg(vnames)

abb <- (is.logical(abbrev) && abbrev) || is.character(abbrev)
if(is.logical(conf.int) && conf.int) conf.int <- c(.7,.9)
if(!is.logical(conf.int) && (length(conf.int)!=length(col.conf)))
  stop("conf.int and col.conf must have same length")
if(is.logical(col.grid) && col.grid) col.grid <- if(under.unix).2 else 5

oldpar <- oPar()  # in Hmisc Misc.s
mgp <- oldpar$mgp
mar <- oldpar$mar
par(mgp=c(mgp[1],lmgp,mgp[3]),mar=c(mar[1],1.1,mar[3],mar[4]))
on.exit(setParNro(oldpar))  ## was par(oldpar) 11Apr02
tck2 <- tck/2

se <- FALSE
if(any(conf.int>0)) {
  se <- TRUE
  zcrit <- qnorm((conf.int+1)/2)
  bar <- function(x, y, zcrit, se, col.conf, nlev=4) {
    y <- rep(seq(y[1], y[2], length=nlev), length.out=length(x))
    for(j in 1:length(x)) {
      xj <- x[j]; yj <- y[j]
      W <- c(0,zcrit)*se[j]
      for(i in 1:length(zcrit)) {
        segments(xj-W[i+1], yj, xj-W[i], yj, col=col.conf[i], lwd=1)
        segments(xj+W[i+1], yj, xj+W[i], yj, col=col.conf[i], lwd=1)
      }
    }
  }
}

nfun <- if(missing(fun)) 0 else if(is.list(fun)) length(fun) else 1
if(nfun>1 && length(funlabel)==1) funlabel <- rep(funlabel, nfun)
if(nfun>0 && is.list(fun) && length(names(fun))) funlabel <- names(fun)

if(!missing(fun.at)) {
  if(!is.list(fun.at)) fun.at <- rep(list(fun.at),nfun)
}
if(!missing(fun.lp.at)) {
  if(!is.list(fun.lp.at)) fun.lp.at <- rep(list(fun.lp.at),nfun)
}
if(!missing(fun.side)) {
  if(!is.list(fun.side)) fun.side <- rep(list(fun.side),nfun)
  if(any(!(unlist(fun.side) %in% c(1,3))))
    stop('fun.side must contain only the numbers 1 and 3')
}
at <- fit$Design
if(!length(at)) at <- getOldDesign(fit)
assume <- at$assume.code
if(any(assume==10))
 warning("does not currently work with matrix factors in model")
name  <- at$name
names(assume) <- name
parms <- at$parms
label <- if(vnames=="labels") at$label else name
if(any(d <- duplicated(name))) stop(paste("duplicated variable names:",
	paste(name[d],collapse=" ")))
label <- name
if(vnames=="labels") {
  label <- at$label
  if(any(d <- duplicated(label))) stop(paste("duplicated variable labels:",
	paste(label[d],collapse=" ")))
}

ia    <- at$interactions

factors <- list(...)
nf <- length(factors)

which <- if(est.all) (1:length(assume))[assume!=8] else 
	(1:length(assume))[assume!=8 & assume!=9]
if(nf>0) {
  jw <- charmatch(names(factors),name,0)
  if(any(jw==0))stop(paste("factor name(s) not in the design:",
		paste(names(factors)[jw==0],collapse=" ")))
  if(!est.all) which <- jw
}

Limval <- Getlim(at, allow.null=TRUE, need.all=FALSE)
values <- Limval$values
lims <- Limval$limits[c(6,2,7),,drop=FALSE]

#Next 4 lines 27Nov99 - keep character variables intact
lims <- oldUnclass(lims)
for(i in 1:length(lims))
   if(is.factor(lims[[i]]))lims[[i]] <- as.character(lims[[i]])
attr(lims, 'class') <- 'data.frame'  # so can subscript later

#Find underlying categorical variables
ucat <- rep(FALSE, length(assume))
names(ucat) <- name
for(i in (1:length(assume))[assume!=5 & assume<8]) {
   ucat[i] <- !is.null(V <- values[[name[i]]]) # did add && is.character(V)
   if(ucat[i]) parms[[name[i]]] <- V
}

discrete <- assume==5 | assume==8 | ucat
names(discrete) <- name

#Number of non-slopes:
nrp <- num.intercepts(fit)
Intercept <- if(nrp>0) fit$coefficients[intercept] else 
             if(!is.null(fit$center)) -fit$center else 0
#non.slopes <- rep(0, nrp)   23Jun95
#non.slopes[intercept] <- 1

intercept.offset <- if(nrp<2) 0 else
   fit$coefficients[intercept]-fit$coefficients[1]
#linear.predictors stored in fit always used first intercept

settings <- list()
for(i in which[assume[which]<9]) {
  ni <- name[i]
  z <- factors[[ni]]
  lz <- length(z)
  if(lz < 2) settings[[ni]] <- value.chk(at, i, NA, -nint, Limval, 
				type.range="full") else
  if(lz > 0 && any(is.na(z))) stop("may not specify NA as a variable value")
  if(lz==1) lims[2,i] <- z else if(lz>1) {
    settings[[ni]] <- z
    if(is.null(lims[[ni]]) || is.na(lims[2,ni])) {
      lims[[ni]] <- c(NA,z[1],NA)
      warning(paste("adjustment values for ",ni,
        " not defined in datadist; taken to be first value specified (", 
        z[1],")" ,sep=""))
    }
  }
}

adj <- lims[2,,drop=FALSE]
if(!missing(adj.to)) for(nn in names(adj.to)) adj[[nn]] <- adj.to[[nn]]
isna <- sapply(adj, is.na)
if(any(isna)) stop(
   paste("adjustment values not defined here or with datadist for",
	paste(name[assume!=9][isna],collapse=" ")))

num.lines <- 0
space.used <- 0
entities <- 0
set <- list()
iset <- 0
start <- len <- NULL
end <- 0

#Sort to do continuous factors first if any interactions present
main.effects <- which[assume[which]<8]  # this logic not handle strata w/intera.
if(any(assume==9)) main.effects <- main.effects[order(10*discrete[main.effects]+
	(name[main.effects] %in% names(interact)))]

#For each predictor, get vector of predictor numbers directly or
#indirectly associated with it
rel <- related.predictors(at)   # Function in Design.Misc.s

already.done <- structure(rep(FALSE,length(name)), names=name)
for(i in main.effects) {
  nam <- name[i]
  if(already.done[nam] || (nam %in% omit)) next
  r <- sort(rel[[i]])
  if(length(r)==0) { #main effect not contained in any interactions
    num.lines <- num.lines + 1
    space.used <- space.used + 1
    entities <- entities+1
    x <- list()
    x[[nam]] <- settings[[nam]]
    iset <- iset+1
    attr(x,'info') <- list(predictor=nam, effect.name=nam,
                           type='main')
    set[[label[i]]] <- x
    start <- c(start, end+1)
    n <- length(settings[[nam]])
    len <- c(len, n)
    end <- end+n
    NULL	#handles bug in S
  }  else {
    namo <- name[r]
    s <- !(name[r] %in% names(interact))
    if(any(s)) {
      if(is.null(interact)) interact <- list()
      for(j in r[s]) {
        nj <- name[j]
        if(discrete[j]) interact[[nj]] <-  parms[[nj]]
	NULL
      }
      s <- !(name[r] %in% names(interact))
    }
    if(any(s)) stop(paste("factors not defined in interact=list(...):",
                          paste(name[r[s]],collapse=",")))
    combo <- expand.grid(interact[namo]) #list[vector] gets sublist
    oldClass(combo) <- NULL   # so combo[[n]] <- as.character will really work
    acombo <- combo
    if(abb) for(n in if(is.character(abbrev))abbrev else names(acombo)) {
      if(discrete[n]) { 
       acombo[[n]] <- abbreviate(parms[[n]],
         minlength=if(minlength==1)4 else minlength)[combo[[n]]]  
       #lucky that abbreviate function names its result
      }
    }
    for(n in names(combo)) if(is.factor(combo[[n]])) {
      combo[[n]] <- as.character(combo[[n]])   #so row insertion will work xadj
      acombo[[n]] <- as.character(acombo[[n]]) #so format() will work
      NULL
    }
    entities <- entities+1
    already.done[namo] <- TRUE
    for(k in 1:length(combo[[1]])) {		# was nrow(combo))
      num.lines <- num.lines+1
      space.used <- space.used + if(k==1) 1 else ia.space
      x <- list()
      x[[nam]] <- settings[[nam]]   #store fastest first
      for(nm in namo) x[[nm]] <- combo[[nm]][k]   #was combo[k,nm] 2Dec94
      iset <- iset+1
      set.name <- paste(nam, " (",sep="")
      for(j in 1:length(acombo)) {		# was ncol
        set.name <- paste(set.name, 
          if(varname.label) paste(namo[j],varname.label.sep,sep="") else "",
          format(acombo[[j]][k]),sep="")   # was acombo[k,j]
        if(j<length(acombo)) set.name <- paste(set.name," ",sep="")
      }
      set.name <- paste(set.name,")",sep="")
     #Make list of all terms needing inclusion in calculation
     #Include interation term names  - interactions.containing in Design.Misc.s
      ia.names <- NULL
      for(j in r) ia.names <- c(ia.names, name[interactions.containing(at, j)])
      ia.names <- unique(ia.names)
      attr(x,'info') <- list(predictor=nam,
                             effect.name=c(nam,namo[assume[namo]!=8],ia.names), 
                             type=if(k==1) "first" else "continuation")
      set[[set.name]] <- x
			##Don't include strata main effects
      start <- c(start, end+1)
      n <- length(settings[[nam]])
      len <- c(len, n)
      end <- end+n
      NULL
    }
    NULL
  }    
}
xadj <- oldUnclass(Design.levels(adj, at))
for(k in 1:length(xadj)) xadj[[k]] <- rep(xadj[[k]], sum(len))

j <- 0
for(S in set) {
  j <- j+1
  ns <- names(S)
  nam <- names(S)
  for(k in 1:length(nam)) {
    xadj[[nam[k]]][start[j]:(start[j]+len[j]-1)] <- S[[k]]
    NULL
  }
}
xadj <- structure(xadj, class='data.frame',
                  row.names=as.character(1:sum(len)))   # 27Nov99

xx <- predictDesign(fit, newdata=xadj, type="terms", center.terms=FALSE, se.fit=FALSE,
	kint=intercept)   #non.slopes=non.slopes) 23Jun95
##20Nov00:
if(any(is.infinite(xx)))
  stop("variable limits and transformations are such that an infinite axis value has resulted.\nRe-run specifying your own limits to variables.")

if(se) xse <- predictDesign(fit, newdata=xadj, se.fit=TRUE, 
	kint=intercept)   #non.slopes=non.slopes) 23Jun95

R <- matrix(NA, nrow=2, ncol=length(main.effects),
            dimnames=list(NULL,name[main.effects]))
R[1,] <- 1e30
R[2,] <- -1e30
# R <- apply(xx, 2, range)  - does not work since some effects are for
# variable combinations that were never used in constructing axes

for(i in 1:num.lines) {
  is <- start[i]; ie <- is+len[i]-1
  s <- set[[i]]
  setinfo <- attr(s,'info')
  nam <- setinfo$effect.name
  xt <- xx[is:ie, nam]
  if(length(nam)>1) xt <- apply(xt, 1, sum)  # add all terms involved
  set[[i]]$Xbeta <- xt
  r <- range(xt)
  pname <- setinfo$predictor
  R[1,pname] <- min(R[1,pname], r[1])
  R[2,pname] <- max(R[2,pname], r[2])
  if(se) {
    set[[i]]$Xbeta.whole <- xse$linear.predictors[is:ie] #note-has right interc.
    set[[i]]$se.fit      <- xse$se.fit[is:ie]
    NULL
  }
  NULL
}


R <- R[,R[1,]<1e30,drop=FALSE]
sc <- maxscale/max(R[2,]-R[1,])
Intercept <- Intercept + sum(R[1,])

xl <- -xfrac*maxscale
if(missing(naxes)) naxes <- 
  if(total.sep.page) max(space.used + 1, nfun + lp + 1) else
                     space.used + 1 + nfun + lp + 1

Format <- function(x) { # like format but does individually
  f <- character(l <- length(x))
  for(i in 1:l) f[i] <- format(x[i])
  f
}

newpage <- function(naxes,xl,maxscale,cex.var,nint,space.used,col.grid,
		cex.axis,tck,tck2,label.every,force.label,
		points=TRUE, points.label='Points',usr) {
  y <- naxes-1
  plot(0,0,xlim=c(xl,maxscale),ylim=c(0,y),
	   type="n",axes=FALSE,xlab="",ylab="")
  if(!missing(usr)) par(usr=usr)
  if(!points) return(y + 1)

  ax <- c(0,maxscale)
  text(xl, y, points.label, adj=0, cex=cex.var)
  x <- pretty(ax, n=nint)
  x2 <- seq((x[1]+x[2])/2, max(x), by=x[2]-x[1])
  if(col.grid>0) {
	segments(x ,y,x, y-space.used,col=col.grid,lwd=1)
	segments(x2,y,x2,y-space.used,col=if(col.grid==1)1 else col.grid/2,lwd=1)
  }  
  axisf(3, at=x,  pos=y, cex=cex.axis, tck=tck, label.every=label.every, 
	   force.label=force.label)
  axisf(3, at=x2, labels=FALSE, pos=y, tck=tck2, cex=cex.axis)
  y
}

y <- newpage(naxes,xl,maxscale,cex.var,nint,space.used,col.grid,
			cex.axis,tck,tck2,label.every=label.every,
			force.label=force.label,points.label=points.label)

i <- 0
ns <- names(set)
Abbrev <- list()
for(S in set) {
  i <- i+1
#  type <- attr(S,"type")  24Sep00
  setinfo <- attr(S,'info')
  type <- setinfo$type
  y <- y - (if(type=="continuation") ia.space else 1)
  if(y < -.05) {
	y <- newpage(naxes,xl,maxscale,cex.var,nint,space.used,col.grid,
			cex.axis,tck,tck2,
			label.every=label.every,force.label=force.label,
			points.label=points.label) -
			  (if(type=="continuation") ia.space else 1)
  }

  text(xl, y, ns[[i]], adj=0, cex=cex.var)
  x <- S[[1]]
  nam <- names(S)[1]  #stored with fastest first
  fx <- if(is.character(x)) x else 
          sedit(Format(x)," ","") #axis not like bl   - was translate()
  if(abb && discrete[nam] && (is.logical(abbrev) || nam %in% abbrev)) {
    old.text <- fx
    fx <- if(abb && minlength==1)letters[1:length(fx)] else 
          abbreviate(fx, minlength=minlength)
    Abbrev[[nam]] <- list(abbrev=fx, full=old.text)
  }

  j <- match(nam, name, 0)
  if(any(j==0)) stop("program logic error 1")
#  if(!discrete[nam] && label.every>1) {
#    sq <- seq(along=x, by=label.every)
#    fx[-sq] <- ""}
  is <- start[i]; ie <- is+len[i]-1
  xt <- (S$Xbeta - R[1,nam])*sc
  set[[i]]$points <- xt
 #Find flat pieces and combine their labels
  r <- rle(xt)
  if(any(r$length>1)) {
    is <- 1
    for(j in r$length) {
      ie <- is+j-1
      if(j>1) {
        fx[ie] <- if(discrete[nam] || ie < length(xt))
		  paste(fx[is], "-", fx[ie],sep="") else
		paste(fx[is], '+', sep='')

        fx[is:(ie-1)] <- ""
	xt[is:(ie-1)] <- NA
      }
      is <- ie+1
    }
  fx <- fx[!is.na(xt)]
  xt <- xt[!is.na(xt)]
  }
 #Find direction changes
  ch <- if(length(xt)>2) c(FALSE,FALSE,diff(diff(xt)>0)!=0) else rep(FALSE, length(xt))
  if(discrete[nam] && length(xt)>1) { # categorical - alternate adjacent levels
    j <- order(xt)
    lines(range(xt),rep(y,2))   # make sure all ticks are connected
    for(k in 1:2) {
      is <- j[seq(k, length(j), by=2)]
	  axisf(1+2*(k==2), at=xt[is], labels=fx[is], pos=y, cex=cex.axis,
		   tck=tck, force.label=force.label || (
		   abb && minlength==1 && (is.logical(abbrev) || nam %in% abbrev)),
		   disc=TRUE)
      if(se) bar(xt[is], if(k==1) y-conf.space-.32 else y+conf.space+.32, 
				 zcrit, sc*S$se.fit[is], col.conf)
    }
  } else if(!any(ch)) {
    axisf(1, at=xt, labels=fx, pos=y, cex=cex.axis, tck=tck,
		 label.every=label.every, force.label=force.label, disc=discrete[nam])
    if(se)bar(xt, y+conf.space, zcrit, sc*S$se.fit, col.conf)
  }
  else {
    lines(range(xt), rep(y,2))  # make sure all ticks are connected
    j <- (1:length(ch))[ch]
    if(max(j)<length(ch)) j <- c(j, length(ch)+1)
    side <- 1
    is <- 1
    for(k in j) {
      ie <- k-1
      axisf(side, at=xt[is:ie], labels=fx[is:ie], pos=y, cex=cex.axis, tck=tck,
		   label.every=label.every, force.label=force.label, 
		   disc=discrete[nam])
      if(se) bar(xt[is:ie], if(side==1) y-conf.space-.32 else y+conf.space+.32, 
                 zcrit, sc*S$se.fit[is:ie], col.conf)
      side <- if(side==3)1 else 3
      is <- ie+1
    }
  }
}

if(missing(lp.at)) {
  xb <- fit$linear.predictors
  if(is.null(xb)) xb <- fit$fitted
  if(is.null(xb)) stop("lp.at not given and fit did not store linear.predictors or fitted.values")
  if(nrp>1) xb <- xb + intercept.offset
  lp.at <- pretty(range(xb), n=nint)
}

sum.max <- if(entities==1) maxscale else max(maxscale,sc*max(lp.at-Intercept))
x <- pretty(c(0, sum.max), n=nint)

new.max <- max(x)
xl.old <- xl
xl <- -xfrac*new.max
u <- par()$usr
if(!missing(total.fun)) total.fun()
usr <- c(xl*u[1]/xl.old, new.max*u[2]/maxscale, u[3:4])
par(usr=usr)

x.double <- seq(x[1], new.max, by=(x[2]-x[1])/2)

y <- y-1

if(y < -.05 || total.sep.page) 
  y <- newpage(naxes,xl,maxscale,cex.var,nint,space.used,col.grid,
			cex.axis,tck,tck2,
			label.every=label.every,force.label=force.label,
			points=FALSE,usr=usr) - 1

text(xl, y, total.points.label, adj=0, cex=cex.var)
axisf(1, at=x, pos=y, cex=cex.axis, tck=tck, label.every=label.every,
	 force.label=force.label)
axisf(1, at=x.double, labels=FALSE, pos=y, tck=tck2, cex=cex.axis)
set$total.points <- list(x=x, y=y)

if(lp) {
  x2 <- seq(lp.at[1], max(lp.at), by=(lp.at[2]-lp.at[1])/2)
  scaled.x <- (lp.at-Intercept)*sc
  scaled.x2 <- (x2-Intercept)*sc
  y <- y-1
  if(y < -.05) 
	y <- newpage(naxes,xl,maxscale,cex.var,nint,space.used,col.grid,
			cex.axis,tck,tck2,
			label.every=label.every,force.label=force.label,
			points=FALSE,usr=usr) - 1

  text(xl, y, lplabel, adj=0, cex=cex.var)
  axisf(1, at=scaled.x,  labels=Format(lp.at), pos=y, cex=cex.axis, tck=tck,
	   label.every=label.every, force.label=force.label)
  axisf(1, at=scaled.x2, labels=FALSE, tck=tck2,   pos=y, cex=cex.axis)
  set$lp <- list(x=scaled.x, y=y, x.real=lp.at)
  if(se && conf.lp!="none") {
    xxb <- NULL
    xse <- NULL
    for(S in set) { xxb <- c(xxb, S$Xbeta.whole); xse <- c(xse, S$se.fit) }
    i <- order(xxb)
    if(length(xxb)<16 | conf.lp=="representative") 
	{nlev <- 4; w <- 1} else {nlev <- 8; w <- 2}
    if(conf.lp=="representative") {
      deciles <- cut2(xxb[i], g=10)
      mean.xxb <- tapply(xxb[i], deciles, mean)
      median.se <- tapply(xse[i], deciles, median)
      bar((mean.xxb-Intercept)*sc, 
	  y+c(conf.space[1],conf.space[1]+w*diff(conf.space)),
          zcrit, sc*median.se, col.conf, nlev=nlev)
    } else
    bar((xxb[i]-Intercept)*sc, y+c(conf.space[1],
	conf.space[1]+w*diff(conf.space)), 
        zcrit, sc*xse[i], col.conf, nlev=nlev)
  }
}

if(nfun>0) {
  if(!is.list(fun)) fun <- list(fun)
  i <- 0
  for(func in fun) {
    i <- i+1
   #Now get good approximation to inverse of fun evaluated at fat
   #Unless inverse function given explicitly
    if(!missing(fun.lp.at)) {
      xseq <- fun.lp.at[[i]]
      fat <- func(xseq)
      w <- xseq
    } else {
      if(missing(fun.at)) fat <- pretty(func(range(lp.at)), n=nint)
      else fat <- fun.at[[i]]
      if(verbose) {
        cat('Function',i,'values at which to place tick marks:\n')
        print(fat)
      }
      xseq <- seq(min(lp.at),max(lp.at),length=1000)
      fu <- func(xseq)
      s <- !is.na(fu)
      w <- approx(fu[s], xseq[s], fat)$y
      if(verbose) {
        cat('Estimated inverse function values (lp):\n')
        print(w)
      }
    }
    s <- !(is.na(w)|is.na(fat))
    w <- w[s]
    fat <- fat[s]
    fat.orig <- fat
    fat <- if(is.category(fat)) as.character(fat) else Format(fat)
    scaled <- (w-Intercept)*sc
    y <- y-1
	if(y < -.05) 
	  y <- newpage(naxes,xl,maxscale,cex.var,nint,space.used,col.grid,
			cex.axis,tck,tck2,
			label.every=label.every,force.label=force.label,
			points=FALSE,usr=usr) - 1

    text(xl, y, funlabel[i], adj=0, cex=cex.var)
    sides <- if(missing(fun.side)) rep(1, length(fat)) else (fun.side[[i]])[s]
    if(length(sides)!=length(fat)) 
      stop('fun.side vector not same length as fun.at or fun.lp.at')
    for(jj in 1:length(fat)) 
		if(.R.)axis(sides[jj], at=scaled[jj], label=fat[jj],
                    pos=y, cex.axis=cex.axis, tck=tck) else
		       axis(sides[jj], at=scaled[jj], label=fat[jj],
                    pos=y, cex=cex.axis, tck=tck)
    lines(range(scaled),rep(y,2))  #make sure all ticks are connected
    set[[funlabel[i]]] <- list(x=scaled, y=y, x.real=fat.orig)
  }
}
set$abbrev <- Abbrev
oldClass(set) <- "nomogram"
invisible(set)
}

print.nomogram <- function(x, dec=0, ...) {

  obj <- x
  
  w <- diff(range(obj$lp$x))/diff(range(obj$lp$x.real))
  cat('Points per unit of linear predictor:',format(w),
	  '\nLinear predictor units per point   :',format(1/w),'\n\n')

fun <- FALSE
for(x in names(obj)) {

  k <- x=='total.points' || x=='lp' || x=='abbrev'
  if(k) { fun <- TRUE; next }
  y <- obj[[x]]
  if(fun) {
    z <- cbind(round(y[[1]],dec), y$x.real)
    dimnames(z) <- list(rep('',nrow(z)), c('Total Points',x))
  } else {
    if(.R.) {
      z <- cbind(format(y[[1]]), format(round(y$points,dec)))
      dimnames(z)  <- list(rep('',length(y$points)), c(x, 'Points'))
    } else {
      z <- list(x=y[[1]], Points=round(y$points,dec))
      names(z) <- c(x,'Points')
      attr(z,'row.names') <- rep('',length(y$points))
      attr(z,'class') <- 'data.frame'
    }
    ## didn't use data.frame since wanted blank row names
  }
  cat('\n')
  if(.R.) print(z, quote=FALSE) else print(z)
  cat('\n')
}
invisible()
}

legend.nomabbrev <- function(object, which, x, y=NULL, ncol=3, ...) {
abb <- object$abbrev[[which]]
if(length(abb)==0) stop(paste('no abbreviation information for',which))
if(max(nchar(abb$abbrev))==1) if(length(y)) legend(x, y, abb$full, ncol=ncol,
                                     pch=paste(abb$abbrev,collapse=''), ...)
	else legend(x, abb$full, ncol=ncol, pch=paste(abb$abbrev,collapse=''),
		 ...)
else if(length(y)) legend(x, y, paste(format(abb$abbrev),':',abb$full,sep=''), 
		ncol=ncol, ...) else
	legend(x, paste(format(abb$abbrev),':',abb$full,sep=''), ncol=ncol, 
            ...)
invisible()
}


##Version of axis allowing tick mark labels to be forced, or to
##label every 'label.every' tick marks

axisf <- function(side, at, labels=TRUE, pos, cex, tck, 
				  label.every=1, force.label=FALSE, disc=FALSE) {

  ax <- if(.R.) function(..., cex) axis(..., cex.axis=cex) else
		        function(..., cex) axis(..., cex=cex)
		
  ax(side, at, labels=FALSE, pos=pos, cex=cex, tck=tck)

  if(is.logical(labels) && !labels) return(invisible())

  if(label.every>1 && !disc) {
	sq <- seq(along=at, by=label.every)
	at[-sq] <- NA
  }
  if(is.logical(labels)) labels <- format(at)

  if(force.label) {
    for(i in 1:length(labels))
	if(!is.na(at[i])) ax(side, at[i], labels[i], pos=pos, cex=cex, tck=tck)
    }
  else ax(side, at[!is.na(at)], labels[!is.na(at)], 
          pos=pos, cex=cex, tck=tck)

  invisible()
}
