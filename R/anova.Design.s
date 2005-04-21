#main.effect=F to suppress printing main effects when the factor in
#question is involved in any interaction.

anova.Design <- function(object,...,main.effect=FALSE, tol=1e-9, 
						 test=c('F','Chisq'), ss=TRUE)	{

ava <- function(idx,coef,cov,tol) {
	chisq <- coef[idx] %*% solvet(cov[idx,idx], coef[idx], tol=tol)
	c(chisq, length(idx))
  }

obj.name <- as.character(sys.call())[2]
itype <- 1	#Wald stats. Later sense score stats from object$est
misstest <- missing(test)   ## R updates missing  8Apr02
test <- match.arg(test)
is.ols <- inherits(object,'ols') ||
 (length(object$fitFunction) && any(object$fitFunction=='ols')) ##14Nov00 22May01
if(misstest) test <- if(is.ols) 'F' else 'Chisq'
if(!is.ols && test=='F') stop('F-test not allowed for this type of model')
if(!is.ols) ss <- FALSE

at <- object$Design
if(!length(at)) at <- getOldDesign(object)
assign <- object$assign
name <- at$name
nama <- names(assign)[1]
asso <- 1*(nama=="(Intercept)" | nama=="Intercept")
names(assign)[-asso] <- name

ia <- at$interactions
if(!length(ia))nia <- 0 else nia <- ncol(ia)
assume <- at$assume.code
#if(is.null(assume))stop("fit does not have Design information")
parms <- at$parms
f <- length(assume)
dotlist <- if(!.SV4. && !.R.) (((sys.frame())[["..."]])[[1]])[-1] else
{  ## 12Nov00
  ncall <- names(sys.call())[-(1:2)]
  other.arg <- as.character(sys.call())[-(1:2)]
  if(length(other.arg) && length(ncall))
    other.arg <- other.arg[ncall=='']
  other.arg
}
#  (as.character(sys.call()[-1])[-1])[if(length(ncall))ncall=='' else TRUE]

if(length(dotlist)==0) which <- 1:f	else		{
  if(!.SV4. && !.R.) {
	alist <- NULL
	for(i in 1:length(dotlist))
		alist <- c(alist, deparse(dotlist[[i]]))
  } else alist <- dotlist ## 12Nov00
	jw <- charmatch(alist,name,0)
	if(any(jw==0))	stop(paste("factor names not in design: ",
			paste(alist[jw==0],collapse=" ")))
	which <- jw	}

if(length(object$est) && !length(object$u))
	stop("est in fit indicates score statistics but no u in fit")

if(itype==1)	{
	if(!length(object$coefficients))
	   stop("estimates not available for Wald statistics")
	coef <- object$coefficients	}	else
	{
	if(!length(object$u)) stop("score statistics not available")
#	if(pr)cat("\n                     Score Statistics\n\n")
	coef <- object$u	}
np <- length(coef)

#Compute # intercepts to skip in testing
nrp <- num.intercepts(object)
if(itype==2 & nrp!=0)stop("fit score statistics and x are incompatible")

nc <- length(coef)

cov <- Varcov(object, regcoef.only=TRUE)  #Omit row/col for scale parameters

stats <- NULL
lab <- NULL
W <- list()
s <- 0
all.slopes <- rep(FALSE, nc)
all.ia <- rep(FALSE, nc)
all.nonlin <- rep(FALSE, nc)
num.ia <- 0
num.nonlin <- 0
issue.warn <- FALSE

for(i in which)	{
	j <- assume[i]
	parmi <- parms[[name[i]]]
	if(j!=9) low.fact <- i
	else low.fact <- (parmi[,1])[parmi[,1]>0]
	if(!length(names(at$nonlinear))) nl <- at$nonlinear[[i]]
	else nl <- at$nonlinear[[name[i]]]
	if(!length(nl)) nl <- rep(FALSE,length(assign[[name[i]]]))
#Factor no. according to model matrix is 1 + number of non-strata factors
#before this factor
       if(j!=8)	{				#ignore strata
	if(i==1) jfact <- 1
	else jfact <- 1 + sum(assume[1:(i-1)]!=8)
	main.index <- assign[[jfact+asso]]
	nonlin.ia.index <- NULL	#Should not have to be here. Bug in S?
	all.slopes[main.index] <- TRUE
	if(nia==0)ni <- 0 else ni <- sum(ia==i)
	if(nia==0)ni <- 0
	else for(k in 1:ncol(ia)) ni <- ni + !any(is.na(match(low.fact,ia[,k])))
	if(ni==0 | main.effect)			{
		w <- ava(main.index,coef,cov,tol=tol)
        s <- s+1; W[[s]] <- main.index
		stats <- rbind(stats,w)
		lab <- c(lab, name[i])		}
	#If term is involved in any higher order effect, get pooled test
	#by adding in all high-order effects containing this term
	#For 2nd order interaction, look for 3rd order interactions
	#containing both factors
#	nonlin.ia.index <- NULL	#Used to be here.  Bug in S?
	if(ni>0)	{
		ia.index <- NULL
		mm <- (1:f)[assume==9]
		mm <- mm[mm!=i]
		for(k in mm)		{

			parmk <- parms[[name[k]]]
			hi.fact <- parmk[,1]
			m <- match(low.fact, hi.fact)
			if(!any(is.na(m)))		{
				if(k==1)kfact <- 1
				else kfact <- 1 + sum(assume[1:(k-1)]!=8)
				idx <- assign[[kfact+asso]]
				ia.index <- c(ia.index,idx)
				if(ncol(parmk)>1)for(jj in 1:length(m))	{
				  nonlin.ia.index <- c(nonlin.ia.index,
				     idx[parmk[m[jj],-1]==1])		}
				  nonlin.ia.index <- if(length(nonlin.ia.index))
					unique(nonlin.ia.index) else NULL
				#Highest order can be counted twice
#				c(nonlin.ia.index, added 17 Sep 91
							}
									 }
			idx <- c(main.index,ia.index)
			all.slopes[idx] <- TRUE
			w <- ava(idx,coef,cov,tol=tol)
            s <- s+1; W[[s]] <- idx
			stats <- rbind(stats,w)
			lab <- c(lab, paste(name[i], 
			   " (Factor+Higher Order Factors)"))
			#If factor i in >1 interaction, print summary
			#Otherwise, will be printed later
			if(j!=9 & ni>1)	{
				w <- ava(ia.index,coef,cov,tol=tol)
                s <- s+1; W[[s]] <- ia.index
				stats<-rbind(stats,w)
				lab <- c(lab, " All Interactions")
									}
					    				 }
#	if((any(nl) & j!=9) | (j==9 && parmi[3,i]==1))	{
	if(any(nl))					{
		# Tests of adequacy of linear relationship
		idx <- c(main.index[nl], nonlin.ia.index)
		num.nonlin <- num.nonlin+1
		all.nonlin[idx] <- TRUE
		w <- ava(idx,coef,cov,tol=tol)
        s <- s+1; W[[s]] <- idx
		stats <- rbind(stats,w)
		lab <- c(lab, if(!length(nonlin.ia.index))" Nonlinear"
			else " Nonlinear (Factor+Higher Order Factors)")	
							} 
		#If interaction factor involves a non-linear term from an
		#expanded polynomial, lspline, rcspline, or scored factor,
		#do tests to see if a simplification (linear interaction) is
		#adequate.  Do for second order only.
		if(j==9)	{
			num.ia <- num.ia+1
			all.ia[main.index] <- TRUE
			if(parmi[3,1]>0) issue.warn <- TRUE
			if(parmi[3,1]==0 && ncol(parmi)>1)		{
			nonlin.x <- as.logical(parmi[1,2:ncol(parmi)])
			nonlin.y <- as.logical(parmi[2,2:ncol(parmi)])
			nonlin.xy <- nonlin.x | nonlin.y
			nonlin.xandy <- nonlin.x & nonlin.y
			idx <- main.index[nonlin.xy]
			li <- length(idx)
			if(li>0)	{
			  num.nonlin <- num.nonlin+1
			  all.nonlin[idx] <- TRUE
			  w <- ava(idx,coef,cov,tol=tol)
              s <- s+1; W[[s]] <- idx
			  stats<-rbind(stats,w)
			  lab<-c(lab," Nonlinear Interaction : f(A,B) vs. AB")
			  idx <- main.index[nonlin.xandy]
			  li <- length(idx)
			  if(li>0)	{
			  w <- ava(idx,coef,cov,tol=tol)
              s <- s+1; W[[s]] <- idx
			  stats<-rbind(stats,w)
			  lab<-c(lab," f(A,B) vs. Af(B) + Bg(A)")	}
			  idx <- main.index[nonlin.x]
			  li <- length(idx)
			  if(li>0 & any(nonlin.x!=nonlin.xy)) {
			    w <- ava(idx,coef,cov,tol=tol)
                s <- s+1; W[[s]] <- idx
			    stats<-rbind(stats,w)
			    lab<-c(lab,paste(" Nonlinear Interaction in",
			         name[parmi[1,1]],"vs. Af(B)"))	}
			    idx <- main.index[nonlin.y]
			    li <- length(idx)
			  if(li>0 & any(nonlin.y!=nonlin.xy)) {
			    w <- ava(idx,coef,cov,tol=tol)
                s <- s+1; W[[s]] <- idx
			    stats<-rbind(stats,w)
			    lab<-c(lab,paste(" Nonlinear Interaction in",
				 name[parmi[2,1]],"vs. Bg(A)"))	}
									}}  }
									  } }

#If >1 test of adequacy, print pooled test of all nonlinear effects
if(num.nonlin>1)					{
	idx <- (1:nc)[all.nonlin]
	li <- length(idx)
	w <- ava(idx,coef,cov,tol=tol)
    s <- s+1; W[[s]] <- idx
	stats <- rbind(stats,w)
	lab <- c(lab, "TOTAL NONLINEAR")			}
#If >1 test of interaction, print pooled test of all interactions in list
if(num.ia>1)						{
	idx <- (1:nc)[all.ia]
	li <- length(idx)
	w <- ava(idx,coef,cov,tol=tol)
    s <- s+1; W[[s]] <- idx
	stats <- rbind(stats,w)
	lab <- c(lab,"TOTAL INTERACTION")			}
#If >0 test of adequacy and >0 test of interaction, print pooled test of
#all nonlinear and interaction terms
if(num.nonlin>0 & num.ia>0)	{
	idx <- (1:nc)[all.nonlin | all.ia]
	li <- length(idx)
	w <- ava(idx,coef,cov,tol=tol)
    s <- s+1; W[[s]] <- idx
	stats <- rbind(stats,w)
	lab <- c(lab,"TOTAL NONLINEAR + INTERACTION")
				}
#Get total test for all factors listed
	idx <- (1:nc)[all.slopes | all.ia]
	w <- ava(idx,coef,cov,tol=tol)
s <- s+1; W[[s]] <- idx
stats <- rbind(stats,w)
lab <- c(lab,"TOTAL")

statnam <- c('Chi-Square','d.f.')
if(is.ols) {
  sigma2 <- object$stats['Sigma']^2
  dfe <- object$df.residual
}

if(ss) {
  stats <- cbind(stats[,2], stats[,1]*sigma2, stats[,1]*sigma2/stats[,2], 
				 stats[,1])
  statnam <- c('d.f.','Partial SS','MS','Chi-Square')
  stats <- rbind(stats, Error=c(dfe, sigma2*dfe, sigma2, NA))
  s <- s+1; W[[s]] <- NA
  lab <- c(lab, 'ERROR')
}

j <- statnam=='Chi-Square'
dfreg <- stats[,statnam=='d.f.']
if(test=='F') {
  stats[,j] <- stats[,j] / dfreg
  statnam[j] <- 'F'
  stats <- cbind(stats, P=1-pf(stats[,j], dfreg, dfe))
  attr(stats,'df.residual') <- dfe
} else stats <- cbind(stats,1-pchisq(stats[,j], dfreg))

statnam <- c(statnam, 'P')
dimnames(stats) <- list(lab, statnam)
## attr(stats,"formula") <- formula(object$terms) 30may02
attr(stats,'formula') <- formula(object)
## was attr(object$terms,"formula") 17Apr02
attr(stats,"obj.name") <- obj.name
attr(stats,"class") <- if(.SV4.)'anova.Design' else c("anova.Design","matrix")
names(W) <- lab
attr(stats,"which") <- W
attr(stats,"coef.names") <- names(coef)
attr(stats,"non.slopes") <- nrp
if(issue.warn) 
 warning("tests of nonlinear interaction with respect to single component \nvariables ignore 3-way interactions")
stats
}


print.anova.Design <-
  function(x, which=c('none','subscripts','names','dots'), ...) {
  stats <- x
  digits <- c('Chi-Square'=2, F=2, 'd.f.'=0, 'Partial SS'=15, MS=15, P=4)
  cstats <- matrix('', nrow=nrow(stats), ncol=ncol(stats), 
                   dimnames=dimnames(stats))

which <- match.arg(which)

do.which <- which!='none' && length(W <- attr(stats,'which'))
if(do.which) {
  if(which=='subscripts') simplifyr <- function(x) {
    x <- sort(unique(x))
    n <- length(x)
    ranges <- character(n)
    m <- 0
    s <- x
    while(length(s) > 0) {
      j <- s == s[1] + (1:length(s))-1
      m <- m+1
      ranges[m] <- if(sum(j)>1) paste(range(s[j]),collapse='-') else s[1]
      s <- s[!j]
    }
    ranges[1:m]
  }
  
  k <- length(W)
  w <- character(k)
  coef.names <- attr(stats,'coef.names')
  nrp <- attr(stats,'non.slopes')
  for(i in 1:k) {
    z <- W[[i]]
    if(all(is.na(z))) w[i] <- '' else {
      z <- sort(z)
      w[i] <- switch(which,
                     subscripts=paste(simplifyr(z - nrp), collapse=','),
                     names=paste(coef.names[z],collapse=','),
                     dots={
                       dots <- rep(' ',length(coef.names)-nrp)
                       dots[z - nrp] <- '.'
                       paste(dots,collapse='')} )
    }
  }
}
sn <- dimnames(cstats)[[2]]
for(j in 1:ncol(cstats)) cstats[,j] <- format(round(stats[,j], digits[sn[j]]))

cstats[is.na(stats)] <- ''
j <- sn=='P'
cstats[stats[,j] < 0.00005,j] <- '<.0001'
cstats <- cbind(dimnames(stats)[[1]], cstats)

#cstats<-cbind(dimnames(stats)[[1]],format(round(stats[,1],2)),
#	format(stats[,2]),format(round(stats[,3],4)))
dimnames(cstats) <- list(rep("",nrow(stats)),
	c("Factor    ",dimnames(stats)[[2]]))

heading <- paste("                ",
				 if(any(dimnames(stats)[[2]]=='F'))"Analysis of Variance" else
				 "Wald Statistics", "          Response: ", 
		as.character(attr(stats, "formula")[2]), sep = "")
cat(heading,"\n\n")
if(any(sn=='MS')) cstats[cstats[,1]=='TOTAL',1] <- 'REGRESSION'
if(do.which) cstats <- cbind(cstats, Tested=w)
print(cstats,quote=FALSE)
if(do.which && which!='names') {
  cat('\nSubscripts correspond to:\n')
  print(if(nrp > 0)coef.names[-(1:nrp)] else coef.names, quote=FALSE)
}
if(!any(sn=='MS') && length(dfe <- attr(stats,'df.residual'))) 
  cat('\nError d.f.:', dfe, '\n')
invisible()
}

latex.anova.Design <-
 function(object,
          title=if(under.unix)
           paste('anova',attr(object,'obj.name'),sep='.') else
           paste("ano",substring(first.word(attr(object,"obj.name")),
                                 1,5),sep=""), 
           psmall=TRUE,
           dec.chisq=2, dec.F=2, dec.ss=NA, dec.ms=NA, dec.P=4, ...) {

## expr in first.word 18Nov00 removed 25May01

rowl <- dimnames(object)[[1]]
#Translate interaction symbol (*) to times symbol
rowl <- sedit(rowl, "*", "$\\times$", wild.literal=TRUE)
#Put TOTAL rows in boldface
rowl <- ifelse(substring(rowl,1,5) %in% c("TOTAL","ERROR"), paste("{\\bf",rowl,"}"),rowl)
rowl <- ifelse(substring(rowl,1,1)==" ",
	paste("~~{\\it ",substring(rowl,2),"}",sep=""), rowl) # preserve leading blank
P <- object[,3]

dstats <- as.data.frame(object)
attr(dstats, 'row.names') <- rowl

if(psmall) {
  psml <- !is.na(dstats$P) & dstats$P < 0.00005
  if(any(psml)) dstats$P <- ifelse(is.na(dstats$P),'',ifelse(psml, 
				"$<0.0001$",
				paste("~",format(round(dstats$P,dec.P)),sep="")))
}

digits <- c('Chi-Square'=dec.chisq, F=dec.F, 'd.f.'=0,
            'Partial SS'=dec.ss, MS=dec.ms, P=dec.P)

sn <- dimnames(object)[[2]]
dig <- digits[sn]
sn[sn=='Chi-Square'] <- '\\chi^2'
names(dstats) <- paste('$',sn,'$',sep='')

#Make LaTeX preserve spaces in heading
head <- paste(if(any(sn=='F'))"Analysis of Variance" else "Wald Statistics", "for {\\tt",
	as.character(attr(object,"formula")[2]),"}")
latex(dstats, cdec=dig, title=title, caption=head, rowlabel="",
	col.just=rep('r',length(sn)), ...)
}

text.anova.Design <- function(x, at, cex=.5, font=2, ...)	{

#Note: a bug in text() prevents writing long character strings
ltext <- function(z, line, label, cex = 0.5, font=2, adj = 0) {
  zz <- z
  zz$y <- z$y - ((line - 1) * 1.2 * cex * par("csi") * (
	par("usr")[4] - par("usr")[3]))/(par("fin")[2])
  text(zz, label, cex = cex, adj = adj, font=font)
}

fi <- tempfile()
sink(fi)
print.anova.Design(x)
sink()

k <- if(.R.) scan(fi, list(z=""), sep="\n", quiet=TRUE)$z else
             scan(fi, list(z=""), sep="\n")$z
if(!.R. && existsFunction('unlink')) unlink(fi)

for(l in 1:length(k)) ltext(at, l, k[l], font=font, cex=cex)

invisible(k)
}

plot.anova.Design <- function(x,
    what=c("chisqminusdf","chisq","aic","P","partial R2","remaining R2",
      "proportion R2"),
	xlab=NULL,
	pch=16, rm.totals=TRUE, rm.ia=FALSE, rm.other=NULL, newnames,
	sort=c("descending","ascending","none"), pl=TRUE, ...) {

what <- match.arg(what)
sort <- match.arg(sort)

if(!length(xlab)) xlab <-
  switch(what, chisq=if(.R.)expression(chi^2) else "Chi-square", 
         chisqminusdf=if(.R.)expression(chi^2~-~df) else
         "Chi-Square Minus Degrees of Freedom", 
         aic="Akaike Information Criterion",
         P="P-value",
         "partial R2"=if(.R.)expression(paste("Partial",~R^2)) else "Partial R^2",
         "remaining R2"=if(.R.)expression(paste("Remaining~",R^2,
             "~After Removing Variable")) else
         "Remaining R^2 After Removing Variable",
         "proportion R2"=
         if(.R.)expression(paste("Proportion of Overall",~R^2))
         else "Proportion of Overall R^2")

if(.SV4.) x <- matrix(oldUnclass(x), nrow=nrow(x),
                          dimnames=dimnames(x))  ##14Nov00

rm <- c(if(rm.totals) c("TOTAL NONLINEAR","TOTAL NONLINEAR + INTERACTION",
	"TOTAL INTERACTION","TOTAL"), 
        " Nonlinear"," All Interactions", "ERROR", rm.other)
rn <- dimnames(x)[[1]]
rm <- c(rm, rn[substring(rn,2,10)=="Nonlinear"])
k <- !(rn %in% rm)
if(rm.ia) k[grep("\\*", rn)] <- FALSE
an <- x[k,,drop=FALSE]

dof <- an[,'d.f.']
P <- an[,'P']
chisq <- if(any(dimnames(an)[[2]]=='F')) an[,'F']*dof else an[,'Chi-Square']

if(what %in% c("partial R2","remaining R2","proportion R2")) {
  if("Partial SS" %nin% dimnames(x)[[2]])
    stop('to plot R2 you must have an ols model and must not have specified ss=F to anova')
  sse <- x['ERROR','Partial SS']
  ssr <- x['TOTAL','Partial SS']
  sst <- sse + ssr
}

an <- switch(what,
             chisq=chisq,
             chisqminusdf=chisq-dof,
             aic=chisq-2*dof,
             P=P,
             "partial R2" = an[,"Partial SS"]/sst,
             "remaining R2" = (ssr - an[,"Partial SS"]) / sst,
             "proportion R2" = an[,"Partial SS"] / ssr)

if(missing(newnames)) newnames <- sedit(names(an),
	"  (Factor+Higher Order Factors)", "")
names(an) <- newnames
an <- switch(sort, descending=-sort(-an), ascending=sort(an), none=an)

if(pl) dotchart2(an, xlab=xlab, pch=pch, ...)
invisible(an)
}

