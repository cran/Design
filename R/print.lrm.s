print.lrm <- function(x, digits=4, strata.coefs=FALSE, ...) {
sg <- function(x,d) 	{
#  .Options$digits <- d  14Sep00
  oldopt <- options(digits=d)
  on.exit(options(oldopt))
	format(x)	}
rn <- function(x,d) format(round(as.single(x),d))

cat("\n")
if(x$fail)	{
	cat("Model Did Not Converge\n")
	return()
		}

cat("Logistic Regression Model\n\n")
dput(x$call)
cat("\n\nFrequencies of Responses\n")
print(x$freq)
if(length(x$sumwty)) {
  cat('\n\nSum of Weights by Response Category\n')
  print(x$sumwty)
}
cat("\n")
if(!is.null(x$nmiss))	{  #for backward compatibility
   cat("Frequencies of Missing Values Due to Each Variable\n")
   print(x$nmiss)
   cat("\n")		}
else if(!is.null(x$na.action)) naprint(x$na.action)

ns <- x$non.slopes
nstrata <- x$nstrata
if(!length(nstrata)) nstrata <- 1

pm <- x$penalty.matrix
if(length(pm)) {
   psc <- if(length(pm)==1) sqrt(pm) else
	sqrt(diag(pm))
   penalty.scale <- c(rep(0,ns),psc)
   cof <- matrix(x$coef[-(1:ns)], ncol=1)
   cat("Penalty factors:\n\n"); print(as.data.frame(x$penalty, row.names=''))
   cat("\nFinal penalty on -2 log L:",
	rn(t(cof) %*% pm %*% cof,2),"\n\n")
}

#est.exp <- 1:ns
#if(length(f$est)) est.exp <- c(est.exp, ns+f$est[f$est+ns <= length(f$coef)])
vv <- diag(x$var)
cof <- x$coef
if(strata.coefs) {
  cof <- c(cof, x$strata.coef)
  vv  <- c(vv,  x$Varcov(x,which='strata.var.diag'))
  if(length(pm)) penalty.scale <- c(penalty.scale,rep(NA,x$nstrat-1))
}
score.there <- nstrata==1 && (length(x$est) < length(x$coef)-ns)
stats <- x$stats
stats[2] <- signif(stats[2],1)
stats[3] <- round(stats[3],2)
stats[4] <- round(stats[4],2)
stats[5] <- round(stats[5],4)
stats[6] <- round(stats[6],3)
stats[7] <- round(stats[7],3)
if(nstrata==1) { ##17Dec97
  stats[8] <- round(stats[8],3)   ##21Aug97
  stats[9] <- round(stats[9],3)
  stats[10] <- round(stats[10],3)
  if(length(stats)>10) {
    stats[11] <- round(stats[11],3)
    if(length(x$weights)) stats[12] <- round(stats[12],3)
  }
} else stats <- c(stats,Strata=x$nstrat)

if(.R.) {   ## 8Apr02
  nst <- length(stats)
  cstats <- character(nst)
  names(cstats) <- names(stats)
  for(i in 1:nst) cstats[i] <- format(stats[i])
  print(cstats, quote=FALSE)
} else if(!score.there) print(stats)	else	{
	print(stats[1:10])
	cat("\n")
	st <- stats[11:13]
	st[1] <- round(st[1],2)
	st[3] <- round(st[3],4)
	print(st)			}
cat("\n")

##if(length(f$var)==0) vv <- NULL	#doesn't bother with this for x=NULL
z <- cof/sqrt(vv)
stats <- cbind(sg(cof,digits), sg(sqrt(vv),digits), 
	rn(cof/sqrt(vv),2))
stats <- cbind(stats, rn(1-pchisq(z^2,1),4))
dimnames(stats) <- list(names(cof),
	c("Coef","S.E.","Wald Z","P"))
if(length(pm)) stats <- cbind(stats, "Penalty Scale"=sg(penalty.scale,digits))
print(stats,quote=FALSE)
cat("\n")


if(score.there)	{
	q <- (1:length(cof))[-est.exp]
	if(length(q)==1) vv <- x$var[q,q] else vv <- diag(x$var[q,q])
	z <- x$u[q]/sqrt(vv)
	stats <- cbind(rn(z,2), rn(1-pchisq(z^2,1),4))
	dimnames(stats) <- list(names(cof[q]),c("Score Z","P"))
	print(stats,quote=FALSE)
        cat("\n")
}
invisible()

}

