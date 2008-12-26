robcov <- function(fit, cluster, method=c('huber','efron')) {

  method <- match.arg(method)

var <- fit$var
vname <- dimnames(var)[[1]]

if(inherits(fit, "ols") ||
   (length(fit$fitFunction) && any(fit$fitFunction=='ols'))) {
   ##14Nov00 22May01
   var <- fit$df.residual * var/sum(fit$residuals^2)  #back to X'X
##   warning("printing the fit object from robcov (with print.ols) will not print the adjusted std. err.\nUse sqrt(diag(fit$var)) to get adjusted std. err.,\nwhere fit is the result of robcov.") 22Dec01
 } else if(method=='efron') stop('method="efron" only works for ols fits')

X <- as.matrix(residuals(fit, type=if(method=='huber')"score" else "hscore"))

n <- nrow(X)
if(missing(cluster)) cluster <- 1:n
else if(any(is.na(cluster))) stop("cluster contains NAs")
if(length(cluster)!=n)
 stop("length of cluster does not match number of observations used in fit")
cluster <- as.factor(cluster)

p <- ncol(var)
j <- is.na(X %*% rep(1, p))
if(any(j))	{
  X <- X[!j,,drop=FALSE]    # 12Apr02
  cluster <- cluster[!j]
  n <- length(cluster)
		}

j <- order(cluster)
X <- X[j,,drop=FALSE]
clus.size <- table(cluster)
clus.start <- c(1,1+cumsum(clus.size))
nc <- length(levels(cluster))
clus.start <- clus.start[-(nc+1)]
storage.mode(clus.start) <- "integer"

#dyn.load("/suserlib/robcovf.o")
W <- matrix(if(.R.)
            .Fortran("robcovf", n, p, nc, clus.start, clus.size, X, 
                     double(p), double(p*p), w=double(p*p),
                     PACKAGE="Design")$w else
            .Fortran("robcovf", n, p, nc, clus.start, clus.size, X, 
                     double(p), double(p*p), w=double(p*p),
                     PACKAGE="Design")$w,
            nrow=p)
#The following has a small bug but comes close to reproducing what robcovf does
#W <- tapply(X,list(cluster[row(X)],col(X)),sum)
#W <- t(W) %*% W
#The following logic will also do it, also at great cost in CPU time
#for(j in levels(cluster))		{
#   s <- cluster==j
#   if(sum(s)==1) sx <- X[s,,drop=F]
#   else {sx <- apply(X[s,,drop=F], 2, sum); dim(sx) <- c(1,p)}
#
#   sc <- sc + t(sx) %*% sx
#
#					}

adjvar <- var %*% W %*% var

#var.new <- diag(adjvar)
#deff <- var.new/var.orig; names(deff) <- vname
#eff.n <- n/exp(mean(log(deff)))

#if(pr)			{
#   v <- cbind(var.orig, var.new, deff)
#   dimnames(v) <- list(vname, c("Original Variance","Adjusted Variance",
#			"Design Effect"))
#   .Options$digits <- 4
#   cat("\n\nEffect of Adjustment for Cluster Sampling on Variances of Parameter #Estimates\n\n")
#   print(v)
#   cat("\nEffective sample size:",format(round(eff.n,1)),"\n\n")
#   nn <- n^2/sum(clus.size^2)
#   cat("\nN^2/[sum of Ni^2]    :",format(round(nn,1)),"\n\n")
#			}

fit$orig.var <- fit$var
fit$var <- adjvar
#fit$design.effects <- deff
#fit$effective.n <- eff.n
#oldClass(fit) <- c("robcov",oldClass(fit))  ##14Nov00

fit

								}
