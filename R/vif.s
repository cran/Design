vif <- function(fit) {

v <- Varcov(fit, regcoef.only=TRUE)
nam <- dimnames(v)[[1]]
ns <- num.intercepts(fit)
#  v <- solve(v)  this should never have been there   06oct03
if(ns>0) {v <- v[-(1:ns),-(1:ns),drop=FALSE]; nam <- nam[-(1:ns)]}
d <- diag(v)^.5
v <- diag(solve(v/(d %o% d)))
names(v) <- nam
v

}
