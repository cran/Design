#This is Terry Therneau's old print.coxreg with conf.int default to F
#Add Nagelkerke R2 9Jun92
#Remove printing hazard ratios 17Jun92
#Removed stats 23Jun92 since in print.cph now, delete print x$n if 1 stratum
print.cph.fit <- 
function(x, table = TRUE, coef = TRUE, conf.int = FALSE, scale = 1,
         digits = NULL, ...)
{
	if(table && !is.null(x$n) && is.matrix(x$n))
		print(x$n)
	if(is.null(digits))
		digits <- 3
	savedig <- options(digits = digits)
	on.exit(options(savedig))
	beta <- x$coef
	se <- sqrt(diag(x$var))
	if(is.null(beta) | is.null(se))
		stop("Input is not valid")
	if(coef) {
		tmp <- cbind(beta, se, beta/se, 1 - pchisq((beta/
			se)^2, 1))
		dimnames(tmp) <- list(names(beta), c("coef", 
			"se(coef)", "z", "p"))
		cat("\n")
		prmatrix(tmp)
	}
	if(conf.int) {
		z <- qnorm((1 + conf.int)/2, 0, 1)
		beta <- beta * scale
		se <- se * scale
		tmp <- cbind(exp(beta), exp( - beta), exp(beta - z * se), exp(
			beta + z * se))
		dimnames(tmp) <- list(names(beta), c("exp(coef)", "exp(-coef)",
			paste("lower .", round(100 * conf.int, 2), sep = ""),
			paste("upper .", round(100 * conf.int, 2), sep = "")))
		cat("\n")
		prmatrix(tmp)
	}
	invisible(x)
}
