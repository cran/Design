print.cph <- function(x, long=FALSE, digits=3, conf.int=FALSE,
                      table=TRUE,  ...) { 

cat("\n")
if(x$fail)	{
	cat("Model Did Not Converge\n")
	return()
		}

cat("Cox Proportional Hazards Model\n\n")
dput(x$call)
cat("\n")
if(!is.null(z <- x$na.action)) naprint(z)
if(!is.null(x$coef))						{
   stats <- x$stats
   stats[3] <- round(stats[3],2)
   stats[5] <- round(stats[5],4)
   stats[6] <- round(stats[6],2)
   stats[7] <- round(stats[7],4)
   stats[8] <- round(stats[8],3)
   if(.R.) print(format.sep(stats), quote=FALSE) else print(stats)
   cat("\n")
   print.cph.fit(x, digits=digits, conf.int=conf.int, table=table, ...)
   if(long)cat("Centering constant:",format(x$center),"\n")
 }
else if(table) print(x$n)
invisible()
									}

