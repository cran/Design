latex.pphsm <- function(object, title,
    file=paste(first.word(deparse(substitute(object))),".tex",sep=""),
    append=FALSE, which, varnames, 
    columns=65, inline=FALSE, 
    before=if(inline)"" else "& &",pretrans=TRUE, caption=NULL, ...) {

w <- if(length(caption)) paste('\\begin{center} \\bf',caption,'\\end{center}')

sc <- exp(f$parms)
at <- f$Design
if(!length(at)) at <- getOldDesign(f)

if(missing(which) & !inline)			{
  dist <- paste("\\exp\\{-t^{",format(1/sc),"} \\exp(X\\hat{\\beta})\\}")
  w <- c(w,paste("\\[{\\rm Prob}\\{T\\geq t\\} = ",dist,
	"{\\rm \\ \\ where} \\\\ \\]",sep=""))
						}				
if(missing(which)) which <- 1:length(at$name)
if(missing(varnames)) varnames <- at$name[at$assume.code!=9]
cat(w, file=file, sep=if(length(w))"\n" else "", append=append)
latexDesign(f, file=file, append=TRUE, which=which, varnames=varnames, 
	columns=columns, 
	before=before, prefix=if(missing(which))"X\\hat{\\beta}" else NULL, 
	inline=inline,pretrans=pretrans)  ## 4Dec00
}


