latex.psm <- function(object,  title,
   file=paste(first.word(deparse(substitute(object))),".tex",sep=""),
   append=FALSE, which, varnames, 
   columns=65, inline=FALSE, 
   before=if(inline)"" else "& &",pretrans=TRUE, caption=NULL, ...) {

  f <- object
  
w <- if(length(caption)) paste('\\begin{center} \\bf',caption,'\\end{center}')

if(missing(which) & !inline)			{
  if(.SV4. || .R.) {
    dist <- f$dist
    w <- c(w, paste("\\[{\\rm Prob}\\{T\\geq t\\} = ",
	survreg.auxinfo[[dist]]$latex(f$scale),
	"{\\rm \\ \\ where} \\\\ \\]",sep=""))
  } else {
  fam <- f$family[1:2]
  dist <- fam[1]
  transform <- fam[2]
  w <- c(w,paste("\\[{\\rm Prob}\\{T\\geq t\\} = ",
	survreg.auxinfo[[dist]]$latex(f$parms, transform),
	"{\\rm \\ \\ where} \\\\ \\]",sep=""))
}
}
atr <- f$Design
if(!length(atr)) atr <- getOldDesign(f)

if(missing(which)) which <- 1:length(atr$name)
if(missing(varnames)) varnames <- atr$name[atr$assume.code!=9]

cat(w, sep=if(length(w)) "\n" else "", file=file, append=append)
latexDesign(f, file=file, append=TRUE, which=which,
            varnames=varnames, columns=columns, 
            before=before,
            prefix=if(missing(which))"X\\hat{\\beta}" else NULL, 
            inline=inline,pretrans=pretrans)  ## 4Dec00
}


