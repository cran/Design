# $Id: AFirst.lib.s,v 1.2 2004/07/11 14:39:08 harrelfe Exp $
#under.unix   <- !(version$os=='Microsoft Windows' ||
#                  version$os=='Win32' || version$os=='mingw32')
.R.          <- TRUE
.SV4.        <- FALSE
.newSurvival. <- TRUE

.noGenenerics <- TRUE  # faster loading as new methods not used

.First.lib <- function(lib, pkg, ...) {
  verbose <- .Options$Hverbose
  if(!length(verbose) || verbose)
    cat("Design library by Frank E Harrell Jr\n\n",
        "Type library(help='Design'), ?Overview, or ?Design.Overview')\n",
        "to see overall documentation.\n\n",
        sep='')
  library.dynam("Design", pkg, lib)
  require('Hmisc') || stop('Hmisc package not available')
  invisible()
}

#  requirePos('survival',pos=3)
#  p <- .packages()
#  if(match('Design',p) > match('survival',p))
#    warning('By not specifying library(Design,pos=2), functions in the\nsurvival package such as survfit will override those in Design.')
