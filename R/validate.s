validate <-
  function(fit,  method="boot", B=40,
           bw=FALSE, rule="aic", type="residual", sls=0.05, aics=0, 
           pr=FALSE,...) UseMethod("validate")

