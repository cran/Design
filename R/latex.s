if(FALSE) latex <-  function(x, ...)   # duplicates what's in print.display
{
  if(is.null(oldClass(x)))
    oldClass(x) <- data.class(x)
  UseMethod("latex")
}
