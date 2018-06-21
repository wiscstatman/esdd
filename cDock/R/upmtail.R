upmtail <-
function(svec, alpha=0.01, n0, x0)
{
  ytemp <- max(svec, na.rm=TRUE)
  ntemp <- sum(!is.na(svec))
  
  thetalpha <- x0 * exp(-log(alpha)/n0)
  if(ytemp >= thetalpha)
  {
    return(1)
  }else{
    tailtemp <- (ytemp/thetalpha)^(n0+ntemp)
    return(tailtemp)
  }
}
