inbin <-
function(u)
{
  xb <- apply(u, 1, mean, na.rm=TRUE)
  ss <- apply(u, 1, sd, na.rm=TRUE)
  
  yard <- xb + 2*ss
  yard[yard>1] <- 1
  y <- 1* ((u - yard)>0)
  
  return(list(y=y, thresh=yard))
}
