thredec <-
function(u, informer, newinh)
{
  source("inbin.R")
  thresh <- inbin(u)$thresh
  M1 <- glm(thresh ~ u[,informer], family="binomial")
  
  ## generating newthreshold
  newthresh <-sum( coef(M1) * c(1,newinh))
  newthresh <- exp(newthresh)/(exp(newthresh)+1)
  return(newthresh)
  
}
