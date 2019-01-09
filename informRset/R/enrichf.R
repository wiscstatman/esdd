enrichf <-
function(x, y, quantile)
{
  ## assign min value to NA element 
  x[is.na(x)] <- min(x, na.rm=TRUE)
  # x: theta value
  # y: label
  # quantile: quantile paramenter for enrichment factor 
  #  ok <- !is.na(x)
  # x <- x[ok]
  # y <- y[ok]
  
  if(length(x)!= length(y)) {
    stop(paste("The length of theta should be equal to the length of labels"))
  }
  N <- length(y) ## num of ligands
  n <- sum(y) ## num of actives
  
  ord <- order(x, decreasing = TRUE, na.last=TRUE) ## rank theta by decreasing order
  yard <- cumsum(y[ord])  
  pos <- round(N * quantile)
  
  if(n==0){
    stop(paste("number of active ligands should be positive"))
  }else{  
    ef <- yard[pos]/(n*quantile)
  }
  
  return(ef)
}
