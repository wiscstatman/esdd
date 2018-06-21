upmrval <-
function(scores, alpha=c(seq(0.0001, 0.0099, length.out = 80), seq(0.01, 0.99, length.out=20)), n0=8, x0=max(scores, na.rm=TRUE)){
  source("upmtail.R")
  source("Vupmtail.R")
  ## generate summary data
  xtemp <- apply(scores, 1, max, na.rm=TRUE)
  ntemp <- apply(scores, 1, function(x){sum(!is.na(x))})
  
  ## slow but eay to code
  post <- apply(quants, 1, Vtailvol, alpha=alpha, n=n0, x0=x0)
  
  ## generate rvalue
  yard <- apply(post, 1, quantile, probs=1-alpha)
  yard <- diag(yard)
  
  tyard <- post - yard
  
  ## function to get intersection of twe curves
  intersect <- function(x)
  {
    m <- length(x)
    
    base <- x[1]
    if(base>=0){
      rank <- 1
      return(rank)
    }
    
    for(jj in 2:m){
      compare <- x[jj]
      if(base<0 && compare>=0){
        rank <- jj
        return(rank)
      }else{
        base <- compare
      }
    }
    rank <- length(x)
    return(rank)
  }
  
  
  rank <- apply(tyard, 2, intersect)
  theta <- alpha[rank]
  return(theta)
}
