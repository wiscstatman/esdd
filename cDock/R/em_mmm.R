em_mmm <-
function( quants, pi0=.98 , deltaA=1.76, deltaB = 1)
{
  ## quants should contain no NA value
  testna <- sum(is.na(quants))

  cmax <- 13 ## EM max iteration
  
  X <- as.matrix(quants)
  n <- nrow(X) # number of ligands
  m <- ncol(X)
  
  xb <- apply( quants, 1, mean, na.rm=TRUE )
  muA <- apply(quants, 2, mean, na.rm=TRUE) - 1  ## shifted a bit from average; gets updated
  vxA <- var(quants, na.rm=TRUE)/(0.8*deltaB) ## shifted bit from covariance matrix; gets updated
  
  ## interpolate NA value using row means
  if(testna > 0){
    quants <- t(apply(quants, 1, function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); return(x)}))
    warning("NA values found in input; Imputation is performed.")
  }
  
  notdone <- TRUE
  count <- 0
  while( notdone )
  {
    ldenA <- dmvnorm(quants, mean=muA, sigma=vxA, log=TRUE)
    ldenB <- dmvnorm(quants, mean=muA+deltaA, sigma = deltaB*vxA, log=TRUE)
    
    pp <- pmax(ldenA,ldenB)
    tmp0 <-  pi0*exp(ldenA - pp ) 
    tmp1 <- (1-pi0)*exp(ldenB - pp) 
    pi0.update <- max( mean( tmp0/(tmp0+tmp1) ), 0.001 )  ## put in a floor
    
    count <- count + 1
    notdone <- count < cmax
    pi0 <- pi0.update
    ## update
    post1 <- tmp1/(tmp0+tmp1)
    tmp <- quants*(1-post1) + (quants-deltaA)*post1 /deltaB
    muA <- apply(tmp,2,sum)/(sum((1-post1)+post1/deltaB))
    
    tmp <-t((1- post1)*t((t(quants)-muA))) %*% t(t(quants)-muA) +
      t(post1/deltaB*t(t(quants)-muA-deltaA)) %*% t(t(quants)-muA-deltaA)
    vxA <- tmp/n
  }
  params <- list(pi0, muA, vxA, deltaA, deltaB)
  names( params ) <- c("pi0", "muA", "vxA", "deltaA", "deltaB")
  print(params)
  return(list(post1=post1, params=params, xb=xb) )
}
