em_mvm <-
function( quants,  deltaA=1.6, deltaB=0.4 , pi0=.98 )
{
  cmax <- 13 ## EM max iteration
  
  ## both fixed
  ## deltaA = average estimated difference between active and decoy xbar
  ## deltaB = average estimated difference between active and decoy sd
  
  xb <- apply( quants, 1, mean, na.rm=TRUE )
  ss <- apply( quants, 1, sd, na.rm=TRUE ) 
  ss[is.na(ss)] <- min(min(ss, na.rm=TRUE)/2, 0.05) ## modifying ss for compounds with 1 observatins
  # initialize
  muA <- mean(xb) - 1  ## shifted a bit from average; gets updated
  sdA <- median(ss)  ## gets updated
  shapeB <- ( mean(ss) )^2/var(ss)  ## assume fixed in both components
  rateB <- mean(ss)/var(ss)  ##  to be adjusted
  
  notdone <- TRUE
  count <- 0
  while( notdone )
  {
    # on the means
    lden0A <- dnorm( xb, mean=muA, sd=sdA , log=TRUE )
    lden1A <- dnorm( xb, mean=muA+deltaA, sd=sdA , log=TRUE )
    
    # on the sd's
    lden0B <- dgamma( ss, shape=shapeB, rate=rateB  , log=TRUE )
    rr <- 1/( 1/rateB + deltaB/shapeB )  ## keeps constant shape and fixed diff on sd
    lden1B <- dgamma( ss, shape=shapeB, rate=rr, log=TRUE )
    
    lden0 <- lden0A+lden0B
    lden1 <- lden1A+lden1B
    
    pp <- pmax(lden0,lden1)
    
    tmp0 <-  pi0*exp(lden0 -pp ) 
    tmp1 <- (1-pi0)*exp(lden1-pp) 
    pi0.update <- max( mean( tmp0/(tmp0+tmp1) ), 0.001 )  ## put in a floor
    count <- count + 1
    notdone <- count < cmax
    pi0 <- pi0.update
    ## update 
    post1 <- tmp1/(tmp0+tmp1)
    
    tmp <- xb*(1-post1) + (xb-deltaA)*post1 
    muA <- mean(tmp)
    
    # update sdA
    tmp <- (1-post1)*( xb - muA )^2 + post1 * ( xb -muA-deltaA)^2
    sdA <- sqrt( mean(tmp) )
    
    # update rateB (just use null cases for simplicity)
    tmp <- (ss)*(1-post1) 
    mss0  <-  sum( tmp )/sum( 1-post1 ) # mean ss on null
    
    tmp <- (ss)*(1-post1)  + (ss-deltaB)*post1
    rateB <- shapeB/mss0
    
  }
  params <- c(pi0, muA, sdA, rateB, deltaA, deltaB,  shapeB )
  names( params ) <- c("pi0", "muA", "sdA", "rateB", "deltaA", "deltaB", "shapeB" )
  return(list(post1=post1, params=params, xb=xb, ss=ss) )
}
