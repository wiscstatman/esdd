codingMC <- function( data, nbig=1000, clusternumbers=c( 5, 10, 15, 20), lambda=5, isize=16 )
 {
 # Do a Monte Carlo search to score good informer compounds [coding selection]

 n <- ncol(data) # number of compounds
 m <- nrow(data) # number of targets

 # threshold, by row 2 SD rule

 ss <- apply( data, 1, sd )
 xb <- apply( data, 1, mean )

 y <- 1*(data > xb+2*ss )
 ###

 ## vectors to hold compound information 
 Asum <- rep(0,n)
 Anum <- rep(0,n)

 # initiate MC search
 a <- sample(1:n, size=isize ) # initial informer set
 A <- rep(FALSE, n)
 A[a] <- TRUE
 im <- imeasure(A, y, nclust=clusternumbers[1] ) # score set

 Anum[a] <- 1

 best <- vector( mode="list" )
 best$A <- A
 best$score <- im$sumd - lambda*length(im$mapping) 
 best$clust <- im$cluster
 best$map <- im$mapping

 Asum[a] <- best$score

 ## Now do further MC searching...
 for( nc in clusternumbers )
  { 
  print( paste("number of clusters", nc ) ) 
  for( ii in 2:nbig )
   {
    a <- sample(1:n, size=isize )
    A <- rep(FALSE, n)
    A[a] <- TRUE
  
    im <- imeasure(A, y, nclust=nc )
    score <- im$sumd - lambda*length( im$mapping )  ## longer mapping's better
    Asum[a] <- Asum[a] + score 
    Anum[a] <- Anum[a] +1
  
    if( score < best$score )
    {
      best$score <- score
      best$A <- A
      best$clust <- im$cluster
      best$map <- im$mapping
    }
   }
 }

 compoundscore <- Asum/Anum
 names(compoundscore) <- dimnames(data)[[2]]
 return( list(compoundscore=compoundscore, best=best)  )
}
