imeasure <- function(A, y, nclust)
{
  # measure the informer-set-ness of a putative informer set A
  # small values may give good informer sets
  
  # find a code...a good one will not entail much extrapolation
  # for a new target [i.e., it will cover a lot of bases and put
  # observed activity profiles readily into clusters]
  
  # it will be a function from informer-set activities into
  # clusters [so based on codewords of the transformation]
  
  uu <- apply( y[,A], 1, paste, collapse="" )
  uuu <- unique(uu)  ## the distinct codewords; functions can't 
  ##distinguish individual codewords
  ncodewords <- length(uuu)
  
  # now make nclust by merging these...[get a random partition]
  
  if(ncodewords < nclust)
  {
    hh <- c(1:ncodewords)
  }else{
    
    ## random block codewords to form cluster
    bpoint <- sample(c(2:ncodewords), size=nclust-1, replace=FALSE)
    bpoint <- bpoint[order(bpoint)]
    hh <- rep(1, ncodewords)
    
    ## asign codewords to cluster
    for(ii in bpoint)
    {
      hh[-c(1:(ii-1))] <- hh[-c(1:(ii-1))] + 1
    }
    hh <- sample(hh)
    
  }
  
  names(hh) <- uuu
  gg <- as.numeric(hh[uu])  ##  cluster id's 
  
  z <- y[,!A]  ## non-informer
  # get with cluster sum of pairwise distances
  
  sumd <- 0
  for( b in 1:nclust )
  {
    zz <- z[gg==b,]
    dd <- dist( zz, method="binary" )
    if(sum(gg==b)==1){ sumd <- sumd + 1}
    if(sum(gg==b)>1)
    {
      sumd <- sumd + sum(dd)
    }
  }
  return( list(sumd=sumd,cluster=gg,mapping=hh) )
}
