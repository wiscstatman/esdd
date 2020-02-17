targetcluster <- function(u, maxclust=41, epsilon=0.04, seed=123)
{
  set.seed(seed)
  
  # define vector for storing clustering results 
  u.ss <- numeric(maxclust-1)
  names(u.ss) <- c(2:maxclust)


  n <- ncol(u) # of compounds
  m <- nrow(u) # of targets
  
  ## define matrix to store clustring labels
  u.label <- matrix(NA, nrow=m, ncol=maxclust-1)
  
  for( ii in 2:maxclust ) 
  {
    u.kclu <- kmeans(u, centers = ii, iter.max = 500, nstart=20)
    
    ## assign within cluster distance
    u.ss[ii-1] <- u.kclu$tot.withinss
    ##u.ss[ii-1] <- u.kclu$tot.withinss/ii
    
    ## assign clustering results/labels
    u.label[,ii-1] <- u.kclu$cluster
  }
  
  
  ## initialize number of clusters
  nk<- 2
  
  ## choose number of clusters
  while(abs(1-u.ss[nk]/u.ss[nk-1]) >  epsilon){
    nk <- nk+1
  }
  
  return(list(label=u.label[,(nk-1)], withinss=u.ss[nk-1], K=nk))
}
