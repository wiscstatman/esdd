adpstep <-
function(u, cinformer) {
  
  ## split informer and non-informer part
  informer <- u[,c(cinformer)]
  noninformer <- u[,-c(cinformer)]  ## inhibition

  centroidsum.informer <- apply(informer, 1, sum)
  newcentroid <- (noninformer + centroidsum.informer)/(length(cinformer)+1)
  
  ## didnt remove the compound being evluated but easy to code
  distance <- (noninformer - newcentroid)^2
  distance <- colSums(distance)
  add.point <- names(distance)[distance==max(distance)]
  add.id <- which(colnames(u)==add.point)
  out.informer <- c(cinformer, add.id)
  return(out.informer)
}
