oor <-
function(scorelist){
  source("oor_prob.R")
  
  tmp.result <- oor_prob(scorelist)$probs
  tmp.score <- sapply(tmp.result, function(x){apply(x, 1, mean, na.rm=TRUE); return(x)})
  
  return(list(oor.cs=tmp.score))
  

}
