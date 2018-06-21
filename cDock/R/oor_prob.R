oor_prob <-
function(scorelist){
  
  n <- length(scorelist) ## number of score sets
  
  tmp.result <- vector("list", length=n) 
  names(tmp.result) <- names(scorelist)
  modellist <- vector("list", length=n)
  names(modellist) <- names(scorelist)
  
  for(ii in 1:n)
  {
    train.tmp <- scorelist[ii]
    modellist[[ii]] <- glm(label~., data=train.tmp, family="binamial")
  }
  
  for(ii in 1:n)
  {
    test.data <- scorelist[[ii]]
    result.test <- matrix(NA, nrow=nrow(test.data), ncol=n)
    
    for(jj in 1:n)
    {
      model.tmp <- modellist[[jj]]
      test.pred <- predict(model.tmp, newdata=test.data, type="response")
      result.test[,jj] <- test.pred
    }
    result.test <- result.test[,-ii]
    tmp.result[[ii]] <- result.test
    
    
  }
  
  return(list(models=modellist, probs=tmp.result))
}
