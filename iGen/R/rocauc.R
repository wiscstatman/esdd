rocauc <-
function(obs, pred) {
  
  ## assign min value to NA element
  pred[is.na(pred)] <- min(pred, na.rm=TRUE)
  # ok <- !is.na(pred)
  # pred <- pred[ok]
  # obs <- obs[ok]
  
  pred <- prediction(pred, obs)
  auc  <- performance(pred, "auc")@y.values[[1]]
  return(auc)
}
