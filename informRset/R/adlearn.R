adlearn <-
function(u, label, nset=16){
  require(glmnet)
  # find initial informers as nset/2 compounds predictive of input clustering
  fit.cur <- glmnet(as.matrix(u), as.factor(label), family="multinomial", alpha=1, dfmax=20, type.multinomial="grouped")
  
  lambd <- fit.cur$lambda
  jj <- length(lambd)
  lamd.cur <- lambd[jj]
  tmp_coeffs <- coef(fit.cur, s=lamd.cur) 
  infor.og <- tmp_coeffs[[1]]@i[-1]
  
  while(length(infor.og) > round(nset/2) ){
    jj <- jj-1
    lamd.cur <- lambd[jj]
    tmp_coeffs <- coef(fit.cur, s = lamd.cur  ) ## 86 104 125 137 253 313 349 357
    infor.og <- tmp_coeffs[[1]]@i[-1]
  }
  
  ## build remaining informer set to predict non-informers adaptively
  infor.tmp <- infor.og
  while(length(infor.tmp) < nset)
  {
    infor.tmp <- adpstep(u, infor.tmp)
  } 
  
  return(infor.tmp)
}
