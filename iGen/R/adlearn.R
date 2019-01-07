adlearn <-
function(u, label, nset=16){
  ##library(glmnet)
  ##source("adpstep.R")
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
  
  ## adaptive step
  infor.tmp <- infor.og
  while(length(infor.tmp) < nset)
  {
    infor.tmp <- adpstep(u, infor.tmp)
  } 
  
  return(infor.tmp)
}
