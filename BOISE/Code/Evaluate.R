### Evaluation of selected informer set given the DPMM sampling
## Input: DPMM cluster sampling "cl_sample", informer set "inform",
##        measurement criteria "measure", test set, train set, nA, nT, 
##        MCMC sample size iter, priors a and b, prior mass alpha
## Output: Evaluated value under given criteria

Evaluate <- function(cl_sample, inform, measure,
                     test, train, nA, nT, iter, a, b ,alpha){
  source("pel2.R")
  Score = rep(0, ncol(train))
  xA = test[inform]
  for(k in 1:iter){
    tmp_cl = list(K = cl_sample$KK[k], N = cl_sample$NN[k,], C = cl_sample$CC[k,])
    post_theta = pel2_beta(tmp_cl, x0 = train, xA, nA, A = inform, nT, a, b, alpha)
    Score = Score + post_theta
  }
  Score[inform[which(xA==1)]] = rep(max(Score) + 1, sum(xA))
  Score[inform[which(xA==0)]] = rep(min(Score) - 1, nA - sum(xA))
  test = as.vector(test)
  Score = as.vector(Score)
  if(measure == "nef"){
    top = order(Score,decreasing = T)[1:nT]
    pred_hit = sum(test[top])
    hit = sum(test)
    maxhit = min(hit,nT)
    result = ((pred_hit/nT - hit/ncol(train)) / (maxhit/nT - hit/ncol(train)) + 1)/2
  } else if(measure == "rocauc"){
    rocobj = pROC::roc(test,Score)
    result = rocobj$auc
  } else if(measure %in% c("mat", "f")){
    pred.obj = ROCR::prediction(Score, test)
    perform.obj = ROCR::performance(pred.obj, measure)
    result = max(unlist(perform.obj@y.values),na.rm = T)
  } else{
    print("Criteria is not supported.")
    result = 0
  }
  return(result)
}
