cvadlearn <-
function(u, label, nset=16){
  source("adlearn.R")
  source("rocauc.R")
  source("inbin.R")
  library(pROC)
  
  m <- nrow(u) ## number of targets
  n <- ncol(u) ## number of compounds
  
  ## variables to store statistics
  auc.result <- numeric(m)
  names(auc.result) <- rownames(u)
  
  nef.result <- auc.result
  
  y <- inbin(u)$y
  thresh <- inbin(u)$thresh
  
  neval <- round(0.1*n)
  
  
  for(ii in 1:m)
  {
    print(ii)
    u.out <- u[ii,]
    y.out <- y[ii,]
    
    cur.df <- u[-ii,]
    cur.y.df <- y[-ii,]
    
    cur.inform <- adlearn(cur.df, label[-ii], nset)
    cur.sum <- cpdrank(u.out[cur.inform], cur.inform, cur.df)$sum
    cur.fdr <- cpdrank(u.out[cur.inform], cur.inform, cur.df)$fdr

    
    if(sum(y.out)==0){
      auc.result[ii] <- NA
      nef.result[ii] <- NA
    }else{
      if(sum(y.out[cur.inform])>0){
        tmp.response <- c(rep(1,sum(y.out[cur.inform])), y.out[-cur.inform])
        tmp.predict <- c(rep(max(cur.sum),sum(y.out[cur.inform])), cur.sum )
        auc.result[ii] <- roc(tmp.response, tmp.predict)$auc
        
        ef.max <- (max(sum(y.out),neval+sum(y.out[cur.inform]))/sum(y.out))*(n/(neval+sum(y.out[cur.inform])))
        ef.result <- ((sum(y.out[cur.inform])+ sum(head(y.out[names(cur.fdr)], neval)))/sum(y.out))*(n/(neval+sum(y.out[cur.inform])))
        nef.result[ii] <- (1+(ef.result-1)/(ef.max-1))/2
      }else{
        tmp.response <- c( y.out[-cur.inform])
        tmp.predict <- c( cur.sum )
        auc.result[ii] <- roc(tmp.response, tmp.predict)$auc
        
        ef.max <- (max(sum(y.out),neval)/sum(y.out))*(n/neval)
        ef.result <- ((sum(head(y.out[names(cur.fdr)], neval)))/sum(y.out))*(n/(neval+sum(y.out[cur.inform])))
        nef.result[ii] <- (1+(ef.result-1)/(ef.max-1))/2
      }

    }
    
    
  }
  return(list(rocauc=auc.result, nef=nef.result))
  
}
