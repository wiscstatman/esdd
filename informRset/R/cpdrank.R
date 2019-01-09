cpdrank <-
function(newinh, A, u)
{
  ##source("thredec.R") ## load threhhold deciding function
  ##source("inbin.R")

  y <- inbin(u)$y
  
  threhnew <- thredec(u, A, newinh)
  newlabel <- 1 * (newinh > threhnew)
  atomnew <- paste(newlabel, collapse = "")
  
  
  atomset <- apply(y[,A], 1, paste, collapse="")
  
  if(sum(atomset==atomnew) > 1 ){
    part <- y[atomset==atomnew,]
    part <- part[,-A]
    natom <- dim(part)[1]
    inhibitor.sum <- apply(part,2,sum)
    kinase.sum <- apply(part,1,sum)
  }else if(sum(atomset==atomnew) == 1){
    part <- y[atomset==atomnew,]
    part <- part[-A]
    natom <- 1
    inhibitor.sum <- part
    kinase.sum <- sum(part)
  }else{
    atomdif <- colSums(t(y[,A]-newlabel))
    
    part <- y[atomdif==min(atomdif),]
    part <- part[, -A]
    
    if(is.vector(part)==TRUE)
    {
      natom <- 1
      inhibitor.sum <- part
      kinase.sum <- sum(part)
    }else{
      natom <- dim(part)[1]
      inhibitor.sum <- apply(part,2,sum)
      kinase.sum <- apply(part,1,sum)
    }
  }
  
  sum.order <- inhibitor.sum[order(-inhibitor.sum)]
  inhibitor.cum <- cumsum(sum.order)
  fdr.inhibitor <- 1-inhibitor.cum/(natom*c(1:length(inhibitor.cum)))
  rank.newton <- names(fdr.inhibitor)
  
  kinase.tb <- table(kinase.sum)
  kinase.tb <- data.matrix(kinase.tb)
  
  expected.nactive <- round(sum(kinase.tb * as.numeric(rownames(kinase.tb)))/sum(kinase.tb))
  
  return(list(expected=expected.nactive, fdr=fdr.inhibitor, sum=inhibitor.sum))
  
}
