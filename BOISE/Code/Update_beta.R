### Dirichlet process updating for one iteration (based on Neal's 2000 Paper, Algorithm 3)
## Input: prior vectors a, b, divergence alpha, data, old cl list of {K, N, C}
## Output a new cl list {K, N, C}

Update_beta <- function(cl, dat, a, b, alpha = 2){
  n = dim(dat)[1]
  m = dim(dat)[2]
  new_K = cl$K
  new_N = cl$N
  new_C = cl$C
  
  ## Impute missing values
  dat = unlist(sapply(1:m, function(j){
    imp = dat[,j]
    miss_label = which(is.na(imp))
    tmp= unlist(sapply(miss_label,function(i){
      target = which(new_C == new_C[i])
      at= sum(dat[target,j], na.rm = T) + a[j] 
      bt = length(target) - 1 - at + a[j] + b[j]
      p = rbeta(1,at,bt)
      return(rbinom(1,1,p))
    }))
    imp[miss_label] = tmp
    return(imp)
  }))
  
  for (i in 1:n) {
    if(new_N[new_C[i]] == 1){
      p = rep(0, new_K)
      for (k in 1:new_K) {
        if(new_N[k] == 0){
          next
        } else if(k == new_C[i]){
          p[k] = exp(log(alpha) + sum(dat[i,] * log(a / (a + b))) + 
                       sum((1 - dat[i,]) * log(b / (a + b))))
        } else{
          targets = which(new_C == k)
          tmp = a
          for (t in targets) {
            tmp = tmp + dat[t, ]
          }
          q = tmp / (a + b + new_N[k])
          p[k] = exp(log(new_N[k]) + sum(dat[i, ] * log(q)) + 
                       sum((1 - dat[i,]) * log(1 - q)))
        }
      }
      p = p / sum(p)
      class = which(rmultinom(1,1,p) == 1)
      if(class != new_C[i]){
        new_K = new_K - 1
        new_N[new_C[i]] = 0
        new_N[class] = new_N[class] + 1
        new_C[i] = class
      }
    } else{
      new_N[new_C[i]] = new_N[new_C[i]] - 1
      new_C[i] = 0
      p = rep(0, new_K + 1)
      p[new_K + 1] = exp(log(alpha) + sum(dat[i,] * log(a / (a + b))) + 
                           sum((1 - dat[i,]) * log(b / (a + b))))
      for (k in 1:new_K) {
        if(new_N[k] == 0){
          next
        } else{
          targets = which(new_C == k)
          tmp = a
          for (t in targets) {
            tmp = tmp + dat[t, ]
          }
          q = tmp / (a + b + new_N[k])
          p[k] = exp(log(new_N[k]) + sum(dat[i, ] * log(q)) + 
                       sum((1 - dat[i,]) * log(1 - q)))
        }
      }
      p = p / sum(p)
      class = which(rmultinom(1,1,p) == 1)
      new_C[i] = class
      new_N[class] = new_N[class] + 1
      if(class == new_K + 1){new_K = new_K + 1}
    }
  }
  
  new_cl = list("K" = 0, "N" = rep(0, 2 * n), "C" = rep(0,n))
  new_cl$K = length(which(new_N > 0))
  rank = order(new_N, decreasing = T)
  for (i in 1:new_cl$K) {
    r = rank[i]
    new_cl$N[i] = new_N[r]
    new_cl$C[which(new_C == r)] = i
  }
  
  return(new_cl)
}


# #Test
# load("clustering.RData")
# load("pkis1.rda")
# u <- foo$scaled.x
# dat <- 1*(u > .5 )
# rm(u)
# rm(foo)
# dat <- t(apply(pkis1, 1, function(x){
#   thres = mean(x) + 2 * sd(x)
#   return(as.numeric(x>thres))
# }))
# a = rep(3*mean(dat,na.rm = T), dim(dat)[2])
# b = 3-a
# cl = Initial_beta(dat,a,b,alpha = 30)
# cl = Update_beta(cl, dat,a,b, alpha = 30)
# for (i in 1:500) {
#   cl = Update_beta(cl, dat,a,b, alpha = 30)
# }
# tmp = dat[which(cl$C==1),]
# apply(tmp,2,function(x){return(mean(x,na.rm = T))})
# tmp = dat[which(cl$C==2),]
# apply(tmp,2,function(x){return(mean(x,na.rm = T))})
# tmp = dat[which(cl$C==3),]
# apply(tmp,2,function(x){return(mean(x,na.rm = T))})
# sum(cl$xi)
# which(cl$xi == 1)
