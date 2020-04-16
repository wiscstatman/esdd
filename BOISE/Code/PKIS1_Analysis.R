source("dpmm_beta.R")
source("pel1_beta.R")
source("pel2.R")
source("inform_beta_v1.R")
source("inform_beta_v2.R")
source("Evaluate.R")
#load data
load("pkis1.rda")
## Use 2sd criteria to create binary matrix
dat <- t(apply(pkis1, 1, function(x){
  thres = mean(x) + 2 * sd(x)
  return(as.numeric(x>thres))
}))
rm(pkis1)

#choose iterations, warm up step, step length
warm = 500
iter = 100
step = 10
#find an empirical best prior mass alpha
alpha = 15

# One experiment in leave-one-out cross validation for BOISE framework
# First create cl_sample objects for ith training set
value <- commandArgs(trailingOnly=TRUE)
i <- as.numeric(value) + 1
test = dat[i, ]
train = dat[-i, ]
a = rep(mean(train),ncol(train))
b = 1 - a
cl_sample = dpmm_beta(a,b,x0 = train, warm, iter, step, alpha)

nA = 3
nT = 36

### Posterior sample of x_i*
size = 1000
n = nrow(train)
m = ncol(train)
cl_sample$XX = rep(0, iter*size*m)
dim(cl_sample$XX) = c(iter, size, m)

## Create a matrix summarize information of each clustering.
## K by (n+1) matrix. Each row for one cluster, last column as number of targets in cluster
P = list(rep(0,iter))
n = nrow(train)
m = ncol(train)
for (j in 1:iter) {
  K = cl_sample$KK[j]
  N = cl_sample$NN[j,1:K]
  tmp_cl = matrix(0,K,(m+1))
  tmp_cl[,m+1] = N
  tmp_cl[,1:m] = t(sapply(1:K,function(k){
    target = which(cl_sample$CC[j,] == k)
    if(length(target)==1){
      return((a + train[target,])/2)
    }else{
      return((a + apply(train[target,], 2, sum))/(N[k] + 1))
    }
  }))
  P[[j]] = tmp_cl
}
## Sample for x_i*
for (j in 1:iter) {
  K = cl_sample$KK[j]
  p = rep(0, K + 1)
  p[K + 1] = alpha / (n + alpha)
  p[1:K] = P[[j]][,m+1] / (n+alpha)
  cl_sample$XX[j, , ] = t(as.matrix(sapply(1:size, function(s){
    classi = which(rmultinom(1, 1, p) == 1)
    if(classi == K+1){
      post_theta = a / (a + b)
    } else{
      post_theta = P[[j]][classi,1:m]
    }
    new_xi = sapply(1:m, function(x){
      return(as.numeric(rbinom(1, 1, p = post_theta[x])))
    })
    return(new_xi)
  })))
}

inform = inform_beta1(cl_sample,iter,size, nA = nA, nT = nT,a,b, x0 = train,alpha)
#tmp = read.table("Test6_result.txt")
#pre_inform = as.numeric(unlist(strsplit(as.character(tmp$V2[i]), split = ' ')))
# inform = inform_beta2(cl_sample,iter,size,nT,a,b,x0 = train,alpha,
#                      inform = pre_inform, nAdd = 1)
## To maintain consistency, we would keep the initial clustering samples 
## and use the same assignments in all the computations.
nef.result = Evaluate(cl_sample, inform, measure = "nef",test,train,
                      nA,nT,iter,a,b,alpha)
auc.result = Evaluate(cl_sample, inform, measure = "rocauc",test,train,
                      nA,nT,iter,a,b,alpha)
mcc.result = Evaluate(cl_sample, inform, measure = "mat",test,train,
                      nA,nT,iter,a,b,alpha)
f1.result = Evaluate(cl_sample, inform, measure = "f",test,train,
                     nA,nT,iter,a,b,alpha)
