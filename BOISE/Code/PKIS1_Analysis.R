### An example of BOISE computation pipeline using PKIS1 data.

source("Initial_beta.R")
source("Update_beta.R")
source("dpmm_beta.R")
source("clust_sum.R")
source("npel1_beta.R")
source("npel2.R")
source("Boise.R")
source("Boise_Aug.R")
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
iter = 10
step = 10
#find an empirical best prior mass alpha
alpha = 15

# One experiment in leave-one-out cross validation for BOISE framework
# First create cl_sample objects for ith training set
i = 1
test = dat[i, ]
train = dat[-i, ]
a = rep(mean(train),ncol(train))
b = 1 - a
cl_sample = dpmm_beta(a,b,x0 = train, warm, iter, step, alpha)

nA = 2
nT = 36
size = 1000
inform = Boise(cl_sample, iter,size, nA = nA, nT = nT,a,b, x0 = train,alpha)
## Add 1 more informer
inform = Boise_Aug(cl_sample,iter,size,nT,a,b,x0 = train,alpha,
                    inform, nAdd = 1)
## To maintain consistency, we would keep the initial clustering samples 
## and use the same assignments in all the computations.
nef.result = Evaluate(cl_sample,inform, measure = "nef",test,train,
                      nT,iter,a,b,alpha)
auc.result = Evaluate(cl_sample, inform, measure = "rocauc",test,train,
                      nT,iter,a,b,alpha)
mcc.result = Evaluate(cl_sample, inform, measure = "mat",test,train,
                      nT,iter,a,b,alpha)
f1.result = Evaluate(cl_sample, inform, measure = "f",test,train,
                     nT,iter,a,b,alpha)
