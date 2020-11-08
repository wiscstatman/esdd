\name{BOISE-package}
\alias{BOISE-package}
\alias{BOISE}
\docType{package}
\title{
\packageTitle{BOISE}
}
\description{
\packageDescription{BOISE}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{BOISE}
\packageIndices{BOISE}

}
\author{
\packageAuthor{BOISE}

Maintainer: \packageMaintainer{BOISE}
}
\references{
Yu, P., Ericksen, S.S., Gitter, A. and Newton, M.A. (2020) \emph{Bayes Optimal Informer Sets for Drug Discovery}. University of Wisconsin, Technical Report, 2020.
}

\keyword{ package }
\seealso{
\url{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006813}
}
\examples{
library(BOISE)
#load data
data(pkis1)
## Use 2sd criteria to create binary matrix
dat <- t(apply(pkis1, 1, function(x){
  thres = mean(x) + 2 * sd(x)
  return(as.numeric(x>thres))
}))
#choose iterations, warm up step, step length, prior mass
warm = 500
iter = 10
step = 10
alpha = 15
# One experiment in leave-one-out cross validation for BOISE framework
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
}
