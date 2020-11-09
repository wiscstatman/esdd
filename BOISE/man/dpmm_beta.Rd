\name{dpmm_beta}
\alias{dpmm_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  function to do MCMC sampling for conjugate prior (Beta-binomial) clustering assignments
}
\description{
Output: A list cl_sample, including clustering assignments of size "iter"
}
\usage{
dpmm_beta(a, b, x0, warm = 200, iter = 10, step = 5, alpha = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
  Prior hyperparameters alpha0
}
  \item{b}{
  Prior hypereparameters beta0
}
  \item{x0}{
  Original interaction data x0
}
  \item{warm}{
  Burn-in parameter. How many samples are dropped in the beginning.
}
  \item{iter}{
  Sample size of DPMM
}
  \item{step}{
  Thinning parameters. Gap between consecutive samples.
}
  \item{alpha}{
  Chinese restaurant prior mass m0.
}
}
\details{
Typically, we choose hyperparameters a and b with the mean of actives/inactives in the data set. We choose alpha with several tests that aim to minimize the gap of cluster numbers between Chinese restaurant process prior and posterior. 
}
\value{
cl_sample, is a list of clustering assignment samples. cl_sample = {KK, NN, CC}, where KK is an array  of size "iter" with each element as the number of clusters; NN is a "iter" * "row(x0)" matrix with each row as the cluster assignment of the targets; CC is summary of every cluster in each assignment. 
}
\references{
Yu, P., Ericksen, S.S., Gitter, A. and Newton, M.A. (2020) \emph{Bayes Optimal Informer Sets for Drug Discovery}. University of Wisconsin, Technical Report, 2020.
}
\author{
Peng Yu, Michael A. Newton
}
\note{
It may take hours for sampling, please try small sample size and burn-in/thinning step first. Only works for binary data.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\url{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006813}
}
\examples{
data(pkis1)
## Use 2sd criteria to create binary matrix
dat <- t(apply(pkis1, 1, function(x){
  thres = mean(x) + 2 * sd(x)
  return(as.numeric(x>thres))
}))
a = rep(mean(dat,na.rm = T),ncol(dat))
b = 1-a
cl_sample = dpmm_beta(a,b,x0 = dat, warm=500, iter=10, step=10, alpha=15)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
