\name{Boise}
\alias{Boise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  function to identify a Bayes optimal informer set (corresponding to Algorithm 2 in Yu et al.)
}
\description{
Input includes a binary matrix of bioactivity data, Dirichlet-process-mixture-model samples  of
clustering over targes (with dpmm_beta.R to get cl_sample object). This function 
computes and minimizes posterior expected loss and returns an informer set A of size nA compounds.
}
\usage{
Boise(cl_sample, iter, size, nA, nT, a, b, x0, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl_sample}{
Posterior sampling of clusterings: A list of clustering assignments. cl_sample = {KK, NN, CC}.
}
  \item{iter}{
Sample size of DPMM clustering
}
  \item{size}{
Sample size of possible intermediate outcomes xA
}
  \item{nA}{
Length of informer set
}
  \item{nT}{
Length of top set, used to define loss function in BOISE.
}
  \item{a}{
Prior hyperparameters alpha0
}
  \item{b}{
Prior hyperparameters beta0
}
  \item{x0}{
Original interaction data.
}
  \item{alpha}{
Chinese restaurant prior mass m0
}
}
\details{
x0 is the original data used throughout the process. a, b, alpha and iter are the same parameters for \code{dpmm_beta}. size is a new parameter for sample size of possible outcomes on the new target.  
}
\value{
Bayes optimal informer set A with size nA.
}
\references{
Yu, P., Ericksen, S.S., Gitter, A. and Newton, M.A. (2020) \emph{Bayes Optimal Informer Sets for Drug Discovery}. University of Wisconsin, Technical Report, 2020.
}
\author{
Peng Yu, Michael A. Newton
}
\note{
  Please test the function for small nA. When nA >= 8 we recommend using Boise_Aug to augment informer set step by step, otherwise it may take days to complete the selection. 
}


\seealso{
\url{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006813}
}
\examples{
#suppose you have the dat and cl_sample as in dpmm_beta
data(pkis1)
dat <- t(apply(pkis1, 1, function(x){
  thres = mean(x) + 2 * sd(x)
  return(as.numeric(x>thres))
}))
a = rep(mean(dat,na.rm = T),ncol(dat))
b = 1-a
cl_sample = dpmm_beta(a,b,x0 = dat, warm=500, iter=10, step=10, alpha=15)
size = 1000
A = Boise(cl_sample, iter = 10, size = 1000, nA = 2, nT = 10, a, b, x0 = dat, alpha = 15)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
