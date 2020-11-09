\name{Initial_beta}
\alias{Initial_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  intermediate function to initial Chinese Restaurant Process for conjugate prior DPMM
}
\description{
Output a list of number of clusters, number of targets in each cluster, and cluster assignment.
}
\usage{
Initial_beta(dat, a, b, alpha = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
  The original interaction data x0
}
  \item{a}{
  Beta prior hyperparameters alpha0 
}
  \item{b}{
  Beta prior hyperparameters beta0
}
  \item{alpha}{
  Chinese restaurant process prior mass m0.
}
}
\details{
a, b and alpha are hyperparameters determined by users. 
}
\value{
cl_sample, a list of {K, N, C}, where K = number of clusters, N = #of targets in each cluster, C = cluster assignment.
}
\references{
Yu, P., Ericksen, S.S., Gitter, A. and Newton, M.A. (2020) \emph{Bayes Optimal Informer Sets for Drug Discovery}. University of Wisconsin, Technical Report, 2020.
}
\author{
Peng Yu, Michael A. Newton
}
\note{
This is just an internal function for \code{dpmm_beta}
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
rm(pkis1)
a = rep(mean(dat),ncol(dat))
b = 1-a
cl = Initial_beta(dat,a,b,10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
