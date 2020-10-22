\name{Boise_Aug}
\alias{Boise_Aug}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  function to augment an existing informer set, usually used in parallel computing.
}
\description{
Pre-step: Sample DPMM samples of cl_sample. 
          A previous selected smaller informer set inform.
Input: list cl_sample, corresponding parameters a,b,iter,size,alpha; nT; x0;
Pre-selected informer set inform, number of informers to add nAdd
Output: Informer set A with size n+nAdd (via advanced adaptive selection)
}
\usage{
Boise_Aug(cl_sample, iter, size, nT, a, b, x0, alpha, inform, nAdd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl_sample}{
A list of {K, N, C}. K = number of clusters, N = #of targets in each cluster
}
  \item{iter}{
%%     ~~Describe \code{iter} here~~
}
  \item{size}{
%%     ~~Describe \code{size} here~~
}
  \item{nT}{
%%     ~~Describe \code{nT} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{x0}{
%%     ~~Describe \code{x0} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{inform}{
%%     ~~Describe \code{inform} here~~
}
  \item{nAdd}{
%%     ~~Describe \code{nAdd} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cl_sample, iter, size, nT, a, b, x0, alpha, inform, 
    nAdd) 
{
    if (!require("parallel")) {
        stop("The package parallel was not installed")
    }
    P = clust_sum(cl_sample, x0, iter, a, b)
    n = nrow(x0)
    m = ncol(x0)
    cl = cl_sample
    cl$XX = rep(0, iter * size * m)
    dim(cl$XX) = c(iter, size, m)
    for (j in 1:iter) {
        K = cl$KK[j]
        p = rep(0, K + 1)
        p[K + 1] = alpha/(n + alpha)
        p[1:K] = P[[j]][, m + 1]/(n + alpha)
        cl$XX[j, , ] = t(as.matrix(sapply(1:size, function(s) {
            classi = which(rmultinom(1, 1, p) == 1)
            if (classi == K + 1) {
                post_theta = a/(a + b)
            }
            else {
                post_theta = P[[j]][classi, 1:m]
            }
            new_xi = sapply(1:m, function(x) {
                return(as.numeric(rbinom(1, 1, p = post_theta[x])))
            })
            return(new_xi)
        })))
    }
    step = 1
    n = length(inform)
    candidate = (1:dim(x0)[2])[-inform]
    while (step <= nAdd) {
        pel = rep(0, length(candidate))
        pel = unlist(mclapply(candidate, function(x) {
            return(pel1_beta(cl, P, iter, size, A = c(inform, 
                x), nA = n + step, nT, a, b, x0, alpha))
        }, mc.cores = detectCores()))
        tmp = candidate[order(pel)[1]]
        inform = c(inform, tmp)
        step = step + 1
        candidate = candidate[-which(candidate == tmp)]
    }
    return(inform)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
