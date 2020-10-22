\name{Initial_beta}
\alias{Initial_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  intermediate function to initial Chinese Restaurant Process for conjugate prior DPMM
}
\description{
Input: number of targets to be clustered n, divergency parameter alpha
Output: A list of {K, N, C,xi}. K = number of clusters, N = #of targets in each cluster, C = cluster assignment (including new target), x_i* = sampled outcome on new target.
}
\usage{
Initial_beta(dat, a, b, alpha = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
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
function (dat, a, b, alpha = 2) 
{
    n = dim(dat)[1]
    m = dim(dat)[2]
    N = rep(0, 2 * n)
    C = rep(0, n)
    K = 1
    N[1] = 1
    C[1] = 1
    for (i in 2:n) {
        p = rep(0, K + 1)
        p[K + 1] = alpha/(i - 1 + alpha)
        p[1:K] = sapply(1:K, function(x) {
            return(N[x]/(i - 1 + alpha))
        })
        class = which(rmultinom(1, 1, p) == 1)
        C[i] = class
        N[class] = N[class] + 1
        if (class == K + 1) {
            K = K + 1
        }
    }
    cl = list(K = K, N = N, C = C)
    return(cl)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
