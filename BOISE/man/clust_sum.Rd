\name{clust_sum}
\alias{clust_sum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  An intermediate function that summarizes the cl_sample to P.
}
\description{
  P is a list of matrices, with P[[i]] a matrix contain all information of ith cluster assignment. P[[i]] is K by (n+1) matrix. Each row for one cluster, last column as number of targets in the cluster.
}
\usage{
clust_sum(cl_sample, dat, iter, a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl_sample}{
%%     ~~Describe \code{cl_sample} here~~
}
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{iter}{
%%     ~~Describe \code{iter} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
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
function (cl_sample, dat, iter, a, b) 
{
    n = nrow(dat)
    m = ncol(dat)
    P = list(rep(0, iter))
    for (j in 1:iter) {
        K = cl_sample$KK[j]
        N = cl_sample$NN[j, 1:K]
        tmp_cl = matrix(0, K, (m + 1))
        tmp_cl[, m + 1] = N
        tmp_cl[, 1:m] = t(sapply(1:K, function(k) {
            target = which(cl_sample$CC[j, ] == k)
            if (length(target) == 1) {
                missing = which(is.na(dat[target, ]))
                ak = a + dat[target, ]
                bk = b + (1 - dat[target, ])
                ak[missing] = a[missing]
                bk[missing] = b[missing]
                return(ak/(ak + bk))
            }
            else {
                ak = a + colSums(dat[target, ], na.rm = T)
                bk = b + colSums(1 - dat[target, ], na.rm = T)
                return(ak/(ak + bk))
            }
        }))
        P[[j]] = tmp_cl
    }
    return(P)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
