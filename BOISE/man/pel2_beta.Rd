\name{pel2_beta}
\alias{pel2_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  function to calculate PEL-2 loss based on intermediate data xA
}
\description{
Input: Matrix P summarize the cluster assignment information (For computation)
        Informer set A, Data matrix x0, Intermediate data xA, 
        Informer size nA, Top size nT, 
        prior a,b, divergence alpha
Output: PEL2 value (given x0,xA,Cl) [posterior expected theta E(theta|C,x0,xA)]
}
\usage{
pel2_beta(P, x0, xA, nA, A, nT = 10, a, b, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{
%%     ~~Describe \code{P} here~~
}
  \item{x0}{
%%     ~~Describe \code{x0} here~~
}
  \item{xA}{
%%     ~~Describe \code{xA} here~~
}
  \item{nA}{
%%     ~~Describe \code{nA} here~~
}
  \item{A}{
%%     ~~Describe \code{A} here~~
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
function (P, x0, xA, nA, A, nT = 10, a, b, alpha) 
{
    K = nrow(P)
    m = ncol(x0)
    n = nrow(x0)
    p = rep(0, K + 1)
    theta = matrix(0, K + 1, m)
    theta[K + 1, ] = a/(a + b)
    theta[K + 1, A] = (a[A] + xA)/(a[A] + b[A] + 1)
    p[K + 1] = exp(sum(xA * log(a[A]/(a[A] + b[A]))) + sum((1 - 
        xA) * log(b[A]/(a[A] + b[A]))))
    p[1:K] = apply(P[, 1:m], 1, function(q) {
        return(exp(sum(xA * log(q[A])) + sum((1 - xA) * log(1 - 
            q[A]))))
    })
    theta[1:K, ] = P[, 1:m]
    if (length(A) == 1) {
        theta[1:K, A] = (P[, A] * (P[, m + 1] + rep(a[1] + b[1], 
            K)) + xA)/(P[, m + 1] + rep(a[1] + b[1], K) + 1)
    }
    else {
        theta[1:K, A] = t(as.matrix(apply(P[, c(A, m + 1)], 1, 
            function(q) {
                l = length(A)
                counts = q[1:l] * (q[l + 1] + a[1] + b[1])
                counts = counts + xA
                return(counts/(q[l + 1] + a[1] + b[1] + 1))
            })))
    }
    theta[1:K, A] = t(matrix(rep(xA, K), nrow = nA, ncol = K))
    weight = c(P[, m + 1], alpha)
    weight = weight/(sum(weight))
    post_prob = sum(weight * p)
    w = c(P[, m + 1], alpha)
    w = w * p
    w = as.matrix(w/sum(w))
    post_theta = t(w) \%*\% theta
    postls <- list(post_theta = post_theta, post_prob = post_prob)
    return(postls)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
