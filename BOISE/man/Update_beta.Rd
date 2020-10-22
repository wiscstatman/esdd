\name{Update_beta}
\alias{Update_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Dirichlet process updating for one iteration (based on Neal's 2000 Paper, Algorithm 3)
}
\description{
Input: prior vectors a, b, divergence alpha, data, old cl list of {K, N, C}
Output a new cl list {K, N, C}
}
\usage{
Update_beta(cl, dat, a, b, alpha = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl}{
%%     ~~Describe \code{cl} here~~
}
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
function (cl, dat, a, b, alpha = 2) 
{
    n = dim(dat)[1]
    m = dim(dat)[2]
    new_K = cl$K
    new_N = cl$N
    new_C = cl$C
    for (i in 1:n) {
        if (new_N[new_C[i]] == 1) {
            p = rep(0, new_K)
            for (k in 1:new_K) {
                if (new_N[k] == 0) {
                  next
                }
                else if (k == new_C[i]) {
                  p[k] = exp(log(alpha) + sum(dat[i, ] * log(a/(a + 
                    b)), na.rm = T) + sum((1 - dat[i, ]) * log(b/(a + 
                    b)), na.rm = T))
                }
                else {
                  subdat = dat[which(new_C == k), ]
                  if (is.vector(subdat)) {
                    missing = which(is.na(subdat))
                    ak = a + subdat
                    ak[missing] = a[missing]
                    bk = b + 1 - subdat
                    bk[missing] = b[missing]
                  }
                  else {
                    ak = a + colSums(subdat, na.rm = T)
                    bk = b + colSums(1 - subdat, na.rm = T)
                  }
                  q = ak/(ak + bk)
                  p[k] = exp(log(new_N[k]) + sum(dat[i, ] * log(q), 
                    na.rm = T) + sum((1 - dat[i, ]) * log(1 - 
                    q), na.rm = T))
                }
            }
            p = p/sum(p)
            class = which(rmultinom(1, 1, p) == 1)
            if (class != new_C[i]) {
                new_K = new_K - 1
                new_N[new_C[i]] = 0
                new_N[class] = new_N[class] + 1
                new_C[i] = class
            }
        }
        else {
            new_N[new_C[i]] = new_N[new_C[i]] - 1
            new_C[i] = 0
            p = rep(0, new_K + 1)
            p[new_K + 1] = exp(log(alpha) + sum(dat[i, ] * log(a/(a + 
                b)), na.rm = T) + sum((1 - dat[i, ]) * log(b/(a + 
                b)), na.rm = T))
            for (k in 1:new_K) {
                if (new_N[k] == 0) {
                  next
                }
                else {
                  subdat = dat[which(new_C == k), ]
                  if (is.vector(subdat)) {
                    missing = which(is.na(subdat))
                    ak = a + subdat
                    ak[missing] = a[missing]
                    bk = b + 1 - subdat
                    bk[missing] = b[missing]
                  }
                  else {
                    ak = a + colSums(subdat, na.rm = T)
                    bk = b + colSums(1 - subdat, na.rm = T)
                  }
                  q = ak/(ak + bk)
                  p[k] = exp(log(new_N[k]) + sum(dat[i, ] * log(q), 
                    na.rm = T) + sum((1 - dat[i, ]) * log(1 - 
                    q), na.rm = T))
                }
            }
            p = p/sum(p)
            class = which(rmultinom(1, 1, p) == 1)
            new_C[i] = class
            new_N[class] = new_N[class] + 1
            if (class == new_K + 1) {
                new_K = new_K + 1
            }
        }
    }
    new_cl = list(K = 0, N = rep(0, 2 * n), C = rep(0, n))
    new_cl$K = length(which(new_N > 0))
    rank = order(new_N, decreasing = T)
    for (i in 1:new_cl$K) {
        r = rank[i]
        new_cl$N[i] = new_N[r]
        new_cl$C[which(new_C == r)] = i
    }
    return(new_cl)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
