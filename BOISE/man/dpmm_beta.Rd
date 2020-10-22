\name{dpmm_beta}
\alias{dpmm_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dpmm_beta(a, b, x0, warm = 200, iter = 10, step = 5, alpha = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{x0}{
%%     ~~Describe \code{x0} here~~
}
  \item{warm}{
%%     ~~Describe \code{warm} here~~
}
  \item{iter}{
%%     ~~Describe \code{iter} here~~
}
  \item{step}{
%%     ~~Describe \code{step} here~~
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
function (a, b, x0, warm = 200, iter = 10, step = 5, alpha = 2) 
{
    n = dim(x0)[1]
    m = dim(x0)[2]
    KK = rep(0, iter)
    NN = matrix(0, iter, 2 * n)
    CC = matrix(0, iter, n)
    cl = Initial_beta(x0, a, b, alpha)
    for (i in 1:warm) {
        cl = Update_beta(cl, dat = x0, a, b, alpha)
    }
    for (i in 1:iter) {
        tmp = 0
        while (tmp < step) {
            cl = Update_beta(cl, dat = x0, a, b, alpha)
            tmp = tmp + 1
        }
        KK[i] = cl$K
        NN[i, ] = cl$N
        CC[i, ] = cl$C
    }
    cl_sample = list(KK = KK, NN = NN, CC = CC)
    return(cl_sample)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
