\name{codingMC}
\alias{codingMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
codingMC(data, nbig = 1000, clusternumbers = c(5, 10, 15, 20), lambda = 5, isize = 16)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{nbig}{
%%     ~~Describe \code{nbig} here~~
}
  \item{clusternumbers}{
%%     ~~Describe \code{clusternumbers} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{isize}{
%%     ~~Describe \code{isize} here~~
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
function (data, nbig = 1000, clusternumbers = c(5, 10, 15, 20), 
    lambda = 5, isize = 16) 
{
    n <- ncol(data)
    m <- nrow(data)
    ss <- apply(data, 1, sd)
    xb <- apply(data, 1, mean)
    y <- 1 * (data > xb + 2 * ss)
    Asum <- rep(0, n)
    Anum <- rep(0, n)
    a <- sample(1:n, size = isize)
    A <- rep(FALSE, n)
    A[a] <- TRUE
    im <- imeasure(A, y, nclust = clusternumbers[1])
    Anum[a] <- 1
    best <- vector(mode = "list")
    best$A <- A
    best$score <- im$sumd - lambda * length(im$mapping)
    best$clust <- im$cluster
    best$map <- im$mapping
    Asum[a] <- best$score
    for (nc in clusternumbers) {
        print(paste("cluster size", nc))
        for (ii in 2:nbig) {
            a <- sample(1:n, size = isize)
            A <- rep(FALSE, n)
            A[a] <- TRUE
            im <- imeasure(A, y, nclust = nc)
            score <- im$sumd - lambda * length(im$mapping)
            Asum[a] <- Asum[a] + score
            Anum[a] <- Anum[a] + 1
            if (score < best$score) {
                best$score <- score
                best$A <- A
                best$clust <- im$cluster
                best$map <- im$mapping
            }
        }
    }
    compoundscore <- Asum/Anum
    return(list(compoundscore = compoundscore, best = best))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
