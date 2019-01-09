\name{imeasure}
\alias{imeasure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
imeasure(A, y, nclust)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{nclust}{
%%     ~~Describe \code{nclust} here~~
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
function (A, y, nclust) 
{
    uu <- apply(y[, A], 1, paste, collapse = "")
    uuu <- unique(uu)
    ncodewords <- length(uuu)
    if (ncodewords < nclust) {
        hh <- c(1:ncodewords)
    }
    else {
        bpoint <- sample(c(2:ncodewords), size = nclust - 1, 
            replace = FALSE)
        bpoint <- bpoint[order(bpoint)]
        hh <- rep(1, ncodewords)
        for (ii in bpoint) {
            hh[-c(1:(ii - 1))] <- hh[-c(1:(ii - 1))] + 1
        }
        hh <- sample(hh)
    }
    names(hh) <- uuu
    gg <- as.numeric(hh[uu])
    z <- y[, !A]
    sumd <- 0
    for (b in 1:nclust) {
        zz <- z[gg == b, ]
        dd <- dist(zz, method = "binary")
        if (sum(gg == b) == 1) {
            sumd <- sumd + 1
        }
        if (sum(gg == b) > 1) {
            sumd <- sumd + sum(dd)
        }
    }
    return(list(sumd = sumd, cluster = gg, mapping = hh))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
