\name{enrichf}
\alias{enrichf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
enrichf(x, y, quantile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{quantile}{
%%     ~~Describe \code{quantile} here~~
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
function (x, y, quantile) 
{
    x[is.na(x)] <- min(x, na.rm = TRUE)
    if (length(x) != length(y)) {
        stop(paste("The length of theta should be equal to the length of labels"))
    }
    N <- length(y)
    n <- sum(y)
    ord <- order(x, decreasing = TRUE, na.last = TRUE)
    yard <- cumsum(y[ord])
    pos <- round(N * quantile)
    if (n == 0) {
        stop(paste("number of active ligands should be positive"))
    }
    else {
        ef <- yard[pos]/(n * quantile)
    }
    return(ef)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
