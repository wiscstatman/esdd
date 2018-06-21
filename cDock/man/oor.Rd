\name{oor}
\alias{oor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
One-Out Regression model result
%%  ~~function to do ... ~~
}
\description{
One-Out Regression model result
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
oor(scorelist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scorelist}{list of data.frames or matrices
%%     ~~Describe \code{scorelist} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
list of oor results
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Ericksen, S. S., Wu, H., Zhang, H., Michael, L. A., Newton, M. A., Hoffmann, F. M., & Wildman, S. A. (2017). Machine learning consensus scoring improves performance across targets in structure-based virtual screening. Journal of chemical information and modeling, 57(7), 1579-1590.

%% ~put references to the literature/web site here ~
}
\author{
Huikun Zhang
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
function (scorelist) 
{
    source("oor_prob.R")
    tmp.result <- oor_prob(scorelist)$tmp.result
    tmp.score <- sapply(tmp.result, function(x) {
        apply(x, 1, mean, na.rm = TRUE)
        return(x)
    })
    return(list(oor.cs = tmp.score))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
