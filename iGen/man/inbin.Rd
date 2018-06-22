\name{inbin}
\alias{inbin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to decide thresholds of binary activity label on inhibition rates for a library of targets
%%  ~~function to do ... ~~
}
\description{
function to decide thresholds of binary activity label on inhibition rates for a library of targets
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
inbin(u)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{matrix of inhibition rates, rows are targets; columns are compounds
%%     ~~Describe \code{u} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{y }{thresholded binary activity matrix}
\item{thresh }{vectors of threholds on inhibition rates for each target}
%% ...
}
\references{
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
data(pkis1)
y <- inbin(pkis1)$y
y_thresh <- inbin(pkis1)$thresh
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
