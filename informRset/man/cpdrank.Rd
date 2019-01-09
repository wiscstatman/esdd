\name{cpdrank}
\alias{cpdrank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
given the assay results on some new target, function to rank non-inform compounds against this target 
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cpdrank(newinh, A, u)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{newinh}{vector of inhibition rate of inform compounds (assay result on new target)
%%     ~~Describe \code{newinh} here~~
}
  \item{A}{ index of inform compounds
%%     ~~Describe \code{A} here~~
}
  \item{u}{ training matrix of inhibition rates
%%     ~~Describe \code{u} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{expected}{expected number of active compounds on this new target}
\item{fdr}{ranked false discovery rate of all non-inform compounds against this new target}
\item{sum}{statistics of all non-inform compounds against this new target for ranking}
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
A <- c(1:16)
newinh <- seq(0.1,0.5, length=16)

cpd.rkresult <- cpdrank(newinh, A, pkis1)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
