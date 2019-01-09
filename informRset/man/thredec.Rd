\name{thredec}
\alias{thredec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to decide binary activity threshold of inhibition on new target
%%  ~~function to do ... ~~
}
\description{
function to decide binary activity threshold of inhibition on new target
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
thredec(u, informer, newinh)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{ inhibition matrix 
%%     ~~Describe \code{u} here~~
}
  \item{informer}{ index of inform compounds
%%     ~~Describe \code{informer} here~~
}
  \item{newinh}{ inhibition data of inform compounds on some new target
%%     ~~Describe \code{newinh} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
threshold of inhibition rate on some new target for binary activity label
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
inform <- c(1:16)
newinh <- seq(0.1,0.5, length=16)

newthresh <- thredec(pkis1, inform, newinh)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
