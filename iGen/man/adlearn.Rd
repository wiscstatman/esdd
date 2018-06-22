\name{adlearn}
\alias{adlearn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to implement adaptive selection
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
adlearn(u, label, nset = 16)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{training matrix of inhibition rates, rows are targets, columns are compounds
%%     ~~Describe \code{u} here~~
}
  \item{label}{clustering labels of targets, used in selection an initial set of inform compounds
%%     ~~Describe \code{label} here~~
}
  \item{nset}{size of set of inform compounds
%%     ~~Describe \code{nset} here~~
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
index of inform compounds
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
label <- sample(c(1:8), nrow(pkis1), replace=TRUE)

inform.index <- adlearn(u=pkis1, label=label, nset=16)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
