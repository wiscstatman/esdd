\name{adpstep}
\alias{adpstep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to implement adaptive step in adaptive selection
%%  ~~function to do ... ~~
}
\description{
function to implement adaptive step in adaptive selection
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
adpstep(u, cinformer)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{ matrix of inhibition rates, rows are targets, columns are compounds
%%     ~~Describe \code{u} here~~
}
  \item{cinformer}{ index of initial inform compounds 
%%     ~~Describe \code{cinformer} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
index set of all inform compounds
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
cinform <- c(1:8)
final.inform <- adpstep(pkis1, cinform)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
