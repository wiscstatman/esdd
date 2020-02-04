\name{imeasure}
\alias{imeasure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
a function to score a potential informer set of compounds
}
\description{
Given bioactivity data, the function computes the coding-selection score
of an input informer set of compounds.   It scores them in conjunction with
a random partition of the induced code words.
}
\usage{
imeasure(A, y, nclust)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
  a logical vector; TRUE where compound is in informer set
}
  \item{y}{
binary matrix, targets by compounds, holding bioactivity calls
}
  \item{nclust}{
number of clusters to be entertained on the target space
}
}
\details{
The code first randomly partitions the targets by blocking targets according
to code-words induced by the input informer set.   Then it scores the resulting
mapping for within-cluster purity on the non-informers
}
\value{
A list containing the coding-selection score and details of the target-space clustering/mapping
}
\references{
Zhang et al. 2019
}
\author{
Michael A. Newton
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
used as core function by `codingMC`
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
