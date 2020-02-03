\name{targetcluster}
\alias{targetcluster}
\title{
function to do an initial clustering of targets [for Adaptive Selection]
%%  ~~function to do ... ~~
}
\description{
applies k-means on the rows of the initial data matrix
}
\usage{
targetcluster(u, maxclust=41, epsilon=0.04, seed=123 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{training matrix of inhibition rates, rows are targets, columns are compounds}
  \item{maxclust}{a maximal number of target clusters}
  \item{epsilon}{a threshold to do with selecting number of clusters}
  \item{seed}{input for set.seed}
}
\details{
applies k-means on the protein targets for various numbers of clusters
}
\value{
A list holding cluster labels and statistics for best number of clusters
}
\references{
index of inform compounds
%% ~put references to the literature/web site here ~
}
\author{
Huikun Zhang
}
\note{
this function is used to initiate adaptive selection
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(pkis1)

tc <- targetcluster(u=pkis1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
