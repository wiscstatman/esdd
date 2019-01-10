\name{codingMC}
\alias{codingMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to deploy the coding selection method for informer based ranking 
}
\description{
The method takes a targets by compounds chemogenomic matrix, and returns a 
vector of scores (over compounds) for how useful the compounds are as members
of an informer set. (Lower scores are better.)
Such compounds will be tested experimentally against new
compounds to enable predictions on non-informer compound bio-activity.
}
\usage{
codingMC(data, nbig = 1000, clusternumbers = c(5, 10, 15, 20), lambda = 5, isize = 16)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
 a targets by compounds matrix of bioactivity levels (e.g. percent inhibition)
}
  \item{nbig}{
 size of the Monte Carlo search
}
  \item{clusternumbers}{
 entertained numbers of clusters in CS method
}
  \item{lambda}{
penalty parameter
}
  \item{isize}{
size of the informer set
}
}
\details{
The code converts input data to binary records of bioactivity of each compound
against each target.   It deploys a Monte Carlo search strategy, randomly 
sampling informer sets and induced code-word partitions of the target space,
and scores each with a penalized score indicating purity within clusters across
non-informer compounds.
}
\value{
  \item{compoundscore}{A vector over compounds holding average marginal
score.   Lowest scores constitute the methods informer set}}
  \item{best}{The single best identified set, by joint analysis}
}
\references{
Zhang et al, 2019
}
\author{
Michael A. Newton
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
adlearn for informer set via the adaptive selection method
}
\examples{


data(pkis1)
target <- 7  ## drop one out for a cross validation example

f <- codingMC(data=pkis1[-target,] , nbig=1000 )

      ## In paper, we use condor and do `nbig=10^7`

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
