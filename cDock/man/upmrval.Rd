\name{upmrval}
\alias{upmrval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rvalue calculation for Uniform-Pareto Mixture model
%%  ~~function to do ... ~~
}
\description{
Rvalue calculation for Uniform-Pareto Mixture model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
upmrval(scores, alpha , n0 = 8, x0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scores}{ data.frame or matrix of docking scores
%%     ~~Describe \code{scores} here~~
}
  \item{alpha}{ vector or scalar of percentile values for R-value calculation
%%     ~~Describe \code{alpha} here~~
}
  \item{n0}{ hyperparameter in Pareto prior
%%     ~~Describe \code{n0} here~~
}
  \item{x0}{hyperparameter in Pareto prior
%%     ~~Describe \code{x0} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
R-values
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Henderson, N. C., & Newton, M. A. (2016). Making the cut: improved ranking and selection for large‐scale inference. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 78(4), 781-804.
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
data(ace)
score <- ace[,-1] ## first score vector
alpha.grid <-  seq(0.01, 0.1, 0.01)
upmval(score, alpha = alpha.grid, n0=8, x0= max(score, na.rm=TRUE))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
