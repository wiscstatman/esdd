\name{Vupmtail}
\alias{Vupmtail}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate tail probability vector on multiple alpha of Uniform-Pareto Model
%%  ~~function to do ... ~~
}
\description{
Generate tail probability vector on multiple alpha of Uniform-Pareto Model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Vupmtail(svec, alpha = 0.01, n0, x0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{svec}{ vector of observations
%%     ~~Describe \code{svec} here~~
}
  \item{alpha}{ the percentile of tail, P(theta>theta_alpha) = alpha 
%%     ~~Describe \code{alpha} here~~
}
  \item{n0}{ hyperparameter in Pareto prior
%%     ~~Describe \code{n0} here~~
}
  \item{x0}{ hyperparameter in Pareto prior
%%     ~~Describe \code{x0} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
tail probability P(theta > theta_alpha | data) for multiple alpha values
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
valpha <- seq(0.01, 0.1, 0.01)
score <- ace[1,-1] ## first score vector
Vupmtail(score, alpha = 0.01, n0=8, x0= max(score, na.rm=TRUE))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
