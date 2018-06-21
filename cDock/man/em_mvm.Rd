\name{em_mvm}
\alias{em_mvm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EM algorithm for fitting Mean-Variance Mixture Model
%%  ~~function to do ... ~~
}
\description{
EM algorithm for fitting free parameters in  Mean-Variance Mixture Model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
em_mvm(quants, deltaA = 1.6, deltaB = 0.4, pi0 = 0.98)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{quants}{
    data.frame or matrix of docking score vectors
%%     ~~Describe \code{quants} here~~
}
  \item{deltaA}{
fixed distance of score average between active and in-active components
%%     ~~Describe \code{deltaA} here~~
}
  \item{deltaB}{
  fixed distance of standar deviation between active and in-active components
%%     ~~Describe \code{deltaB} here~~
}
  \item{pi0}{
  initial proportion of in-active component 
%%     ~~Describe \code{pi0} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{post1}{posterior probability to be in the active component}
 \item{params}{free parameters that are fitted via EM algorithm}
 \item{xb}{average docking score of each compound for reference purpose}
 \item{ss}{score standard deviation of each compound for reference purpose  }
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Huikun Zhang, Michael A. Newton
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
scores <- ace[,-1]
em_mvm(quants = scores, deltaA = 1.6, deltaB = 0.6, pi0 = 0.98 )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
