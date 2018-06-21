\name{em_mmm}
\alias{em_mmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EM algorithm for fitting Multivariate Mixture Model
%%  ~~function to do ... ~~
}
\description{
EM algorithm for fitting free parameters in  Multivariate Mixture Model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
em_mmm(quants, pi0 = 0.98, deltaA = 1.76, deltaB = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{quants}{
  data.frame or matrix of docking score vectors
%%     ~~Describe \code{quants} here~~
}
  \item{pi0}{
  initial proportion of in-active component 
%%     ~~Describe \code{pi0} here~~
}
  \item{deltaA}{
  fixed distance of mean vectors between active and in-active components
%%     ~~Describe \code{deltaA} here~~
}
  \item{deltaB}{
  fixed ratio of variance between active and in-active components
%%     ~~Describe \code{deltaB} here~~
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
data(ace)
scores <- ace[,-1]
em_mmm(quants = scores, pi0 = 0.98, deltaA = 1.76, deltaB = 1 )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
