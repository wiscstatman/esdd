\name{oor_prob}
\alias{oor_prob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predicting results for One-Out Regression
%%  ~~function to do ... ~~
}
\description{
Predicting results for One-Out Regression
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
oor_prob(scorelist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scorelist}{ list of data.frames or matrices
%%     ~~Describe \code{scorelist} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{models }{list of logistic regression models}
\item{probs }{list of predicted probability matrices using One-Out model}
%% ...
}
\references{
Ericksen, S. S., Wu, H., Zhang, H., Michael, L. A., Newton, M. A., Hoffmann, F. M., & Wildman, S. A. (2017). Machine learning consensus scoring improves performance across targets in structure-based virtual screening. Journal of chemical information and modeling, 57(7), 1579-1590.

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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
