\name{rocauc}
\alias{rocauc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
internal function to compute area under ROC curve 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rocauc(obs, pred)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{
%%     ~~Describe \code{obs} here~~
}
  \item{pred}{
%%     ~~Describe \code{pred} here~~
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
%% ~put references to the literature/web site here ~
}
\author{
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (obs, pred) 
{
    pred[is.na(pred)] <- min(pred, na.rm = TRUE)
    pred <- prediction(pred, obs)
    auc <- performance(pred, "auc")@y.values[[1]]
    return(auc)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
