\name{Evaluate}
\alias{Evaluate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  function to evaluate selected informer set given the DPMM sampling.
}
\description{
Input: DPMM cluster sampling "cl_sample", informer set "inform", measurement criteria "measure", test set, train set, nA, nT, MCMC sample size iter, priors a and b, prior mass alpha
Output: Evaluated value under given criteria
}
\usage{
Evaluate(cl_sample, inform, measure, test, train, nT, iter, a, b, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl_sample}{
%%     ~~Describe \code{cl_sample} here~~
}
  \item{inform}{
%%     ~~Describe \code{inform} here~~
}
  \item{measure}{
%%     ~~Describe \code{measure} here~~
}
  \item{test}{
%%     ~~Describe \code{test} here~~
}
  \item{train}{
%%     ~~Describe \code{train} here~~
}
  \item{nT}{
%%     ~~Describe \code{nT} here~~
}
  \item{iter}{
%%     ~~Describe \code{iter} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
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
function (cl_sample, inform, measure, test, train, nT, iter, 
    a, b, alpha) 
{
    P = clust_sum(cl_sample, train, iter, a, b)
    Score = rep(0, ncol(train))
    xA = test[inform]
    nA = length(inform)
    m = ncol(train)
    post_probs = matrix(0, 1, iter)
    post_thetas = matrix(0, iter, m)
    for (k in 1:iter) {
        postls = pel2_beta(P[[k]], x0 = train, xA, nA, A = inform, 
            nT = 36, a, b, alpha = 15)
        post_probs[k] = postls$post_prob
        post_thetas[k, ] = postls$post_theta
    }
    post_probs = post_probs/(sum(post_probs))
    Score = post_probs \%*\% post_thetas
    Score[inform[which(xA == 1)]] = rep(max(Score) + 1, sum(xA))
    Score[inform[which(xA == 0)]] = rep(min(Score) - 1, nA - 
        sum(xA))
    test = as.vector(test)
    Score = as.vector(Score)
    if (measure == "nef") {
        top = order(Score, decreasing = T)[1:nT]
        pred_hit = sum(test[top])
        hit = sum(test)
        maxhit = min(hit, nT)
        result = ((pred_hit/nT - hit/ncol(train))/(maxhit/nT - 
            hit/ncol(train)) + 1)/2
    }
    else if (measure == "rocauc") {
        rocobj = pROC::roc(test, Score)
        result = rocobj$auc
    }
    else if (measure \%in\% c("mat", "f")) {
        pred.obj = ROCR::prediction(Score, test)
        perform.obj = ROCR::performance(pred.obj, measure)
        result = max(unlist(perform.obj@y.values), na.rm = T)
    }
    else {
        print("Criteria is not supported.")
        result = 0
    }
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
