\name{pel1_beta}
\alias{pel1_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pel1_beta(cl_sample, P, iter, size, A, nA, nT = 10, a, b, x0, alpha = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl_sample}{
%%     ~~Describe \code{cl_sample} here~~
}
  \item{P}{
%%     ~~Describe \code{P} here~~
}
  \item{iter}{
%%     ~~Describe \code{iter} here~~
}
  \item{size}{
%%     ~~Describe \code{size} here~~
}
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{nA}{
%%     ~~Describe \code{nA} here~~
}
  \item{nT}{
%%     ~~Describe \code{nT} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{x0}{
%%     ~~Describe \code{x0} here~~
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
function (cl_sample, P, iter, size, A, nA, nT = 10, a, b, x0, 
    alpha = 2) 
{
    XA = matrix(0, iter * size, nA)
    for (i in 1:iter) {
        XA[(size * i - size + 1):(size * i), ] = cl_sample$XX[i, 
            , A]
    }
    XA = apply(XA, 1, function(x) {
        return(paste(x, collapse = ""))
    })
    tab = table(XA)
    YA = names(tab)
    l = length(tab)
    tmp_pel2 = sapply(1:l, function(x) {
        xA = as.numeric(strsplit(YA[x], "")[[1]])
        post_probs = matrix(0, 1, iter)
        post_thetas = matrix(0, iter, ncol(x0))
        for (k in 1:iter) {
            postls = pel2_beta(P[[k]], x0, xA, nA, A, nT, a, 
                b, alpha)
            post_probs[k] = postls$post_prob
            post_thetas[k, ] = postls$post_theta
        }
        post_probs = post_probs/(sum(post_probs))
        Score = post_probs \%*\% post_thetas
        return(sum(sort(1 - Score)[1:nT]))
    })
    pel1 = (tmp_pel2 * tab)/iter
    return(sum(pel1)/size)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
