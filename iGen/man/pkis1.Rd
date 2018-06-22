\name{pkis1}
\alias{pkis1}
\docType{data}
\title{
PKIS1: chemogenomic data of inhibition rates between human kinases and inhibitors
%%   ~~ data name/kind ... ~~
}
\description{Assay data matrix of inhibition rates between 224 human kinases and 366 inhibitors
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("pkis1")}
\format{
  The format is:
 num [1:224, 1:366] 0.1397 0.0568 0.1041 0.1502 0.0924 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:224] "ABL1" "ABL1_E255K" "ABL1_H396P" "ABL1_M351T" ...
  ..$ : chr [1:366] "6539592" "6539600" "6539599" "6539597" ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
H Drewry, D., M Willson, T., & J Zuercher, W. (2014). Seeding collaborations to advance kinase science with the GSK Published Kinase Inhibitor Set (PKIS). Current topics in medicinal chemistry, 14(3), 340-342.
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(pkis1)
## maybe str(pkis1) ; plot(pkis1) ...
}
\keyword{datasets}
