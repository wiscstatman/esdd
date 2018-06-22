\name{pkis2}
\alias{pkis2}
\docType{data}
\title{
PKIS2: chemogenomic data of inhibition rates between human kinases and inhibitors
%%   ~~ data name/kind ... ~~
}
\description{
Assay data matrix of inhibition rates between 406 human kinases and 415 inhibitors
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("pkis2")}
\format{
  The format is:
 num [1:406, 1:415] 0.14 0.28 0.2 0.05 0 0 0 0 0.09 0.11 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:406] "AAK1" "ABL1.nonphosphorylated" "ABL1.phosphorylated" "ABL2" ...
  ..$ : chr [1:415] "CID070683337" "CID053466951" "CID044529279" "CID002787267" ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Drewry, D. H., Wells, C. I., Andrews, D. M., Angell, R., Al-Ali, H., Axtman, A. D., ... & Gileadi, O. (2017). Progress towards a public chemogenomic set for protein kinases and a call for contributions. PloS one, 12(8), e0181585.
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(pkis2)
## maybe str(pkis2) ; plot(pkis2) ...
}
\keyword{datasets}
