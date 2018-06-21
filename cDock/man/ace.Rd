\name{ace}
\alias{ace}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
A data.frame of 8 docking scores and binary label between a library of DUD-E compounds and protein target ace
}
\usage{data("ace")}
\format{
  A data frame with 16413 observations on the following 9 variables.
  \describe{
    \item{\code{label}}{binary label: 0 for inactive; 1 for active}
    \item{\code{dock6}}{a numeric vector: docking score via using docking program DOCK v6.7}
    \item{\code{rdock}}{a numeric vector: docking score via using docking program rDock v2013.1}
    \item{\code{surf}}{a numeric vector: docking score via using docking program Surflex-Dock v3.040}
    \item{\code{hybrid}}{a numeric vector: docking score via using docking program HYBRID v3.0.1}
    \item{\code{ad4}}{a numeric vector: docking score via using docking program AutoDock v4.2.6}
    \item{\code{fred}}{a numeric vector: docking score via using docking program FRED v3.0.1}
    \item{\code{plants}}{a numeric vector: docking score via using docking program PLANTS v1.2}
    \item{\code{smina}}{a numeric vector: docking score via using docking program Smina v1.1.2}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(ace)
## maybe str(ace) ; plot(ace) ...
}
\keyword{datasets}
