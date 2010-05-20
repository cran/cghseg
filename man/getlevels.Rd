\name{getlevels}
\alias{getlevels}
\alias{getlevels-methods}
\title{Levels frequency}
\description{Gives the calls status at every position}
\usage{getlevels(.Object)}
\arguments{
  \item{.Object}{an object of the CGHresults class}
}
\value{
  \item{levelslist}{list of dataframes containing calls for each position}
}

\examples{
#set.seed(1)
#simul          = simulprofiles(M=5,n=100,k.mean=2,SNR=5,lambda=1)
#CGHd           = new("CGHdata",Y=simul$Y)
#CGHo           = new("CGHoptions")
#calling(CGHo)  = TRUE
#CGHr           = multiseg(CGHd,CGHo)
#levelslist     = getlevels(CGHr)
}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
\keyword{methods}

