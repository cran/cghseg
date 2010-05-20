\name{getsegprofiles}
\alias{getsegprofiles}
\alias{getsegprofiles-methods}
\title{getsegprofiles}
\description{Extract segmented profiles from a CGHresults object}
\usage{getsegprofiles(.Object)}
\arguments{
  \item{.Object}{an object of the CGHresults class}
}
\value{
  \item{segprofiles}{list of segmented profiles}
}
\examples{
#set.seed(1)
#simul          = simulprofiles(M=5,n=100,k.mean=2,SNR=5,lambda=1)
#CGHd           = new("CGHdata",Y=simul$Y)
#CGHo           = new("CGHoptions")
#calling(CGHo)  = FALSE
#CGHr           = multiseg(CGHd,CGHo)
#segprofiles    = getsegprofiles(CGHr)

#calling(CGHo)  = TRUE
#CGHr           = multiseg(CGHd,CGHo)
#segprofiles    = getsegprofiles(CGHr)
}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
\keyword{methods}

