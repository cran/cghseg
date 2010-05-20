\name{wavenorm<-}
\alias{wavenorm<-}
\alias{wavenorm<--methods}
\title{Replace slot wavenorm}
\description{changes slot \code{wavenorm} of a CGHoptions object}
\arguments{
  \item{.Object}{an object of the CGHoptions class}
}
\details{
  \code{CGHo["wavenorm"]} determines the method to remove the
  wave-effect: \code{none}, \code{position}, \code{spline}.
}
\examples{
#set.seed(1)
#simul          = simulprofiles(M=5,n=100,k.mean=2,SNR=5,lambda=1)
#CGHd           = new("CGHdata",Y=simul$Y)
#CGHo           = new("CGHoptions")
#wavenorm(CGHo) = "spline"
#CGHr           = multiseg(CGHd,CGHo)
}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
\keyword{ts}
\keyword{methods}
