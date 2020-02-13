\name{multiseg}
\alias{multiseg}
\alias{multiseg-methods}
\title{Joint segmentation of multiple CGH profiles}
\description{Joint segmentation of multiple CGH profiles}
\usage{multiseg(.Object,CGHo,uniKmax=NULL,multiKmax=NULL)}
\arguments{
  \item{.Object}{an object of the CGHdata class}
  \item{CGHo}{an object of the CGHoption class}
  \item{uniKmax}{list containing the maximum number of segments per
    profile. If !is.null(uniKmax) in the entry of the method, the method
    checks its validity.}
  \item{multiKmax}{total number of segments for joint segmentation}
}
\value{
  \item{CGHr}{object of class CGHresult}
}
\examples{
set.seed(1)
simul          = simulprofiles(M=5,n=100,k.mean=2,SNR=5,lambda=1)
CGHd           = new("CGHdata",Y=simul$Y)
CGHo           = new("CGHoptions")
CGHr           = multiseg(CGHd,CGHo)
}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
\keyword{methods}

