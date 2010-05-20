\name{uniseg}
\alias{uniseg}
\alias{uniseg-methods}
\title{Segmentation of univariate CGH profiles}
\description{Segmentation of univariate CGH profiles}
\usage{uniseg(.Object,CGHo,uniKmax=NULL)}
\arguments{
  \item{.Object}{an object of the CGHdata class}
  \item{CGHo}{an object of the CGHoption class}
  \item{uniKmax}{list containing the maximum number of segments per
    profile. If !is.null(uniKmax) in the entry of the method, the method
    checks its validity.}
}
\value{
  \item{CGHr}{object of class CGHresult}
}
\examples{
#set.seed(1)
#simul = simulprofiles(M=5,n=100,k.mean=2,SNR=5,lambda=1)
#CGHd  = new("CGHdata",Y=simul$Y)
#CGHo  = new("CGHoptions")
#CGHr  = uniseg(CGHd,CGHo)
}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
\keyword{methods}
