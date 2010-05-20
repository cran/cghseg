\name{getbgoutliers}
\alias{getbgoutliers}
\alias{getbgoutliers-methods}
\title{getbgoutliers}
\description{Identification of outliers in the background intensity}
\usage{getbgoutliers(.Object,CGHr.smooth,fdr)}
\arguments{
  \item{.Object}{an object of the CGHresults class for which \code{CGHo["wavenorm"]="position"}} 
  \item{CGHr.smooth}{an object of the CGHresults class for which
    \code{CGHo["wavenorm"]} is \code{"spline"}} 
  \item{fdr}{fdr level to control for multiple testing}
}
\value{
  \item{bgoutliers}{dataframe containing background outliers}
}
\examples{
#set.seed(1)
#simul          = simulprofiles(M=5,n=100,k.mean=2,SNR=5,lambda=1)
#CGHd           = new("CGHdata",Y=simul$Y)
#CGHo           = new("CGHoptions")
#calling(CGHo)  = FALSE
#wavenorm(CGHo) = "position"
#CGHr.pos       = multiseg(CGHd,CGHo)
#wavenorm(CGHo) = "spline"
#CGHr.spl       = multiseg(CGHd,CGHo)
#getbgoutliers(CGHr.pos,CGHr.spl,fdr=0.05)
}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
\keyword{methods}

