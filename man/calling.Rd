\name{calling<-}
\alias{calling<-}
\alias{calling<--methods}
\title{Replace slot calling}
\description{changes slot \code{calling} of a CGHoptions object}
\arguments{
  \item{.Object}{an object of the CGHoptions class}
}
\details{
  \code{CGHo["calling"]} determines if calling should be done during
  segmentation (logical). The number of levels is fixed and set in \code{CGHo["nblevels"]} 
}
\examples{
#set.seed(1)
#simul         = simulprofiles(M=5,n=100,k.mean=2,SNR=5,lambda=1)
#CGHd          = new("CGHdata",Y=simul$Y)
#CGHo          = new("CGHoptions")
#calling(CGHo) = TRUE
#CGHr          = multiseg(CGHd,CGHo)
}
\seealso{nblevels}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
\keyword{ts}
\keyword{methods}
