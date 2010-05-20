\name{itermax<-}
\alias{itermax<-}
\alias{itermax<--methods}
\title{Replace slot itermax}
\description{changes slot \code{itermax} of a CGHoptions object}
\arguments{
  \item{.Object}{an object of the CGHoptions class}
}
\examples{
set.seed(1)
simul         = simulprofiles(M=5,n=100,k.mean=2,SNR=5,lambda=1)
CGHd          = new("CGHdata",Y=simul$Y)
CGHo          = new("CGHoptions")
itermax(CGHo) = 10
}
\details{
  \code{CGHo["itermax"]} sets the number of iterations when using wave correction
}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
\keyword{ts}
\keyword{methods}

