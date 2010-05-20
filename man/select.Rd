\name{select<-}
\alias{select<-}
\alias{select<--methods}
\title{Replace slot select}
\description{changes slot \code{select} of a CGHoptions object}
\arguments{
  \item{.Object}{an object of the CGHoptions class}
}
\details{
  \code{CGHo["select"]} sets the model selection criterion to be used in the
  procedure. The default criterion is based on the mBIC proposed by
  Zhang & Siegmund (2007). If  \code{CGHo["select"]=="none"} then the
  algorithm is run for the maximum number of segments.
}
\examples{
#set.seed(1)
#simul        = simulprofiles(M=5,n=100,k.mean=2,SNR=5,lambda=1)
#CGHd         = new("CGHdata",Y=simul$Y)
#CGHo         = new("CGHoptions")
#select(CGHo) = "mBIC"
#CGHr         = multiseg(CGHd,CGHo)

#select(CGHo) = "none"
#CGHr         = multiseg(CGHd,CGHo)
}
\seealso{alpha, beta, getmultiKmax, getuniKmax}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
\keyword{ts}
\keyword{methods}
