\name{beta<-}
\alias{beta<-}
\alias{beta<--methods}
\title{Replace slot beta}
\description{changes slot \code{beta} of a CGHoptions object}
\arguments{
  \item{.Object}{an object of the CGHoptions class}
}
\details{
  \code{CGHo["beta"]} sets the proportion of each Kmax in \code{multiKmax}. If \code{uniKmax[[i]]} is
  the Kmax of profile i then
  \deqn{\code{multiKmax} = \beta \times \sum_i \code{uniKmax[[i]]}}
  MultiKmax must be greater than the number of profiles, but should not
  be too small since it should be greater than the selected number of
  segments (unknown).
}
\examples{
#set.seed(1)
#simul       = simulprofiles(M=5,n=100,k.mean=2,SNR=5,lambda=1)
#CGHd        = new("CGHdata",Y=simul$Y)
#CGHo        = new("CGHoptions")
#beta(CGHo) = 0.5
}
\seealso{alpha, getmultiKmax}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
\keyword{ts}
\keyword{methods}

