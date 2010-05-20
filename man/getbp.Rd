\name{getbp}
\alias{getbp}
\alias{getbp-methods}
\title{Breakpoints frequencies}
\description{Gives the frequency of breakpoints}
\usage{getbp(.Object)}
\arguments{
  \item{.Object}{object of class CGHresults}
}
\value{
  \item{bplist}{list of dataframes containing 0/1 for the position of the breaks} 
}
\examples{
#set.seed(1)
#simul          = simulprofiles(M=5,n=100,k.mean=2,SNR=5,lambda=1)
#CGHd           = new("CGHdata",Y=simul$Y)
#CGHo           = new("CGHoptions")
#CGHr           = multiseg(CGHd,CGHo)
#bplist         = getbp(CGHr)
}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
\keyword{ts}
\keyword{methods}
