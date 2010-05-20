\name{getmultiKmax}
\alias{getmultiKmax}
\alias{getmultiKmax-methods}
\title{Determining the maximum number of segments for joint segmentation}
\description{Determining the total maximum number of segments according to the maximum number of segments for each profile
}
\usage{getmultiKmax(.Object,CGHo,uniKmax=NULL,multiKmax=NULL)}
\arguments{
  \item{.Object}{an object of the CGHdata class}
  \item{CGHo}{an object of the CGHoptions class}
  \item{uniKmax}{NULL if no value is proposed, list of Kmax for each profile}
  \item{multiKmax}{NULL or a proposed value for multiKmax, to check validity}
}
\value{
  \item{multiKmax}{total number of segments for joint segmentation}
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
#simul   = simulprofiles(M=5,n=100,k.mean=2,SNR=5,lambda=1)
#CGHd    = new("CGHdata",Y=simul$Y)
#CGHo    = new("CGHoptions")

## Default initialization

#multiKmax = getmultiKmax(CGHd,CGHo)

## Modification of all fields:
## to set multiKmax to 50\% of the sum of Kmax of each profile:

#beta(CGHo) = 0.5
#multiKmax = getmultiKmax(CGHd,CGHo)

## Modification of uniKmax

#uniKmax      = getuniKmax(CGHd,CGHo)
#uniKmax[[1]] = 80
#multiKmax    = getmultiKmax(CGHd,CGHo,uniKmax)

## Modification of uniKmax and multiKmax
## and check that the value is correct wrt options

#uniKmax[[1]] = 80
#multiKmax    = getmultiKmax(CGHd,CGHo,uniKmax,multiKmax = 10)
}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
\keyword{ts}
\keyword{methods}

