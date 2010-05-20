\name{getuniKmax}
\alias{getuniKmax}
\alias{getuniKmax-methods}
\title{Maximum number of segments per profile}
\description{Maximum number of segments per profile according 
to a given proportion of the length of the data}
\usage{getuniKmax(.Object,CGHo,uniKmax=NULL)}
\arguments{
  \item{.Object}{an object of the CGHdata class}
  \item{CGHo}{an object of the CGHoption class}
  \item{uniKmax}{NULL if no value is proposed, list of Kmax for each profile}
}
\value{
  \item{uniKmax}{list containing the maximum number of segments per
    profile. If !is.null(uniKmax) in the entry of the method, the method
    checks its validity.}
}
\details{
  \code{CGHo["alpha"]} sets the proportion of length of each profile that
  will be used for \code{uniKmax}. If ni is the length of profile i then
  \deqn{\code{uniKmax[[i]]} = \alpha n_i} 
  The speed of execution of the method decreases with Kmax.
  Each Kmax must fulfill some conditions wrt \code{CGHo}: 
  \deqn{\code{uniKmax[[i]]} > \code{nblevels}}
  \deqn{\code{uniKmax[[i]]}< n_i}

}
\examples{
#set.seed(1)
#simul = simulprofiles(M=5,n=100,k.mean=2,SNR=5,lambda=1)
#CGHd  = new("CGHdata",Y=simul$Y)
#CGHo  = new("CGHoptions")

## Default initialization

#uniKmax = getuniKmax(CGHd,CGHo)

## Modification of all fields:
## to set Kmax to 50\% of the length of each profile:

#alpha(CGHo) = 0.5
#uniKmax     = getuniKmax(CGHd,CGHo)

## Modification of one field for one profile
## and check that the value is correct wrt options

#uniKmax[[1]] = 80
#uniKmax      = getuniKmax(CGHd,CGHo,uniKmax)
}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
\keyword{ts}
\keyword{methods}
