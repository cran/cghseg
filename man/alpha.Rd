\name{alpha<-}
\alias{alpha<-}
\alias{alpha<--methods}
\title{Replace slot alpha}
\description{changes slot \code{alpha} of a CGHoptions object}
\arguments{
  \item{.Object}{an object of the CGHoptions class}
}
\examples{
set.seed(1)
simul       = simulprofiles(M=5,n=100,k.mean=2,SNR=5,lambda=1)
CGHd        = new("CGHdata",Y=simul$Y)
CGHo        = new("CGHoptions")
alpha(CGHo) = 0.5
}
\details{
   \code{CGHo["alpha"]} sets the proportion of length of each profile that
  will be used for . If ni is the length of profile i then
  \deqn{\code{uniKmax[[i]]} = \alpha n_i} 
  The speed of execution of the method decreases with Kmax.
  Each Kmax must fulfill some conditions wrt \code{CGHo}: 
  \deqn{\code{uniKmax[[i]]} > \code{nblevels}}
  \deqn{\code{uniKmax[[i]]}< n_i} 
}
\seealso{beta, getuniKmax}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
\keyword{ts}
\keyword{methods}

