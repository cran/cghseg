\name{nbprocs<-}
\alias{nbprocs<-}
\alias{nbprocs<--methods}
\title{Replace slot nbprocs}
\description{changes slot \code{nbprocs} of a CGHoptions object}
\arguments{
  \item{.Object}{an object of the CGHoptions class}
}
\details{
  \code{CGHo["nbprocs"]} sets the number of computer processes to use for parallel computations. 
}
\examples{
#CGHo         = new("CGHoptions")
#nbprocs(CGHo) = 2 # if 2 processes are available
}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
\keyword{ts}
\keyword{methods}
