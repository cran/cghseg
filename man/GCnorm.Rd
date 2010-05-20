\name{GCnorm<-}
\alias{GCnorm<-}
\alias{GCnorm<--methods}
\title{Replace slot GCnorm}
\description{changes slot \code{GCnorm} of a CGHoptions object}
\arguments{
  \item{.Object}{an object of the CGHoptions class}
}
\details{ \code{CGHo["GCnorm"]} determines if the bias due to GC content
should be removed from the data: \code{"linear"} for linear regression
(degree 2) \code{"none"} for no correction}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
\keyword{ts}
\keyword{methods}
