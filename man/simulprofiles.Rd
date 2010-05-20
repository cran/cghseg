\name{simulprofiles}
\alias{simulprofiles}
\title{simulprofiles}
\description{Simulates multiple CGH profiles}
\usage{
simulprofiles(M, n,k.mean,SNR,lambda)
}
\arguments{
  \item{M}{Number of profiles}
  \item{n}{length of profiles}
  \item{k.mean}{average number of segments per profile}
  \item{SNR}{Signal to Noise Ratio between the mean and the measurement noise. The mean of segments can take 3 values -1;0;0.58.}
  \item{lambda}{Variance ratio between measurement noise and bakground noise}
}
\value{
  \item{simul}{list containing Y (dataframe of simulated CGH signals), Ki (the true number of segments per profile), mu0 (the true signal), theta0 (the true background intensity)}
}
\examples{
#set.seed(1)
#simul          = simulprofiles(M=5,n=100,k.mean=2,SNR=5,lambda=1)
#CGHd           = new("CGHdata",Y=simul$Y)
}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
