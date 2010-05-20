\name{CGHresults-class}
\docType{class}
\alias{CGHresults-class}
\alias{[,CGHresults-method}
\alias{getbp,CGHresults-method}
\alias{initialize,CGHresults-method}
\alias{getprobebias,CGHresults-method}
\alias{getsegprofiles,CGHresults-method}
\alias{getlevels,CGHresults-method}
\alias{getbgoutliers,CGHresults-method}
\alias{summary,CGHresults-method}

\title{Class "CGHresults"}
\description{The CGHresults class is built to store the results of \code{cghseg}.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("CGHresults", CGHd, CGHo)}.
}
\section{Slots}{
  \describe{
    \item{\code{mu}:}{list of dataframes containing the position of breakpoints and estimated means (list)}
    \item{\code{theta}:}{List of effects used for normalization :
      waveffect, GCeffect (listOrNULL)}
    \item{\code{loglik}:}{Log-likelihood of the selected model (numeric)}
    \item{\code{nbiter}:}{number of iterations (numeric)}
    \item{\code{from}:}{Name of the function that produced the results (character)}
    \item{\code{probeID}:}{Contains the ID of each probe (factor)}     
    \item{\code{genomic.position}:}{Contains the genomic position of each probe (numeric)}
    \item{\code{options}:}{list of options (of class CGHoptions)}
  }
}
\section{Methods}{
  \describe{
    \item{[}{Accessor to slots of the class}
    \item{initialize}{to construct the class}
    \item{summary}{summary of CGHr}        
    \item{getbp}{Extract the frequency of breakpoints}
    \item{getsegprofiles}{Extract model predictions}
    \item{getlevels}{Extract levels of segments}
    \item{getbgoutliers}{Identification of outliers in the background intensity}
  }
}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}\examples{
#set.seed(1)
#simul          = simulprofiles(M=5,n=100,k.mean=2,SNR=5,lambda=1)
#CGHd           = new("CGHdata",Y=simul$Y)
#CGHo           = new("CGHoptions")
#CGHr           = multiseg(CGHd,CGHo)
#bplist         = getbp(CGHr)
#segprofiles    = getsegprofiles(CGHr)

#calling(CGHo)  = TRUE
#CGHr           = multiseg(CGHd,CGHo)
#bplist         = getbp(CGHr)
#levelslist     = getlevels(CGHr)
#segprofiles    = getsegprofiles(CGHr)
}
\keyword{classes}
