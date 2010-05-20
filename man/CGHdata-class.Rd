\name{CGHdata-class}
\docType{class}
\alias{CGHdata-class}
\alias{[,CGHdata-method}
\alias{initialize,CGHdata-method}
\alias{getuniKmax,CGHdata-method}
\alias{getmultiKmax,CGHdata-method}
\alias{multiseg,CGHdata-method}
\alias{uniseg,CGHdata-method}
\alias{print,CGHdata-method}
\alias{show,CGHdata-method}
\alias{summary,CGHdata-method}

\title{Class "CGHdata"}
\description{The CGHdata class is built to store the dataset}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("CGHdata", Y)}. 
  
  \describe{
    \item{\code{Y}:}{dataframe containing recorded values. One column
      per profile (of the same length), with patient IDs as column names. If the data contains
      probeID, then one column of Y must be named \code{probeID}. If the data contains the genomic position, then one
      column of Y must be named \code{genomic.position}. If the data
      contains GC content of the probe then one column of Y must be named
      \code{GCcontent}.}
  }
    
  The dataset is checked while constructing the class. If colnames are
  not \code{GCcontent}, \code{genomic.position}, \code{probeID}, all columns will be
  considered as microarray signals. 
}
\section{Slots}{
  \describe{
    \item{\code{Y}:}{Contains the data one field per patient (list)}
    \item{\code{probeID}:}{Contains the ID of each probe (factor)}     
    \item{\code{genomic.position}:}{Contains the genomic position of each probe (numeric)}
    \item{\code{GCcontent}:}{Contains the GC content of each probe (numeric)}
  }
}
\section{Methods}{
  \describe{
    \item{[}{Accessor to slots of the class}
    \item{initialize}{to construct the class}
    \item{getuniKmax}{Assesses the maximum number of segments per profile}
    \item{getmultiKmax}{Assesses the maximum total number of segments}
    \item{multiseg}{Joint segmentation of multiple CGH profiles}
    \item{uniseg}{Separate segmentation single CGH profiles}
    \item{print}{print class object}
    \item{show}{show class object}
    \item{summary}{summary of class object}
  }
}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
\examples{
#set.seed(1)
#simul          = simulprofiles(M=5,n=100,k.mean=2,SNR=5,lambda=1)
#CGHd           = new("CGHdata",Y=simul$Y)
#show(CGHd)
#summary(CGHd)
}
\keyword{classes}
