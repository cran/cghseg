\name{CGHoptions-class}
\docType{class}
\alias{CGHoptions-class}
\alias{[,CGHoptions-method}
\alias{calling<-,CGHoptions-method}
\alias{alpha<-,CGHoptions-method}
\alias{beta<-,CGHoptions-method}
\alias{nblevels<-,CGHoptions-method}
\alias{wavenorm<-,CGHoptions-method}
\alias{GCnorm<-,CGHoptions-method}
\alias{itermax<-,CGHoptions-method}
\alias{print,CGHoptions-method}
\alias{select<-,CGHoptions-method}
\alias{nbprocs<-,CGHoptions-method}
\alias{show,CGHoptions-method}

\title{Class "CGHoptions"}
\description{The CGHoptions class is built to store informations
  regarding the methods used in the \code{cghseg} package}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CGHoptions")}.
}
\section{Slots}{
	 \describe{
    \item{\code{select}:}{Model selection to select the number of
      segments; \code{"mBIC"},\code{"none"} (character)}
    \item{\code{calling}:}{calling of segments (logical)}
    \item{\code{wavenorm}:}{Wave normalization; \code{"none"},
      \code{"position"}, \code{"spline"} (character)}
    \item{\code{GCnorm}:}{GCcontent normalization;
      \code{"none"},\code{"linear"} (character)}
    \item{\code{nblevels}:}{Number of levels when calling (numeric)}
    \item{\code{alpha}:}{Proportion of the individual signal length that gives uniKmax.}
    \item{\code{beta}:}{Proportion of the sum of uniKmax that gives multiKmax.}
    \item{\code{itermax}:}{Maximal number of iterations}
    \item{\code{nbprocs}:}{Number of available computer processes}
  }
}
\section{Prototype}{
\code{CGHo             = new('CGHoptions')}\cr
\code{CGHo@select      = "mBIC"}\cr
\code{CGHo@calling     = FALSE}\cr
\code{CGHo@wavenorm    = "position"}\cr
\code{CGHo@GCnorm      = "none"}\cr
\code{CGHo@nblevels    = 3}\cr
\code{CGHo@alpha       = 0.2}\cr
\code{CGHo@beta        = 0.5}\cr
\code{CGHo@itermax     = Inf}\cr
\code{CGHo@nbprocs     = 1}\cr
}
\section{Methods}{
  \describe{
    \item{[}{Accessor to slots of the class}
    \item{initialize}{to construct the class}
    \item{print}{Print CGHo}
    \item{show}{Show CGHo}
    \item{select<-}{Modification of slot select; \code{select(CGHo)<-}}
    \item{calling<-}{Modification of slot calling; \code{calling(CGHo)<-}}
    \item{wavenorm<-}{Modification of slot wavenorm; \code{wavenorm(CGHo)<-}}
    \item{GCnorm<-}{Modification of slot GCnorm; \code{GCnorm(CGHo)<-}}
    \item{nblevels<-}{Modification of slot nblevels; \code{nblevels(CGHo)<-}}
    \item{alpha<-}{Modification of slot alpha; \code{alpha(CGHo)<-}}
    \item{beta<-}{Modification of slot beta; \code{beta(CGHo)<-}}
    \item{itermax<-}{Modification of slot itermax;
      \code{itermax(CGHo)<-}}
    \item{nbprocs<-}{Modification of slot nbproc; \code{nbproc(CGHo)<-}}
  }
}
\references{Joint segmentation, calling and normalization of multiple CGH profiles, Biostatistics (2011)}
\author{F. Picard, E. Lebarbier, M. Hoebeke, G. Rigaill, B. Thiam, S. Robin}
\examples{
#CGHo = new("CGHoptions")
#showClass("CGHoptions")
}
\keyword{classes}
