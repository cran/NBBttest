\name{gtfa}
\alias{gtfa}
\docType{data}
\title{
gtf
}
\description{
gtf is annotation file for genes and exons
}
\usage{data("gtfa")}
\format{
  A data frame with 1461 observations on the following 8 variables.
  \describe{
    \item{\code{V1}}{a factor with levels}
    \item{\code{V2}}{a factor with levels}
    \item{\code{V3}}{a numeric vector}
    \item{\code{V4}}{a numeric vector}
    \item{\code{V5}}{a factor with levels}
    \item{\code{V6}}{a factor with levels }
    \item{\code{V7}}{a numeric vector}
    \item{\code{V8}}{a factor with levels}
  }
}
\details{
V1:chr,V2:exon,V3:start position,V4:end position,V5:strand,V6:gene id,V7: exon,V8:gene
}
\examples{
data(gtfa)
}
\keyword{datasets}
