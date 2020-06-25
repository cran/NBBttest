\name{simSplicing}
\alias{simSplicing}
\docType{data}
\title{
Simulated alternative splicing
}
\description{
This alternative splicing count data were created 
by simulating dorsal and ventral RNA-sequence 
count data.
}
\usage{data("simSplicing")}
\format{
A data frame with 5000 observations on the following 
9 variables.
  \describe{
    \item{\code{Isoform}}{a factor with levels}
    \item{\code{geneid}}{a factor with levels}
    \item{\code{label}}{a numeric vector}
    \item{\code{D1}}{a numeric vector}
    \item{\code{D2}}{a numeric vector}
    \item{\code{D3}}{a numeric vector}
    \item{\code{V1}}{a numeric vector}
    \item{\code{V2}}{a numeric vector}
    \item{\code{V3}}{a numeric vector}
  }
}
\details{
D1-D3 are samples from mouse dorsal tissue and V1-V3 
are samples from mouse ventral tissue. This dataset 
was simulated with negative binomial ditribution with 
10 percent of isoforms of being differentially spliced  
and differential effect of 500U where U is uniform variable.
}

\examples{
data(simSplicing)
}
\keyword{datasets}
