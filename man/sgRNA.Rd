\name{sgRNA}
\alias{sgRNA}
\docType{data}
\title{
sgRNA dataset
}
\description{
This dataset was created by simulating single 
guide RNAs to edit genes. 
}
\usage{data("sgRNA")}
\format{
A data frame with 1000 observations on the 
following 11 variables.
  \describe{
    \item{\code{gene}}{a factor with levels }
    \item{\code{sgRNA}}{a factor with levels}
    \item{\code{class}}{a factor with levels}
    \item{\code{L1}}{a numeric vector}
    \item{\code{L2}}{a numeric vector}
    \item{\code{L3}}{a numeric vector}
    \item{\code{L4}}{a numeric vector}
    \item{\code{H1}}{a numeric vector}
    \item{\code{H2}}{a numeric vector}
    \item{\code{H3}}{a numeric vector}
    \item{\code{H4}}{a numeric vector}
  }
}
\details{
The dataset is CRISPR screening data. Each gene was 
edited by 10 sgRNAs that contains a targeting sequence 
(crRNA sequence) and a Cas9 nuclease-recruiting sequence 
(tracrRNA). L1-L4 samples were lowly targeted by sgRNAs 
and H1-H4 were highly targeted by sgRNAs.
}

\examples{
data(sgRNA)
}
\keyword{datasets}
