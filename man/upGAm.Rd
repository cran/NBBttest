\name{upGAm}
\alias{upGAm}
\docType{data}
\title{
Count data of group A treated breast cancer in mice
}
\description{
Data of RNA-seq reads were obtained from drug 
treated breast cancer in mice and mapped by using 
STAR onto mm10 to create count data. Dataset 
upGAm is output of NBBttest.
}
\usage{data("upGAm")}
\format{                                                                                              
A data frame with 263 observations on the following 
19 variables.
  \describe{
    \item{\code{gene}}{a factor with 263 DE genes}
    \item{\code{A.2_S35S}}{a numeric vector}
    \item{\code{A.4_S36S}}{a numeric vector}
    \item{\code{A.42_S45S}}{a numeric vector}
    \item{\code{A.39_S41S}}{a numeric vector}
    \item{\code{A.9_S44R}}{a numeric vector}
    \item{\code{A.12_S40R}}{a numeric vector}
    \item{\code{A.18_S37R}}{a numeric vector}
    \item{\code{A.29_S39R}}{a numeric vector}
    \item{\code{A.31_S38R}}{a numeric vector}
    \item{\code{A.38_S34R}}{a numeric vector}
    \item{\code{tvalue}}{a numeric vector}
    \item{\code{rho}}{a numeric vector}
    \item{\code{pvalue}}{a numeric vector}
    \item{\code{adjp}}{a numeric vector}
    \item{\code{w}}{a numeric vector}
    \item{\code{order}}{a numeric vector}
    \item{\code{FDR}}{a numeric vector}
    \item{\code{significance}}{a numeric vector}
  }
}
\details{
This dataset is a demo dataset having 263 DE genes 
identifiedby NBBttest. The dataset was created from 
types of two tomur cells: sensitive to drug(S) and 
resistant to drug (R). 4 sensitive cells and 6 
resistant cells are available for differential 
analysis. Tvalue is t-statistic, \eqn{\rho} is 
gene-wise variable, pvalue is p-value for t-test, 
\eqn{\omega} is threshold for \eqn{\rho}, FDR is 
false discovery rate, significance = 1 if pvalue < FDR, 
0, otherwise. 
}
\examples{
data(upGAm)
}
\keyword{datasets}
