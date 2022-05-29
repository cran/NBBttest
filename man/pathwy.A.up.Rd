\name{pathwy.A.up}
\alias{pathwy.A.up}
\docType{data}
\title{
Pathway or function data
}
\description{
pathway data are functions or pathways of up-regulated 
genes given by gene annotation tool such as DAVID 
gene function annotion.
}
\usage{data("pathwy.A.up")}
\format{
A data frame with 39 observations on the following 159 
variables.
  \describe{ 39 observations, 159 variables             
}
}
\details{
The 39 observations are 39 functions/pathways 
listed each being a set up-regulated genes found 
by NBBttest in knockdowned DDX39 cell line. Log10 
of pvalue of gene in differential expression 
detection is used to weight conpoment of this 
gene playing role in this function. Then heatmap 
of pathways or functions are made using myheatmap. 
V1-V159 are gene sets of function or pathways   
}
\examples{
data(pathwy.A.up)
}
\keyword{datasets}
