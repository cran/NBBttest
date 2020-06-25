\name{DDX39_100}
\alias{DDX39_100}
\docType{data}
\title{
DDX39 exon data with 100 exons  
}
\description{
DDX39_100 dataset is an example for implementing annotation of genes within there are some exons with gtf.
}
\usage{data("DDX39_100")}
\format{
  A data frame with 100 observations on the following 14 variables.
  \describe{
    \item{\code{X}}{a numeric vector}
    \item{\code{gene_id}}{a factor with ENSG levels}
    \item{\code{CTL_1_S3}}{a numeric vector}
    \item{\code{CTL_2_S7}}{a numeric vector}
    \item{\code{CTL_3_S11}}{a numeric vector}
    \item{\code{DDX39_1_S1}}{a numeric vector}
    \item{\code{DDX39_2_S5}}{a numeric vector}
    \item{\code{DDX3_S9}}{a numeric vector}
    \item{\code{tvalue}}{a numeric vector}
    \item{\code{rho}}{a numeric vector}
    \item{\code{df}}{a numeric vector}
    \item{\code{pvalue}}{a numeric vector}
    \item{\code{adjp}}{a numeric vector}
    \item{\code{w}}{a numeric vector}
  }
}
\details{
DDX39_100 is dataset that was RNA-seq data from human cell breast cells with knockout DDX39 gene. DDX39 RNA-seq reads were mapped to human reference annotated with GRCh38 by HTseq with DEXseq annotation and generated count data of RNA reads of exons within genes. DDX39_100 is output of NBBttest where X is row order column; gene_id is ENSG ID; CTL is control cells with replicates 1, 2,and 3; DDX39 is treatment cells with knockout gene DDX39; tvalue is NBBttest t-statistic; rho is exon-wise variable of NBBttest; df is degree of freedom; pvalue is p-value for t-statistics and w is threshold of rho.    
}
\examples{
data(DDX39_100)
}
\keyword{datasets}
