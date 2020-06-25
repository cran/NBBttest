\name{prime3_PRP8_50}
\alias{prime3_PRP8_50}
\docType{data}
\title{
3'UTR splicing data of 50 genes detected in the 
knockdowned PRP8 cell line
}
\description{
Splicing events occurring in 3'UTR of 50 genes 
were detected in the knockdowned PRP8 cell line.
}
\usage{data("prime3_PRP8_50")}
\format{
  A data frame with 50 observations on the following 24 variables.
  \describe{
    \item{\code{chr}}{ a factor with levels}
    \item{\code{strand}}{a factor with levels}
    \item{\code{isoform}}{a factor with levels}
    \item{\code{Gene}}{a factor with levels}
    \item{\code{exon_const_start}}{a numeric vector}
    \item{\code{exon_const_end}}{a numeric vector}
    \item{\code{exon_alt1_start}}{a numeric vector}
    \item{\code{exon_alt1_end}}{a numeric vector}
    \item{\code{exon_alt2_start}}{a numeric vector}
    \item{\code{exon_alt2_end}}{a numeric vector}
    \item{\code{Siha_Ctl_I_S3}}{a numeric vector}
    \item{\code{Siha_Ctl_II_S7}}{a numeric vector}
    \item{\code{Siha_Ctl_III_S11}}{a numeric vector}
    \item{\code{PRP8_I_S2}}{a numeric vector}
    \item{\code{PRP8_II_S6}}{a numeric vector}
    \item{\code{PRP8_III_S10}}{a numeric vector}
    \item{\code{tvalue}}{a numeric vector}
    \item{\code{rho}}{a numeric vector}
    \item{\code{pvalue}}{a numeric vector}
    \item{\code{adjp}}{a numeric vector}
    \item{\code{w}}{a numeric vector}
    \item{\code{X}}{a numeric vector}
    \item{\code{X.1}}{a numeric vector}
    \item{\code{X.2}}{a numeric vector}
  }
}
\details{
This dataset is an example for annotation of genes within 
which splicing events occurring in 3'UTR were detected by 
spladder-NBBttest in knockdowned PRP8 cell line. Gene 
information columns are assigned to "chr","Gene", "isoform", 
"exon_const_start" and "exon_const_end" that are start 
and end positions of constituted exons and "exon_alt1_start" 
and "exon_alt1_end" that are start and end positions of 
alternative exons. Count data of RNA reads are assigned 
to three replicate control columns (Siha_Ctl_I_S3,II-S7,
III_S11) and three replicate knokdowned PRP8 columns 
(PRP8_I_S2,II_S6 and III_S10). The tvalue is t-statistic, 
rho is gene-wise variable, pvalue is p-value for t-statistic, 
w is \eqn{\omega},a threshold for rho, FDR (X.1)is false 
discovery rate with \eqn{\alpha}= 0.05, selection (X.2)= 1 
if pvalue < FDR, 0, otherwise.    
}

\examples{
data(prime3_PRP8_50)
}
\keyword{datasets}
