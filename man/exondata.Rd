\name{exondata}
\alias{exondata}
\docType{data}
\title{
Exon data for NBBplot
}
\description{
This dataset is an example for NBBplot to show 
differential exons within genes between two 
conditions.
}
\usage{data("exondata")}
\format{
  A data frame with 16 observations on the 
  following 22 variables.
  \describe{
    \item{\code{ID}}{a factor with levels 
    	\code{ENSMUSG00000015597:E1} 
    	\code{ENSMUSG00000015597:E10} 
    	\code{ENSMUSG00000015597:E2} 
    	\code{ENSMUSG00000015597:E3} 
    	\code{ENSMUSG00000015597:E4} 
    	\code{ENSMUSG00000015597:E5} 
    	\code{ENSMUSG00000015597:E6} 
    	\code{ENSMUSG00000015597:E7} 
    	\code{ENSMUSG00000015597:E8} 
    	\code{ENSMUSG00000015597:E9} 
    	\code{ENSMUSG00000079547:E1} 
    	\code{ENSMUSG00000079547:E2} 
    	\code{ENSMUSG00000079547:E3} 
    	\code{ENSMUSG00000079547:E4} 
    	\code{ENSMUSG00000079547:E5} 
    	\code{ENSMUSG00000079547:E6}}
    \item{\code{chr}}{a numeric vector}
    \item{\code{element}}{a factor with levels \code{exon}}
    \item{\code{start}}{a numeric vector}
    \item{\code{end}}{a numeric vector}
    \item{\code{strand}}{a factor with levels \code{+}}
    \item{\code{gene_id}}{a factor with levels 
    	\code{ENSMUSG00000015597} 
    	\code{ENSMUSG00000079547}}
    \item{\code{exon}}{a numeric vector}
    \item{\code{gene}}{a factor with levels 
    	\code{H2-DMb1} 
    	\code{Zfp318}}
    \item{\code{WT_28}}{a numeric vector}
    \item{\code{WT_43}}{a numeric vector}
    \item{\code{A25_f0}}{a numeric vector}
    \item{\code{cKO_30}}{a numeric vector}
    \item{\code{A26_cKO}}{a numeric vector}
    \item{\code{A6_f0}}{a numeric vector}
    \item{\code{tv}}{a numeric vector}
    \item{\code{rho}}{a numeric vector}
    \item{\code{pvalue}}{a numeric vector}
    \item{\code{w}}{a numeric vector}
    \item{\code{order_number}}{a numeric vector}
    \item{\code{FDR0.05}}{a numeric vector}
    \item{\code{significance}}{a numeric vector}
  }
}
\details{
This dataset is an object of NBBttest and annoted 
with a R function. It contains columns start and 
end, gene,isoform, tv and pvalue, data, significance 
that NBBplut uses. WT and A25 are count data from 
wild type cell lines 28,43 and fo and A and cKO are 
also count data from knockout a gene, tv is t-statistic, 
rho is exon-wise variable and w is threshold of rho. 
Column "significance" is 1 if pvalue < FDR0.05, 0, 
otherwise. Gene Zfp318 has 10 exons and H2-DMb1 has 
6 exons.

}

\examples{
data(exondata)
}
\keyword{datasets}
