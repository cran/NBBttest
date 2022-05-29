\name{gbetattest}
\alias{gbetattest}
\title{
Beta t-tests within groups
}
\description{
Beta t-tests are conducted within groups,genes,or libraries. 
}
\usage{
gbetattest(xx, W, nci, na, nb, level, padjust_methods,C=1.222, side)
}

\arguments{
  \item{xx}{
a datasheet consisting of \code{n} columns and \code{m}
rows. Columns contain information and count data columns 
\code{n} must be 1 or more and \code{m} must be over 100.
}
  \item{W}{
numeric value. It is omega estimated from null 
simulation.
}
  \item{nci}{
int numeric value indicating number of  
information columns.
}
  \item{na}{
int numeric value indicating number of replicates 
in condition a.
}
  \item{nb}{
int numeric value indicating number of replicates 
in condition b.
}
  \item{level}{
string value. It has 6 options: "isoform", "sgRNA", 
"RNA", "polyA.gene", "CRISPR.gene" and "splicing.gene".
}
  \item{padjust_methods}{
 padjust.methods can choose one of "holm", "hochberg", 
 "hommel", "bonferroni", "BH", "BY",  "fdr", TX, and "none" 
 where "fdr" 	= "BH", "TX" is Tan and Xu's method (2015) 
 with C=1.222 for adjusting p-value.
	
}
  \item{C}{
float numeric value for specifying a multiple procedure. 
C=0 tells mbetattest to perform single tests, C=1.222 
tells mbetattest to perform BH correction of pvalues, 
C>1000 tells mbetattest to perform Bonferroni correction 
of pvalues.
}
  \item{side}{
string value for specifying one-tail test or two-tail test: 
side="up" for left-tail test, side="down" for right-tail 
test and side="both" for two-tail tests. 
}
}
\details{
Beta t-test will be conducted within a specified group or at
a specified level. If level="RNA", then beta t-tests will be 
conducted within a whole library or the whole data. If level=
"isoform", then data will be sparated in two parts: 
single-isoform and multi-isoform datasets. Single-isoform RNA
indicates that there is only one RNA isform within a gene, 
while multi-isoform RNAs indicate that there are at least two 
RNA isoformswithin a gene. For single-isoforms, data are as a 
group and beta t-tests will be performed in the group. 
For the multi-isoforms, t-test will be performed within genes. 
If level="polyA.gene" or "CRISPR.gene", then t-test will be 
performed at gene level. If level="splicing.gene", then t-values 
and p-values will be selected from gene groups with the least 
p-values.
}
\value{
return a list containing dataset, t-values, corrected p-values, 
rhos and w.
}
\references{
Baggerly KA, Deng L, Morris JS, Aldaz CM (2003) 
Differential expression in SAGE: accounting for 
normal between-library variation. 
\cr
Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery 
rate: a practical and powerful approach to multiple testing. 
\emph{Journal of the Royal Statistical Society Series}, \bold{B}, \bold{57}5, 289-300. 
\url{doi:10.1111/j.2517-6161.1995.tb02031.x},
\url{https://www.jstor.org/stable/2346101}.	
\cr
Benjamini, Y., and Yekutieli, D. (2001). The control of the false 
discovery rate in multiple testing under dependency. 
\emph{Annals of Statistics}, \bold{29}, 1165 -1188. 
\url{doi:10.1214/aos/1013699998}.
\cr
Tan YD, Xu H. A general method for accurate estimation of false discovery 
rates in identification of differentially expressed genes. 
\emph{Bioinformatics}. 2014 Jul 15;30(14):2018-25. 
\url{doi:10.1093/bioinformatics/btu124}. Epub 2014 Mar 14. PMID: 24632499.
}
\author{
Yuan-De Tan
\email{tanyuande@gmail.com}
}

\seealso{
\code{\link{betattest}}

}
\examples{
data(jkttcell)
colnames(jkttcell)[3]<-"Gene"
res.isfo<-gbetattest(xx=jkttcell[1:100,], W=1, nci=7,
na=3, nb=3, level="isoform", padjust_methods="fdr",C=0,side="both")
	
}
\keyword{group }
\keyword{beta t-test }
