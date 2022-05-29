\name{mbetattest}
\alias{mbetattest}
\title{
Performance of multiple beta t-test on count data
}
\description{
This function is used to peform multiple beta t-test 
method on real count data. The result lists "geneid" or 
"isoformid", gene name, the other information, t-value, 
p-value,  rho, and w. 
}
\usage{
mbetattest(X, nci, na, nb, alpha=0.05, norm="no", 
side="both", level="sgRNA",padjust_methods,C=1.222)
}
\arguments{
  \item{X}{
count data of RNA sequence reads with na replicates 
in condition A and nb replicates in condition B.
}
\item{nci}{
nonnegative int value: number of columns for data information, 
such as geneID, isoformID, gene name etc.	
}
  \item{na}{
nonnegative int value: number of replicate libraries in condition A.
}
  \item{nb}{
int numeric value: number of replicate libraries in condition B.
}
  \item{alpha}{
float numeric value, a probabilistic threshold. The value must 
be in [0,1]. User can set alpha=0.05 or 0.01 or the other values. 
Defalt value is 0.05
}
  \item{norm}{
logistic value:"yes" or "no". If norm="yes", the count 
data will be normalized and mbetattest will work on the 
normalized data, if norm="no", then mbetattest will 
work on the unnormalized data.
}

  \item{side}{
string for specifying tail(s) of t-distribution. 
If side="up", then p-value is given with t-test in the 
left tail. If side="down", p-value is given with 
t-test in right tail. If side ="both", p-value is given 
with t-test in both sides.
 
}
  \item{level}{
string for specifying which level mbetattest work on. 
In the current version, level has 6 options: "isoform", 
"sgRNA", "RNA", "splicing.gene","polyA.gene", and 
"CRISPR.gene". 
}
  \item{padjust_methods}{
string for specifying a method for a multiple procedure. 
padjust_methods can choose one of "holm", "hochberg", "hommel", 
"bonferroni", "BH", "BY",  "fdr", TX, and "none" where "fdr" 	
= "BH", "TX" is Tan and Xu's method (2015) with C=1.222 
for adjusting p-value.
}
  \item{C}{
float numeric value for specifying a multiple procedure. 
C=0 tells mbetattest to perform single tests, C=1.222 
tells mbetattest to perform BH correction of pvalues, 
C>1000 tells mbetattest to perform Bonferroni correction 
of pvalues. 
}
}
\details{
see MBttest2-manual. 
}
\value{
return a data and result list: data columns, t-values, rho. 
}

\references{
Baggerly KA, Deng L, Morris JS, Aldaz CM (2003) 
Differential expression in SAGE:  accounting for 
normal between-library variation. 
\emph{Bioinformatics} \bold{19}: 1477-1483.	
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
 \code{\link{smbetattest}}, \code{\link{mtpvadjust}}, 
 \code{\link{normalized}},\code{\link{omega}}.
}
\examples{
data(jkttcell) 
res<-mbetattest(X=jkttcell[1:70, ], nci=7, na=3, 
nb=3, alpha=0.05, norm="yes", side="both", 
level="isoform",padjust_methods="fdr",C=0)
}
\keyword{ beta and negative binomial }
\keyword{ t-tests }