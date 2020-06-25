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
C=0, side="both", level="sgRNA")
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
  \item{C}{
real numeric value for specifying a multiple procedure. 
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
}
\details{
see MBttest2-manual. 
}
\value{
return a data and result list: data columns, t-values, rho. 
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
nb=3, alpha=0.05, norm="yes", C=0,side="both", 
level="isoform")
}
\keyword{ beta and negative binomial }
\keyword{ t-tests }