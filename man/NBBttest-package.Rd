\name{NBBttest-package}
\alias{NBBtatest-package}
\alias{NBBtatest}
\docType{package}
\title{
Negative Beta Binormail t-Test Package
}
\description{
This package consists of 20 functions, of which 
betaparametab.R, betaparametVP.R, betaparametw.R, 
gbetattest.R,betattest.R, mbetattest.R, normalize.R, 
oddratio.R, pratio.R, omega.R, simulat.R, and 
smbetattest.R are used to estimate beta, alpha, weight, 
t-statistics,rho and omega, p-value,and multiple 
test at gene level or RNA isoform level. NBBplot.R is
used to visualize count data of exons within a specified 
gene in given conditions. QC.R is used to check quality 
of data, and myheatmap.R,myheatmap2.R, and pathwayheatmap.R 
are used to show heatmap of differential expressions of 
DE genes,DE RNA isoform or pathways. Run of mbetattest.R 
would output beta t-test results including geneid or 
isoformid, gene name, the other information, t-value, 
p-value, rho, and omega(W). 
}
\details{
\tabular{ll}{
Package: \tab NBBttest\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2018-01-11\cr
License: \tab GPL-3\cr
}

}
\author{
Yuan-De Tan

Maintainer: Yuan-De Tan 
\email{tanyuande@gmail.com}

}
\references{
Baggerly KA, Deng L, Morris JS, Aldaz CM (2003) 
Differential expression in SAGE: accounting for 
normal between-library variation. \emph{Bioinformatics}, 
\bold{19}: 1477-1483.
\cr
Yuan-De Tan Anita M. Chandler, Arindam Chaudhury, 
and Joel R. Neilson(2015) A Powerful Statistical Approach 
for Large-scale Differential Transcription Analysis.
\emph{Plos One},10.1371/journal.pone.0123658. 

}

\keyword{ package }
\seealso{
\code{\link{betaparametab}}, \code{\link{betaparametVP}}, 
\code{\link{betaparametw}}, \code{\link{gbetattest}},
\code{\link{betattest}}, \code{\link{mbetattest}},
\code{\link{myheatmap}}, \code{\link{pathwayHeatmap}}
\code{\link{myheatmap2}}, \code{\link{normalized}}, 
\code{\link{oddratio}}, \code{\link{NBBplot}},
\code{\link{omega}},\code{\link{pratio}}, 
 \code{\link{QC}},\code{\link{simulat}}, 
\code{\link{smbetattest}}
}
\examples{
data(jkttcell) 	
res<-mbetattest(X=jkttcell[1:100, ], nci=7, na=3, 
nb=3, alpha=0.05, norm="yes", side="both", 
level="isoform",padjust_methods="fdr",C=0)
}
