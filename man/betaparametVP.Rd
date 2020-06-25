\name{betaparametVP}
\alias{betaparametVP}
\title{
Estimation of parameters V and P in count data of RNA reads
}
\description{
This function is used to estimate parameters P and V 
by optimalizing estimates of parameters: alpha and beta. 
}
\usage{
betaparametVP(X, NX)
}

\arguments{
  \item{X}{
count dataset derived from \code{m} replicate libraries 
in one condition.
}
  \item{NX}{
vector of m library sizes. Library size is sum of 
counts over the whole library. 
}
}
\details{
Count data of \emph{RNA} sequence reads are assumed 
to follow binomial distribution with parameters (\code{P}) 
and (\code{n}) or negative binomial distribution with 
parameters (\code{P}) and (\code{r}) , while \code{P} 
is frequency of a gene or an isoform in RNA sequence 
population and assumed to follow beta distribution with 
parameters \code{alpha} (\code{a}) and \code{beta}(\code{b}). 
Parameters  \code{P and V} are estimated by using optimal 
estimation of parameters a and b. The optimal method is an 
iteration algorithm drived by weighting proportion of gene 
or isoform in each replicate library. This is a large-scale 
method for estimating these parameters. Estimation of 
parameters \code{P and V} is core of the multiple beta 
t-test method because \code{P} and \code{V} will be 
used to calculate t-value. 
 
}
\value{
return a list:
\item{P}{N proportions estimated.}
\item{V}{N variances estimated.}

}
\references{
Baggerly KA, Deng L, Morris JS, Aldaz CM (2003) 
Differential expression in SAGE: accounting for normal 
between-library variation. \emph{Bioinformatics}, 
\bold{19}: 1477-1483.
\cr
Yuan-De Tan, Anita M. Chandler, Arindam Chaudhury, 
and Joel R. Neilson(2015) A Powerful Statistical 
Approach for Large-scale Differential Transcription Analysis.
\emph{Plos One},10.1371/journal.pone.0123658. 
}
\author{
Yuan-DE Tan
\email{tanyuande@gmail.com}
}
\note{
betaparametVP requres functions betaparametab and betaparametw.
}

\seealso{
 \code{\link{betaparametab}}, \code{\link{betaparametw}}	
}
\examples{
data(jkttcell) 
X<-jkttcell[1:100,]
na<-3
nb<-3
cn<-length(X[1,])
rn<-length(X[,1])
XC<-X[,1:(cn-na-nb)]
XX<-X[,(cn-na-nb+1):cn]
n<-na+nb
XA<-XX[,1:na]
SA<-apply(XA,2,sum)
PA<-betaparametVP(XA,SA)
}
\keyword{ proportion }
\keyword{ variance }
