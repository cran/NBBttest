\name{betaparametw}
\alias{betaparametw}
\title{
Estimation of proportion weights
}
\description{
Function betaparametw is used to calculate weights.
}
\usage{
betaparametw(xn, a, b)
}
\arguments{
  \item{xn}{
a vector of \code{m} library sizes. Library size 
is sum of counts over the whole library. 
}
  \item{a}{
parameter alpha(\eqn{\alpha}) in beta distribution 
derived from output of function betaparametab.
}
  \item{b}{
parameter beta (\eqn{\beta})in beta distribution 
derived from output of function betaparametab.
}
}
\details{
alpha and beta (\eqn{\alpha},\eqn{\beta}) are 
used to calculate weight. Then weight is in turn 
used to correct bias of estimation of alpha and 
beta in betaparametab function. 
}
\value{
return weight(W).
}
\references{
Baggerly KA, Deng L, Morris JS, Aldaz CM (2003) 
Differential expression in SAGE: accounting for 
normal between-library variation. 
\emph{Bioinformatics}, \bold{19}: 1477-1483.
\cr
Yuan-De Tan, Anita M. Chandler, Arindam Chaudhury, 
and Joel R. Neilson(2015) A Powerful Statistical 
Approach for Large-scale Differential Transcription Analysis. 
\emph{Plos One}. 2015 DOI: 10.1371/journal.pone.0123658. 
}
\author{
Yuan-De Tan
\email{tanyuande@gmail.com}
}

\seealso{
 \code{\link{betaparametab}} and \code{\link{betaparametVP}}.
}
\examples{
XX<-c(2000,2000,2000)
a<-1.1458
b<-6.4932
betaparametw(xn=XX,a=a,b=b)
}
\keyword{ weight}
