\name{betaparametab}
\alias{betaparametab}
\title{
 Estimation of parameters alpha (\eqn{\alpha}) 
 and beta (\eqn{\beta}) of beta distribution
}
\description{
Parameters alpha(a) and beta (b) in betat distribution 
are estimated by using an iteration algorithm. 
}
\usage{
betaparametab(xn, w, P, V)
}
\arguments{
  \item{xn}{
column vector, a set of library sizes.
}
  \item{w}{
column vector, a set of weights.
}
  \item{P}{
proportion of counts of a gene or an isoform.
}
  \item{V}{
variance of proportions of counts of a gene or an isoform over m 
replicate libraries in a condition.
}
}

\value{
 return parameters a and b.
}
\references{
Baggerly KA, Deng L, Morris JS, Aldaz CM (2003) 
Differential expression in SAGE:  accounting for 
normal between-library variation. 
\emph{Bioinformatics} \bold{19}: 1477-1483.
}
\author{
Yuan-De Tan
\email{tanyuande@gmail.com}
}

\seealso{
 \code{\link{betaparametVP}}, \code{\link{betaparametw}}	
}
\examples{
XX<-c(2000,2000,2000)
p<-0.15
V=0.004
w<-c(0.3,0.3,0.3)
betaparametab(xn=XX,w=w,P=p,V=V)
}
\keyword{alpha}
\keyword{beta}
