\name{smbetattest}
\alias{smbetattest}

\title{
Performance of multiple beta t-test on simulated data
}
\description{
This function is to peform beta t-test with \eqn{\rho}
=1 and \eqn{\omega}=1 on simulated data. The result 
lists differentially expressed genes or isoforms and 
their \eqn{\rho} values. The \eqn{\rho} values are used 
to calculate \eqn{\omega} value for performance of beta 
t-tests on the real data.
}
\usage{
smbetattest(X, na, nb, alpha)
}

\arguments{
  \item{X}{
simulated count data with N genes or isoforms.

}
  \item{na}{
number of replicate libraries in condition A.
}
  \item{nb}{
number of replicate libraries in condition B.

}
  \item{alpha}{
statistical probabilistic threshold, default is 0.05.

}

}
\details{
Before performing NBBttest on real data, user needs 
\eqn{\omega} value for the threshold of \eqn{\rho}. 
To determine \eqn{\omega} value, user is requred to 
generate a set of null data having the same gene or 
isoform number and the same numbers of replicate 
libraries in two conditions and then performs beta 
t-test on the null datasets by setting \eqn{\rho} =1 
and \eqn{\omega} =1.  In current package, NBBttest 
can automatically perform the simulation of null data,
multiple beta t-test to estimate \eqn{\omega}. 
}
\value{
Return a set of null \eqn{\rho} values.

}
\references{
Yuan-De Tan Anita M. Chandler, Arindam Chaudhury, and 
Joel R. Neilson(2015) A Powerful Statistical Approach 
for Large-scale Differential Transcription Analysis. 
\emph{Plos One}. DOI: 10.1371/journal.pone.0123658. 
}
\author{
Yuan-De Tan
\email{tanyuande@gmail.com}
}

\seealso{
See Also as \code{\link{mbetattest}}
}
\examples{
data(skjt) 
nrho<-smbetattest(X=skjt[1:60,],na=3,nb=3,alpha=0.05)

}
\keyword{ t-test }
\keyword{ simulation}
