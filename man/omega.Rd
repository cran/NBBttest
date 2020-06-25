\name{omega}
\alias{omega}
\title{
Omega calcularion
}
\description{
Omega function is a function that 
is used to estimate omega using simulate  
null data from negative bionomial distribution. 
is a null rho that is used as 
a threshold for real rho. Simulation is 
dependent on the original data. 
}
\usage{
omega(XX, nci, r1, r2, sn, alpha = 0.05)
}
\arguments{
  \item{XX}{
the real dataset.
}
  \item{nci}{
number of columns for information of data, like gene id, 
isoform id, gene name, etc.
}
  \item{r1}{
size of sample 1 or number of replicates in condition 1.
}
  \item{r2}{
size of sample 2 or number of replicates in condition 2.
}
  \item{sn}{
number of simulations specified.
}
  \item{alpha}{
significance level of test. Default value is 0.05.
}
}
\details{
This function is to use  null data to calculate 
omega value with rho  = 1. 
}
\value{
return a numeric value.
}
\references{
Yuan-De Tan Anita M. Chandler, Arindam Chaudhury, and 
Joel R. Neilson(2015)  A Powerful Statistical Approach 
for Large-scale Differential Transcription Analysis. 
\emph{Plos One}. 2015 DOI: 10.1371/journal.pone.0123658. 
}
\author{
Yuan-De Tan
\email{tanyuande@gmail.com}
}

\seealso{
 \code{\link{pratio}} and \code{\link{oddratio}}

}
\examples{
data(jkttcell)	

w<-omega(XX=jkttcell[1:100,],nci=7,r1=3,r2=3,sn=2,alpha=0.05)	
}
\keyword{ rho }
\keyword{omega}
