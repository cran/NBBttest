\name{mtpvadjust}
\alias{mtpvadjust}

\title{
Adjust p-values for multiple comparisons
}
\description{
Given a set of p-values and chosen a C-value, 
returns a set of adjusted p-values
}
\usage{
mtpvadjust(pv, C)
}

\arguments{
  \item{pv}{
numeric vector of p-values (possibly with NAs). 
Any other R object is coerced by as.numeric.
}
  \item{C}{
real numeric value for specifying a multiple procedure. 
}
}
\details{
C=0 indicates that p-values are not adjusted,  C=1.22 
indicates that p-values are adjusted with Benjamini and 
Hochberg (1995) ("BH"). The adjusted p-values are 
called "fdr". When \eqn{C\ge 1000},p-values are adjusted 
with the Bonferroni method. \eqn{C<1.22} indicates that 
p-values are adjusted by a relaxed BH method while 
C>1.22, p-values are adjusted by a more strict BH method.
}
\value{
A numeric vector of corrected p-values (of the same 
length as p, with names copied from p)
}
\references{
Benjamini, Y., and Hochberg, Y. (1995). Controlling the 
false discovery rate: a practical and powerful approach to 
multiple testing. Journal of the Royal Statistical Society 
Series B, 57, 289 - 300. 
\cr
Yuan-De Tan, Hongyan Xu; A general method for 
accurate estimation of false discovery rates in identification 
of differentially expressed genes, Bioinformatics, 
Volume 30, Issue 14, 15 July 2014, Pages 2018 - 2025.
}
\author{
Yuan-De Tan
\cr
\email{tanyuande@gmaail.com}
}
\seealso{
\code{\link[stats]{p.adjust}}
}
\examples{
set.seed(123)
x <- rnorm(50, mean = c(rep(0, 25), rep(3, 25)))
p <- 2*pnorm(sort(-abs(x)))
p.1.22<-mtpvadjust(pv=p, C=1.22)
}
\keyword{ multiple test }
\keyword{ adjust pvalue}
