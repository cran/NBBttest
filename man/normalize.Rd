\name{normalized}
\alias{normalized}
\title{
Normalization of data
}
\description{
Function normalize makes all libraries in dataset 
have the same library size.
}
\usage{
normalized(dat, nci, m=0, lg2="no")
}

\arguments{
  \item{dat}{
count data of RNA reads.
}
  \item{nci}{
number of columns for the information of genes or 
isoforms in dataset.
}
  \item{m}{
numeric value for choosing genes or isoforms. If user 
wants to discard genes or isoforms with mean < 5, 
then m = 5. The default value is 0.
}
  \item{lg2}{
logistic value. lg2="yes"indicates that data are 
transformed in logarithm of 2.
}
}
\details{
Due to difference in RNA abstraction between libraries 
or cell samples or tissues, PCR amounts of RNA libraries 
would have difference that is not due to biological 
effects. To correctly compare differential expressisons 
of genes between conditions or samples, one must should 
give the same RNA abstraction in all given samples. 
This is impossible. To address this problem, only one 
way is to normalize these count data across all given 
samples so that all experimental samples (libraries) 
have the same total counts. 
}
\value{
output a standard datasheet.
}
\references{
Yuan-De Tan Anita M. Chandler, Arindam Chaudhury, and 
Joel R. Neilson(2015) A Powerful Statistical Approach 
for Large-scale Differential Transcription Analysis. 
\emph{Plos One}, 10.1371/journal.pone.0123658.  
\cr
Anders S, Huber W (2010) Differential expression 
analysis for sequence count data. \emph{Genome Biol}, 
11: R106.
}
\author{
Yuan-De Tan
\email{tanyuande@gmail.com}
}

\examples{
data(jkttcell)	
njkttcell<-normalized(dat=jkttcell[1:50,],nci=7)
}

\keyword{ normalize }
\keyword{ library size }
