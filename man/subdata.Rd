\name{subdata}
\alias{subdata}
\title{
Splite data into two subsets
}
\description{
Data are splited into two subsets: gene 
single-isoform  data and gene multi-isoform data.
}
\usage{
subdata(xx, sg)
}
\arguments{
  \item{xx}{
real data containing single-isoforms and 
multi-isoforms of genes in rows.
}
  \item{sg}{
int value.
}
}
\details{
For the RNA count data, some genes have only 
one isoform, some genes have multiple isoforms.
so data are divided by subdata.R into two datasets: 
single-isoform data and multi-isoform data.
}
\value{
return dataset by setting sg=1 or sg=2.
}

\author{
Yuan-De Tan
\email{tanyuande@gmail.com}
}

\examples{
data(jkttcell)	
colnames(jkttcell)[3]<-"Gene"
jk.mg<-subdata(xx=jkttcell, sg=2)
}
\keyword{data}
\keyword{splite}
