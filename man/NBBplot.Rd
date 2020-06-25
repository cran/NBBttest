\name{NBBplot}
\alias{NBBplot}
\title{
Plot differential expression of exons within a 
specified gene using result outputed by NBBttest.
}
\description{
After performing NBBttest, NBBplot can be used to 
show differential expression of exons within a 
specified gene in na and nb replicates between
conditions A and B.
}
\usage{
NBBplot(res, gene, nci, na, nb, C1, C2)
}

\arguments{
  \item{res}{
object of NBBttest containing information of genes 
including gene name,strand, chromosome,exons, introns,
data, and t-value, p-value, significance/selection etc.
}
  \item{gene}{
gene name or symbol specified by user.
}
  \item{nci}{
number of columns for gene information.
}
  \item{na}{
replicate number in condition A
}
  \item{nb}{
replicate number in condition B
}
  \item{C1}{
name for condition 1(A)
}
  \item{C2}{
name for condition 2(B)
}
}
\details{
NBBplot consists of two parts: top is expression 
value of each exon in each replicate in two 
conditions marked in red and blue and bottom is 
boxes for exon and solid lines for introns. 
Differential expression of an exon is marked in red.   
}
\value{
output NBBplot figure for given gene.
}

\author{
Yuan-De Tan
\email{tanyuande@gmail.com}
}

\seealso{
\code{\link[DEXSeq]{plotDEXSeq}}
}
\examples{
data(exondata)

NBBplot(res=exondata, gene="H2-DMb1", nci=9, na=3, nb=3, C1="WT", C2="KO")

}
\keyword{ differentil expression }
\keyword{ exons }
