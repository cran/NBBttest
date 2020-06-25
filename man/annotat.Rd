\name{annotat}
\alias{annotat}

\title{
Annotation of genes within which alternative 
splicing occurs
}
\description{
Alternative splicing is detected in any element 
of 3'UTR, 5'UTR, exons and introns within a gene 
using RNA-seq data where RNA reads are mapped to 
a reference genome. As an example for annotation, 
the RNA-seq reads derived from human samples can 
be mapped onto human genome reference (GRCh38) 
using different methods, for example, HTSeq, spladder, 
rMAT, cufflinks, etc. These methods can detect 
alternative splicing sites within genes. However, 
none of these methods does gene annotation for users. 
Our NBBttest offers a R function for annotating genes 
with exons or isoforms.      
}
\usage{
annotat(infile, mfile, type="gene", colunmset = "NLL")
}

\arguments{
  \item{infile}{
input data file with ENSGid column for annotation
}
  \item{mfile}{
reference genome file with ENSG id column.
}
  \item{type}{
has three options: "gene", "isoform" or ""isof" and 
"exon", see details. The default is "gene".
}
  \item{colunmset}{
EGNSid column set. If the RNA count data are made by 
using HTSeq and DEXSeq annotation file, then some of 
genes have many different ENSGids. For example, ENSG00000285476+ENSG00000182230+ENSG00000251623 has 
three ENSG ids ENSG00000285476,ENSG00000182230, and 
ENSG00000251623 that share one gene, so it is splited 
into three columns (2,3,4)in excel. The default is 
2(column 2). 
 }
}
\details{
If type = "gene", then count data of RNA reads are 
obtained at gene level, annotation would be executed 
at gene level or if type = "isof" or "isoform", then 
RNA reads were mapped onto elements (for example, 3'UTR,
5'UTR, exon or cassette, intron) within genes and  
annotation would be executed at isoform level or if 
type = "exon", then RNA count data were obtained by 
mapping RNA reads onto exome by DEXseq and annotation 
would be done at exon level defined by DEXSeq. Note that 
GRCh38 is too big so it was removed from data. User may 
request to get it from yuande/github. 
}
\value{
return original data with an additional column for gene. 
}
\author{
Yuan-De Tan
\email{tanyuande@gmail.com}
}

\examples{
data(DDX39_100)
data(gtfa)
DDX39_30<-annotat(infile=DDX39_100[1:30,],mfile=gtfa,type="gene")
}
\keyword{ annotation }
\keyword{ gene }
