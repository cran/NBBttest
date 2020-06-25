\name{pathwayHeatmap}
\alias{pathwayHeatmap}

\title{
Heatmap for pathways found by gene ontology analysis
}
\description{
This function is used to show differential expressions 
of pathways or functions between conditions. These 
pathways or functions were detected by function 
annotation or gene ontology methods such as David function 
analysis tools or Ingenuity pathway analysis. Pathway score 
or pathway value is a weighted expression value across 
genes in a pathway or a function. The weigths of genes are 
given by p-values of enrichment or hit in function analysis. 
}

\usage{
pathwayHeatmap(dat, pathway, nci, r1, r2, colclass, rowclass, 
colrs, maptitle)
}

\arguments{
  \item{dat}{
count dataset that contains a column for gene name 
and count data columns in two conditions.
}
  \item{pathway}{
a list that lists a function colomn and at least two 
gene columns. Row is function name.
}
  \item{nci}{
number of columns for gene information such as gene name, 
strand, chromosome, etc in dataset dat.
}
  \item{r1}{
number of replicates in condition 1 in dataset dat.
}
  \item{r2}{
number of replicates in condition 2 in dataset dat.
}
  \item{colclass}{
colonm class, equivalent to replicates in two conditions, 
such as c(1,1,1,1,2,2,2,2),meaning that condition A has 4 
replicates and condition B has 4 replicates. Default =NULL.
}
  \item{rowclass}{
row class, pathway class, or pathway module. Default = NULL.
}
  \item{colrs}{
heatmap colors. User has 8 options: "redgreen", 
"greenred", "redblue", "bluered", "cm.colors",
"terrain.colors","topo.colors", and"heat.colors".  
}
  \item{maptitle}{
heatmap title, default = "".
}
}
\details{
This function uses sigificance to choose gene function/
pathwaydata and then to normalize the selected data by 
using z-scale. This function has multiple options to 
select map color, distance, cluster and x- and y-lab 
angles. 
} 
\value{
return a graph for heatmap.	
}

\author{
Yuan-De Tan
\email{tanyuande@gmail.com}
}
\note{
requres gplots and grDevices.
}

\seealso{
    \code{\link[gplots]{heatmap.2}},\code{\link{myheatmap}},
    \code{\link{myheatmap2}}
}

\examples{
data(upGAm)
data(pathwy.A.up)
pathwayup<-pathwy.A.up
colclass=c("1","1","1","1","2","2","2","2","2","2")
oldpar <- par(no.readonly =TRUE)
on.exit(par(oldpar))
par(mar=c(7.5,5.5,3.5,3))
par(oma=c(3,1,1,10))
pathwayHeatmap(dat=upGAm,pathway=pathwayup,nci=1,
r1=4,r2=6,colclass=colclass,rowclass=NULL,
colrs="greenred",maptitle="pathway up-expression 
in Group A")
}
\keyword{ heatmap }
\keyword{ pathway }
