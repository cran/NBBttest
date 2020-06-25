\name{QC}
\alias{QC}
\title{
Count data quality check }
\description{
QC is used to check quality of count data in 
two ways: scatter plot of two replicated samples 
and correlation heatmap of all samples.
}
\usage{
QC(dat, nci, S1="NULL", S2="NULL", method=
"plot", colrs="greenred", rwcex=1.8, clcex=1.8, 
x=10, tree="none", log="none", col="blue", 
pch=19, labsize=1.5, axis=1.5)
}


\arguments{
  \item{dat}{
count matric dataset.
}
  \item{nci}{
number of columns containing data information suach 
as gene id, library id, target id, gene name, 
strand etc.
}
  \item{S1}{
numeric int value, indicating which column in data 
matrix is specified in x-axis. S1 > nci. 
}
  \item{S2}{
numeric int value, indicating which column in data 
matrix is specified in y-axis. S2> nci. 
}
  \item{method}{
string. Here two methods are given for choice: 
"plot" and "heatmap". Default is "plot". 
}
  \item{colrs}{
string. 8 color sets are given for choice in this 
function: "redgreen", "heat.colors", "redblue", 
"greenred", "bluered", "cm.colors", "terrain.colors", 
"topo.colors". The default color set is "redgreen".
}
  \item{rwcex}{
positive numbers, used as cex.axis for the row axis 
labeling. The default value is 1.8. 
}
  \item{clcex}{
positive numbers, used as cex.axis for the column 
axis labeling. The default value is 1.8. 
}
  \item{x}{
numeric argument: positive number, used as argument 
of cm.colors(x), terrain.colors(x) and topo.colors(x), 
the default value is 10.	
}
  \item{tree}{
tree(s) drwan on row or column or both. User has 
four options: "both" for drawing trees on both row 
and column,"row" for drawing treeonly on row,"column" 
for drawing tree on only column, and "none" for no 
tree specified. Default is "none".
}
  \item{log}{
string. Two options are given for choice: "none" and 
"log". log="log" indicates that data value would be 
transfomed with log2. Default is "none". 
}
  \item{col}{
string, used to specify scatter plot dot color. 
}
  \item{pch}{
numeric value, used to specify dot type.
}
\item{labsize}{
numeric value for size of xlabel and ylabel	
}
\item{axis}{
numeric value for axis scale.	
}
}
\details{
S1, S2 and nci must be given numeric int values for 
plot and heatmap. However, when method is chosen to 
be "heatmap", then S1 and S2  are not specified.  
Columns of information should be left of matrix and 
count data should be after columns of information. 
}

\value{
not return values but create scatter plot or heatmap plot.
}
%%}
\author{
Yuan-De Tan
\email{tanyuande@gmail.com} 
}
\note{
requires gplots and grDevices.
}

\seealso{
\code{\link[gplots]{heatmap.2}} and \code{\link{grDevices}}
}
\examples{
data(jkttcell)
QC(dat=jkttcell, nci=7, S1=8, S2=9, method = "plot", 
log = "log", col = "blue", pch = 19)
QC(dat=jkttcell, nci=7, S1=8, S2=9, method = "plot", 
log = "log", col = "blue", pch = 19)
QC(dat=jkttcell, nci=7, method = "heatmap", log = "log")
}
\keyword{QC }
\keyword{ scatter plot}
