\name{myheatmap2}
\alias{myheatmap2}

\title{
Heatmap2
}
\description{
This function is used to display heatmap of 
differential expressions of genes or isoforms  
detected by NBBttest in the real count data.
 }
\usage{
myheatmap2(dat, IDcol=1, nci=NULL, r, colrs=
"greenred", rwcex=2.8, clcex=2.8, x=10, tree=
"both", method="euclidean", ky=1.5, rowBarColor
=NULL,  colBarColor=NULL, labrow="yes", 
labcol="yes", rsort="yes", adjrow=c(0.2, 0.0 ), 
adjcol = c(1, 1) , rwangle=30, clangle=30, 
maptitle="",keyvalue)
}
\arguments{
\item{dat}{
object of NBBttest, including information 
columns, data columns, tvalue column, pvalue and 
selection column.
}
\item{nci}{
int value, number of columns for information 
of genes or isoforms, like "geneid", ""tagetid"", 
"straind", annotation etc.	
}
\item{IDcol}{
int value, indicating which column is ID column 
where IDs will be shown in heatmap rows.
}
  \item{r}{
int value, indicating column number of count data 
including empty columns if there are multiple datasets 
with the same row names.
}

\item{colrs}{
heatmap colors. User has 8 options: "redgreen", 
"greenred", "redblue", "bluered", "cm.colors",
"terrain.colors","topo.colors", and"heat.colors". 
Default colrs is "redgreen".   
}
\item{rwcex}{
numeric argument: nonnegative number, used as cex.axis 
for the row axis labeling. Default value is 1.8.
}
\item{clcex}{
numeric argument: nonnegative number, used as cex.axis 
for the column axis labeling. Default value is 1.8.
}
\item{x}{
numeric argument: nonnegative number, used as argument 
of cm.colors(x), terrain.colors(x) and topo.colors(x), 
Default value is 10.	
}
\item{tree}{
tree(s) drawn on row or column or both. User has four 
options: "both" for drawing trees in both row and 
column,"row" for drawing tree only in row,"column" 
for drawing tree in only column, and "none" for no 
tree specified on rows and columns. If tree =both, 
then colomns and rows are sorted by trees. If tree = 
"row", the colmns are not sorted.  If tree ="column", 
then rows are not sorted. Default is "both".
}
\item{method}{
method to be chosen to calculate distance between 
columns or rows. It has four options: "euclidean", 
"pearson","spearman" and "kendall". The latter three 
are d=1-cc where  cc is correlation coefficients. 
Default is "euclidean". 	
}
\item{ky}{
numeric argument: nonnegative number, is used to 
determine key size. The defalut =1.5.	
}
\item{rowBarColor}{
(optional) character vector for RowSideColors and
colRow. Length of rowBarColor equals to gene or 
isoform number or row number.  rowBarColor contains 
the color names of classes or typpes or groups of 
row names and show row color side bars and color 
row names. If row names, for example, genes in row 
are not classified, then we suggest user use its 
default value: rowBarColor=NULL.
}
\item{colBarColor}{
(optional) character vector for colSideColors and 
colCol. Length of colBarColor equals to sample number.  
colBarColor contains the color names of classes or 
typpes or groups of  coloumn names and show column 
color side bars and color column names. If column 
names, for example, samples in column are just two 
types, then we suggest user use its default value: 
colBarColor=NULL.	
}
\item{labrow}{
a string, logistic value. Rows on heatmap are labeled 
with genes or targets if labrow="yes", otherwise, 
the rows are not labeled.
}
\item{labcol}{
a string, logistical value. Columns on heatmap are labeled 
with samples or treatments if labcol="yes", otherwise, the 
columns are not labeled. 	
}

\item{adjrow}{
two numeric values. The first value used to adjust left or 
right position (x-axis) of labels and  the second value is 
used to adjust up or down position (y-axis) of labels. The 
fadault values are c(0.3, 0.0). 	
}
\item{adjcol}{
two numeric values. The first value used to adjust left or 
right position (x-axis) of labels, and  the second value is 
used to adjust up or down position (y-axis) of labels. 
The fadault values are c(1, 1 ). 	
}
\item{rsort}{
logistical value. If choose tree="both" or "row", then 
rsort does not work. However, if tree ="none" or "columnn", 
then rsort ="yes" will force rows to be sorted in descent 
way. The default is "yes"	
}
\item{rwangle}{
heatmap object: angle of xlab. Thed default value is 30. 
}
\item{clangle}{
heatmap object: angle of ylab. The default value is 30. 
}
\item{maptitle}{
string argument for giving heatmap title. The default 
value is set to be empty string ' '.
}
\item{keyvalue}{
sting argument for name of key x-axis.	
}
}
\details{
This function uses selection to choose genes or isoforms 
in the data and then to normalize the selected data by 
using n-scale. Different from z-score, n-score does not 
follow standard normal distribution with mean = 0 and 
variance =1 for all rows but it has the same largest 
count in all rows and shows multiple colors for numeric 
difference between two conditions. This function has 
multiple options to select map color, distance, cluster 
and x-lab and y-lab angles. This function can be able to 
display multiple datasets in two ways: if multiple 
datasets have the same row names or features, the these 
datasets are put onto the different columns separated 
with empty column named with dataset names.If multiple 
datasets have the same column names of the datasets, 
then put them on different rows saparated with empty 
rows named with dataset names or whatever names user 
specifies.  
}
\value{
not return value but create a heatmap.
}

\author{
\email{tanyuande@gmail.com} 
}
\note{
requires gplots and grDevices. If the data for heatmap 
are multiple datasets, then tree="none" and sort="no", 
otherwise, myheatmap2 will get error. So before 
performing myheatmap2, user should sort the data in excel. 
}

\seealso{
\code{\link[gplots]{heatmap.2}}, \code{\link{grDevices}}, 
and \code{\link{myheatmap}}.
}
\examples{

data(result)
colclass=c("1","1","1","2","2","2")
oldpar <- par(no.readonly =TRUE)
on.exit(par(oldpar))
par(mar=c(7.5,5.5,3.5,1.2))
par(oma=c(3,1,1,3))
myheatmap2(dat=result, IDcol=1, nci=7, r=6, 
colrs="greenred", rwcex=1.8, clcex=1.8, x=10, 
tree="both", method="euclidean", ky=1.5, 
rowBarColor=NULL,  colBarColor=colclass, 
labrow="no", labcol="yes", adjrow=c(0.2, 0.0 ), 
adjcol = c(1, 1) , rwangle=0, clangle=30,
maptitle="My heatmap2",keyvalue="count")

}

\keyword{heatmap}
\keyword{differential expressions }
