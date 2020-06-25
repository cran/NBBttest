\name{myheatmap}
\alias{myheatmap}
\title{
Heatmap 
}
\description{
This function is used to display heatmap of differential 
expressions of genes or isoforms  detected by NBBttest 
in the real count data.
}
\usage{
myheatmap(dat, IDcol, nci, r, r1, r2, 
colrs="greenred", rwcex=2.0,clcex=2.0,x=10, tree="both",
method="euclidean", ky=1.5, rowBarColor=NULL,  
colBarColor=NULL,  labrow="yes", labcol="yes", rsort="yes",
adjrow=c(0.3, 0.0 ), adjcol = c(1, 1) , rwangle=30,
clangle=30,  maptitle)
}

\arguments{
\item{dat}{
is NBBttest object that is outputted by mbetattest, 
includes information columns, data columns, t-value 
column, pvalue column,selection column.
}
\item{IDcol}{
nonnegative int value, indicating which column is 
ID column where IDs will be shown in heatmap rows,
is required.
}
\item{nci}{
nonnegative int value,indicating column number of 
informaction (gene or isoform annotation) columns 
in the data, is required. 	
}
\item{r}{
nonnegative int value, indicating column number of count 
data including empty columns if there are multiple 
datasets with the same row names.
}
\item{r1}{
nonnegative int value, indicating replicate number 
in condition 1.
}
\item{r2}{
nonnegative int value, indicating replicate number 
in cocndition 2.
}
 \item{colrs}{
heatmap colors. User has 8 options: "redgreen", 
"greenred", "redblue", "bluered", "cm.colors",
"terrain.colors","topo.colors", and"heat.colors". 
}
\item{rwcex}{
numeric argument: nonnegative value, used as 
cex.axis in for the row axis labeling.Default=1.8.
}
\item{clcex}{
numeric argument: nonnegative value, used as 
cex.axis for the column axis labeling. Default=1.8.	
}
\item{x}{
numeric argument: nonnegative value, used as 
argument of cm.colors(x), terrain.colors(x) and 
topo.colors(x), the default value is 10.	
}
\item{tree}{
tree(s) drawn on row or colum or both. User has four 
options: "both" for drawing trees on both row and 
column,"row" for drawing tree only on row,"column" 
for drawing tree only on column, and "none" for no 
tree specified on rows and columns. Default is "both".
}

\item{method}{
method to be chosen to calculate distance between 
columns or rows. It has four options: "euclidean", 
"pearson","spearman" and "kendall". The latter three 
are d=1-cc where cc is correlation coefficients. 
Default is "euclidean". 	
}
\item{ky}{
numeric argument: nonnegative value is used to 
determine key size. The defalut =1.5	
}
\item{rowBarColor}{(optional) character vector for 
RowSideColors and colRow. Length of rowBarColor 
equals to gene or isoform number or row number.  
rowBarColor contains the color names of classes or 
typpes or groups of row names and show row color 
side bars and color row names. If row names, for 
example, genes in row are not classified, then we 
suggest user use its default value: rowBarColor 
is NULL.
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
string, logistic value. Rows on heatmap are labeled 
with genes or targets if labrow="yes", otherwise, the 
rows are not labeled. 
}
\item{labcol}{
string, logistical value. columns on heatmap are labeled 
with samples or treatments if labcol="yes", otherwise, 
the columns are not labeled. 
}
\item{rsort}{
logistical value. If choose tree="both" or "row", then 
rsort does not work. However, if tree ="none" or "columnn", 
then rsort ="yes" will force rows to be sorted in descent 
way. The default is "yes"	
}

\item{adjrow}{
two numeric values. The first value is used to adjust left 
or right position (x-axis) of labels, and  the second value 
is used to adjust up or down position (y-axis) of labels. 
The fadault values are c(0.4, 1). 	
}
\item{adjcol}{
two numeric values. The first value is used to adjust left 
or right position (x-axis) of labels, and  the second value 
is used to adjust up or down position (y-axis) of labels. 
The fadault values are c(1, 1 ). 	
}
 \item{rwangle}{
angle of xlab under heatmap. The Default value is 30. 
}
 \item{clangle}{
angle of ylab. The default value is 30
}
\item{maptitle}{
string argument for giving heatmap title, default is 
set to be empty string ' '.
}
}

\details{
This function uses sigificance to choose data and then 
to normalize the selected data by using z-scale. This 
function has multiple options to select map color, 
distance, cluster and x- and y-lab angles.  
If r>(r1+r2), then data dat are multiple datasets. 
}        
                                      
\value{
no return value but create a heatmap.
}

\author{
Yuan-De Tan
\email{tanyuande@gmail.com}
}
\note{
requres gplots and grDevices.
}

\seealso{
       \code{\link[gplots]{heatmap.2}},
       \code{\link{myheatmap2}}
}
\examples{
 
data(result)
colclass=c("1","1","1","2","2","2")
oldpar <- par(no.readonly =TRUE)
on.exit(par(oldpar))
par(mar=c(7.5,5.5,3.5,1.2))
par(oma=c(3,1,1,3))
myheatmap (dat=result, IDcol=1, nci=7, r=6, r1=3,
r2=3, colrs="terrain.colors", rowBarColor=NULL,  
colBarColor=colclass, labrow="no", labcol="yes", 
rsort="yes", adjrow=c(0.3, 0.0 ), 
adjcol = c(1, 1) , maptitle="My heatmap")

}
\keyword{heatmap}
