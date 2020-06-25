library(NBBttest)
#require(gplots, warn.conflicts = FALSE)


#########################################################
# function QC for quality chek
########################################################
data(jkttcell)
QC(dat=jkttcell, nci=7, S1=8, S2=9, method = "plot", log = "log", col = "blue", pch = 19)
QC(dat=jkttcell, nci=7, S1=8, S2=9, method = "plot", log = "log", col = "blue", pch = 19)
QC(dat=jkttcell, nci=7, method = "heatmap", log = "log")
##############################################################
# perform mbetattest
##############################################################
data(jkttcell) 

res<-mbetattest(X=jkttcell[1:200,], nci=7, na=3, nb=3, alpha=0.05, norm="yes",
C=0,side="both", level="isoform")

#data(jkttcell)
#head(jkttcell)
#res.gene<-mbetattest(X=jkttcell[1:200, ], nci=7, na=3, nb=3, alpha=0.05, norm="yes",
#C=0,side="both", level="polyA.gene")

##############################################################
# perform mbetattest  on the CRISPR data at gene level
##############################################################
data(sgRNA)
res.gene<-mbetattest(X=sgRNA[1:200,], nci=3, na=4, nb=4, alpha=0.05, norm="no",
C=1.22,side="both", level="CRISPR.gene")
#######################################################
# perform heatmap for differential expression
#######################################################
# both trees in row and column
data(result)
# # "terrain.colors
myheatmap (dat=result, IDcol=1, nci=7, r=6, r1=3, r2=3,  colrs="terrain.colors",
rowBarColor=NULL,  colBarColor=c("1","1","1","2","2","2"), labrow="no", labcol="yes",
rsort="yes", adjrow=c(0.3, 0.0 ), adjcol = c(1, 1) , maptitle="My heatmap")
myheatmap2(dat=result, IDcol=1, nci=7, r=6,  colrs="terrain.colors",
rowBarColor=NULL,  colBarColor=c("1","1","1","2","2","2"), labrow="no", labcol="yes",
rsort="yes", adjrow=c(0.3, 0.0 ), adjcol = c(1, 1) , maptitle="My heatmap")
##"topo.colors"
myheatmap (dat=result, IDcol=1, nci=7, r=6, r1=3, r2=3,  colrs="topo.colors",
x=10, rowBarColor=NULL,  colBarColor=c("1","1","1","2","2","2"), labrow="no", labcol="yes",
rsort="yes", adjrow=c(0.3, 0.0 ), adjcol = c(1, 1) , maptitle="My heatmap")
myheatmap2(dat=result, IDcol=1, nci=7, r=6,  colrs="topo.colors", x=10,
rowBarColor=NULL,  colBarColor=c("1","1","1","2","2","2"), labrow="no", labcol="yes",
rsort="yes", adjrow=c(0.3, 0.0 ), adjcol = c(1, 1) , maptitle="My heatmap")
#none
#myheatmap (dat=result, IDcol=1, nci=7, r=6, r1=3, r2=3, colrs="topo.colors",
#x=10, tree="none",rowBarColor=NULL,  colBarColor=c("1","1","1","2","2","2"), labrow="no",
#labcol="yes", rsort="yes", adjrow=c(0.3, 0.0 ), adjcol = c(1, 1) , maptitle="My heatmap")
#myheatmap2(dat=result, IDcol=1, nci=7, r=6,  colrs="topo.colors", x=10, tree="none",
#rowBarColor=NULL,  colBarColor=c("1","1","1","2","2","2"), labrow="no", labcol="yes",
#rsort="yes", adjrow=c(0.3, 0.0 ), adjcol = c(1, 1) , maptitle="My heatmap")

