pathwayHeatmap<-function(dat,pathway,nci,r1,r2,colclass,rowclass,colrs,maptitle){

dat<-as.data.frame(dat)	
#pathway<-as.data.frame(pathway)	
r=r1+r2
pathdata<-function(dat,geneset,nci,r){
gene<-dat$gene
gene<-toupper(gene)
#print(gene[1:20])
geneset<-na.omit(geneset)
geneset<-t(geneset)
ind<-is.element(gene,geneset)
sbdat<-dat[ind,]	
#print(sbdat[1:5,])
pvalue<-sbdat$pvalue
np<-length(pvalue)
for(i in seq(np)){
	if(pvalue[i]==0){
		pvalue[i]<-1e-08
	}
}
weight<-(-log10(pvalue)/sum(-log10(pvalue)))
sbdat1<-sbdat[,c((nci+1):(nci+r))]*weight
sbdat2<-apply(sbdat1,2,sum)
return(sbdat2)
}

nr=nrow(pathway)

pathwaydat<-matrix(NA,nrow=nr,ncol=r)

pathwayname<-pathway[,1]
colnames(pathwaydat)<-colnames(dat)[c((nci+1):(nci+r))]

for(i in seq(nr)){
pathwaydat[i,]<-pathdata(dat=dat,geneset=pathway[i,],nci=nci,r=r)	
   	
}

pathwaydat<-cbind(as.data.frame(pathwayname),as.data.frame(pathwaydat))
#print(pathwaydat[1:5,])

#par(mar=c(5,4,2,20),oma=c(4,1,2,19))

myheatmap(dat=pathwaydat, IDcol=1, nci=1, r=r, r1=r1, r2=r2,
colrs=colrs, rwcex=1.5,clcex=2.0,x=10, tree="both",
method="euclidean", ky=1.5, rowBarColor=rowclass,
colBarColor=colclass,  labrow="yes", labcol="yes", 
rsort="yes", adjrow=c(0.0, 0.0 ), adjcol = c(1, 1) ,
rwangle=0, clangle=30,  maptitle=maptitle)
	
}
