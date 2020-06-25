NBBplot<-function(res,gene,nci,na,nb,C1,C2){
    #res is result file from differential analysis
    #gene is specified for NBBplot
    #nci is number of columns of information in result files such as gene,
    # exon, start,end, strand etc
    #na,nb are sample sizes or replicate numbers in conditions a and b
    #C1 and C2 are string, names of conditions 1(a) and 2(b), such as "control"
    # "treatment" or "unknockout","knockout", etc.
clname<-colnames(res)
genename<-c("gene","Gene","geneid","Geneid","gene_name")

indx<-is.element(clname, genename)

jj<-which(indx==TRUE)
try(if(jj==0)
  stop("result file does not have gene column")
)
colnames(res)[jj]<-"Gene"

res<-as.data.frame(res)
Genename<-res$Gene

res.gene<-subset(res,res$Gene==gene)

exonname<-c("exon","Exon","isoform","Isoform","isof","Isof","featureID")
indx<-is.element(clname, exonname)
jj<-which(indx==TRUE)
try(if(jj==0)
    stop("result file does not have exon or isoform column")
)

colnames(res.gene)[jj]<-"exon"
starname<-c("start","Start","posit1","posit.start","pos.start")
indx<-is.element(clname,starname)
jj<-which(indx==TRUE)
try(if(jj==0) stop("result file does not have start column"))
colnames(res.gene)[jj]<-"start"
endname<-c("end","End","posit2","posit.end","pos.end")
indx<-is.element(clname,endname)
jj<-which(indx==TRUE)
try(if(jj==0) stop("result file does not end column"))
colnames(res.gene)[jj]<-"end"
res.gene.sort<-res.gene[order(res.gene$start),]
res.gene.sort<-res.gene.sort[order(res.gene.sort$end),]
dat<-res.gene.sort[,c(c(nci+1):c(nci+na+nb))]
dat<-apply(dat,2,as.numeric)

ddat<-matrix(NA,nrow=2*nrow(dat),ncol=ncol(dat))
for(i in seq(nrow(dat))){
	i1<-i*2-1
	i2=i*2
	for(j in seq(ncol(dat))){
	ddat[i1,j]<-dat[i,j]	
	ddat[i2,j]<-dat[i,j]
	}
}

en<-nrow(ddat)
x<-seq(en)
start<-res.gene.sort$start
end<-res.gene.sort$end

tddat<-t(ddat)
#print(tddat)

exon<-res.gene.sort$exon
lty<-rep(2,en)
for(i in seq(nrow(dat))){
	j=i*2-1
#	j2<-i*2
	lty[j]<-1
}

oldpar <- par(no.readonly =TRUE)
on.exit(par(oldpar))
par(mar=c(7,7,2,1),oma=c(3,2,1,1))
par(mfrow=c(2,1))

#op <- par(bg = "thistle")
plot(x,tddat[1,],type="l",lty=2,lwd=1, ylim=c(0,150+max(max(dat))), col="blue", xaxt="n",bty="n",xlab="",ylab="expression",cex.lab=2.8,cex.axis=1.8, main=gene,cex.main=3.0)
for(i in 1:na){
	lines(x,tddat[i,],type="l",lty=2,lwd=1,col="blue")
}
for(i in c(na+1):c(na+nb)){
	lines(x,tddat[i,],type="l",lty=2,lwd=1,col="red")
}

for(i in 1:na){
	for(j in 1:nrow(dat)){
		j1<-j*2-1
		j2<-j*2
	lines(x[j1:j2],tddat[i,j1:j2],type="l",lwd=3,col="blue")
}
}
for(i in c(na+1):c(na+nb)){
	for(j in 1:nrow(dat)){
		j1<-j*2-1
		j2<-j*2
	lines(x[j1:j2],tddat[i,j1:j2],type="l",lwd=3,col="red")
}
}

xx<-seq(nrow(ddat))
posit<-rep(NA,nrow(ddat))
hposit<-rep(20,nrow(ddat))
Exon<-rep(" ", nrow(ddat))
for(i in seq(nrow(dat))){
	j1=i*2-1
	j2<-i*2
	posit[j1]<-start[i]
	posit[j2]<-end[i]
	Exon[j1]<-as.character(exon[i])
	Exon[j2]<-as.character(exon[i])
}
#axis(1,at=xx,labels=Start,las=2)
#axis(1,at=xx,labels=End,las=2)
if(length(Exon)/2<10){
axis(1,at=xx,labels=Exon,las=0,cex.axis=2)
}else if(length(Exon)/2>=10&length(Exon)/2<20){
axis(1,at=xx,labels=Exon,las=0,cex.axis=1.5)
}else{
axis(1,at=xx,labels=Exon,las=0,cex.axis=1.0)
}

nn=length(xx)
if(nn==4){
    k1=1
    k2=4
}else if(nn==6){
    k1=2
    k2=5
}else if(nn==8){
  k1=(nn/2)-(nn/4)+1
  k2=(nn/2)+(nn/4)-1
}else if(nn>=10&nn<=14){
    k1=(nn/2)-(nn/4)+2
    k2=(nn/2)+(nn/4)-2
}else{
    k1=(nn/2)-(nn/4)+4
    k2=(nn/2)+(nn/4)-4
}
mtext(C1,side=1, line=4,col="blue",at=k1,cex=3.0)
mtext(C2,side=1, line=4,col="red",at=k2,cex=3.0)
#legend(x1,y,lty=1,lwd=3,c(C1,C2),col=c("blue","red"),cex=2)
#legend(x1,y2,lty=1,lwd=3,C2,col="red",cex=2)

lab<-res.gene.sort[,ncol(res.gene.sort)]

boxcolor<-rep("darkseagreen3",length(lab))
for(i in 1:length(lab)){
	if(lab[i]==1){
	boxcolor[i]<-"red"	
	}
}

plot(posit,hposit/2,bty="n", type="l", lwd=6, yaxt="n",xaxt="n",ylim=c(0,30),ylab="",xlab="",cex.axis=3)

axis(1,at=posit,labels=posit,las=2)
Start<-start
End<-end

for(i in 1:c(length(start)-1)){
    if(end[i]-start[i+1]>100){
    Start[i+1]<-start[i+1]-20
}
if(end[i]-start[i+1]<10){
Start[i+1]<-start[i+1]+20
}
if(end[i+1]-start[i+1]<20){
End[i+1]<-end[i+1]+20
}

}

for(i in 1:length(start)){
	rect(Start[i],0,End[i],hposit,lwd=0.5,col=boxcolor[i])

}
for(i in 1:length(start)){
text(c(start[i]+end[i])/2,c(2+hposit),exon[i])
}
mtext("exon map(bp)",side=1, line=6,cex=3.0)

}
