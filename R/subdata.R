subdata<-
     function(xx,sg){
m<-nrow(xx)
data<-xx[with(xx,order(Gene)),]	
#data<-xx[order(xx$Gene),]	
gm<-table(data$Gene)
gm<-gm[gm>0]
gn<-length(gm)

gg<-rep(0,m)
a<-0
for(i in seq(gn)){
	for(j in seq(gm[i])){
		a<-a+1
		gg[a]<-gm[i]
		}
	}

data<-cbind(as.data.frame(data),gg)	

if(sg==1){
x<-subset(data,data$gg==1)
}else{
x<-subset(data,data$gg>1)	
}
nc<-ncol(x)
x<-x[,-nc]
return(x)
	
	}

