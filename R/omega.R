omega<-
function(XX,nci,r1,r2,sn,alpha=0.05){	
rho<-function(dat,nci,na,nb,alpha){
simdat<-simulat(yy=dat,nci=nci,r1=na,r2=nb,q=0)
simdata<-simdat[,(nci+2):ncol(simdat)]
#print(head(simdata))
#simdata<-apply(simdata,2,as.numeric)
simdat<-as.data.frame(simdat)
#print(c(length(simdat$geneid), nrow(simdata)))
geneid<-simdat$geneid
#simid<-geneid

simdatb<-cbind(geneid, geneid,simdata)
simdatb<-as.data.frame(simdatb)
colnames(simdatb)[seq(2)]<-c("geneid","simid")

rho<-smbetattest(X=simdatb, na=r1,nb=r2)
 #symb<-sim$symb
# rho<-subset(sim$rho,symb=="+")
K<-length(rho)
if(K==0){
rhov<-max(rho)
#rhov<-max(sim$rho)
}else if(K<=5){
	
rhov<-min(rho)

}else if(K>5){

rho.sort<-sort(rho, decreasing = FALSE)

i<-seq(K)/K
rho.i<-cbind(rho.sort, i)
rho.85<-(subset(rho.i, i>=0.85))
rhov<-mean(rho)
	}
#print(rhov)
return(rhov)
}

rhov<-rep(NA, sn)

for(i in seq(sn)){
rhov[i]<-rho(dat=XX,  nci=nci, na=r1, nb=r2, alpha=alpha)
}

w<-round(mean(na.omit(rhov)), 2)
return(w)

}
