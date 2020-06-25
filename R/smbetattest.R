smbetattest <-
function(X, na, nb, alpha=0.05){
 ###############################################	
 #	smbetattest is equivalent to mbetattest but used to simulated data
 # to calculate W values
 ################################################
 
	cn<-length(X[1,])
	rn<-length(X[,1])
	XC<-X[, seq(cn-na-nb)]
	Xa<-X[,(cn-na-nb+1):cn]
    XX<-na.omit(Xa)
	size<-max(apply(XX, 2, sum))
	pvalue<-rep(1,rn)
	betattest<-betattest(XX,na=na,nb=nb, level="RNA")
	prat<-pratio(xx=XX,na=na,nb=nb)
	odrat<-oddratio(XX=XX, na=na,nb=nb)
	t_value<-betattest[,1]
 
	rho<-(prat*odrat)^0.5
	beta_t<-t_value#*rho
	df<-betattest[,2]
#	for(i in 1:rn){
#		if(!(is.na(beta_t[i]))){
#			if(abs(beta_t[i])>0){
#			pvalue[i]<-2*(1-pt(abs(beta_t[i]),df=df[i]))
             pvalue <-2*(1-pt(abs(beta_t ),df=df))
#			}
#			}

#		}

#	XD<-cbind(XC, beta_t, pvalue,rho)
	XD<-cbind(beta_t, pvalue,rho)
	XD<-as.data.frame(XD)
	XD1<-subset(XD, pvalue<alpha)
#	print(length(XD1$rho)==0)
 if(length(XD1$rho)==0){
   XD2<-subset(XD, pvalue<0.1)
	#	print(rho[1:50])
	Rho<-XD2$rho
}else{
 Rho<-	XD1$rho
 	}
	return(Rho)

	
	}
