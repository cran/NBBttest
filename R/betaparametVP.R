betaparametVP <-
    function(X, NX){
  
rn=length(X[,1])
cn=length(X[1,])
MX=apply(X,1,sum, na.rm=TRUE)
 P<-rep(0,rn) 	
 V<-rep(0,rn)
VX2<-rep(0,rn)
 VN<-rep(0,rn)
 VA<-rep(0,rn)
AA<-matrix(0,nrow=rn,ncol=100)
BB<-matrix(0,nrow=rn,ncol=100)
VVX2<-matrix(0,nrow=rn,ncol=100)	
WW<-array(0,c(rn,100,cn))
WW2<-array(0,c(rn,100,cn))	
VV2<-array(0,c(rn,100,cn))	
SWW2<-matrix(0,nrow=rn,ncol=100)
SW<-rep(0,rn)
PP<-matrix(0,nrow=rn,ncol=100)
PW<-array(0,c(rn,100,cn))
VV<-matrix(0,nrow=rn,ncol=100)
D<-matrix(NA,nrow=rn,ncol=100)
Q<-X/NX# matrix(0,nrow=rn,ncol=cn)
w<-NX/sum(NX)	
W2<-w^2
for(i in seq(rn)){		   	
   if(MX[i]>0){		 	
      for( j in seq(cn)){
 P[i]<-P[i]+w[j]*Q[i,j]
           }	 
          for(j in seq(cn)){
 VX2[i]<-VX2[i]+(w[j]*Q[i,j])^2
         }			 	
V[i]<-(V[i]-sum(W2)*P[i])^2/(1-sum(W2))
ab<-betaparametab(NX,w,P[i],V[i])			
AA[i,1]<-ab[1]		
BB[i,1]<-ab[2]						
VVX2[i,1]<-VX2[i]
for(j in seq(cn)){
  WW2[i,1,j]<-W2[j]
   WW[i,1,j]<-w[j]	
  }
           
PP[i,1]<-P[i]			
VV[i,1]<-V[i]
 D[i,1]<-10
 k=1
 while(D[i,k]>0.01){
 k<-k+1	
AB<-betaparametab(xn=NX,w=WW[i,c(k-1),],
P=PP[i,(k-1)],V=VV[i,c(k-1)])
                 
AA[i,k]<-AB[1]
 BB[i,k]<-AB[2]				
 EW<-betaparametw(xn=NX,a=AA[i,k],b=BB[i,k])		
WW[i,k,]<-EW
 SWW2[i,k]<-0
                    
for(j in seq(cn)){
PP[i,k]<-PP[i,k]+WW[i,k,j]*Q[i,j]
SWW2[i,k]<-SWW2[i,k]+WW[i,k,j]^2
 VVX2[i,k]<-VVX2[i,k]+(WW[i,k,j]*Q[i,j])^2			    	
     }		 			     
for(j in seq(cn)){
VV[i,k]<-VV[i,k]+(WW[i,k,j]*Q[i,j]
                  -PP[i,k]/cn)^2/(1-SWW2[i,k])
   }		
                  
   if(k<100){
    D[i,k]<-abs(BB[i,c(k-1)]-BB[i,k])
     }else{
    D[i,k]<-0
      }
    if(is.na(D[i,k])){
        D[i,k]<-0
    }
 }
                
  VN[i]<-VV[i,k]
  SW[i]<-SWW2[i,k]				
   P[i]<-PP[i,k]	
VA[i]<-((1+MX[i])/sum(NX))*(1-(1+
                   MX[i])/sum(NX))/sum(NX)
 if(is.na(VN[i])==TRUE){
    VN[i]<-0
     }                  
   V[i]<-max(VN[i],VA[i])	    		
    if(is.na(P[i])){
     P[i]=0	
     for(j in seq(cn)){
      P[i]<-P[i]+w[j]*Q[i,j]
       }	  	
    }					               	    	
 }		
}     
paramet<-cbind(P,V)
return(paramet)
    }
