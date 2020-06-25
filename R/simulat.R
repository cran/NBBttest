simulat <-
    function(yy,nci,r1,r2,p=0,q=0, A=0){
        m<-nrow(yy)
        cn<-ncol(yy)
        xxc<-yy[,seq(nci)]
        #        print(head(yy))
        #       print(dim(yy))
        xxa<-yy[,c(nci+1):c(nci+r1+r2)]
        r<-max(r1,r2)
        colsm1<-rep(0,r)
        for(i in seq(r)){
            
            for(j in 2:m){
                colsm1[i]<-colsm1[i]+xxa[j,i]
            }
        }
      
        v<-matrix(0,nrow=m,ncol=3)
        pp<-rep(1,m)
        sv<-rep(0,m)
        sm<-rep(0,m)
        y1<-matrix(0,nrow=m,ncol=r1)
        y2<-matrix(0,nrow=m,ncol=r2)
        sn<-rep(0,m)
        # condition effect and label of tag expression
        effect<-rep(1,m)
        label<-rep("-", m)
        Prob<-rep(0,m)
        
        for(i in seq(m)){
            if(runif(1)<p){
                effect[i]<-runif(1)*A
                label[i]<-"+"		
            }
        }
         
        sm<-apply(xxa,1,mean) 
        sv<-apply(xxa,1,var)
        #      print(cbind(sm,sv))
        for(i in seq(m)){	
            Prob[i]<-sm[i]/mean(colsm1)				
            rn<-runif(1)
  #          for(j in 1:r){		
       #     	m<-sample(r, 1)
                if(is.na(sm[i])|is.na(sv[i])|sv[i]==0|sm[i]==0){	  		
                if(rn<0.5){
                    y1[i,]<-round((rnbinom(r1,mu=10,
                                            size=10)+effect[i]),digits=0)
                    y2[i,]<-rnbinom(r2,mu=10,10)
                } else{		
                    y1[i,]<-rnbinom(r1,mu=10,size=10)
                    y2[i,]<-round((rnbinom(r2,mu=10,
                                            size=10)+effect[i]),digits=0)	   
                }
                }else{ 		
                if(rn<0.5){
                    y1[i,]<-round((rnbinom(r1,mu=sm[i],
                     size=sv[i]^0.5)+effect[i]),digits=0)
                    y2[i,]<-rnbinom(r2,mu=sm[i],size=sv[i]^0.5)
                } else{
    
                    y1[i,]<-rnbinom(r1,mu=sm[i],size=sv[i]^0.5)
                    y2[i,]<-round((rnbinom(r2,mu=sm[i],
                     size=sv[i]^0.5)+effect[i]),digits=0)	   
                }	

                }	
       #     }		
          
            # artificial effect
        }	
     
        y1<-y1[, seq(r1)]
        y2<-y2[, seq(r2)]		
  
        y<-cbind(y1,y2)
        nr<-nrow(y)
     geneid<-paste(seq(nr), label[seq(nr)],sep="")	

        xx<-cbind(xxc,y)
        colnames(xx)<-colnames(yy)
        xx<-cbind(geneid,xx)
       
        return(xx)
    }
