pratio <-
function(xx,na,nb){
        xx1<-abs(xx[,seq(na)])
        xx2<-abs(xx[,(na+1):(na+nb)])	
        XM<-apply(xx1,1,max)
        Xm<-apply(xx1,1,min)
        YM<-apply(xx2,1,max)
        Ym<-apply(xx2,1,min)	
        Ng<-nrow(xx)
        #pratio<-rep(0,Ng)	
        #for(i in 1:Ng){
        A<-Ym/(XM+1)
        B<-Xm/(YM+1)
        AB<-cbind(A,B)
        pratio<-apply(AB, 1, max)	
        #}		
        return(pratio)		
    }
