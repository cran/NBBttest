betaparametw <-
    function(xn,a,b){	
        n=length(xn)
        w<-rep(0,n)
        W<-rep(0,n)
        for(i in seq(n)){
            w[i]<-(a+b)*xn[i]/(a+b+xn[i])
        }
        for(i in seq(n)){ 
            if(!is.na(sum(w))){
                W[i]<-w[i]/sum(w) 
            }
        } 
        return(W)	
    }
