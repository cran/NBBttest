betattest <-
    function(X, na, nb, NX=100, level){
     
        N=length(X[,1])
        n<-na+nb
        XA<-X[, seq(na)]
        XB<-X[, (na+1):n]	
        sa<-apply(XA, 2, sum, na.rm=TRUE)
        sb<-apply(XB, 2, sum, na.rm=TRUE)
        ma<-apply(XA, 1, mean, na.rm=TRUE)
        mb<-apply(XB, 1, mean, na.rm=TRUE)
        SM<-(max(sa,sb))	
        #	print(SM)
        if(level=="sgRNA"|level=="isoform"){
            SA<-rep(SM, na)
            SB<-rep(SM, nb)
        }else{
            SA<-sa
            SB<-sb
        }
        PA<-betaparametVP(XA, SA)
        Pa<-PA[,1]
        Va<-PA[,2]
        PB<-betaparametVP(XB, SB)
        Pb<-PB[,1]
        Vb<-PB[,2]
  
        tv<-rep(0,N)
        df<-rep(0,N)
        aa<-rep(NA,N)
        for(i in seq(N)){
            if(Va[i]+Vb[i]>0){			
                tv[i]<-(Pa[i]-Pb[i])/(Va[i]+Vb[i])^0.5	
            }
            if(ma[i]<=1|is.na(ma[i])){
                ma[i]<-2
            }
            if(mb[i]<=1|is.na(mb[i])){
                mb[i]<-2
            }	
            if(level=="sgRNA"){
                if(abs(ma[i]-mb[i])==0){
                    aa[i]<-1
                }else{ 
                 aa[i]<-sqrt(abs(ma[i]-mb[i])/(1+(ma[i]+mb[i])))
                    #  aa[i]<-abs(ma[i]-mb[i])/15
              }
                df[i]<-aa[i]*(Va[i]+Vb[i])^2/(Va[i]^2/(sum(sa)-1)
                                        +Vb[i]^2/(sum(sb)-1))
            }else if(level=="isoform"){
                 if(abs(ma[i]-mb[i])==0){
                    aa[i]<-1
                }else{ 
                 aa[i]<-sqrt(abs(ma[i]-mb[i])/(1+(ma[i]+mb[i])))
                    #  aa[i]<-abs(ma[i]-mb[i])/15
              }
              
                df[i]<-aa[i]*((Va[i]+Vb[i])^2/(Va[i]^2/(ma[i]-1)+Vb[i]^2/(mb[i]-1)))
                #                          df[i]<-(Va[i]+Vb[i])^2/(Va[i]^2/(na-1)
                #                       +Vb[i]^2/(nb-1))
            }else if(level=="RNA"){
                #                 df[i]<-sqrt((Va[i]+Vb[i])^2/(Va[i]^2/(mean(sa)-1)
                if(abs(ma[i]-mb[i])==0){
                    aa[i]<-1
                }else{
                    aa[i]<-sqrt(abs(ma[i]-mb[i])/(1+(ma[i]+mb[i])))
                    #  aa[i]<-abs(ma[i]-mb[i])/15
                }
             df[i]<-aa[i]*(Va[i]+Vb[i])^2/(Va[i]^2/(sqrt((ma[i]))-1)+Vb[i]^2/(sqrt((mb[i]))-1))
        
            }
        }
        for(i in seq(N)){
            if(is.na(df[i])|df[i]==0){
                df[i]<-2
            }
        }
        
        #  tv<-tv^2
        
        ttest<-cbind(tv,df)
        return(ttest)
        
    }
