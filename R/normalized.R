normalized<-
    function(dat, nci, m=0,lg2="no"){
        if (nci==0) {
            dat<-as.data.frame(dat)	
            rowname<-row.names(dat)
            dat1<-subset(dat, rowname!="N_unmapped") 
            dat1<-subset(dat1, rowname!="N_noFeature") 
            dat1<-subset(dat1, rowname!="N_multimapping") 
            dat1<-subset(dat1, rowname!="MIR6087") 
            dat1<-subset(dat1, rowname!="N_ambiguous") 
            
            dat2<-na.omit(dat1)
            mean<-apply(dat2, 1, mean)
            dat3<-cbind(dat2, mean)
            dat4<-subset(dat3, mean>m)
            dat5<-dat4[,-ncol(dat4)]
            csm<-apply(dat5,2,sum)
            
            meanc<-mean(csm)
            
            dat8<-matrix(NA,nrow(dat5),ncol(dat5))
            
            for(i in seq(nrow(dat5))){
                for(j in seq(ncol(dat5))){
                    dat8[i,j]<-dat5[i,j]*meanc/csm[j]	
                }
            }
            
            row.names(dat8)<-row.names(dat5)
            colnames(dat8)<-colnames(dat5)
        }else if(nci==1){
            dat1<-dat[,-1]
            rowname<-dat[,1]
            dat1<-as.data.frame(dat1)
            #print(dat1)
            row.names(dat1)<-rowname
            dat2<-subset(dat1, rowname!="N_unmapped") 
            dat2<-subset(dat2, rowname!="N_noFeature") 
            dat2<-subset(dat2, rowname!="N_multimapping") 
            dat2<-subset(dat2, rowname!="MIR6087") 
            dat2<-subset(dat2, rowname!="N_ambiguous") 
            dat3<-na.omit(dat2)
            dat4<-apply(dat3, 2, as.numeric)
            row.names(dat4)<-row.names(dat3)
            mean<-apply(dat4,1,mean)
            dat5<-cbind(dat4,mean)
            dat6<-subset(dat5, mean>m)
            dat7<-dat6[,-ncol(dat6)]
            csm<-apply(dat7, 2, sum)
            meanc<-mean(csm)
            dat8<-matrix(NA, nrow(dat7), ncol(dat7))
            for(i in seq(nrow(dat7))){
                for(j in seq(ncol(dat7))){
                    dat8[i,j]<-dat7[i,j]*meanc/csm[j]	
                }
            }
            row.names(dat8)<-row.names(dat7)
            colnames(dat8)<-colnames(dat3)
        }else {
            nc<-ncol(dat)	
            dat1<-dat[, c((nci+1):nc)]
            dat2<-na.omit(dat1)
            dat3<-apply(dat2, 2, as.numeric)
            mean<-apply(dat3,1,mean)
            dat4<-cbind(as.data.frame(dat[, seq(nci)]), dat3, mean)
            # print(head(dat4))
            dat5<-subset(dat4, mean>m)
            dat6<-dat5[, c((1+nci):nc)]
            csm<-apply(dat6, 2, sum)
            meanc<-mean(csm)
            dat7<-matrix(NA,nrow(dat6), ncol(dat6))
            for(i in seq(nrow(dat6))){
                for(j in seq(ncol(dat6))){
                	
                    dat7[i, j]<-dat6[i, j]*meanc/csm[j]	
                }
            }
            if(lg2=="yes"){
            	dat7<-log2(1+dat7)
            }
            colnames(dat7)<-colnames(dat6)	
            dat8<-cbind(dat5[,seq(nci)],dat7)
            #print(csm<-apply(dat8, 2, sum))
        }
        return(dat8)
    }


