mbetattest<-
    function(X, nci, na, nb, alpha=0.05, norm="no",
                   side="both", level="sgRNA",padjust_methods,C=1.222){
        
  #     X<-as.data.frame(X)
  cn<-ncol(X)
  rn<-nrow(X)
  
        #	X<- X[with(X, order(Gene)), ]
        colname<-colnames(X)
        
        genename<-c("gene","genes", "Genes","Gene",
                    "gene name","Gene name","gene_name",
                    "Gene_name", "name","Name")
                    
       geneid<-c("geneid","gene_id","Geneid","gene id","Gene ID","Gene_id",
                   "Gene_ID", "gene ID","gene ID")
        indx<-is.element(colname, genename)
        jj<-which(indx==TRUE)
        
        if(length(jj)==0){
            Gene<-rep("gene", nrow(X))	
            X<-cbind(as.data.frame(Gene), X)
            nci<-nci+1
            
        }else if(length(jj)==1){
            colnames(X)[jj]<-"Gene"
        }else if(length(jj)>1){
            colnames(X)[jj[2]]<-"Gene"	
        }
        
        indx<-is.element(colname, geneid)
        jj<-which(indx==TRUE)
        
        if(length(jj)==0){
            Gene_id<-rep("Gene_id", nrow(X))
            X<-cbind(as.data.frame(Gene_id), X)
            nci<-nci+1
            
        }else if(length(jj)==1){
            colnames(X)[jj]<-"Gene_id"
        }else if(length(jj)>1){
            colnames(X)[jj[2]]<-"Gene_id"
        }
        
         sm<-apply(X[, c((nci+1) : ncol(X))],1, sum)
         X<-subset(X, sm>0&!is.na(sm))
         X<-na.omit(X)
        if(norm=="lg2"){	
            NX<-normalized(dat=X, nci=nci, m=0,lg2="yes")
        }else if(norm=="yes"){
            NX<-normalized(dat=X, nci=nci, m=0)
            if (nci==1){
            NX<-cbind(as.data.frame(row.names(NX)),NX)
            }
        }else{	
            NX<-X	
        }
   
    #   NX<-cbind(as.data.frame(row.names(NX)),NX)
    #  print(head(NX))
        
        rn<-nrow(NX)

        if(rn>=10000){
            m<-sample(rn, 10000)
            w<-omega(XX=NX[m,], nci=nci, 
                     r1=na, r2=nb, sn=4, alpha=alpha)
        }else if(rn<10000&rn>=5000){
            w<-omega(XX=NX, nci=nci, r1=na, 
                     r2=nb, sn=8, alpha=alpha)
        }else if(rn<5000&rn>=1000){
            w<-omega(XX=NX, nci=nci, r1=na, 
                     r2=nb, sn=15, alpha=alpha)		
        }else if(rn<1000&rn>=500){
            w<-omega(XX=NX, nci=nci, r1=na, 
                     r2=nb, sn=20, alpha=alpha)			
        }else{
            w<-omega(XX=NX, nci=nci, r1=na, 
                     r2=nb, sn=30, alpha=alpha)			
        }
        
        Betattest<-gbetattest(xx=NX, W=w, nci=nci,
         na=na, nb=nb,level=level, padjust_methods,C=1.222,side=side)

        res<-cbind(Betattest, w)
        res<-as.data.frame(res)	
        if(length(jj)==0){
            res<-res[,-1]	
        }
        
        return(res)
        
    }
