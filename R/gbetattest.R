gbetattest<-function(xx, W, nci, na, nb, level, padjust_methods,C=1.222,side){
    
    cn<-ncol(xx)	
    rn<-nrow(xx)
    
    size<-max(apply(xx[, c((nci+1):cn)], 2, sum))
    
    genes<-unique(xx$Gene)
    
    Gene<-xx$Gene
    
    colname<-colnames(xx)
    
    labname<-c("Class","class","essential", "Essential",
               "lab","Lab", "label","Label")
    
    indx<-is.element(colname, labname)
    jj<-which(indx==TRUE)
    
    if(length(jj)==0){
        geneclass<-Gene
    }else{
        lab<-xx[,jj]
        geneclass<-paste(Gene, lab, sep="_")
    }
    gn<-length(genes)
    
    ggroup<-function(dat, na, nb,W, NX=100, level="RNA"){
        cn<-ncol(dat)	
        ssx<-dat[, c((cn-na-nb+1):cn)]
        prat<-pratio(xx=ssx, na=na, nb=nb)
        odrat<-oddratio(XX=ssx, na=na,nb=nb)
        rho<-sqrt(prat*odrat)/W
        
        grho<-rep(mean(rho),length(rho))
        
        ttest<-betattest(X=ssx, na=na, nb=nb,NX=NX, level=level)
        ttest<-cbind(dat, ttest, rho, grho)	
        return(ttest)
    }
    mgroup<-function(xx,nci, na,nb, gn, genes){
        rhov<-function(xx, na,nb, cn)	{
            ssx<-xx[, c((cn-na-nb+1):cn)]
            prat<-pratio(xx=ssx, na=na, nb=nb)
            odrat<-oddratio(XX=ssx, na=na,nb=nb)
            rho<-sqrt(prat*odrat)/W
            return(rho)
        }
        cn<-ncol(xx)	
        grho<-rep(NA,gn)
        grho[1]<-mean(rhov(subset(xx, xx$Gene==genes[1]),
                           na=na, nb=nb,cn=cn))
        dat<-apply((subset(xx[, c((cn-na-nb+1):cn)], 
                           xx$Gene==genes[1])), 2, mean)
        
        for(i in 2:gn){
            grho[i]<-mean(rhov(subset(xx, xx$Gene==genes[i]),
                               na=na, nb=nb,cn=cn))
            dat1<-apply((subset(xx[, c((cn-na-nb+1):cn)], 
                                xx$Gene==genes[i])), 2, mean)	
            dat<-rbind(dat, dat1)
        }
        ttest<-betattest(X=dat, na=na, nb=nb, level="gene")
        ttest<-cbind(dat, grho, ttest)
        return(ttest)
    }
    
    if(level=="RNA"){
        ttest1<-ggroup(dat=xx,na=na, nb=nb,W=W)	
    }else if(level=="isoform"){	
        data.sg<-subdata(xx=xx, sg=1)
        
        data.mg<-subdata(xx=xx,sg=2)
        #       print(nrow(data.mg))
        sum.mg<-apply(data.mg[, c((cn-na-nb+1):cn)],  
                      1, sum, na.rm=TRUE)	
        N.mg=max(sum.mg)
        ttest<-ggroup(dat=data.sg,na=na, nb=nb,W=W)	
        genes<-unique(data.mg$Gene)
        ng<-length(genes)
        for (i in seq(ng)){
            sbdata.mg<-(subset(data.mg, data.mg$Gene==genes[i]))	
            ttest1<-ggroup(dat=sbdata.mg, na=na, nb=nb,
                           W=W,  NX=N.mg, level="isoform")	 
            ttest<-rbind(ttest, ttest1)	
        }
        ttest2<-ttest		   
    }else if(level=="sgRNA"){
        genes<-unique(xx$Gene)
        ng<-length(genes)	
        ttest<-ggroup(dat=subset(xx, xx$Gene==genes[1]),
                      na=na, nb=nb,W=W, level="sgRNA")	
        for (i in 2:ng){	
            sbdata<-(subset(xx, xx$Gene==genes[i]))	
            ttest1<-ggroup(dat=sbdata, na=na, nb=nb,
                           W=W, level="sgRNA")	
            ttest<-rbind(ttest, ttest1)	
        }
        ttest2<-ttest
        
    }else if(level=="polyA.gene"||level=="PA.gene"|level=="CRISPR.gene"){
    
        ttest3<-mgroup(xx, nci=nci, na=na,nb=nb, 
                       gn=gn, genes=genes)
                    
        ttest3<-as.data.frame(ttest3)	
    }
    
    if(level=="RNA"){
        xc<-ttest1[, seq(nci+na+nb)]
        tv<-ttest1$tv
        df<-ttest1$df
        rho<-ttest1$rho
        
        tvalue<-tv*rho
        
        if(side=="both"){
            pvalue<-2*(1-pt(abs(tvalue), df=df))
        }else if(side=="up"){
            pvalue<-(1-pt(tvalue, df=df))	
        }else if(side=="down"){
            pvalue<-pt(tvalue, df=df)		
        }
        pv<-pvalue
        result<-cbind(xc, tvalue, rho, pvalue)
    }else if(level=="sgRNA"){
        
        xc<-ttest2[, seq(nci+na+nb)]
        
        tvalue<-ttest2$tv
        df<-ttest2$df
        rho<-ttest2$grho
        
        tv<-tvalue*rho
        
        if(side=="both"){
            
            pvalue<-2*(1-pt(abs(tv), df=df))
        }else if(side=="up"){
            pvalue<-(1-pt(tv, df=df))		
        }else if(side=="down"){
            pvalue<-pt(tv, df=df)		
        }
        pv<-pvalue
        result<-cbind(xc, tv, rho, pvalue)	
        
    }else if(level=="isoform"){
        xc<-ttest2[, seq(nci+na+nb)]
        
        tvalue<-ttest2$tv
        df<-ttest2$df
        rho<-ttest2$rho
        #tv[i]<-tvalue[i]*rho[i]
        tv<-tvalue*rho
        
        if(side=="both"){
            pvalue<-2*(1-pt(abs(tv), df=df))
        }else if(side=="up"){
            pvalue<-(1-pt(tv, df=df))		
        }else if(side=="down"){
            pvalue<-pt(tv, df=df)		
        }
        
        pv<-pvalue
        result<-cbind(xc, tv, rho, pvalue)	
        
    }else if(level=="CRISPR.gene"){
        
        
        grho<-ttest3$grho
        
        gene.Treatment<-unique(geneclass)
        
        gt<-ttest3$tv
        gdf<-ttest3$df
        
        gtvalue<-gt*grho
        
        if(side=="both"){
            gpvalue<-2*(1-pt(abs(gtvalue),df=gdf))
        }else if(side=="up"){
            gpvalue<-(1-pt(gtvalue, df=gdf))		
        }else if(side=="down"){
            gpvalue<-pt(gtvalue, df=gdf)
        }
        
        Gene<-as.character(gene.Treatment)
        HG<-unlist(strsplit(Gene,split="_"))
        n<-length(Gene)
        i<-seq(n)
        k1<-i*2-1
        k2<-i*2
        gene<-HG[k1]
        class<-(HG[k2])
        ttest3<-ttest3[,-c(ncol(ttest3)-1,ncol(ttest3))]
        result<-cbind( gene, class, ttest3,  gtvalue, gpvalue)
        pv<-gpvalue
        
    }else if(level=="polyA.gene"|level=="PA.gene"){
        grho<-ttest3$grho
        gene.Treatment<-unique(Gene)
        
        gt<-ttest3$tv
        gdf<-ttest3$df
        
        gtvalue<-gt*grho
        
        if(side=="both"){
            gpvalue<-2*(1-pt(abs(gtvalue),df=gdf))
        }else if(side=="up"){
            gpvalue<-(1-pt(gtvalue, df=gdf))		
        }else if(side=="down"){
            gpvalue<-pt(gtvalue, df=gdf)
        }
        
        gene<-as.character(gene.Treatment)
        ttest3<-ttest3[,-c(ncol(ttest3)-1,ncol(ttest3))]
        result<-cbind( gene, ttest3,  gtvalue, gpvalue)
        pv<-gpvalue
        
    }else if(level=="splicing.gene"){
        
        data.sg<-subdata(xx=xx, sg=1)
        
        data.mg<-subdata(xx=xx,sg=2)	
        sum.mg<-apply(data.mg[, c((cn-na-nb+1):cn)], 1, 
                      sum, na.rm=TRUE)	
        N.mg=max(sum.mg)
        ttest3<-ggroup(dat=data.sg,na=na, nb=nb,W=W)	
        genes.mg<-unique(data.mg$Gene)
        ng<-length(genes.mg)
        for (i in seq(ng)){
            sbdata.mg<-subset(data.mg, 
                              data.mg$Gene==genes.mg[i])
            ttest1<-ggroup(dat=sbdata.mg, na=na, nb=nb,
                           W=W,  NX=N.mg, level="isoform")	 
            maxt<-max(abs(ttest1$tv))
            sm2<-apply(ttest1[, c((cn-na-nb+1):cn)], 2, 
                       sum, na.rm=TRUE)	
            ttest2<-subset(ttest1, abs(ttest1$tv)==maxt)
            ttest2[1, c((cn-na-nb+1):cn)]<-sm2
            ttest3<-rbind(ttest3, ttest2[1, ])	
        }
        
        grho<-ttest3$grho
        
        #grho<-unique(rhoa)
        gt<-ttest3$tv
        gdf<-ttest3$df
        #grho<-ttest2$rho
        #gtv[j]<-gt[j]*grho[j]
        
        gtvalue<-gt*grho
        
        if(side=="both"){
            gpvalue<-2*(1-pt(abs(gtvalue),df=gdf))
        }else if(side=="up"){
            gpvalue<-(1-pt(gtvalue, df=gdf))		
        }else if(side=="down"){
            gpvalue<-pt(gtvalue, df=gdf)
        }
        nc<-ncol(ttest3)
        # print(ttest3[1:4,])
        result<-cbind( ttest3[,c(1:(nc-4)),nc],  gtvalue, gpvalue)
        pv<-gpvalue
        
    }
    
    res<-cbind(as.data.frame(result), pv)
    res<-as.data.frame(res)
    res<-res[order(-pv),]
    pv<-res$pv
    if(padjust_methods=="TX"){
    	if(C==0){
    	adjp<-pv	
    	}else{
     padj<-mtpvadjust(pv=pv,C=C)
     }
     #res<-res[,-c(ncol(res)-1)]
   }else{
     padj<-p.adjust(p=pv, method = padjust_methods, n = length(pv))
    }

   res<-res[,-ncol(res)]
   res<-cbind(res,padj)
    return(res)
    
}
