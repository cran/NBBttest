QC<-
    function(dat, nci, S1="NULL",S2="NULL", method="plot", 
             colrs="greenred", rwcex=1.8, clcex=1.8, x=10, tree="none", 
             log="none", col="blue", pch=19, labsize=1.5, axis=1.5){
        dat<-as.data.frame(dat)	
        nc<-ncol(dat)
        X<-dat[, -seq(nci)]
        
        opar<-par(mar=c(5,5,1,1),oma=c(1,2,1,1), lwd=3)
        on.exit(par(opar))
        if(method=="plot"){
            if(is.null(S1)||is.null(S2)){
                stop("S1 and S2 must be numeric constants for plot")	
            }else if((S1>ncol(dat)||S1<nci)||(S2>ncol(dat)||S2<nci)){
                stop("S1 or S2 is out of data")	
            }	
            X1<-dat[, S1]
            X2<-dat[, S2]
            sn1<-colnames(dat)[S1]
            sn2<-colnames(dat[S2])
            if(log=="log"){
                plot(log2(X1+1),log2(X2+1), col=col,pch=pch, xlab=sn1,
                     yla=sn2, cex.lab=labsize, cex.axis=axis)
                abline(lm(log2(X1+1)~log2(X2+1)), col="red", lwd=3)
                y<-max(log2(X2+1))
                x<-min(log(X1+1))+3
                text(x,y, paste("Person r =", 
                                round(cor(log2(X1+1),log2(X2+1)),4)),cex=1.5)
                box(lwd=3)
            }else{
                plot(X1,X2, col=col, pch=pch, 
                     xlab=sn1,yla=sn2, cex.lab=labsize, cex.axis=axis)		
                abline(lm(X1~X2), col="red", lwd=3)	
                box(lwd=3)
                y<-max(X2)
                x<-max(X1)/10
                text(x,y, paste("Person r =", 
                                round(cor(X1,X2),4)),cex=1.5)	
            }
            
        }else if(method=="heatmap"){
            if(log=="none"){	
                z<-cor(X)
            }else if(log=="log"){
                z<-cor(log2(X+1))	
            }
            
            if (colrs=="redgreen"){
                #cc <- redgreen(ncol(z))
                pltt<-redgreen
            } else if (colrs=="heat.colors"){
                #cc <- heat.colors(ncol(z))
                pltt<-heat.colors
            }else if (colrs=="redblue"){
                #cc <- redblue(ncol(z))
                pltt<-redblue
            }else if(colrs=="greenred"){
                # cc <- greenred(ncol(z))
                pltt<-greenred
            }else if(colrs=="bluered"){
                #cc <- bluered(ncol(z)) 
                pltt<-bluered
            } else if(colrs=="cm.colors"){
                pltt<-cm.colors(x)
            }else if(colrs=="terrain.colors"){
                pltt<-terrain.colors(x)
            }else if(colrs=="topo.colors"){
                pltt<-topo.colors(x)
            }
            
            
            Rowv=TRUE
            Colv=TRUE	
            par(mar=c(5,2,1,7),oma=c(1,1,1,5))
            heatmap.2(z,Rowv=Rowv, Colv=Colv, 
                      col=pltt, trace="none",
                      distfun=function(z)
                      {dist(z,method="euclidean")},
                      key.title=NA, key.xlab=NA, 
                      key.ylab=NA,
                      cexRow=rwcex, cexCol=clcex, 
                      scale ="none", dendrogram=tree, 
                      keysize=1.5, 
                      srtRow=30, srtCol=30, 
                      offsetRow=0, offsetCol=0)	
        }
    }
