myheatmap2 <-
    function(dat,
             IDcol = 1,
             nci = NULL,
             r,
             colrs = "greenred",
             rwcex = 2.8,
             clcex = 2.8,
             x = 10,
             tree = "both",
             method = "euclidean",
             ky = 1.5,
             rowBarColor = NULL,
             colBarColor = NULL,
             labrow = "yes",
             labcol = "yes",
             rsort = "yes",
             adjrow = c(0.2, 0.0),
             adjcol = c(1, 1) ,
             rwangle = 30,
             clangle = 30,
             maptitle = "",
             keyvalue
             ) {
        # import data
        #require(gplots)
        dat <- as.data.frame(dat)
        try (if (is.null(dim(dat)))
            stop("No data for making heatmap")
            )
            
          colname<-colnames(dat) 
          tname<-c("tv","tvalue","statistic","t_value","t.value") 
          signame<-c("signi","significant","signif","significance","selection",
          "select")
          indt<-is.element(colname,tname)
          jt<-which(indt==TRUE)
            indp<-is.element(colname,signame)
          jp<-which(indp==TRUE)
               sig <- dat[,jp]
       # tvalue <- dat$tvalue
        if (length(sig) == 0) {
            dat1 <- dat
            #dat2<-dat1[, c((nci+1):(nci+r))]
        } else {
            dat1 <- subset(dat, sig==1)
            try (if (is.null(dim(dat1)))
                stop("No differential data for making heatmap")
                )
            if (rsort == "yes") {
        	if(length(jt)!=0){
            dat1 <- dat1[order(dat1[, jt]), ]
            }else{
         dat1 <- dat1[order(dat1[, jp]), ]    	
            } 	
        }
     }
    
        dat2 <- dat1[, (nci + 1):(nci + r)]
        #print(head(dat2))
        mxr <- apply(dat2, 1, max, na.rm = TRUE)
        mx <- max(mxr)
        nr1 <- nrow(dat2)
        nc1 <- ncol(dat2)
        z <- matrix(NA, nrow = nr1, ncol = nc1)
        for (i in seq(nr1)) {
            if (mxr[i] != Inf) {
                for (j in seq(nc1)) {
                    if (!is.na(dat2[i, j])) {
                        z[i, j] <- dat2[i, j] * mx / mxr[i]
                    }
                }
            }
        }
        colnames(z) <- colnames(dat2)
        
        if (is.null(IDcol)) {
            row.names(z) <- row.names(dat1)
        } else{
            row.names(z) <- dat1[, IDcol]
        }
        rc <- redgreen(nrow(z))
        #print(z)
        if (colrs == "redgreen") {
            #cc <- redgreen(ncol(z))
            pltt <- redgreen(x)
        } else if (colrs == "heat.colors") {
            #cc <- heat.colors(ncol(z))
            pltt <- function(x)rev(heat.colors(x))
        } else if (colrs == "redblue") {
            #cc <- redblue(ncol(z))
            pltt <- redblue(x)
        } else if (colrs == "greenred") {
            # cc <- greenred(ncol(z))
            pltt <- greenred(x)
        } else if (colrs == "bluered") {
            #cc <- bluered(ncol(z))
            pltt <- bluered
        } else if (colrs == "cm.colors") {
            pltt <- cm.colors(x)
        } else if (colrs == "terrain.colors") {
            pltt <- terrain.colors(x)
        } else if (colrs == "topo.colors") {
            pltt <- topo.colors(x)
        } else if (colrs == "cyanred") {
            pltt = colorpanel(256,
                              low = "cyan",
                              mid = "black",
                              high = "red")
        }
        if (tree != "none") {
            if (method == "pearson") {
                hr <-
                  hclust(as.dist(1 - cor(t(z),
                      method = "pearson")), method = "complete")
                hc <-
                    hclust(as.dist(1 - cor(z, 
                    method = "pearson")), method = "complete")
            } else if (method == "spearman") {
                hr <-
                    hclust(as.dist(1 - cor(t(z), 
                     method = "spearman")), method = "complete")
                hc <-
                    hclust(as.dist(1 - cor(z, 
                    method = "spearman")), method = "complete")
            } else if (method == "kendall") {
                hr <-
                    hclust(as.dist(1 - cor(t(z), 
                    method = "kendall")), method = "complete")
                hc <-
                    hclust(as.dist(1 - cor(z, 
                   method = "kendall")), method = "complete")
            } else if (method == "euclidean") {
                hr <- hclust(dist(t(z), method = "euclidean"))
                hc <- hclust(dist(z, method = "euclidean"))
            }
            if (tree == "row") {
                if (method == "euclidean") {
                    Rowv = TRUE
                }
                else{
                    Rowv = as.dendrogram(hr)
                }
                Colv = FALSE
            } else if (tree == "column") {
                Rowv = FALSE
                if (method == "euclidean") {
                    Colv = TRUE
                } else{
                    Colv = as.dendrogram(hc)
                }
            } else if (tree == "both") {
                if (method == "euclidean") {
                    Rowv = TRUE
                    Colv = TRUE
                } else{
                    Rowv = as.dendrogram(hr)
                    Colv = as.dendrogram(hc)
                }
            }
        }
        if (labrow == "yes") {
            rwnames = rownames(z)
        } else{
            rwnames = ""
        }
        if (labcol == "yes") {
            clnames = colnames(z)
        } else{
            clnames = ""
        }
        
      
        if (tree == "none") {
            if (is.null(rowBarColor) & is.null(colBarColor)) {
                heatmap.2(
                    z,
                    Rowv = FALSE,
                    Colv = FALSE,
                    col = pltt,
                    trace = "none",
                    scale = "none",
                    dendrogram = "none",
                    keysize = ky,
                    key.title = NA,
                    key.xlab = keyvalue,
                    key.ylab = NA,
                    main = maptitle,
                    cex.main = 5.0,
                    adjRow = adjrow,
                    adjCol = adjcol,
                    labRow = rwnames,
                    labCol = clnames,
                    cexRow = rwcex,
                    cexCol = clcex,
                    srtRow = rwangle,
                    srtCol = clangle,
                    offsetRow = 1.5,
                    offsetCol = 0,
                     key.par=list(mar=c(4,0,1,2), cex=1.0, cex.lab=1.5, cex.axis=1.0),
                     lwid = c(2, 4)
                     #  lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(1.5, 4, 2 )
                )
                
            } else if (is.null(rowBarColor) & !is.null(colBarColor)) {
                heatmap.2(
                    z,
                    Rowv = FALSE,
                    Colv = FALSE,
                    col = pltt,
                    trace = "none",
                    scale = "none",
                    dendrogram = "none",
                    keysize = ky,
                    key.title =NA,
                    key.xlab = keyvalue,
                    key.ylab = NA,
                    colCol = colBarColor,
                    ColSideColors = colBarColor,
                    main = maptitle,
                    cex.main = 5.0,
                    adjRow = adjrow,
                    adjCol = adjcol,
                    labRow = rwnames,
                    labCol = clnames,
                    cexRow = rwcex,
                    cexCol = clcex,
                    srtRow = rwangle,
                    srtCol = clangle,
                    offsetRow = 1.5,
                    offsetCol = 0,
                    key.par=list(mar=c(4,0,1,2), cex=1.0, cex.lab=1.3, cex.axis=1.5),
                    lwid = c(3, 4)
                    #lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(1.5, 4, 2 )
                )
            } else if (!is.null(rowBarColor) & is.null(colBarColor)) {
                heatmap.2(
                    z,
                    Rowv = FALSE,
                    Colv = FALSE,
                    col = pltt,
                    trace = "none",
                    scale = "none",
                    dendrogram = "none",
                    keysize = ky,
                    key.title = NA,
                    key.xlab = keyvalue,
                    key.ylab = NA,
                    colRow = rowBarColor,
                    RowSideColors = rowBarColor,
                    main = maptitle,
                    cex.main = 5.0,
                    adjRow = adjrow,
                    adjCol = adjcol,
                    labRow = rwnames,
                    labCol = clnames,
                    cexRow = rwcex,
                    cexCol = clcex,
                    srtRow = rwangle,
                    srtCol = clangle,
                    offsetRow = 1.5,
                    offsetCol = 0,
                   key.par=list(mar=c(4,0,1,2), cex=1.0, cex.lab=1.3, cex.axis=1.5),
                   lwid =c(3, 4)
                   #   lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(1.5, 4, 2 )
                )
                
            } else{
                heatmap.2(
                    z,
                    Rowv = FALSE,
                    Colv = FALSE,
                    col = pltt,
                    trace = "none",
                    scale = "none",
                    dendrogram = "none",
                    keysize = ky,
                    key.title = NA,
                    key.xlab = keyvalue,
                    key.ylab = NA,
                    colRow = rowBarColor,
                    RowSideColors = rowBarColor,
                    colCol = colBarColor,
                    ColSideColors = colBarColor,
                    main = maptitle,
                    cex.main = 5.0,
                    adjRow = adjrow,
                    adjCol = adjcol,
                    labRow = rwnames,
                    labCol = clnames,
                    cexRow = rwcex,
                    cexCol = clcex,
                    srtRow = rwangle,
                    srtCol = clangle,
                    offsetRow = 1.5,
                    offsetCol = 0,
                   key.par=list(mar=c(4,0,1,2), cex=1.0, cex.lab=1.30, cex.axis=1.5),
                   lwid=c(3, 4)
                   #      lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(1.5, 4, 2 )
                )
                
            }
            
        } else{
            if (is.null(rowBarColor) & is.null(colBarColor)) {
                heatmap.2(
                    z,
                    Rowv = Rowv,
                    Colv = Colv,
                    col = pltt,
                    trace = "none",
                    distfun = function(z)
                    {
                        if (method == "euclidean") {
                       dist(z <-
                       ifelse(!is.na(z), z, 0), method = "euclidean")
                        } else if (method == "pearson") {
                        as.dist(1 - cor(t(z), method = "pearson"))
                        } else if (method == "spearman") {
                         as.dist(1 - cor(z, method = "spearman"))
                        } else if (method == "kendall") {
                         as.dist(1 - cor(z, method = "kendall"))
                        }
                    },
                    scale = "none",
                    dendrogram = tree,
                    keysize = ky,
                    key.title = NA,
                    key.xlab = keyvalue,
                    key.ylab = NA,
                    main = maptitle,
                    cex.main = 5.0,
                    adjRow = adjrow,
                    adjCol = adjcol,
                    labRow = rwnames,
                    labCol = clnames,
                    cexRow = rwcex,
                    cexCol = clcex,
                    srtRow = rwangle,
                    srtCol = clangle,
                    offsetRow = 1.5,
                    offsetCol = 0,
                     key.par=list(mar=c(4,0,1,2), cex=1.0, cex.lab=1.30, cex.axis=1.5),
                     lwid = c(3, 4)
       #                    lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(1.5, 4, 2 )
                )
            } else if (!is.null(rowBarColor) 
                       & is.null(colBarColor)) {
                heatmap.2(
                    z,
                    Rowv = Rowv,
                    Colv = Colv,
                    col = pltt,
                    trace = "none",
                    distfun = function(z)
                    {
                        if (method == "euclidean") {
                        dist(z <-
                        ifelse(!is.na(z), z, 0), method = "euclidean")
                        } else if (method == "pearson") {
                        as.dist(1 - cor(t(z), method = "pearson"))
                        } else if (method == "spearman") {
                         as.dist(1 - cor(z, method = "spearman"))
                        } else if (method == "kendall") {
                        as.dist(1 - cor(z, method = "kendall"))
                        }
                    },
                    scale = "none",
                    dendrogram = tree,
                    keysize = ky,
                    key.title = NA,
                    key.xlab = keyvalue,
                    key.ylab = NA,
                    colRow = rowBarColor,
                    RowSideColors = rowBarColor,
                    main = maptitle,
                    cex.main = 5.0,
                    adjRow = adjrow,
                    adjCol = adjcol,
                    labRow = rwnames,
                    labCol = clnames,
                    cexRow = rwcex,
                    cexCol = clcex,
                    srtRow = rwangle,
                    srtCol = clangle,
                    offsetRow = 1.5,
                    offsetCol = 0,
                     key.par=list(mar=c(4,0,1,2), cex=1.0, cex.lab=1.30, cex.axis=1.5),
                     lwid = c(3, 4)
                     #                  lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(1.5, 4, 2 )
   
                )
            } else if (is.null(rowBarColor)
                       & !is.null(colBarColor)) {
                heatmap.2(
                    z,
                    Rowv = Rowv,
                    Colv = Colv,
                    col = pltt,
                    trace = "none",
                    distfun = function(z)
                    {
                        if (method == "euclidean") {
                        dist(z <-
                        ifelse(!is.na(z), z, 0), method = "euclidean")
                        } else if (method == "pearson") {
                        as.dist(1 - cor(t(z), method = "pearson"))
                        } else if (method == "spearman") {
                         as.dist(1 - cor(z, method = "spearman"))
                        } else if (method == "kendall") {
                        as.dist(1 - cor(z, method = "kendall"))
                        }
                    },
                    scale = "none",
                    dendrogram = tree,
                    keysize = ky,
                    key.title = NA,
                    key.xlab = keyvalue,
                    key.ylab = NA,
                    colCol = colBarColor,
                    ColSideColors = colBarColor,
                    main = maptitle,
                    cex.main = 5.0,
                    adjRow = adjrow,
                    adjCol = adjcol,
                    labRow = rwnames,
                    labCol = clnames,
                    cexRow = rwcex,
                    cexCol = clcex,
                    srtRow = rwangle,
                    srtCol = clangle,
                    offsetRow = 1.5,
                    offsetCol = 0,
                    key.par=list(mar=c(4,0,1,2), cex=1.0, cex.lab=1.30, cex.axis=1.5),
                    lwid = c(3, 10)
                    #                   lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(1.5, 4, 2 )
                )
            } else {
                heatmap.2(
                    z,
                    Rowv = Rowv,
                    Colv = Colv,
                    col = pltt,
                    trace = "none",
                    distfun = function(z)
                    {
                        if (method == "euclidean") {
                       dist(z <-
                       ifelse(!is.na(z), z, 0), method = "euclidean")
                        } else if (method == "pearson") {
                         as.dist(1 - cor(t(z), method = "pearson"))
                        } else if (method == "spearman") {
                        as.dist(1 - cor(z, method = "spearman"))
                        } else if (method == "kendall") {
                        as.dist(1 - cor(z, method = "kendall"))
                        }
                    },
                    scale = "none",
                    dendrogram = tree,
                    keysize = ky,
                    key.title = NA,
                    key.xlab = keyvalue,
                    key.ylab = NA,
                    colRow = rowBarColor,
                    RowSideColors = as.character(rowBarColor),
                    colCol = colBarColor,
                    ColSideColors = colBarColor,
                    main = maptitle,
                    cex.main = 5.0,
                    adjRow = adjrow,
                    adjCol = adjcol,
                    labRow = rwnames,
                    labCol = clnames,
                    cexRow = rwcex,
                    cexCol = clcex,
                    srtRow = rwangle,
                    srtCol = clangle,
                    offsetRow = 1.5,
                    offsetCol = 0,
                    key.par=list(mar=c(4,0,1,2), cex=1.0, cex.lab=1.30, cex.axis=1.5),
                    lwid = c(3, 4)
                    #                    lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(1.5, 4, 2 )
                 
                   )
                }
            }

    }
