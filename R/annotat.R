annotat<-
          function (infile, mfile, type="gene", colunmset = "NLL")
{
   
    annotat_ENSG <- function(infile, mfile) {
    
        chrmatch <- function(inf, mf) {
            colname <- colnames(inf)
            genename <- c("gene", "gene_id", "Gene", "Gene_id")
            indx <- is.element(colname, genename)
            jj <- which(indx == TRUE)
            colnames(inf)[jj] <- "ENSG"
            ENSGn <- inf$ENSG
            ENSGm <- mf$ENSG
            ind <- is.element(ENSGm, ENSGn)
            chr <- mf[ind, c(6, 8)]
            chrb <- unique(chr)
            chrc <- chrb[order(chrb[, 1]), ]
            infa <- inf[order(inf[, 2]), ]
            inf <- cbind(as.data.frame(chrc), as.data.frame(infa[, 
                -1]))
            return(inf)
        }
        colname <- colnames(infile)
        genename <- c("chr", "contig", "Chr")
        indx <- is.element(colname, genename)
        jj <- which(indx == TRUE)
        colnames(infile)[jj] <- "chr"
        colnames(mfile) <- c("chr", "exon", "start", "end", "strand", 
            "ENSG", "exon#", "Gene")
        chr_dat <- chrmatch(inf = infile, mf = mfile)
        return(chr_dat)
    }
    
    annotat_isof <- function(infile, mfile) {
    
        chrmatch <- function(inf, mf) {
            colname <- colnames(infile)
            genename <- c("gene_id", "Gene", "Gene_id","gene_name")
            indx <- is.element(colname, genename)
            jj <- which(indx == TRUE)
            colnames(inf)[jj] <- "Gene_id"
          
            ENSGn <- inf$Gene_id
            ENSGm <- mf$ENSG
            ind <- is.element(ENSGm, ENSGn)
            chr <- mf[ind, c(6, 8)]
            chrb <- unique(chr)
            gene <- rep(NA, length(ENSGn))
            GS <- as.character(chrb$Gene)
            for (i in seq(length(ENSGn))) {
                for (j in seq(nrow(chrb))) {
                  if (length(is.element(ENSGn[i], chrb[j, 1])) > 
                    0) {
                    if (is.element(ENSGn[i], chrb[j, 1]) == T) {
                      gene[i] <- GS[j]
                    }
                  }
                }
            }
            if (nrow(inf) == 0) {
                inf <- "NA"
            }
            else {
                inf <- cbind(as.data.frame(gene), as.data.frame(inf))
            }
            return(inf)
        }
        colname <- colnames(infile)
        genename <- c("chr", "contig", "Chr")
        indx <- is.element(colname, genename)
        jj <- which(indx == TRUE)
        colnames(infile)[jj] <- "chr"
        list <- as.data.frame(levels(infile$chr))
        colnames(mfile) <- c("chr", "exon", "start", "end", "strand", 
            "ENSG", "exon#", "Gene")
        chr_dat <- matrix(NA,nrow=1, ncol=1+ncol(infile))
        colnames(chr_dat) <- c("gene", colnames(infile))
         colnames(chr_dat)[5]<-"Gene_id"
        for (i in list[, 1]) {
            chri_dat <- chrmatch(inf = subset(infile, infile$chr == 
                i), mf = subset(mfile, mfile$chr == i))
                #           if (length(dim(chri_dat)) == 0) {
                #chri_dat <- matrix(NA, nrow = 1, ncol = 1 + ncol(infile))
                #colnames(chri_dat) <- c("gene", colnames(infile))
                #}
            chr_dat <- rbind(chr_dat, chri_dat)
        }
        chr_dat <- chr_dat[-1,]
        return(chr_dat)
    }
    
    annotat_DEX <- function(dat, annot, colunmset) {
        colnmatch <- function(dat, annot, matchcolnm) {
            chrgroup <- function(dat=dat, annot=annot, ncolmatch,colnm = matchcolnm) {
                annot <- as.data.frame(annot)
                
                geneid <- annot$ENSG
                nexon <- annot$nexon
                exonid <- paste("E", nexon, sep = "")
                gexon1 <- paste(geneid, exonid, sep = ":")
                gid <- dat[, colnm]
                
                clname <- colnames(dat)
                names <- c("isoform", "exon", "Isoform", "Exon", 
                  "EXON", "featureID")
                indx <- is.element(clname, names)
                jj <- which(indx == TRUE)
                isof <- dat[, jj]
                isofid <- paste("E", isof, sep = "")
                gexon2 <- paste(gid, isofid, sep = ":")

                dat <- cbind(as.data.frame(gexon2), as.data.frame(dat))
                annota <- cbind(as.data.frame(gexon1), as.data.frame(annot))
                annota <- as.data.frame(annota)

                gexon1a <- unique(gexon1)
                list <- levels(factor(gexon1))
                dupexon <- function(annot) {
                  if (length(annot$gexon1) == 1) {
                    res <- annot
                  }
                  else {
                    res <- annot[1, ]
                  }
                }
                annotb <- matrix(NA, 1, ncol = ncol(annota))
                colnames(annotb) <- colnames(annota)
               
                for (i in list) {
                  annotb1 <- dupexon(annot = subset(annota, annota$gexon1 == 
                    i))

                  annotb <- rbind(annotb, annotb1)
                }
              
                annotb <- annotb[-1, ]
                gexon3 <- annotb$gexon1
                n2 <- length(gexon2)
                n3 <- length(gexon3)
                common <- intersect(gexon3, gexon2)
                ind3 <- is.element(gexon3, common)
                annotc <- annotb[ind3, ]
          
                ind2 <- is.element(gexon2, common)
                data <- dat[ind2, ]
                annotc <- annotc[order(annotc$gexon1), ]
                data <- data[order(data$gexon2), ]
                newdat <- cbind(annotc, data)
                newdat<-newdat[,-c(10:(9+ncolmatch))]
                nc<-ncol(newdat)
                colnames(newdat)<-seq(nc)
                #               print(head(newdat))
                return(newdat)
            }
            
            colnames(annot) <- c("chr", "element", "start", "end", "strand",
            "ENSG", "nexon", "Gene")
            annotexon <- subset(annot, annot$element == "exon")
            # print(dim(annotexon))
            list <- levels(annotexon$chr)
            cn1 <- ncol(dat)
            cn2 <- ncol(annotexon)
            cn3<-max(colunmset)
            res <- matrix(NA, nrow = 1, ncol = c(2 + cn1 + cn2-cn3))
            #     print(dim(res))
             colnames(res) <- seq(2 + cn1 + cn2-max(colunmset))
            for (i in list) {
                res1 <- chrgroup(dat, annot = subset(annotexon, 
                  annotexon$chr == i),ncolmatch=cn3)
                res <- rbind(res, res1)
            }
            res <- res[-1, ]
            return(res)
        }
        
        res <- colnmatch(dat = dat, annot = annot, matchcolnm = colunmset[1])
        for (i in colunmset[-1]) {
            res1 <- colnmatch(dat, annot, matchcolnm = i)
            res <- rbind(res, res1)
        }
        cn<-max(colunmset)
        colnames(annot) <- c("chr", "element", "start", "end", "strand",
        "ENSG", "nexon", "Gene")
        colnames(res) <- c("ID",colnames(annot),"x1",colnames(dat[,-c(1:cn)]))
        #     print(head(res))
        res<-res[,-10]
        return(res)
    }
    
    if (type == "gene") {
        res <- annotat_ENSG(infile = infile, mfile = mfile)
    }
    else if (type == "isoform" | type == "isof") {
        res <- annotat_isof(infile = infile, mfile = mfile)
    }
    else {
        res <- annotat_DEX(dat = infile, annot = mfile, colunmset = colunmset)
    }
    return(res)
  }
