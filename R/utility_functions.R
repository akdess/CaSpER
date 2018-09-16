
calcROC <- function(chrMat, chrMat2) {
    chrMatPos <- chrMat
    chrMatPos[chrMatPos < 0] <- 0
    
    chrMatNeg <- chrMat
    chrMatNeg[chrMatNeg > 0] <- 0
    
    chrMatPos2 <- chrMat2
    chrMatPos2[chrMatPos2 < 0] <- 0
    
    chrMatNeg2 <- chrMat2
    chrMatNeg2[chrMatNeg2 > 0] <- 0
    
    tp <- sum(chrMat * chrMat2)
    all_pos <- sum(apply(chrMat2, 2, function(x) length(which(!x == 0))))
    recall <- tp/all_pos
    
    all_pred_pos <- sum(apply(chrMat, 2, function(x) length(which(!x == 0))))
    fp <- all_pred_pos - tp
    if (fp < 0) 
        fp <- 0
    
    all_neg <- sum(apply(chrMat2, 2, function(x) length(which(x == 0))))
    fallout <- fp/all_neg
    
    # positive
    tp <- sum(chrMatPos * chrMatPos2)
    all_pos <- sum(apply(chrMatPos2, 2, function(x) length(which(!x == 0))))
    recallAmp <- tp/all_pos
    
    all_pred_pos <- sum(apply(chrMatPos, 2, function(x) length(which(!x == 0))))
    fp <- all_pred_pos - tp
    if (fp < 0) 
        fp <- 0
    
    all_neg <- sum(apply(chrMatPos2, 2, function(x) length(which(x == 0))))
    falloutAmp <- fp/all_neg
    
    tp <- sum(chrMatNeg * chrMatNeg2)
    all_pos <- sum(apply(chrMatNeg2, 2, function(x) length(which(!x == 0))))
    recallDel <- tp/all_pos
    
    all_pred_pos <- sum(apply(chrMatNeg, 2, function(x) length(which(!x == 0))))
    fp <- all_pred_pos - tp
    if (fp < 0) 
        fp <- 0
    
    all_neg <- sum(apply(chrMatNeg2, 2, function(x) length(which(x == 0))))
    falloutDel <- fp/all_neg
    
    list(tpr = recall, fpr = fallout, tprAmp = recallAmp, fprAmp = falloutAmp, tprDel = recallDel, fprDel = falloutDel)
}

plotROC <- function(roc, threshold, cost_of_fp, cost_of_fn) {
    
    norm_vec <- function(v) (v - min(v))/diff(range(v))
    
    idx_threshold = which.min(abs(roc$threshold - threshold))
    
    col_ramp <- colorRampPalette(c("green", "orange", "red", "black"))(100)
    col_by_cost <- col_ramp[ceiling(norm_vec(roc$cost) * 99) + 1]
    p_roc <- ggplot(roc, aes(fpr, tpr)) + geom_line(color = rgb(0, 0, 1, alpha = 0.3)) + coord_fixed() + geom_line(aes(threshold, threshold), 
        color = rgb(0, 0, 1, alpha = 0.5)) + labs(title = sprintf("ROC")) + xlab("FPR") + ylab("TPR") + geom_hline(yintercept = roc[idx_threshold, 
        "tpr"], alpha = 0.5, linetype = "dashed") + geom_vline(xintercept = roc[idx_threshold, "fpr"], alpha = 0.5, linetype = "dashed")
}

extractMUAndCooccurence <- function(finalChrMat, loh, loh.name.mapping)
{
    results <- list()
    chrs<- colnames(finalChrMat)

    for (j in 1:length(names(loh)))
    {
        samples <- as.character(loh.name.mapping[loh.name.mapping$loh.name %in% names(loh)[j], 2])
        finalChrMat.sub<- finalChrMat[rownames(finalChrMat) %in% samples,]

        combin <- expand.grid( paste0(chrs, "del"), paste0(chrs, "del"))
        list.names <- apply(combin, 1, function(x) paste(x[1], x[2], sep="_vs_"))

        filter <- c(paste0(chrs, "del_vs_", chrs, "del"), paste0( chrs[seq(1, 44, by=2)], "del_vs_",chrs[seq(2, 44, by=2)], "del" ), paste0( chrs[seq(2, 44, by=2)], "del_vs_",chrs[seq(1, 44, by=2)], "del" ))
        combin <- combin[!(list.names %in% filter), ]
        list.names <- list.names[!(list.names %in% filter)]

        occurenceDel <- matrix("NA", nrow=nrow(combin), ncol=4)
        rownames(occurenceDel) <- list.names
        colnames(occurenceDel) <- c("bothDelOccurencePval","bothDelOccurenceFDR", 
                              "bothDelMUPval","bothDelMUFDR")
        finalChrMat.sub[finalChrMat.sub >0] <- 0
        finalChrMat.sub[finalChrMat.sub <0] <- 1
        occurenceDel <- data.frame(occurenceDel, stringsAsFactors = F)

        for (i in 1:nrow(combin)){
          a1 <- finalChrMat.sub[, combin[i,1]]
          a2 <- finalChrMat.sub[, combin[i,2]]
          if(dim(table(a1, a2))[1]>1 & dim(table(a1, a2))[2]>1){
            occurenceDel[i, 1] <- fisher.test(table(a1, a2) , alternative="greater")$p.val
            occurenceDel[i, 3] <- fisher.test(table(a1, a2) , alternative="less")$p.val

          }
        }
        occurenceDel[, 2] <-  p.adjust(as.numeric(occurenceDel[, 1]), method="fdr")
        occurenceDel[, 4] <- p.adjust(as.numeric(occurenceDel[, 3]), method="fdr")

        combin <- expand.grid( paste0(chrs, "amp"), paste0(chrs, "amp"))
        list.names <- apply(combin, 1, function(x) paste(x[1], x[2], sep="_vs_"))

        filter <- c(paste0(chrs, "amp_vs_", chrs, "amp"), paste0( chrs[seq(1, 44, by=2)], "amp_vs_", chrs[seq(2, 44, by=2)], "amp"), 
                paste0( chrs[seq(2, 44, by=2)], "amp_vs_",chrs[seq(1, 44, by=2)], "amp" ))
        combin <- combin[!(list.names %in% filter), ]
        list.names <- list.names[!(list.names %in% filter)]

        occurenceAmp <- matrix("NA", nrow=nrow(combin), ncol=4)
        rownames(occurenceAmp) <- list.names
        colnames(occurenceAmp) <- c("bothAmpOccurencePval","bothAmpOccurenceFDR", 
                              "bothAmpMUPval",         "bothAmpMUFDR")

        finalChrMat.sub.2<- finalChrMat[grep(names(loh)[j], rownames(finalChrMat)),]
        finalChrMat.sub.2[finalChrMat.sub.2 <0] <- 0

        for (i in 1:nrow(combin)){
          a1 <- finalChrMat.sub.2[, combin[i,1]]
          a2 <- finalChrMat.sub.2[, combin[i,2]]
          if(dim(table(a1, a2))[1]>1 & dim(table(a1, a2))[2]>1){
            occurenceAmp[i, 1] <- fisher.test(table(a1, a2) , alternative="greater")$p.val
            occurenceAmp[i, 3] <- fisher.test(table(a1, a2) , alternative="less")$p.val
          }
        }
        occurenceAmp[, 2] <-  p.adjust(as.numeric(occurenceAmp[, 1]), method="fdr")
        occurenceAmp[, 4] <- p.adjust(as.numeric(occurenceAmp[, 3]), method="fdr")

        combin <- expand.grid( paste0(chrs, "del"), paste0(chrs, "amp"))
        list.names <- apply(combin, 1, function(x) paste(x[1], x[2], sep="_vs_"))

        filter <- c(paste0(chrs, "del_vs_", chrs, "amp"), 
            paste0( chrs[seq(1, 44, by=2)], "del_vs_",chrs[seq(2, 44, by=2)], "amp" ),
             paste0( chrs[seq(2, 44, by=2)], "del_vs_",chrs[seq(1, 44, by=2)] , "amp"))
        combin <- combin[!(list.names %in% filter), ]
        list.names <- list.names[!(list.names %in% filter)]

        occurenceAmpDel <- matrix("NA", nrow=nrow(combin), ncol=4)
        rownames(occurenceAmpDel) <- list.names
        colnames(occurenceAmpDel) <- c("DelAmpOccurencePval",   "DelAmpOccurenceFDR", 
                              "DelAmpMUPval","DelAmpMUFDR")

        for (i in 1:nrow(combin)){
          a1 <- finalChrMat.sub[, combin[i,1]]
          a2 <- finalChrMat.sub.2[, combin[i,2]]
           if(dim(table(a1, a2))[1]>1 & dim(table(a1, a2))[2]>1){
            occurenceAmpDel[i, 1] <- fisher.test(table(a1, a2) , alternative="greater")$p.val
            occurenceAmpDel[i, 3] <- fisher.test(table(a1, a2) , alternative="less")$p.val
          }
        }
        occurenceAmpDel[, 2] <-  p.adjust(as.numeric(occurenceAmpDel[, 1]), method="fdr")
        occurenceAmpDel[, 4] <- p.adjust(as.numeric(occurenceAmpDel[, 3]), method="fdr")    
            
        combin <- expand.grid( paste0(chrs, "amp"), paste0(chrs, "del"))
        list.names <- apply(combin, 1, function(x) paste(x[1], x[2], sep="_vs_"))

        filter <- c(paste0(chrs, "amp_vs_", chrs, "del"), 
            paste0( chrs[seq(1, 44, by=2)], "amp_vs_",chrs[seq(2, 44, by=2)], "del" ),
             paste0( chrs[seq(2, 44, by=2)], "amp_vs_",chrs[seq(1, 44, by=2)] , "del"))
        combin <- combin[!(list.names %in% filter), ]
        list.names <- list.names[!(list.names %in% filter)]


        occurenceDelAmp <- matrix("NA", nrow=nrow(combin), ncol=4)
        rownames(occurenceDelAmp) <- list.names
        colnames(occurenceDelAmp) <- c("DelAmpOccurencePval",   "DelAmpOccurenceFDR", 
                              "DelAmpMUPval","DelAmpMUFDR")

        for (i in 1:nrow(combin)){
          a1 <- finalChrMat.sub.2[, combin[i,1]]
          a2 <- finalChrMat.sub[, combin[i,2]]
           if(dim(table(a1, a2))[1]>1 & dim(table(a1, a2))[2]>1){
            occurenceDelAmp[i, 1] <- fisher.test(table(a1, a2) , alternative="greater")$p.val
            occurenceDelAmp[i, 3] <- fisher.test(table(a1, a2) , alternative="less")$p.val
          }
        }
        occurenceDelAmp[, 2] <-  p.adjust(as.numeric(occurenceDelAmp[, 1]), method="fdr")
        occurenceDelAmp[, 4] <- p.adjust(as.numeric(occurenceDelAmp[, 3]), method="fdr")    

        mat <- rbind(data.frame(Node1=unlist(lapply(strsplit(rownames(occurenceDel), split="_vs_"), function(x) x[1])), 
            Node2=unlist(lapply(strsplit(rownames(occurenceDel), split="_vs_"), function(x) x[2])),
            Pval=occurenceDel[,1], FDR=occurenceDel[,2], type="occurence"),

            data.frame(Node1=unlist(lapply(strsplit(rownames(occurenceDel), split="_vs_"), function(x) x[1])), 
            Node2=unlist(lapply(strsplit(rownames(occurenceDel), split="_vs_"), function(x) x[2])),
            Pval=occurenceDel[,3], FDR=occurenceDel[,4], type="MU"),

            data.frame(Node1=unlist(lapply(strsplit(rownames(occurenceAmp), split="_vs_"), function(x) x[1])), 
            Node2=unlist(lapply(strsplit(rownames(occurenceAmp), split="_vs_"), function(x) x[2])),
            Pval=occurenceAmp[,1], FDR=occurenceAmp[,2], type="occurence"),

            data.frame(Node1=unlist(lapply(strsplit(rownames(occurenceAmp), split="_vs_"), function(x) x[1])), 
            Node2=unlist(lapply(strsplit(rownames(occurenceAmp), split="_vs_"), function(x) x[2])),
            Pval=occurenceAmp[,3], FDR=occurenceAmp[,4], type="MU"),

            data.frame(Node1=unlist(lapply(strsplit(rownames(occurenceAmpDel), split="_vs_"), function(x) x[1])), 
            Node2=unlist(lapply(strsplit(rownames(occurenceAmpDel), split="_vs_"), function(x) x[2])),
            Pval=occurenceAmpDel[,1], FDR=occurenceAmpDel[,2], type="occurence"),

            data.frame(Node1=unlist(lapply(strsplit(rownames(occurenceAmpDel), split="_vs_"), function(x) x[1])), 
            Node2=unlist(lapply(strsplit(rownames(occurenceAmpDel), split="_vs_"), function(x) x[2])),
            Pval=occurenceAmpDel[,3], FDR=occurenceAmpDel[,4], type="MU"), 

            data.frame(Node1=unlist(lapply(strsplit(rownames(occurenceDelAmp), split="_vs_"), function(x) x[1])), 
            Node2=unlist(lapply(strsplit(rownames(occurenceDelAmp), split="_vs_"), function(x) x[2])),
            Pval=occurenceDelAmp[,1], FDR=occurenceDelAmp[,2], type="occurence"),

            data.frame(Node1=unlist(lapply(strsplit(rownames(occurenceDelAmp), split="_vs_"), function(x) x[1])), 
            Node2=unlist(lapply(strsplit(rownames(occurenceDelAmp), split="_vs_"), function(x) x[2])),
            Pval=occurenceDelAmp[,3], FDR=occurenceDelAmp[,4], type="MU")
            )

        mat[, 3] <- as.numeric(as.character(mat[, 3]))
        mat[,1] <- as.character(mat[,1])
        mat[,2] <- as.character(mat[,2])
        results[[j]] <- mat
    }

    names(results) <- names(loh)
    return(results)
}