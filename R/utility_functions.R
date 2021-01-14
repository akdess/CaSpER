#' @title calcROC()
#'
#' @description Calculates tpr and fpr values using genotyping array as gold standard
#'
#' @param chrMat large scale event matrix generated using CaSpER
#' 
#' @param chrMat2 large scale event matrix generated using genotyping array
#'
#' @return accuracy measures 
#'
#' @export
#'
#' 
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

#' @title extractMUAndCooccurence()
#'
#' @description calculates significant mutually exclusive and co-occurent events
#'
#' @param finalChrMat large scale event matrix generated using CaSpER
#' 
#' @param loh  original baf signal 
#'
#' @param loh.name.mapping  contains the cell (sample) name and the matching baf signal sample name
#'
#' @return list of mutually exclusive and co-occurent events
#'
#' @export
#'
#' 
extractMUAndCooccurence <- function(finalChrMat, loh, loh.name.mapping) {
    results <- list()
    chrs <- colnames(finalChrMat)
    
    for (j in 1:length(names(loh))) {
        samples <- as.character(loh.name.mapping[loh.name.mapping$loh.name %in% names(loh)[j], 2])
        finalChrMat.sub <- finalChrMat[rownames(finalChrMat) %in% samples, ]
        
        combin <- expand.grid(paste0(chrs, "del"), paste0(chrs, "del"))
        list.names <- apply(combin, 1, function(x) paste(x[1], x[2], sep = "_vs_"))
        
        filter <- c(paste0(chrs, "del_vs_", chrs, "del"), paste0(chrs[seq(1, 44, by = 2)], "del_vs_", chrs[seq(2, 44, by = 2)], 
            "del"), paste0(chrs[seq(2, 44, by = 2)], "del_vs_", chrs[seq(1, 44, by = 2)], "del"))
        combin <- combin[!(list.names %in% filter), ]
        list.names <- list.names[!(list.names %in% filter)]
        
        occurenceDel <- matrix("NA", nrow = nrow(combin), ncol = 4)
        rownames(occurenceDel) <- list.names
        colnames(occurenceDel) <- c("bothDelOccurencePval", "bothDelOccurenceFDR", "bothDelMUPval", "bothDelMUFDR")
        finalChrMat.sub[finalChrMat.sub > 0] <- 0
        finalChrMat.sub[finalChrMat.sub < 0] <- 1
        occurenceDel <- data.frame(occurenceDel, stringsAsFactors = F)
        
        for (i in 1:nrow(combin)) {
            a1 <- finalChrMat.sub[, combin[i, 1]]
            a2 <- finalChrMat.sub[, combin[i, 2]]
            if (dim(table(a1, a2))[1] > 1 & dim(table(a1, a2))[2] > 1) {
                occurenceDel[i, 1] <- fisher.test(table(a1, a2), alternative = "greater")$p.val
                occurenceDel[i, 3] <- fisher.test(table(a1, a2), alternative = "less")$p.val
                
            }
        }
        occurenceDel[, 2] <- p.adjust(as.numeric(occurenceDel[, 1]), method = "fdr")
        occurenceDel[, 4] <- p.adjust(as.numeric(occurenceDel[, 3]), method = "fdr")
        
        combin <- expand.grid(paste0(chrs, "amp"), paste0(chrs, "amp"))
        list.names <- apply(combin, 1, function(x) paste(x[1], x[2], sep = "_vs_"))
        
        filter <- c(paste0(chrs, "amp_vs_", chrs, "amp"), paste0(chrs[seq(1, 44, by = 2)], "amp_vs_", chrs[seq(2, 44, by = 2)], 
            "amp"), paste0(chrs[seq(2, 44, by = 2)], "amp_vs_", chrs[seq(1, 44, by = 2)], "amp"))
        combin <- combin[!(list.names %in% filter), ]
        list.names <- list.names[!(list.names %in% filter)]
        
        occurenceAmp <- matrix("NA", nrow = nrow(combin), ncol = 4)
        rownames(occurenceAmp) <- list.names
        colnames(occurenceAmp) <- c("bothAmpOccurencePval", "bothAmpOccurenceFDR", "bothAmpMUPval", "bothAmpMUFDR")
        
        finalChrMat.sub.2 <- finalChrMat[grep(names(loh)[j], rownames(finalChrMat)), ]
        finalChrMat.sub.2[finalChrMat.sub.2 < 0] <- 0
        
        for (i in 1:nrow(combin)) {
            a1 <- finalChrMat.sub.2[, combin[i, 1]]
            a2 <- finalChrMat.sub.2[, combin[i, 2]]
            if (dim(table(a1, a2))[1] > 1 & dim(table(a1, a2))[2] > 1) {
                occurenceAmp[i, 1] <- fisher.test(table(a1, a2), alternative = "greater")$p.val
                occurenceAmp[i, 3] <- fisher.test(table(a1, a2), alternative = "less")$p.val
            }
        }
        occurenceAmp[, 2] <- p.adjust(as.numeric(occurenceAmp[, 1]), method = "fdr")
        occurenceAmp[, 4] <- p.adjust(as.numeric(occurenceAmp[, 3]), method = "fdr")
        
        combin <- expand.grid(paste0(chrs, "del"), paste0(chrs, "amp"))
        list.names <- apply(combin, 1, function(x) paste(x[1], x[2], sep = "_vs_"))
        
        filter <- c(paste0(chrs, "del_vs_", chrs, "amp"), paste0(chrs[seq(1, 44, by = 2)], "del_vs_", chrs[seq(2, 44, by = 2)], 
            "amp"), paste0(chrs[seq(2, 44, by = 2)], "del_vs_", chrs[seq(1, 44, by = 2)], "amp"))
        combin <- combin[!(list.names %in% filter), ]
        list.names <- list.names[!(list.names %in% filter)]
        
        occurenceAmpDel <- matrix("NA", nrow = nrow(combin), ncol = 4)
        rownames(occurenceAmpDel) <- list.names
        colnames(occurenceAmpDel) <- c("DelAmpOccurencePval", "DelAmpOccurenceFDR", "DelAmpMUPval", "DelAmpMUFDR")
        
        for (i in 1:nrow(combin)) {
            a1 <- finalChrMat.sub[, combin[i, 1]]
            a2 <- finalChrMat.sub.2[, combin[i, 2]]
            if (dim(table(a1, a2))[1] > 1 & dim(table(a1, a2))[2] > 1) {
                occurenceAmpDel[i, 1] <- fisher.test(table(a1, a2), alternative = "greater")$p.val
                occurenceAmpDel[i, 3] <- fisher.test(table(a1, a2), alternative = "less")$p.val
            }
        }
        occurenceAmpDel[, 2] <- p.adjust(as.numeric(occurenceAmpDel[, 1]), method = "fdr")
        occurenceAmpDel[, 4] <- p.adjust(as.numeric(occurenceAmpDel[, 3]), method = "fdr")
        
        combin <- expand.grid(paste0(chrs, "amp"), paste0(chrs, "del"))
        list.names <- apply(combin, 1, function(x) paste(x[1], x[2], sep = "_vs_"))
        
        filter <- c(paste0(chrs, "amp_vs_", chrs, "del"), paste0(chrs[seq(1, 44, by = 2)], "amp_vs_", chrs[seq(2, 44, by = 2)], 
            "del"), paste0(chrs[seq(2, 44, by = 2)], "amp_vs_", chrs[seq(1, 44, by = 2)], "del"))
        combin <- combin[!(list.names %in% filter), ]
        list.names <- list.names[!(list.names %in% filter)]
        
        
        occurenceDelAmp <- matrix("NA", nrow = nrow(combin), ncol = 4)
        rownames(occurenceDelAmp) <- list.names
        colnames(occurenceDelAmp) <- c("DelAmpOccurencePval", "DelAmpOccurenceFDR", "DelAmpMUPval", "DelAmpMUFDR")
        
        for (i in 1:nrow(combin)) {
            a1 <- finalChrMat.sub.2[, combin[i, 1]]
            a2 <- finalChrMat.sub[, combin[i, 2]]
            if (dim(table(a1, a2))[1] > 1 & dim(table(a1, a2))[2] > 1) {
                occurenceDelAmp[i, 1] <- fisher.test(table(a1, a2), alternative = "greater")$p.val
                occurenceDelAmp[i, 3] <- fisher.test(table(a1, a2), alternative = "less")$p.val
            }
        }
        occurenceDelAmp[, 2] <- p.adjust(as.numeric(occurenceDelAmp[, 1]), method = "fdr")
        occurenceDelAmp[, 4] <- p.adjust(as.numeric(occurenceDelAmp[, 3]), method = "fdr")
        
        mat <- rbind(data.frame(Node1 = unlist(lapply(strsplit(rownames(occurenceDel), split = "_vs_"), function(x) x[1])), 
            Node2 = unlist(lapply(strsplit(rownames(occurenceDel), split = "_vs_"), function(x) x[2])), Pval = occurenceDel[, 
                1], FDR = occurenceDel[, 2], type = "occurence"), data.frame(Node1 = unlist(lapply(strsplit(rownames(occurenceDel), 
            split = "_vs_"), function(x) x[1])), Node2 = unlist(lapply(strsplit(rownames(occurenceDel), split = "_vs_"), function(x) x[2])), 
            Pval = occurenceDel[, 3], FDR = occurenceDel[, 4], type = "MU"), data.frame(Node1 = unlist(lapply(strsplit(rownames(occurenceAmp), 
            split = "_vs_"), function(x) x[1])), Node2 = unlist(lapply(strsplit(rownames(occurenceAmp), split = "_vs_"), function(x) x[2])), 
            Pval = occurenceAmp[, 1], FDR = occurenceAmp[, 2], type = "occurence"), data.frame(Node1 = unlist(lapply(strsplit(rownames(occurenceAmp), 
            split = "_vs_"), function(x) x[1])), Node2 = unlist(lapply(strsplit(rownames(occurenceAmp), split = "_vs_"), function(x) x[2])), 
            Pval = occurenceAmp[, 3], FDR = occurenceAmp[, 4], type = "MU"), data.frame(Node1 = unlist(lapply(strsplit(rownames(occurenceAmpDel), 
            split = "_vs_"), function(x) x[1])), Node2 = unlist(lapply(strsplit(rownames(occurenceAmpDel), split = "_vs_"), 
            function(x) x[2])), Pval = occurenceAmpDel[, 1], FDR = occurenceAmpDel[, 2], type = "occurence"), data.frame(Node1 = unlist(lapply(strsplit(rownames(occurenceAmpDel), 
            split = "_vs_"), function(x) x[1])), Node2 = unlist(lapply(strsplit(rownames(occurenceAmpDel), split = "_vs_"), 
            function(x) x[2])), Pval = occurenceAmpDel[, 3], FDR = occurenceAmpDel[, 4], type = "MU"), data.frame(Node1 = unlist(lapply(strsplit(rownames(occurenceDelAmp), 
            split = "_vs_"), function(x) x[1])), Node2 = unlist(lapply(strsplit(rownames(occurenceDelAmp), split = "_vs_"), 
            function(x) x[2])), Pval = occurenceDelAmp[, 1], FDR = occurenceDelAmp[, 2], type = "occurence"), data.frame(Node1 = unlist(lapply(strsplit(rownames(occurenceDelAmp), 
            split = "_vs_"), function(x) x[1])), Node2 = unlist(lapply(strsplit(rownames(occurenceDelAmp), split = "_vs_"), 
            function(x) x[2])), Pval = occurenceDelAmp[, 3], FDR = occurenceDelAmp[, 4], type = "MU"))
        
        mat[, 3] <- as.numeric(as.character(mat[, 3]))
        mat[, 1] <- as.character(mat[, 1])
        mat[, 2] <- as.character(mat[, 2])
        results[[j]] <- mat
    }
    
    names(results) <- names(loh)
    return(results)
}

#' @title generateAnnotation()
#'
#' @description retrieves gene chromosomal locations from biomart
#'
#' @param id_type gene list identifier, ensembl_gene_id or hgnc_symbol
#' 
#' @param genes  list of genes
#'
#' @param ishg19  boolean values determining the genome version 
#' 
#' @param centromere centromer regions
#' 
#' @return list of mutually exclusive and co-occurent events
#'
#' @export
#'
#'
generateAnnotation <- function(id_type = "ensembl_gene_id", genes, ishg19, centromere, host="www.ensembl.org") {

    if (ishg19) {
        mart <- useDataset("hsapiens_gene_ensembl", useEnsembl(biomart = "ensembl", GRCh = 37, host = host))
    } else {
        mart <- useDataset("hsapiens_gene_ensembl", useEnsembl(biomart = "ensembl", host = host))
    }
    G_list <- getBM(filters = id_type, attributes = c(id_type, 
        "hgnc_symbol", "chromosome_name", "start_position", 
        "end_position", "band"), values = genes, 
        mart = mart, useCache=F)
    common <- intersect(genes, G_list[, id_type])
    ord <- match(common, G_list[, id_type])
    annotation <- G_list[ord, ]
    annotation <- annotation[order(annotation$start_position), ]
    annotation$cytoband <- paste0(annotation$chromosome_name, substr((annotation$band), 0, 1))
    
    chr.list <- as.character(1:22, "X")
    idx <- unlist(unique(as.vector(sapply(chr.list, function(x) as.vector(unlist(which(as.character(annotation$chromosome_name) == 
        x)))))))
    annotation <- as.data.frame(annotation[idx, ])
    
    colnames(annotation)[c(1:5, 7)] <- c("Gene", "GeneSymbol",  "Chr", "start", "end", "cytoband")
    
    annotation$isCentromer <- rep("no", nrow(annotation))
    
    centromere_snps <- NULL
    for (k in 1:(dim(centromere)[1])) {
        annotation$isCentromer[which(as.character(annotation$Chr) == gsub("chr", "", as.character(centromere$V1[k])) & (as.numeric(as.character(annotation$Position)) >= 
            centromere$V2[k] & as.numeric(as.character(annotation$Position)) <= centromere$V3[k]))] <- "yes"
    }
    
    annotation$Position <- (as.numeric(annotation$start) + as.numeric(annotation$end))/2
    
    annotation$new_positions <- as.vector(unlist(lapply(lapply(split(annotation$cytoband, annotation$cytoband), length)[unique(annotation$cytoband)], 
        function(x) 1:x)))
    return(annotation)
    
}

#' @title goEnrichmentBP()
#'
#' @description GO Term enrichment 
#'
#' @param genes list of genes
#' 
#' @param ontology  ontology  (BP, CC or MF) 
#'
#' @param universe  universe of genes
#' 
#' @param pvalue pvalue cutoff
#' 
#' @param annotation ontology annotation default:org.Hs.eg.db
#' 
#' @return significantly enriched GO Terms
#'
#' @export
#'
#'
goEnrichmentBP<-function (genes, ontology, universe=character(0),  
        pvalue=0.05, annotation='org.Hs.eg.db', conditionalSearch=TRUE, genes2)
{
    
    params = new ("GOHyperGParams", geneIds=unique(genes),  ontology="BP",  
            annotation='org.Hs.eg.db',
            universeGeneIds=unique(universe), pvalueCutoff = pvalue,  
            conditional=TRUE,
            testDirection = "over")
    hgr <- hyperGTest (params)
    tab<-summary(hgr)
    geneIdsByCategory(hgr)
    
    geneSymbol.list<-sapply(tab$GOBPID,function(x){
                entrezIds<-geneIdsByCategory(hgr)[[x]]
                gs <-as.vector(unique(na.omit(unlist(mget (entrezIds,revmap(org.Hs.egALIAS2EG),ifnotfound=NA)))))
                    paste(intersect(genes2, gs),collapse=",")})
    
    sum<-cbind(tab,genes=as.vector(geneSymbol.list))
    sum$FDR<- p.adjust(sum$Pvalue, 'fdr')
    
    return (sum)
}

#' @title getDiffExprGenes()
#'
#' @description get differentially expressed genes between samples having selected specified CNV events
#'
#' @param final.objects list of objects
#' 
#' @param sampleName  sample name
#'
#' @param chrs  selected chromosomes
#' 
#' @param event.type cnv event type
#' 
#' @return differentially expressed genes
#'
#' @export
#'
#'
getDiffExprGenes <- function(final.objects, sampleName, chrs, event.type)
{
    finalChrMat <- extractLargeScaleEvents (final.objects, thr=0.75) 
    group1 <- names(which(finalChrMat[grep(sampleName, rownames(finalChrMat)),chrs[1]] == event.type[1]))
    group2 <- names(which(finalChrMat[grep(sampleName, rownames(finalChrMat)),chrs[2]] == event.type[2]))
    common <- intersect(group1, group2)
    group1 <- group1[!(group1 %in% common)]
    group2 <- group2[!(group2 %in% common)]
    
    data <- final.objects[[1]]@data

    caseind = which(colnames(data) %in% group1)
    controlind = which(colnames(data) %in% group2)
        

    eset = data[, c(caseind, controlind)]
        
    TS = as.factor(c(rep("T", length(caseind)), rep("C", 
                                length(controlind))))
        
    design = model.matrix(~0 + TS)
    colnames(design) = c("C", "T")
    
    fit = lmFit(eset, design)
    cont.matrix = makeContrasts(comp = T - C, levels = design)
    fit2 = contrasts.fit(fit, cont.matrix)
    contrasts.fit = eBayes(fit2)
    
    print("Calculating differential expression...")
        
    lods<-contrasts.fit$lods
    colnames(lods)<-"lods"
    all.summary<-topTable(contrasts.fit,n=Inf,adjust="BH",
            sort.by="none",coef=1)
    all.summary<-all.summary[match(rownames(lods),rownames(all.summary)),]
    all.summary<-data.frame(ID=rownames(all.summary),all.summary,lods)
    results <- merge(all.summary, final.objects[[1]]@annotation.filt, by.x="ID", by.y="Gene", all.x=T, all.y=F)
    return(results)
    
}

#' @title generateEnrichmentSummary()
#'
#' @description generate GO Term enrichment summary 
#'
#' @param results output of getDiffExprGenes() function
#' 
#' @return significantly enriched GO Terms
#'
#' @export
#'
#'
generateEnrichmentSummary <- function(results)
{
    genes <- as.character(results$ID[results$adj.P.Val<0.05])
    entrez.id <-as.vector(unique(na.omit(unlist(mget (genes,org.Hs.egALIAS2EG,ifnotfound=NA)))))
    universe.id <- as.vector(unique(na.omit(unlist(as.list(org.Hs.egALIAS2EG)))))
    go.BP<-goEnrichmentBP(genes=entrez.id , ontology="BP",universe=universe.id,  
            pvalue=0.01, annotation='org.Hs.eg.db', conditionalSearch=TRUE, genes2=genes)
    return(go.BP)
}

#' @title gene.matrix()
#'
#' @description Gene level CNV events represented as matrix where rows represent samples and columns represent samples
#'
#' @param segment CNV segments
#' 
#' @param all.genes  gene names
#'
#' @param all.samples  samp names
#' 
#' @param genes.ann gene symbols within each segments
#' 
#' @return matrix of gene level CNV events
#'
#' @export
gene.matrix <- function(segment, all.genes, all.samples, genes.ann){
  
  CN<-matrix(0,nrow=length(all.samples),ncol=length(all.genes))
  rownames(CN)<-all.samples
  colnames(CN)<-all.genes
  
  for (i in 1:length(rownames(CN))){
    
    ampSegments<-which(rownames(CN)[i]==segment$ID & segment$type=="Gain")
    lossSegments<-which(rownames(CN)[i]==segment$ID & (segment$type=="Loss"))
    
    ampGenes<-unique(unlist(sapply(ampSegments,function(x) unique(genes.ann[[x]]))))
    lossGenes<-unique(unlist(sapply(lossSegments,function(x) unique(genes.ann[[x]]))))
    
     if(length(ampGenes)>0){
      m<-match(ampGenes,colnames(CN))
      CN[i,m]<-1
    }
    if(length(lossGenes)>0){
      m<-match(lossGenes,colnames(CN))
      CN[i,m]<-(-1)
    }
    
  }
  
  CN<-t(CN)
  CN
}
