

library(CaSpER)
## "tcga_gbm" contains the following objects: 
## data: normalized gene expression matrix
## loh.name.mapping: data.frame for mapping loh files to expression files
## annotation
## loh
## control.sample.ids: samples that are used as normal
## samps: sample information
## genoMat: genotyping large scale CNV event summary 1: amplification, -1:deletion, 0: neutral
#Download rdata from: https://www.dropbox.com/s/thx3jt589kps886/tcga_gbm.rda?dl=0
load("tcga_gbm.rda")
## "hg19_cytoband.rda" contains the following objects: 
## cytoband: hg19 cytoband information
## centromere: hg19 centromere information
data("hg19_cytoband")

data <- tcga_gbm$data
loh <-  tcga_gbm$loh
loh.name.mapping <-  tcga_gbm$loh.name.mapping
control.sample.ids <-  tcga_gbm$control.sample.ids
genoMat <-  tcga_gbm$genoMat

## generate annotation data.frame
annotation <- generateAnnotation(id_type="hgnc_symbol", genes=rownames(data), ishg19=T, centromere)
data <- data[match( annotation$Gene,rownames(data)), ]

## create CaSpER object 

object <- CreateCasperObject(raw.data=data, loh.name.mapping=loh.name.mapping, 
    sequencing.type="bulk", cnv.scale=3, loh.scale=3,method="iterative",
              annotation=annotation, loh=loh, matrix.type="normalized", expr.cutoff=1,
              control.sample.ids=control.sample.ids, cytoband=cytoband)

## runCaSpER
## this might take some time
final.objects <- runCaSpER(object, removeCentromere=T, cytoband=cytoband, method="iterative")

## download https://www.dropbox.com/s/348a9qjg3ddyajj/tcga.gbm.final.objects.updated.rda?dl=0             
load("tcga.gbm.final.objects.updated.rda")
## plot median filtered gene expression matrix
obj <- final.objects[[9]]
plotHeatmap(object, fileName="heatmap.png", cnv.scale= 3, cluster_cols = F, cluster_rows = T, show_rownames = T, only_soi = T)

## large scale event summary 
finalChrMat <- extractLargeScaleEvents (final.objects, thr=0.75) 
common <- intersect(rownames(finalChrMat), rownames(genoMat))
finalChrMat <- finalChrMat[match(common, rownames(finalChrMat)), ]
genoMat <- genoMat[match(common, rownames(genoMat)), ]
## calculate TPR and FPR using genotyping array as gold standard
calcROC(chrMat=finalChrMat, chrMat2=genoMat)

## segment based summary    
gamma <- 7
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]

## gene based summary 
all.summary<- rbind(loss.final, gain.final)
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), 
    IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "GeneSymbol")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,2])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)

## genotyping array gene based summary
## download https://www.dropbox.com/s/68rn36gpx6plew7/gbm.gt.matrix.rda?dl=0             
load("gbm.gt.matrix.rda")
                    
## calculate TPR and FPR using genotyping array as gold standard
calcROC(chrMat=rna.matrix, chrMat2=gt.matrix)

#### VISUALIZATION 

## plot large scale events 
plotLargeScaleEvent (object=obj, fileName="large.scale.events.pdf") 
## plot large scale events using event summary matrix 1: amplification, -1:deletion, 0: neutral
plotLargeScaleEvent2 (finalChrMat, fileName="large.scale.events.summarized.pdf") 
## plot BAF deviation for each sample in seperate pages 
plotBAFInSeperatePages (loh =obj@loh.median.filtered.data, folderName="LOHPlotsSeperate") 
## plot gene expression and BAF signal for one sample in one plot
plotGEAndBAFOneSample (object=obj, cnv.scale=3, loh.scale=3, sample= "TCGA-02-0047-01A")
