## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=T, message=FALSE, warning=FALSE--------------------------------
library(CaSpER)
data (yale_meningioma)
yale_meningioma$data[1:5, 1:5]


## ---- eval=T-------------------------------------------------------------
data(hg19_cytoband)
cytoband[1:5, ]

## ---- eval=T-------------------------------------------------------------
data(hg19_cytoband)
centromere[1:5, ]

## ---- eval=F-------------------------------------------------------------
#  annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(yale_meningioma$data), ishg19=T, centromere)
#  annotation[1:5, ]

## ---- eval=F-------------------------------------------------------------
#  loh <- readBAFExtractOutput ( path="./meningioma_baf\\", sequencing.type="bulk")
#  names(loh) <- gsub(".snp", "", names(loh))

## ---- eval=F-------------------------------------------------------------
#  
#  object <- CreateCasperObject(raw.data=data,loh.name.mapping=loh.name.mapping, sequencing.type="bulk",
#    cnv.scale=3, loh.scale=3,
#    annotation=annotation, method="iterative", loh=loh,
#    control.sample.ids=control.sample.ids, cytoband=cytoband)

## ---- eval=T-------------------------------------------------------------
yale_meningioma$loh.name.mapping [1:5, ]
yale_meningioma$control.sample.ids

## ---- eval=F-------------------------------------------------------------
#  final.objects <- runCaSpER(object, removeCentromere=T, cytoband=cytoband, method="iterative")

## ----include=FALSE-------------------------------------------------------
load("C:\\Users\\aharmanci\\Google Drive\\uthealth\\SCell_ProcessedData\\TCGA\\final.objects.yale.mn.rda")

## ---- eval=F-------------------------------------------------------------
#  finalChrMat <- extractLargeScaleEvents (final.objects, thr=0.75)

## ---- eval=F-------------------------------------------------------------
#  gamma <- 6
#  all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
#  segment.summary <- extractSegmentSummary (final.objects)
#  loss <- segment.summary$all.summary.loss
#  gain <- segment.summary$all.summary.gain
#  loss.final <- loss[loss$count>=gamma, ]
#  gain.final <- gain[gain$count>=gamma, ]

## ---- eval=F-------------------------------------------------------------
#  all.summary<- rbind(loss.final, gain.final)
#  colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
#  rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))),
#      IRanges(all.summary$Start, all.summary$End))
#  ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
#  hits <- findOverlaps(geno.rna, ann.gr)
#  genes <- splitByOverlap(ann.gr, geno.rna, "GeneSymbol")
#  genes.ann <- lapply(genes, function(x) x[!(x=="")])
#  all.genes <- unique(final.objects[[1]]@annotation.filt[,2])
#  all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
#  rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)

## ---- fig.show='hold'----------------------------------------------------
obj <- final.objects[[9]]
plotHeatmap(object=obj, fileName="heatmap.png",cnv.scale= 3, cluster_cols = F, cluster_rows = T, show_rownames = T, only_soi = T)

## ---- fig.show='hold'----------------------------------------------------
plotLargeScaleEvent (object=obj, fileName="large.scale.events.pdf") 

## ----pressure, echo=FALSE, fig.cap="A caption", out.width = '100%'-------
knitr::include_graphics("temp.png")

