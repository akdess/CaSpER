## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include=FALSE------------------------------------------------------
library(knitr)

## ----eval=T, message=FALSE, warning=FALSE--------------------------------
library(CaSpER)
data (yale_meningioma)
kable(yale_meningioma$data[1:5, 1:5])


## ---- eval=T-------------------------------------------------------------
data(hg19_cytoband)
kable(cytoband[1:5, ])

## ---- eval=F-------------------------------------------------------------
#  annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(yale_meningioma$data), ishg19=T, centromere)

## ---- include=FALSE------------------------------------------------------
annotation <- yale_meningioma$annotation

## ---- eval=T-------------------------------------------------------------
kable(annotation[1:5, ])

## ---- eval=F-------------------------------------------------------------
#  loh <- readBAFExtractOutput ( path="./meningioma_baf\\", sequencing.type="bulk")
#  names(loh) <- gsub(".snp", "", names(loh))

## ---- eval=F-------------------------------------------------------------
#  object <- CreateCasperObject(raw.data=data,loh.name.mapping=loh.name.mapping, sequencing.type="bulk",
#    cnv.scale=3, loh.scale=3,
#    annotation=annotation, method="iterative", loh=loh,
#    control.sample.ids=control.sample.ids, cytoband=cytoband)

## ---- eval=T-------------------------------------------------------------
kable(yale_meningioma$loh.name.mapping[1:5, ])

## ---- eval=T-------------------------------------------------------------
data("scell_gbm")
kable(scell_gbm$loh.name.mapping[1:5, ])

## ---- eval=F-------------------------------------------------------------
#  final.objects <- runCaSpER(object, removeCentromere=T, cytoband=cytoband, method="iterative")

## ----include=FALSE-------------------------------------------------------
load("final.objects.yale.mn.rda")

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

## ---- eval=F-------------------------------------------------------------
#  obj <- final.objects[[9]]
#  plotHeatmap(object=obj, fileName="heatmap.png",cnv.scale= 3, cluster_cols = F, cluster_rows = T, show_rownames = T, only_soi = T)

## ---- echo=FALSE, out.width = '50%'--------------------------------------
knitr::include_graphics("test.png")

## ---- eval=F-------------------------------------------------------------
#  plotLargeScaleEvent (object=obj, fileName="large.scale.events.png")

## ---- echo=FALSE, out.width = '50%'--------------------------------------
knitr::include_graphics("large.scale.events.png")

## ---- eval=F-------------------------------------------------------------
#  plotGEAndGT (chrMat=finalChrMat, genoMat=genoMat, fileName="RNASeqAndGT.png")

## ---- echo=FALSE, out.width = '50%'--------------------------------------
knitr::include_graphics("RNASeqAndGT.png")

## ---- eval=F-------------------------------------------------------------
#  plotBAFAllSamples (loh = obj@loh.median.filtered.data,  fileName="LOHAllSamples.png")

## ---- echo=FALSE, out.width = '50%'--------------------------------------
knitr::include_graphics("LOHAllSamples.png")

## ---- eval=F-------------------------------------------------------------
#  plotBAFOneSample (object, fileName="LOHPlotsAllScales.pdf")

## ---- echo=FALSE, out.width = '50%'--------------------------------------
knitr::include_graphics("MN-60.png")

## ---- eval=F-------------------------------------------------------------
#  plotBAFInSeperatePages (loh=obj@loh.median.filtered.data, folderName="LOHPlots")

## ---- echo=FALSE, out.width = '50%'--------------------------------------
knitr::include_graphics("MN-5.png")

## ---- eval=F-------------------------------------------------------------
#  plotGEAndBAFOneSample (object=obj, cnv.scale=3, loh.scale=3, sample= "MN-5")

## ---- echo=FALSE,  out.width = '50%'-------------------------------------
knitr::include_graphics("MN-5_GE_BAF.png")

## ---- eval=F-------------------------------------------------------------
#  plotSingleCellLargeScaleEventHeatmap(finalChrMat, sampleName="MGH31", chrs=c("5p", "14q"))

## ---- echo=FALSE,  out.width = '10%'-------------------------------------
knitr::include_graphics("mgh31_ls.png")

## ---- eval=F-------------------------------------------------------------
#  ## calculate significant mutual exclusive and co-occurent events
#  results <- extractMUAndCooccurence (finalChrMat, loh, loh.name.mapping)
#  ## visualize mutual exclusive and co-occurent events
#  plotMUAndCooccurence (results)
#  

## ---- echo=FALSE,  out.width = '50%'-------------------------------------
knitr::include_graphics("MGH30network.png")

