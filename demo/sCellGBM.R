
library(CaSpER)

## "sCell_gbm_data.rda" contains the following objects: 
## data: normalized gene expression matrix
## loh.name.mapping: data.frame for mapping loh files to expression files
## annotation
## loh 
data("scell_gbm")

data <- scell_gbm$data
loh <-  scell_gbm$loh
loh.name.mapping <-  scell_gbm$loh.name.mapping

## generate annotation data.frame
annotation <- generateAnnotation(id_type="hgnc_symbol", genes=rownames(data), ishg19=T, centromere)
data <- data[match(annotation$Gene,rownames(data)), ]

## create CaSpER object
object <- CreateCasperObject(raw.data=data,loh.name.mapping=loh.name.mapping, 
              sequencing.type="single-cell",
               cnv.scale=3, loh.scale=3,
              annotation=annotation, method="iterative", loh=loh, 
              control.sample.ids="REF", cytoband=cytoband)

## runCaSpER
final.objects <- runCaSpER(object, removeCentromere=T, cytoband=cytoband, method="iterative")

## plot median filtered gene expression matrix 
plotHeatmap(object, fileName="heatmap.png", cnv.scale= 3, cluster_cols = F, cluster_rows = T, show_rownames = T, only_soi = T)

## summarize large scale events 
finalChrMat <- extractLargeScaleEvents (final.objects, thr=0.75) 

#### VISUALIZATION 

## plot large scale events using event summary matrix 1: amplification, -1:deletion, 0: neutral
plotLargeScaleEvent2 (finalChrMat, fileName="large.scale.events.summarized.pdf") 
## plot BAF deviation for all samples together in one plot (can be used only with small sample size)
plotBAFAllSamples (loh = final.objects[[1]]@loh.median.filtered.data,  fileName="LOHAllSamples.png") 
## plot BAF signal in different scales for all samples
plotBAFOneSample (object, fileName="LohPlotsAllScales.pdf") 
## plot large scale event summary for selected sample and chromosomes
plotSingleCellLargeScaleEventHeatmap(finalChrMat, sampleName="MGH31", chrs=c("5p", "14q"))
plotSingleCellLargeScaleEventHeatmap(finalChrMat, sampleName="MGH29", chrs=c("10p", "4p"))
plotSingleCellLargeScaleEventHeatmap(finalChrMat, sampleName="MGH30", chrs=c("6p", "7p"))
plotSingleCellLargeScaleEventHeatmap(finalChrMat, sampleName="MGH31", chrs=c("5p", "14q"))
plotSingleCellLargeScaleEventHeatmap(finalChrMat, sampleName="MGH28", chrs=c("10p", "14q"))
plotSingleCellLargeScaleEventHeatmap(finalChrMat, sampleName="MGH28", chrs=c("5q", "19q"))
plotSingleCellLargeScaleEventHeatmap(finalChrMat, sampleName="MGH28", chrs=c("8q", "20p"))
plotSingleCellLargeScaleEventHeatmap(finalChrMat, sampleName="MGH26", chrs=c("5q", "7p"))
plotSingleCellLargeScaleEventHeatmap(finalChrMat, sampleName="MGH26", chrs=c("1q", "22q"))

## calculate significant mutual exclusive and co-occurent events
results <- extractMUAndCooccurence (finalChrMat, loh, loh.name.mapping)
## visualize mutual exclusive and co-occurent events
plotMUAndCooccurence (results)
## visualize CNV pyhlogenetic tree 
## first please download pyhlip package from http://evolution.genetics.washington.edu/phylip/getme-new1.html
#plotSCellCNVTree (finalChrMat, sampleName="MGH31", path="C:\\Users\\aharmanci\\Downloads\\phylip-3.695\\phylip-3.695\\exe", fileName="CNVTree.pdf")


all.summary <- getDiffExprGenes(final.objects,  sampleName="MGH31", chrs=c("5q", "14q"), event.type=c(1, -1))
#enrichment <- generateEnrichmentSummary (results=all.summary)

all.summary <- getDiffExprGenes(final.objects, sampleName="MGH31", chrs=c("1p", "13q"), event.type=c(1, -1))

all.summary <- getDiffExprGenes(final.objects,  sampleName="MGH31", chrs=c("5q", "13q"), event.type=c(1, -1))

all.summary <- getDiffExprGenes(final.objects,  sampleName="MGH28", chrs=c("8q", "20p"), event.type=c(1, -1))
all.summary <- getDiffExprGenes(final.objects,  sampleName="MGH30", chrs=c("7p", "6p"), event.type=c(1, -1))

all.summary <- getDiffExprGenes(final.objects,  sampleName="MGH26", chrs=c("1q", "22q"), event.type=c(1, -1))
genes <- as.character(all.summary$ID[all.summary$adj.P.Val<0.05])

all.summary <- getDiffExprGenes(final.objects,  sampleName="MGH28", chrs=c("5q", "19q"), event.type=c(1, -1))
genes <- as.character(all.summary$ID[all.summary$adj.P.Val<0.05])

all.summary <- getDiffExprGenes(final.objects,  sampleName="MGH29", chrs=c("10p", "4p"), event.type=c(-1, -1))
genes <- as.character(all.summary$ID[all.summary$adj.P.Val<0.05])
