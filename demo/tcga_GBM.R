

library(CaSpER)
library('biomaRt')
## "tcga_gbm_data.rda" contains the following objects: 
## data: normalized gene expression matrix
## loh.name.mapping: data.frame for mapping loh files to expression files
## annotation
## loh
## control.sample.ids: samples that are used as normal
## samps: sample information
## genoMat: genotyping large scale CNV event summary 1: amplification, -1:deletion, 0: neutral
load("./data/tcga_gbm_data.rda")

## "hg19_cytoband.rda" contains the following objects: 
## cytoband: hg19 cytoband information
## centromere: hg19 centromere information
load("./data/hg19_cytoband.rda")

## generate annotation data.frame

annotation <- generateAnnotation(id_type="hgnc_symbol", genes=rownames(data), ishg19=T, centromere)
data <- data[match( annotation$Gene,rownames(data)), ]

## create CaSpER object 

object <- CreateCasperObject(raw.data=data, loh.name.mapping=loh.name.mapping, 
    sequencing.type="bulk", cnv.scale=3, loh.scale=3,method="iterative",
              annotation=annotation, loh=loh, 
              control.sample.ids=control.sample.ids, cytoband=cytoband)

library(scales)

## runCaSpER
final.objects <- runCaSPER(object, removeCentromere=T, cytoband=cytoband, method="iterative")

## plot median filtered gene expression matrix
obj <- final.objects[[9]]
plotHeatmap(object, fileName="heatmap.png", cnv.scale = 3, show_rownames=F)

## summarize large scale events 
finalChrMat <- extractLargeScaleEvents (final.objects, thr=0.75) 
common <- intersect(order.sampleNames, intersect(rownames(finalChrMat), rownames(genoMat)))
finalChrMat <- finalChrMat[match(common, rownames(finalChrMat)), ]
genoMat <- genoMat[match(common, rownames(genoMat)), ]

## calculate TPR and FPR using genotyping array as gold standard 
calcROC(chrMat=finalChrMat, chrMat2=genoMat)

#### VISUALIZATION 

## plot large scale events 
plotLargeScaleEvent (object=obj, fileName="large.scale.events.pdf") 
## plot large scale events using event summary matrix 1: amplification, -1:deletion, 0: neutral
plotLargeScaleEvent2 (finalChrMat, fileName="large.scale.events.summarized.pdf") 
## plot BAF deviation for each sample in seperate pages 
plotBAFInSeperatePages (loh =obj@loh.median.filtered.data, folderName="LOHPlotsSeperate") 
## plot gene expression and BAF signal for one sample in one plot
plotGEAndBAFOneSample (object=obj, cnv.scale=3, loh.scale=3, sample= "TCGA-02-0047-01A")
