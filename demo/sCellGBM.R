
library(CaSpER)
library('biomaRt')

## "sCell_gbm_data.rda" contains the following objects: 
## data: normalized gene expression matrix
## loh.name.mapping: data.frame for mapping loh files to expression files
## annotation
## loh 
load("sCell_gbm_data.rda")

## create CaSpER object
object <- CreateCasperObject(raw.data=data,loh.name.mapping=loh.name.mapping, 
              sequencing.type="single-cell",
               cnv.scale=3, loh.scale=3,
              annotation=annotation, method="iterative", loh=loh, 
              control.sample.ids="REF", cytoband=cytoband)

## runCaSpER
final.objects <- runCaSPER(object, removeCentromere=T, cytoband=cytoband, method="iterative")

## plot median filtered gene expression matrix 
plotHeatmap(object, fileName="heatmap.png", cnv.scale = 3, show_rownames=F)

## summarize large scale events 
finalChrMat <- extractLargeScaleEvents (final.objects, thr=0.75) 

#### VISUALIZATION 

## plot large scale events using event summary matrix 1: amplification, -1:deletion, 0: neutral
plotLargeScaleEvent2 (finalChrMat, fileName="large.scale.events.summarized.pdf") 
## plot BAF deviation for all samples together in one plot (can be used only with small sample size)
plotBAFAllSamples (loh = obj@loh.median.filtered.data,  fileName="LOHAllSamples.png") 
## plot BAF signal in different scales for all samples
plotBAFOneSample (object, fileName="LohPlotsAllScales.pdf") 
## plot large scale event summary for selected sample and chromosomes
plotSingleCellLargeScaleEventHeatmap(finalChrMat, sampleName="MGH31", chrs=c("5p", "14q"))

## calculate significant mutual exclusive and co-occurent events
results <- extractMUAndCooccurence (finalChrMat, loh, loh.name.mapping)
## visualize mutual exclusive and co-occurent events
plotMUAndCooccurence (results)
