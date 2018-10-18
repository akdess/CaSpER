
library(CaSpER)
## "yale_meningioma_data.rda" contains the following objects: 
## data: normalized gene expression matrix
## loh.name.mapping: data.frame for mapping loh files to expression files
## annotation
## control.sample.ids: samples that are used as normal
## samps: sample information
## genoMat: genotyping large scale CNV event summary 1: amplification, -1:deletion, 0: neutral
data("yale_meningioma")

## "hg19_cytoband.rda" contains the following objects: 
## cytoband: hg19 cytoband information
## centromere: hg19 centromere information
data("hg19_cytoband")

data <- yale_meningioma$data
annotation <- yale_meningioma$annotation
loh <-  yale_meningioma$loh
loh.name.mapping <-  yale_meningioma$loh.name.mapping
control.sample.ids <-  yale_meningioma$control.sample.ids
genoMat <-  yale_meningioma$genoMat
samps <-  yale_meningioma$samps

## generate annotation data.frame
#curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" | gunzip -c | grep acen | head
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(data), ishg19=T, centromere)
data <- data[match( annotation$Gene,rownames(data)), ]

## read BAF extract output  
#loh <- readBAFExtractOutput ( path="./meningioma_baf\\", sequencing.type="bulk")
#names(loh) <- gsub(".snp", "", names(loh))

## create CaSpER object 
object <- CreateCasperObject(raw.data=data,loh.name.mapping=loh.name.mapping, sequencing.type="bulk", 
  cnv.scale=3, loh.scale=3,
  annotation=annotation, method="iterative", loh=loh, 
  control.sample.ids=control.sample.ids, cytoband=cytoband)

## runCaSpER
final.objects <- runCaSpER(object, removeCentromere=T, cytoband=cytoband, method="iterative")

## sample plot orders
order.sampleNames <- c("MN-1171",  "MN-60835" ,"MN-1236" , "MN-1237" , "MN-1137" , "MN-1161" , "MN-60" ,   "MN-5" ) 

## plot median filtered gene expression matrix 
obj <- final.objects[[9]]
plotHeatmap(object=obj, fileName="heatmap.png",cnv.scale= 3, cluster_cols = F, cluster_rows = T, show_rownames = T, only_soi = T)

## summarize large scale events 
finalChrMat <- extractLargeScaleEvents (final.objects, thr=0.75) 
common <- intersect(order.sampleNames, intersect(rownames(finalChrMat), rownames(genoMat <- )))
finalChrMat <- finalChrMat[match(common, rownames(finalChrMat)), ]
genoMat <- genoMat[match(common, rownames(genoMat)), ]
    
## calculate TPR and FPR using genotyping array as gold standard
calcROC(chrMat=finalChrMat, chrMat2=genoMat)

#### VISUALIZATION 

## plot large scale events called from genotyping and rna-seq (can be used only with small sample size)
plotGEAndGT (chrMat=finalChrMat, genoMat=genoMat, fileName="RNASeqAndGT.png")
## plot large scale events 
plotLargeScaleEvent (object=obj, fileName="large.scale.events.pdf") 
## plot large scale events using event summary matrix 1: amplification, -1:deletion, 0: neutral
plotLargeScaleEvent2 (finalChrMat, fileName="large.scale.events.summarized.pdf")
## plot BAF deviation for each sample in seperate pages 
plotBAFInSeperatePages (loh =obj@loh.median.filtered.data, folderName="LOHPlotsSeperate") 
## plot BAF deviation for all samples together in one plot (can be used only with small sample size)
plotBAFAllSamples (loh = obj@loh.median.filtered.data,  fileName="LOHAllSamples.png") 
## plot gene expression signal for each sample seperately
plotGEAllSamples (object=obj, fileName="GEAllSamples.pdf", cnv.scale=3) 
## plot gene expression and BAF signal for one sample in one plot
plotGEAndBAFOneSample (object=obj, cnv.scale=3, loh.scale=3, sample= "MN-1237")
## plot BAF signal in different scales for all samples
plotBAFOneSample (object, fileName="LohPlotsAllScales.pdf") 
  
  