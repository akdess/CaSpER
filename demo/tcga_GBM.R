

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
load("tcga_gbm_data.rda")

## "hg19_cytoband.rda" contains the following objects: 
## cytoband: hg19 cytoband information
## centromere: hg19 centromere information
load("hg19_cytoband.rda")

## generate annotation data.frame
mart <- useDataset("hsapiens_gene_ensembl", useEnsembl(biomart="ensembl",GRCh=37))
genes <-  rownames(data)
G_list <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol",'chromosome_name', 'start_position', 'end_position', 'strand'), values=genes, mart= mart)
common <- intersect(rownames(data), G_list$hgnc_symbol)
ord <- match(common, G_list$hgnc_symbol)
annotation <- G_list[ord, ]
annotation <- annotation[order(annotation$start_position), ]
chr.list =as.character(1:22, "X")
idx = unlist(as.vector((sapply(chr.list,function(x) as.vector(unlist(which(as.character(annotation$chromosome_name)==x)))))))
annotation<-as.data.frame(annotation[idx,])
data <- data[match( annotation$hgnc_symbol,rownames(data)), ]

colnames(annotation)[1:4] <- c("Gene", "Chr", "start", "end")
annotation$isCentromer <- rep("no", nrow(annotation))

centromere_snps<-NULL
for (k in 1:(dim(centromere)[1]))
{
  annotation$isCentromer[which(as.character(annotation$Chr)==gsub("chr","",as.character(centromere$V1[k])) & 
                                             (as.numeric(as.character(annotation$Position))>=centromere$V2[k] & as.numeric(as.character(coord$Position))<=centromere$V3[k]))] <- "yes"
}

annotation$Position <- (as.numeric(annotation$start) + as.numeric(annotation$end))/2
annotation$cytoband<-rep("", nrow(annotation))

for (k in 1:(dim(annotation)[1]))
{
	annotation$cytoband[which(as.character(annotation$Chr)==gsub("chr","",as.character(cytoband$V1[k])) & 
		(as.numeric(as.character(annotation$Position))>=cytoband$V2[k] & as.numeric(as.character(annotation$Position))<=cytoband$V3[k]))] <-paste(as.character(cytoband$V1[k]),as.character(cytoband$V4[k]),sep="")
}

annotation$new_positions <- as.vector(unlist(lapply(lapply(split(annotation$cytoband,   annotation$cytoband), length)[unique(annotation$cytoband)], function(x) 1:x)))

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
