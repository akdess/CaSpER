
library(Seurat)
library(CaSpER)
data("hg19_cytoband")

# expression data is downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110499

counts  <- read.delim("GSE110499_GEO_processed_MM_10X_raw_UMI_count_martix.txt", stringsAsFactor=F, header=T)
rownames(counts) <- counts[, 1] 
counts <- counts[, -1]

mm135 <- CreateSeuratObject(counts = counts, project = "mm135", min.cells = 3, min.features = 200)
mm135[["percent.mt"]] <- PercentageFeatureSet(mm135, pattern = "^MT-")
mm135 <- subset(mm135, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
mm135 <- NormalizeData(mm135 , scale.factor = 1e6, normalization.method = "RC")
mm135 <- FindVariableFeatures(mm135, do.plot = T, nfeatures = 1000)
mm135 <- ScaleData(mm135)

mm135 <- RunPCA(mm135, features = VariableFeatures(object = mm135),npcs = 100)
mm135 <- RunTSNE(mm135, dims.use = 1:10)
DimPlot(mm135, reduction = "tsne")
FeaturePlot(mm135, features = c("SDC1", "CD38"))

mm135 <- FindNeighbors(mm135, dims = 1:10)
mm135 <- FindClusters(mm135, resolution = 0.5)
DimPlot(mm135, reduction = "tsne", label=T)

log.ge <- as.matrix(mm135@assays$RNA@data)
control <- names(Idents(mm135) )[Idents(mm135) %in% c(2,7)]

genes <- rownames(log.ge)
annotation <- generateAnnotation(id_type="hgnc_symbol", genes=genes, centromere=centromere, ishg19 = T)
log.ge <- log.ge[match( annotation$Gene,rownames(log.ge)) , ]
rownames(log.ge) <- annotation$Gene
log.ge <- log2(log.ge +1)


load("maf.rda") ## from https://github.com/akdess/CaSpER/blob/master/data/maf.rda
loh<- list()
loh[[1]] <- maf
names(loh) <- "MM135"
loh.name.mapping <- data.frame (loh.name= "MM135" , sample.name=colnames(log.ge))


object <- CreateCasperObject(raw.data=log.ge,loh.name.mapping=loh.name.mapping, sequencing.type="single-cell", 
cnv.scale=3, loh.scale=3, 
expr.cutoff=0.1, filter="median", matrix.type="normalized",
annotation=annotation, method="iterative", loh=loh, 
control.sample.ids=control, cytoband=cytoband)


pdf("MM135.Distrubution.pdf")
plot(density(as.vector(object@control.normalized[[3]])))
plot(density(log2(object@control.normalized.noiseRemoved[[3]]+1)))
dev.off()

## runCaSpER
final.objects <- runCaSpER(object, removeCentromere=T, cytoband=cytoband, method="iterative")

## summarize large scale events 
finalChrMat <- extractLargeScaleEvents (final.objects, thr=0.75) 

obj <- final.objects[[9]]
plotHeatmap10x(object=obj, fileName="heatmap.png",cnv.scale= 3, cluster_cols = F, cluster_rows = T, show_rownames = T, only_soi = T)

#### VISUALIZATION 
chrMat <- finalChrMat
plot.data <- melt(chrMat)
plot.data$value2 <- "neutral"
plot.data$value2[plot.data$value > 0] <- "amplification"
plot.data$value2[plot.data$value < 0] <- "deletion"
plot.data$value2 <- factor(plot.data$value2, levels = c("amplification", 
    "deletion", "neutral"))
plot.data$X2 <- factor(plot.data$X2, levels = colnames(chrMat))
p <- ggplot(aes(x = X2, y = X1, fill = value2), data = plot.data) + 
    geom_tile(colour = "white", size = 0.01) + 
    labs(x = "", 
    y = "") + scale_fill_manual(values = c(amplification = muted("red"), 
    deletion = muted("blue"), neutral = "white")) + theme_grey(base_size = 6) + 
    theme(legend.position = "right", legend.direction = "vertical", 
        legend.title = element_blank(), strip.text.x = element_blank(), 
        legend.text = element_text(colour = "black", size = 7, 
            face = "bold"), legend.key.height = grid::unit(0.8, 
            "cm"), legend.key.width = grid::unit(0.5, "cm"), 
        axis.text.x = element_text(size = 5, colour = "black", 
            angle = -45, hjust = 0), axis.text.y = element_text(size = 6, 
            vjust = 0.2, colour = "black"), axis.ticks = element_line(size = 0.4), 
        plot.title = element_text(colour = "black", hjust = 0, 
            size = 6, face = "bold"))

