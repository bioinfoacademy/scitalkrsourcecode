## We're going to be using the Seurat R package for processing and clustering the scRNA-seq data. 
## Seurat 2.1.0 was used for our analysis.
## First load the relevant packages we'll be using.

library(ggplot2)
library(Seurat)
library(plyr)
library(dplyr)
library(Matrix)
library(cowplot)

## Read in an aggregate of the GFP data containing Sham day-3, Sham day-7, MI-day 3 and MI-day 7
gfp.data <- Read10X("data/Pdgfra-GFP_ShamVsMI_Days3_and_7")

gfp.aggr <- CreateSeuratObject(gfp.data, min.cells = 10, min.genes = 200, project = "GFP aggregate")
remove(gfp.data) #  remove the raw data

## Add some meta-data to the Seurat object. Will add percent of RNA mapped to 
## mitochondrial genes and the batch ID so we can compare conditions later.
mito.genes <- grep("^mt-", rownames(gfp.aggr@data), value = TRUE)
percent.mito <- Matrix::colSums(gfp.aggr@raw.data[mito.genes, ])/Matrix::colSums(gfp.aggr@raw.data)
gfp.aggr <- AddMetaData(gfp.aggr, percent.mito, "percent.mito")

## Add batch ID to meta data
## batches are in the following order:
## 1 - Sham-day 3
## 2 - Sham-day 7
## 3 - MI-day 3
## 4 - MI-day 7
batch.id <- sub(".*-(.*)","\\1", gfp.aggr@cell.names)
batch.id <- replace(batch.id, batch.id=="1", "Sham-day 3")
batch.id <- replace(batch.id, batch.id=="2", "Sham-day 7")
batch.id <- replace(batch.id, batch.id=="3", "MI-day 3")
batch.id <- replace(batch.id, batch.id=="4", "MI-day 7")
table(batch.id)
names(batch.id) = gfp.aggr@cell.names
gfp.aggr <- AddMetaData(gfp.aggr, batch.id, "batch")

#Violin plot of QC metrics
VlnPlot(object = gfp.aggr, features.plot = c("nGene", "nUMI", "percent.mito"), 
        nCol = 3, point.size.use = 0.5)

par(mfrow = c(1, 2))
GenePlot(gfp.aggr, "nUMI", "percent.mito", cex.use=1)
GenePlot(gfp.aggr, "nUMI", "nGene",  cex.use=1)
par(mfrow = c(1, 1))

## Appears to be some cell outliers with unusually large numbers of UMIs or genes - will filter these out.
## Filter out cells likely to be doublets. Set a threshold of 4000
gfp.aggr <- FilterCells(object = gfp.aggr, subset.names = c("nGene", "nUMI", "percent.mito"), 
                        low.thresholds = c(200, 500, -Inf), 
                        high.thresholds = c(5000, 30000, 0.05))

## Normalise the data 
gfp.aggr <- NormalizeData(object = gfp.aggr, normalization.method = "LogNormalize", 
                          scale.factor = 10000, display.progress = F)

## detect highly-variable genes to be used for PCA
gfp.aggr <- FindVariableGenes(object = gfp.aggr, mean.function = ExpMean, 
                              dispersion.function = LogVMR, display.progress = F, 
                              x.low.cutoff = 0.01, x.high.cutoff = 5, 
                              y.cutoff = 0.5, y.high.cutoff = 15)

print(paste0("Number of highly variable genes: ", length(x = gfp.aggr@var.genes)))

## Scale the data, regressing out variation due to number of UMIs
gfp.aggr <- ScaleData(gfp.aggr, vars.to.regress = c("nUMI"))

## PCA analysis. Will run PCA up to the first 60 components
gfp.aggr <- RunPCA(object = gfp.aggr, pc.genes = gfp.aggr@var.genes, do.print = FALSE, pcs.compute=60)

## Check what genes have highest correlation with the top PCs
PCHeatmap(object = gfp.aggr, pc.use = c(1:6), cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

## We use the JackStraw function with 1000 replicates to identify PCs for clustering. 
## This takes a long time to run, but the output we obtained is saved as a pdf file.
gfp.aggr <- JackStraw(object = gfp.aggr, num.pc = 60, num.replicate = 1000, do.print = FALSE)

## The JackStraw plot indiciates that up to the first 53 PCs are significant with P < 0.001
JackStrawPlot(object = gfp.aggr, PCs = 1:60)

## We use the PCs identified by the JackStraw test, but have found that modification of 
## the PCs only causes minor changes in the clustering solutions. 
## Run clustering with a range of resolutions
gfp.aggr <- FindClusters(object = gfp.aggr, reduction.type = "pca", dims.use = 1:53, 
                         resolution = seq(0.6,1.8,0.2), print.output = 0, save.SNN = TRUE)

# Found that a resolution of 0.6 as gives a good higher-level breakup of the data.
# But can increase this to get more clusters
gfp.aggr <- SetAllIdent(gfp.aggr, id="res.0.6")

# t-SNE analysis. Seed was selected for aesthetic purposes. 
gfp.aggr <- RunTSNE(gfp.aggr, dims.use = 1:53, do.fast = T)

clusters = gfp.aggr@ident
cluster.idents = names(table(clusters))
ident.replacements = as.numeric(cluster.idents)+1
names(ident.replacements) <- cluster.idents
clusters=revalue(clusters, ident.replacements)
gfp.aggr@ident = clusters

## Generate a t-SNE plot showing clusters
TSNEPlot(gfp.aggr, pt.size = 0.5)

## Pathway & GO analysis of genes up-regulated in each cluster indicated that cluster 9
## likely corresponds to stressed or apoptotic cells. We therefore remove these from the downstream analysis
remainder.cells = gfp.aggr@cell.names[gfp.aggr@ident != "9"]
gfp.aggr = FilterCells(gfp.aggr, subset.names=NULL, cells.use = remainder.cells)
clusters = gfp.aggr@ident
ident.replacements = c("9", "10", "11")
names(ident.replacements) = c("10", "11", "12")
clusters=revalue(clusters, ident.replacements)
gfp.aggr@ident = clusters

TSNEPlot(gfp.aggr, pt.size = 0.5)

## Finally can do a 'pretty' representation of the data. Based on our analysis we have defined the populations as:
## (1) Fibroblast: Sca1-low (F-SL)
## (2) Myofibroblast (MYO)
## (3) Fibroblast: Activated (F-Act)
## (4) Fibroblast: Sca1-high (F-SH)
## (5) Fibroblast: Transitory (F-Trans)
## (6) Fibroblast: Wnt expressing (F-WntX)
## (7) Fibroblast: Cycling (F-Cyc)
## (8) Fibroblast: IFN stimulated (F-IFNS)
## (9) Macrophage (MAC)
## (10) Epicardial (EPI)
## (11) Endothelial Cell (EC)

current.labels = as.character(1:11)
new.labels = c("F-SL", "MYO", "F-Act", "F-SH", "F-Trans", "F-WntX", "F-Cyc", "F-IFNS", "MAC", "EPI", "EC")

gfp.aggr@ident = plyr::mapvalues(x = gfp.aggr@ident, from = current.labels, to = new.labels)

col.set <- c("#fb8500", "#0099cc", "#f60000", "#fde800","#00b600","#00896e", "#0053c8", "#ff3399", "#b374e0", "#00cc99",  "#fbc9e2")

TSNEPlot(gfp.aggr, colors.use = col.set, pt.size = 0.5)

## Finally can save the object as an R data file
save(gfp.aggr, file = "data/GFP_notebook_Seurat.RData")
