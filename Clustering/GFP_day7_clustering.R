## We're going to be using the Seurat R package for processing and clustering the scRNA-seq data. 
## Seurat 2.1.0 was used for our analysis.
## First load the relevant packages we'll be using.

library(ggplot2)
library(Seurat)
library(plyr)
library(dplyr)
library(Matrix)

## Read in an aggregate of the GFP data containing Sham day-3, Sham day-7, MI-day 3 and MI-day 7
gfp.d7.data <- Read10X("data/Pdgfra-GFP_ShamVsMI_Day7")

gfp.d7 <- CreateSeuratObject(gfp.d7.data, min.cells = 10, min.genes = 200, project = "Day7")

mito.genes <- grep("^mt-", rownames(gfp.d7@data), value = TRUE)
percent.mito <- Matrix::colSums(gfp.d7@raw.data[mito.genes, ])/Matrix::colSums(gfp.d7@raw.data)
gfp.d7 <- AddMetaData(gfp.d7, percent.mito, "percent.mito")

## Add batch ID to meta data
## batches are in the following order:
## 1 - Sham-day 7
## 2 - MI-day 7
batch.id <- sub(".*-(.*)","\\1", gfp.d7@cell.names)
batch.id <- replace(batch.id, batch.id=="1", "Sham-day 7")
batch.id <- replace(batch.id, batch.id=="2", "MI-day 7")
table(batch.id)
names(batch.id) = gfp.d7@cell.names
gfp.d7 <- AddMetaData(gfp.d7, batch.id, "batch")

#Violin plot of QC metrics
VlnPlot(gfp.d7, c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(gfp.d7, "nUMI", "percent.mito", cex.use=1)
GenePlot(gfp.d7, "nUMI", "nGene",  cex.use=1)
par(mfrow = c(1, 1))

## Appears to be some cell outliers with unusually large numbers of UMIs or genes - will filter these out.
gfp.d7 <- FilterCells(object = gfp.d7, subset.names = c("nGene", "nUMI", "percent.mito"), 
                    low.thresholds = c(200, 500, -Inf), high.thresholds = c(5000, 25000, 0.05))

# Normalise data
gfp.d7 <- NormalizeData(object = gfp.d7, normalization.method = "LogNormalize", scale.factor = 10000)

## detect highly-variable genes to be used for PCA
gfp.d7 <- FindVariableGenes(object = gfp.d7, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.01, x.high.cutoff = 5, y.cutoff = 0.5, y.high.cutoff = 15)
length(x = gfp.d7@var.genes)

## Scale the data, regressing out variation due to number of UMIs
gfp.d7 <- ScaleData(gfp.d7, vars.to.regress = c("nUMI"))

## PCA analysis. Will run PCA up to the first 50 components
gfp.d7 <- RunPCA(object = gfp.d7, pc.genes = gfp.d7@var.genes, do.print = TRUE, pcs.compute=50)

## Check what genes have highest correlation with the top PCs
PCHeatmap(object = gfp.d7, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

## We use the JackStraw function with 1000 replicates to identify PCs for clustering. 
## This takes a long time to run, but the number of PCs we obtain as significant is below
gfp.d7 <- JackStraw(object = gfp.d7, num.pc = 50, num.replicate = 1000, do.print = FALSE)
JackStrawPlot(object = gfp.d7, PCs = 1:50)
PCElbowPlot(gfp.d7, num.pc = 50)

## For clustering, we use the PCs identified by the JackStraw test (P<0.001)
gfp.d7 <- FindClusters(object = gfp.d7, reduction.type = "pca", dims.use = 1:44, resolution = seq(0.6,2.0,0.2),
                     print.output = 0, save.SNN = TRUE)

### t-SNE analysis
gfp.d7 <- RunTSNE(object = gfp.d7, dims.use = 1:44, do.fast = T, seed.use=1)
gfp.d7 <- SetAllIdent(gfp.d7, id="res.1.2")

## Update cluster IDs to start from 1
clusters = gfp.d7@ident
cluster.idents = names(table(clusters))
ident.replacements = as.numeric(cluster.idents)+1
names(ident.replacements) <- cluster.idents
clusters=revalue(clusters, ident.replacements)
gfp.d7@ident = clusters
TSNEPlot(gfp.d7)

source("cluster_merging.R")
gfp.d7.merged = generateClusterMerging(gfp.d7, 4, gene.set="full")
gfp.d7 = gfp.d7.merged

## Remove cluster 11, which appears to contain stressed cells and relabel
remainder.cells = gfp.d7@cell.names[gfp.d7@ident != "11"]
gfp.d7 = FilterCells(gfp.d7, subset.names=NULL, cells.use = remainder.cells)
clusters = gfp.d7@ident
ident.replacements = c("11")
names(ident.replacements) = c("12")
clusters=revalue(clusters, ident.replacements)
gfp.d7@ident = clusters

## Finally can do a 'pretty' representation of the data. Based on our analysis we have defined the populations as:
## (1) Fibroblast: Sca1-low (F-SL) 
## (2) Fibroblast: Sca1-high (F-SH)
## (3) Myofibroblast 3 (MYO-1)
## (4) Fibroblast: Activated (F-Act)
## (5) Myofibroblast 3 (MYO-2)
## (6) Fibroblast: Transitory (F-Trans)
## (7) Fibroblast: Wnt expressing (F-WntX)
## (8) Myofibroblast 3 (MYO-3)
## (9) Fibroblast: IFN stimulated (F-IFNS)
## (10) Fibroblast: Cycling (F-Cyc)
## (11) Epicardial (EPI)

current.labels = as.character(1:11)
new.labels = c("F-SL", "F-SH", "MYO-1", "F-Act", "MYO-2", "F-Trans", "F-WntX", "MYO-3", "F-IFNS", "F-Cyc","EPI")
gfp.d7@ident = plyr::mapvalues(x = gfp.d7@ident, from = current.labels, to = new.labels)

col.set <- c("#fb8500", "#fde800", "#0099cc", "#f60000", "#669999","#00b600","#00896e", "#bf00c9", "#ff3399", "#0053c8",  "#00cc99")
TSNEPlot(gfp.d7, colors.use = col.set, pt.size = 0.5, do.label = TRUE)


## Finally can save the object as an R data file
save(gfp.d7, file = "data/GFP_day7_Seurat.RData")

