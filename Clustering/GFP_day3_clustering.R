## We're going to be using the Seurat R package for processing and clustering the scRNA-seq data. 
## Seurat 2.1.0 was used for our analysis.
## First load the relevant packages we'll be using.

library(ggplot2)
library(Seurat)
library(plyr)
library(dplyr)
library(Matrix)
library(cowplot)

gfp.d3.data <- Read10X("data/Pdgfra-GFP_ShamVsMI_Day3")

#filter out genes expressing in under 5 cells, and cells with under 200 genes
gfp.d3 <- CreateSeuratObject(gfp.d3.data, min.cells = 10, min.genes = 200, project = "Day3")
remove(gfp.d3.data)

## Add some meta-data to the Seurat object. Will add percent of RNA mapped to 
## mitochondrial genes and the batch ID so we can compare conditions later.
mito.genes <- grep("^mt-", rownames(gfp.d3@data), value = TRUE)
percent.mito <- Matrix::colSums(gfp.d3@raw.data[mito.genes, ])/Matrix::colSums(gfp.d3@raw.data)
gfp.d3 <- AddMetaData(gfp.d3, percent.mito, "percent.mito")

## Add batch ID to meta data
## batches are in the following order:
## 1 - Sham-day 3
## 2 - MI-day 3
batch.id <- sub(".*-(.*)","\\1", gfp.d3@cell.names)
batch.id <- replace(batch.id, batch.id=="1", "Sham-day 3")
batch.id <- replace(batch.id, batch.id=="2", "MI-day 3")
table(batch.id)
names(batch.id) = gfp.d3@cell.names
gfp.d3 <- AddMetaData(gfp.d3, batch.id, "batch")


#Violin plot of QC metrics
VlnPlot(object = gfp.d3, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(gfp.d3, "nUMI", "percent.mito", cex.use=1)
GenePlot(gfp.d3, "nUMI", "nGene",  cex.use=1)
par(mfrow = c(1, 1))

## Appears to be some cell outliers with unusually large numbers of UMIs or genes - will filter these out.
gfp.d3 <- FilterCells(object = gfp.d3, subset.names = c("nGene", "nUMI", "percent.mito"), 
            low.thresholds = c(200, 500, -Inf), 
            high.thresholds = c(6000, 50000, 0.05))

## Normalise the data 
gfp.d3 <- NormalizeData(object = gfp.d3, normalization.method = "LogNormalize", scale.factor = 10000)

## detect highly-variable genes to be used for PCA
gfp.d3 <- FindVariableGenes(object = gfp.d3, mean.function = ExpMean, dispersion.function = LogVMR, 
                  x.low.cutoff = 0.01, x.high.cutoff = 5, y.cutoff = 0.5, y.high.cutoff = 15)
print(paste0("Number of highly variable genes: ", length(x = gfp.d3@var.genes)))

## Scale the data, regressing out variation due to number of UMIs
gfp.d3 <- ScaleData(gfp.d3, vars.to.regress = c("nUMI"))

## PCA analysis. Will run PCA up to the first 50 components
gfp.d3 <- RunPCA(object = gfp.d3, pc.genes = gfp.d3@var.genes, do.print = TRUE, pcs.compute=50)

## Check what genes have highest correlation with the top PCs
PCHeatmap(object = gfp.d3, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

## We use the JackStraw function with 1000 replicates to identify PCs for clustering. 
## This takes a long time to run, but the number of PCs we obtain as significant is below
gfp.d3 <- JackStraw(object = gfp.d3, num.pc = 50, num.replicate = 1000, do.print = FALSE)
JackStrawPlot(object = gfp.d3, PCs = 1:50)
PCElbowPlot(gfp.d3, num.pc = 50)


## For Sham/MI aggregate, it is clear that the first 41 PCs are statistically significant (P < 0.001).
## Will proceed with the clustering analysis based on PCs 1-41.
gfp.d3 <- FindClusters(object = gfp.d3, reduction.type = "pca", dims.use = 1:41, resolution = 0.6,
                     print.output = 0, save.SNN = TRUE)

### t-SNE analysis
gfp.d3 <- RunTSNE(object = gfp.d3, dims.use = 1:41, do.fast = T)
TSNEPlot(gfp.d3)

## Update cluster numbers to start from 1
clusters = gfp.d3@ident
cluster.idents = names(table(clusters))
ident.replacements = as.numeric(cluster.idents)+1
gfp.d3@ident = plyr::mapvalues(x = gfp.d3@ident, from = cluster.idents, to = ident.replacements)

TSNEPlot(gfp.d3)

## There aren't any sub-populations coming out from the analysis that appear to coorespond to 
## stressed cells. But will remove any detected from the full GFP aggregate clustering before further analysis. 
stressed.cells = scan("GFP_C9_stressed_cell_list.txt", what = "character")
# batches from stressed cells file should be in the following order:
# 1 - Sham-day3
# 2 - Sham-day7
# 3 - MI-day3
# 4 - MI-day7
cellSample <- sub(".*-(.*)","\\1", stressed.cells)
d3Sham.stressed.cells = stressed.cells[cellSample == "1"]
d3MI.stressed.cells = stressed.cells[cellSample == "3"]
d3MI.stressed.cells = paste0(sub("(.*)-.*","\\1", d3MI.stressed.cells), "-2")
d3.stressed.cells = c(d3Sham.stressed.cells, d3MI.stressed.cells)
paste("Removing", length(d3.stressed.cells), "stressed cells from further analysis.")

remainder.cells = setdiff(gfp.d3@cell.names, d3.stressed.cells)
gfp.d3 = FilterCells(gfp.d3, subset.names=NULL, cells.use = remainder.cells)

TSNEPlot(gfp.aggr, pt.size = 0.5)

## Finally can do a 'pretty' representation of the data. Based on our analysis we have defined the populations as:
## (1) Fibroblast: Activated (F-Act)
## (2) Fibroblast: Sca1-low (F-SL)
## (3) Fibroblast: Sca1-high (F-SH)
## (4) Fibroblast: Transitory (F-Trans)
## (5) Fibroblast: Wnt expressing (F-WntX)
## (6) Fibroblast: Cycling intermediate (F-CI)
## (7) Fibroblast: Cycling (F-Cyc)
## (8) Macrophage (MAC)
## (9) Fibroblast: IFN stimulated (F-IFNS)
## (10) Endothelial Cell (EC)

## Labels for cell populations
new.labels = c("F-Act","F-SL", "F-SH", "F-Trans", "F-WntX", "F-CI", "F-Cyc", "MAC", "F-IFNS", "EC")
current.labels = as.character(names(table(gfp.d3@ident)))
gfp.d3@ident = plyr::mapvalues(x = gfp.d3@ident, from = current.labels, to = new.labels)

col.set <- c("#f60000","#fb8500", "#fde800", "#00b600", "#00896e","#00eefd","#0053c8","#b374e0", "#ff3399", "#fbc9e2")
TSNEPlot(gfp.d3, colors.use = col.set, pt.size = 0.5)

## Finally can save the object as an R data file
save(gfp.d3, file = "data/GFP_day3_Seurat.RData")

