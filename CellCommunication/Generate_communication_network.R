library(plyr)
library(network)
library(tidygraph)
library(igraph)
library(ggraph)
library(scales)
library(STRINGdb)
library(Seurat)
library(progress)

source("ligand_receptor_functions.R")

## This script will reproduce the cell communication results presented in Figure 3 of the paper
## It assumes a Seurat object named 'tip.aggr'.
## Later plotting functions assume the name of the clusters are according to naming applied
## in ../Clustering/TIP_aggregate_clustering.R.

## Load up the pre-processed & clustered Seurat object of TIP data
load("../data/TIP_aggregate_Seurat.RData")


col.set <- c("#c10023", "#008e17", "#fb8500", "#f60000", "#fde800", "#bc9000","#4ffc00", "#00bcac", "#0099cc",
             "#D35400", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8",
             "#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")

TSNEPlot(tip.aggr, colors.use = col.set, pt.size = 0.5, do.label = TRUE)

## Define clusters to use in the analysis
clusters.use = names(table(tip.aggr@ident))

## Define a label for output files
out.lab = "TIP"

## Read in mouse to human orthologue mappings
mappings = read.csv("ensembl_mouse_human_orthologues.txt", header=TRUE)
mappings = mappings[mappings[, 4]!='' & mappings[, 5]!='', ]
mappings = mappings[, c(4,5)]
mappings = unique(mappings)

## Read in the ligand-receptor pairs
## These were taken from Ramilowski et al. (2015) Nature Communications
ligand.receptor.pairs = read.csv("All.Pairs-Table 1.csv", header=TRUE, row.names=1)

## Keep pairs that are either literature supported or putative but filter out those annotated as being incorrect
ligand.receptor.pairs = ligand.receptor.pairs[ligand.receptor.pairs[ ,"Pair.Evidence"] %in% c("literature supported", "putative"), ]
receptors = as.character(ligand.receptor.pairs[, 3])
ligands = as.character(ligand.receptor.pairs[, 1])
all.genes = c(receptors, ligands)

## Map the human gene names to mouse
ligand.receptor.mappings = mappings[mappings[, 1] %in% all.genes, ]
mouse.names = rownames(tip.aggr@data)
ligand.receptor.mappings = ligand.receptor.mappings[as.character(ligand.receptor.mappings[, 2]) %in% mouse.names, ]

### Remove non-unique human -> mouse mappings
counts = table(ligand.receptor.mappings[, 1])
counts = counts[counts==1]
ligand.receptor.mappings = ligand.receptor.mappings[as.character(ligand.receptor.mappings[, 1]) %in% names(counts), ]
rownames(ligand.receptor.mappings) = ligand.receptor.mappings[ , 1]

### Overlap pairs with the genes expressed in the data
ligand.receptor.pairs = ligand.receptor.pairs[, c(1, 3)]
ligand.receptor.pairs.expressed = ligand.receptor.pairs[as.character(ligand.receptor.pairs[, 1]) %in% rownames(ligand.receptor.mappings) 
                                                        & as.character(ligand.receptor.pairs[, 2]) %in% rownames(ligand.receptor.mappings), ]

### Produce a table of weighted edges for the following:
### Cluster -> Ligand
### Ligand -> Receptor
### Receptor -> Cluster
populations.use = clusters.use

unique.genes = unique(c(as.character(ligand.receptor.pairs.expressed[,1]), as.character(ligand.receptor.pairs.expressed[, 2])))
unique.mouse.genes = unlist(lapply(unique.genes, function(x) as.character(ligand.receptor.mappings[x, 2])))
de.results = calculateClusterSpecificExpression(unique.mouse.genes, tip.aggr, threshold=0.1, cluster.set = clusters.use)

## Build a table for each ligand-receptor pair as expressed in each cluster
all.de.results = getGenePairExpressionValues(ligand.receptor.pairs.expressed, mapping.table = ligand.receptor.mappings, 
                                             cluster.names = clusters.use, gene.de.results = de.results)
dim(all.de.results)


### Pull out connections for clusters of interest and add weights
indicies = all.de.results$Cluster1 %in% clusters.use & all.de.results$Cluster2 %in% clusters.use

ligand.receptor.edges = unique(all.de.results[indicies , c("Gene1", "Gene2")])
colnames(ligand.receptor.edges) = c("source", "target")
pair.identifiers = as.character(apply(ligand.receptor.edges, 1, function(x) paste0(x[1], ".", x[2])))
rownames(ligand.receptor.edges) = pair.identifiers

## Get weights using the STRING data-base
ligands = as.character(ligand.receptor.edges[, 1])
receptors = as.character(ligand.receptor.edges[, 2])

### Here use the STRING data-base to give mouse-specific scores to ligand-receptor relationships
lr_score_table = get_STRING_table(ligands, receptors)
head(lr_score_table) ## print out some interactions

overlapping.genes = intersect(pair.identifiers, rownames(lr_score_table))
print(paste0(length(overlapping.genes), " overlaps between ligand-receptor map and STRINGdb"))

weights = lr_score_table[overlapping.genes, ]$Combined_score
weights = weights/1000

ligand.receptor.edges.overlap = ligand.receptor.edges[overlapping.genes, ]
ligand.receptor.edges.overlap$weight = weights
ligand.receptor.edges.overlap$relationship = "ligand.receptor"
ligands = as.character(ligand.receptor.edges[, 1])
receptors = as.character(ligand.receptor.edges[, 2])

## Mapping 'source' population to corresponding ligand
cluster.ligand.edges = unique(all.de.results[indicies , c("Cluster1", "Gene1", "Gene1.Log2_fold_change")])
colnames(cluster.ligand.edges) = c("source", "target", "weight")
cluster.ligand.edges$relationship = "cluster.ligand"
cluster.ligand.edges$source = paste0("S:", cluster.ligand.edges$source)

## Mapping 'target' population to corresponding receptor
receptor.cluster.edges = unique(all.de.results[indicies , c("Gene2", "Cluster2", "Gene2.Log2_fold_change")])
colnames(receptor.cluster.edges) = c("source", "target", "weight")
receptor.cluster.edges$relationship = "receptor.cluster"
receptor.cluster.edges$target = paste0("T:", receptor.cluster.edges$target)

## Build a table of edges,
col.order =  c("source", "target", "relationship", "weight")
ligand.receptor.edges.overlap = ligand.receptor.edges.overlap[, col.order]
receptor.cluster.edges = receptor.cluster.edges[, col.order]
cluster.ligand.edges = cluster.ligand.edges[, col.order]

all.edges = rbind(as.matrix(ligand.receptor.edges.overlap), trimws(as.matrix(receptor.cluster.edges)), trimws(as.matrix(cluster.ligand.edges)))

### Write the network edges to file
write.csv(all.edges, file = paste0(out.lab, "_all_ligand_receptor_network_edges.csv"), row.names = FALSE)

## Also write out just the weights based on expression values
expression.edges = rbind(trimws(as.matrix(receptor.cluster.edges)), trimws(as.matrix(cluster.ligand.edges)))
write.csv(expression.edges, file = paste0(out.lab, "_ligand_receptor_weights.csv"),  row.names = FALSE)

##############################################################################
### Given a set of weighted edges, build a network and find weighted paths ###
### between query source nodes and all other nodes                         ###
##############################################################################

edge.score.file = paste0(out.lab, "_all_ligand_receptor_network_edges.csv")
score.table = read.csv(edge.score.file)

all.weights = score.table$weight
all.edges = score.table[, c("source", "target")]

populations.test = names(table(tip.aggr@ident))

## Go through and calculate summed path weights from source to target populations
## Minimum weight of 1.5 to select for paths with some up-regulation (in ligand, receptor or both)
complete.path.table = c()
for (s.pop in populations.test) {
  for (t.pop in populations.test) {
    source.population = paste0("S:", s.pop)
    target.population = paste0("T:", t.pop)
    this.path.table = getWeightedPaths(score.table, source.population = source.population, 
                                       target.population = target.population, min.weight = 1.5)
    complete.path.table = rbind(complete.path.table, this.path.table)
  }
}
dim(complete.path.table)

write.csv(complete.path.table, file = paste0(out.lab, "_network_paths_weight1.5.csv"))

## Generate a background set of paths for permutation testing 
## Set mimimum weight to -100 to capture all paths
background.path.table = c()
for (s.pop in populations.test) {
  for (t.pop in populations.test) {
    source.population = paste0("S:", s.pop)
    target.population = paste0("T:", t.pop)
    this.path.table = getWeightedPaths(score.table, source.population = source.population, 
                                       target.population = target.population, min.weight = -100)
    background.path.table = rbind(background.path.table, this.path.table)
  }
}
dim(background.path.table)

write.csv(background.path.table, file = paste0(out.lab, "_background_paths.csv"))

#############################################################################################
### Here we do permutation testing to identify the most significant cell-cell connections ###
#############################################################################################

## Read in the individual weights for the edges in the network
weights.file = paste0(out.lab, "_all_ligand_receptor_network_edges.csv")
weights.table = read.csv(weights.file)

ligand.receptor.table = weights.table[weights.table$relationship == "ligand.receptor", ]
rownames(ligand.receptor.table) = paste0(ligand.receptor.table$source, "_", ligand.receptor.table$target)

receptor.cluster.table = weights.table[weights.table$relationship == "receptor.cluster", ]
rownames(receptor.cluster.table) = paste0(receptor.cluster.table$source, "_", receptor.cluster.table$target)

cluster.ligand.table = weights.table[weights.table$relationship == "cluster.ligand", ]
rownames(cluster.ligand.table) = paste0(cluster.ligand.table$source, "_", cluster.ligand.table$target)

## Read in the background network file
completeFile = paste0(out.lab, "_background_paths.csv")
background.table = read.csv(completeFile, row.names = 1)
dim(background.table)

## Read in the filtered network file
thisFile = paste0(out.lab, "_network_paths_weight1.5.csv")
complete.path.table = read.csv(thisFile, row.names = 1)

## Permutation testing for determing signifcant cell-cell connections
## Iterate through each combination of populations and do random selections 
## of fold changes for ligands and receptors. Calculate empirical P-value.
pvalue.table = c()
for(s.pop in populations.test) {
  this.source = paste0("S:", s.pop)
  for (t.pop in populations.test) {
    this.target = paste0("T:", t.pop)
    
    ## Add weights together
    s.indicies = which(complete.path.table$Source == this.source)
    t.indicies = which(complete.path.table$Target == this.target)
    sub.table = complete.path.table[intersect(s.indicies, t.indicies), ]
    
    num.paths = nrow(sub.table)
    path.sum = sum(sub.table$Weight)
    
    ## Get number of total paths from background table
    s.indicies = which(background.table$Source == this.source)
    t.indicies = which(background.table$Target == this.target)
    sub.table = background.table[intersect(s.indicies, t.indicies), ]
    
    num_total = nrow(sub.table)
    
    ## Get weights from the background table
    s.indicies.background = which(background.table$Source == this.source)
    t.indicies.background = which(background.table$Target == this.target)
    combined.indicies = intersect(s.indicies.background, t.indicies.background)
    sub.background.table = background.table[combined.indicies, ]
    ligand.receptor.set = paste0(sub.background.table$Ligand, "_", sub.background.table$Receptor)
    
    ## Now get the real weights of ligand-receptor connections
    ligand.receptor.sub.table = ligand.receptor.table[ligand.receptor.set, ]
    ppi.weights = ligand.receptor.sub.table$weight
    
    random.paths = rep(NA, num_total)
    random.weight.sums = rep(NA, num_total)
    for (i in 1:100000) {
      random.weights = getRandomFCWeights(ppi.weights, cluster.ligand.table, receptor.cluster.table)
      random.weights = random.weights[random.weights >= 1.5]
      
      this.random.weight.sum = sum(random.weights)
      this.random.path.count = length(random.weights)
      
      random.paths[i] = this.random.path.count
      random.weight.sums[i] = this.random.weight.sum
      
    }
    print(paste0("Source ", s.pop, " to ", t.pop, ": "))
    p.paths = sum(random.paths >= num.paths)/length(random.paths)
    p.sum = sum(random.weight.sums >= path.sum)/length(random.weight.sums)
    
    print(paste0("Num paths P-value = ", p.paths))
    print(paste0("Path weights P-value = ", p.sum))
    print("------------------------------------")
    
    thisLine = data.frame(Source_population = s.pop, Target_population = t.pop, Num_paths = num.paths,
                          Num_paths_pvalue = p.paths, Sum_path = path.sum, Sum_path_pvalue = p.sum)
    pvalue.table = rbind(pvalue.table, thisLine)
    
  }
}

## Do P-value adjustment for multiple testing
pval.adj = p.adjust(pvalue.table$Sum_path_pvalue, method = "BH")
pvalue.table$Sum_path_padj = pval.adj

## Write test results to file 
out.file = paste0("Permutation_tests_", out.lab, "_network.csv")
write.csv(pvalue.table, file = out.file, row.names = FALSE)

##############################################################################################
### Generate a few plots visualising the interactions.                                     ###
### The table of cell-cell interactions can also be read in to Cytoscape for visualisation ###
##############################################################################################

path.sig.file = paste0("Permutation_tests_", out.lab, "_network.csv")

sig.paths = read.csv(path.sig.file, sep = ",")
sig.paths$Source_population = sub(".:(.*)", "\\1", sig.paths$Source_population)
sig.paths$Target_population = sub(".:(.*)", "\\1", sig.paths$Target_population)
sig.paths = sig.paths[which(sig.paths$Sum_path_padj < 0.01), ]
nrow(sig.paths)
rownames(sig.paths) = paste0(sig.paths$Source_population, "_", sig.paths$Target_population)

## sub-set the collapsed weight table by the edges found to be significant
collapsed.path.table = sig.paths[, c("Source_population", "Target_population", "Sum_path")]
colnames(collapsed.path.table) = c("Source", "Target", "Weight")

## Define population colours for plotting
col.set <- c("#c10023", "#008e17", "#fb8500", "#f60000", "#fde800", "#bc9000","#4ffc00", "#00bcac", "#0099cc",
             "#D35400", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8",
             "#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")

names(col.set) =  c("M1M\u03A6", "EC1", "F-SL", "F-Act", "F-SH", "BC", "M2M\u03A6", "M1Mo", "MYO",
                    "EC3", "DC", "EC2", "TC1-Cd8", "TC2-Cd4", "MAC-TR", "Mural", "Cyc",
                    "MAC6", "MAC-IFNIC", "MAC7", "MAC8", "F-WntX", "NKC", "Glial")

## put edges in format for creating a graph
all.weights = collapsed.path.table$Weight
all.edges = collapsed.path.table[, c("Source", "Target")]
edge.sources = all.edges$Source
edges.list = unlist(lapply(t(all.edges), function(x) c(x[1])))

## Make directed graph and add edge weights
lr.plot <- make_graph(edges.list, directed = TRUE)
E(lr.plot)$weight <- all.weights

#population.cluster.map = names(new.ident)
#names(population.cluster.map) = as.character(new.ident)

populations.use = names(V(lr.plot))
col.set = col.set[populations.use]
cluster.colors = data.frame(cluster=populations.use, color=col.set)
rownames(cluster.colors) = cluster.colors$cluster

## Generate a color pallete for the edges
colfunc <- colorRampPalette(c("grey", "black"))

vertex.colors = c()
for (this.vertex in V(lr.plot)$name) {
  if (this.vertex %in% cluster.colors$cluster) {
    this.col = as.character(cluster.colors[cluster.colors$cluster==this.vertex, ]$color)
    vertex.colors = append(vertex.colors, this.col)
  } else{
    vertex.colors = append(vertex.colors, "#e6e6e6")
  }
}

edge.colours = c()
for (this.source in edge.sources) {
  if (this.source %in% cluster.colors$cluster) {
    this.col = as.character(cluster.colors[cluster.colors$cluster==this.source, ]$color)
    edge.colours = append(edge.colours, this.col)
  } else{
    edge.colours = append(edge.colours, "#e6e6e6")
  }
}

arrow.size = 0.6
arrow.width = 2.0
edge.multi = 0.02
E(lr.plot)$width <- as.numeric(all.weights)*edge.multi

layout = "layout_in_circle"
l <- do.call(layout, list(lr.plot))
par(mar=c(0,0,0,0)+.1)
plot(lr.plot, edge.curved = 0.2, vertex.color=vertex.colors,
     layout=l, vertex.label.cex=1.1, edge.arrow.size=arrow.size, edge.arrow.width=arrow.width, vertex.size=18, 
     vertex.label.font=2, vertex.label.color="black", edge.color = edge.colours)

#################################################################
### Generate a figure showing inbound vs outbound connections

p.labels <- c("M1M\u03A6", "EC1", "F-SL", "F-Act", "F-SH", "BC", "M2M\u03A6", "M1Mo", "MYO",
              "EC3", "DC", "EC2", "TC1-Cd8", "TC2-Cd4", "MAC-TR", "Mural", "Cyc",
              "MAC6", "MAC-IFNIC", "MAC7", "MAC8", "F-WntX", "NKC", "Glial")

col.set <- c("#c10023", "#008e17", "#fb8500", "#f60000", "#fde800", "#bc9000","#4ffc00", "#00bcac", "#0099cc",
             "#D35400", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8",
             "#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")

edge.sources = factor(all.edges$Source, levels = p.labels)
edge.targets = factor(all.edges$Target, levels = p.labels)

source.counts = table(edge.sources)
target.counts = table(edge.targets)

outgoing.weights = c()
incoming.weights = c()
for (this.pop in p.labels) {
  outgoing = sum(sig.paths[which(sig.paths$Source_population == this.pop), ]$Sum_path)
  outgoing.weights = append(outgoing.weights, outgoing)
  
  incoming = sum(sig.paths[which(sig.paths$Target_population == this.pop), ]$Sum_path)
  incoming.weights = append(incoming.weights, incoming)
}

ggData = data.frame(Population = factor(p.labels, levels = p.labels), Outbound = as.numeric(source.counts), 
                    WeightsOut = outgoing.weights, Inbound = as.numeric(target.counts), WeightsIn = incoming.weights)
ggData$Population = factor(ggData$Population, levels = p.labels)
ggplot(ggData, aes(x = WeightsIn, y = WeightsOut, colour = Population)) + geom_point(size = 5) +
  scale_colour_manual(values = col.set) + theme_bw(base_size = 15) + ylab("Total outgoing weights") +
  xlab("Total incoming weights") + geom_text(data=subset(ggData, (WeightsOut > 1000 | WeightsIn > 1000)), 
                                             aes(x=WeightsIn, y = WeightsOut, label=Population), hjust=0.4, vjust=-0.6, colour = "black") +
  theme(legend.position = c(0.7, 0.75)) + guides(color = guide_legend(override.aes = list(size=4), ncol = 3))



