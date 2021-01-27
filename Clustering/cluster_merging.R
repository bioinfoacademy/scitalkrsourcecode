########################################################
### Scripts for iteratively merging similar clusters ###
### Used for merging clusters in GFP day 7 data-set  ###
########################################################

generateClusterMerging <- function(seurat.object, num_iter, gene.set="variable", threshold=0.25, expression.threshold=0) {
  # Get genes expressed in the foreground set
  gene.proportions.foreground = apply(seurat.object@data, 1, function(x)
  {return(sum(x > expression.threshold)/length(x))})
  genes.foreground = names(gene.proportions.foreground[which(gene.proportions.foreground > threshold)])
  
  # Decide the gene set to use
  if (gene.set == "high") {
    geneSet = genes.foreground
  } else if (gene.set == "variable") {
    geneSet = seurat.object@var.genes
  } else{
    geneSet = rownames(seurat.object@data)
  }
  
  aggr.merged = seurat.object
  # Build the classification hierarchy using all genes for interpretation
  aggr.merged <- BuildClusterTree(object = aggr.merged, genes.use = geneSet, 
                                  do.reorder = FALSE, reorder.numeric = FALSE)
  
  for (i in c(1:num_iter)) {
    
    node.scores <- AssessNodes(object = aggr.merged)
    node.scores[order(node.scores$oobe, decreasing = TRUE), ] -> node.scores
    
    # Identify the top node connecting leaf nodes. 
    node.merge <- getTopLeafyNode(aggr.merged, node.scores)
    aggr.merged <- MergeNode(object = aggr.merged, node.use = node.merge$node)
    
    # Rebuild the classification hierarchy using all genes for interpretation (you can also try with variable genes as the default)
    aggr.merged <- BuildClusterTree(object = aggr.merged, do.reorder = TRUE, reorder.numeric = TRUE, 
                                    genes.use = geneSet)
    
    #Generate new TSNE plot
    TSNEPlot(object = aggr.merged, do.label = TRUE)
  }
  
  # Merging is complete - now to re-order the clusters based on size
  cluster.idents = table(aggr.merged@ident)
  idents.ordered = order(cluster.idents, decreasing=T)
  ident.replacements = c(1:length(idents.ordered))
  names(ident.replacements) <- idents.ordered
  
  clusters = aggr.merged@ident
  clusters=revalue(clusters, ident.replacements)
  clusters <- factor(clusters, levels = ident.replacements)
  aggr.merged@ident = clusters
  TSNEPlot(object = aggr.merged)
  return(aggr.merged)
}

################################################################
### Identify nodes representing 'leaves' on the cluster tree ###
################################################################
getTopLeafyNode <- function(seurat.object, node.scores) {
  ct <- seurat.object@cluster.tree
  cluster.idents = names(table(seurat.object@ident))
  # First pull out the edges only connected to clusters
  cluster.edges = ct[[1]]$edge[ct[[1]]$edge[, 2] %in% cluster.idents, ]
  
  # Nodes that are repeated will be connected to two clusters, representing 'leaf' nodes.
  node.counts = plyr::count(cluster.edges[, 1])
  leafy.connecting.nodes = node.counts[plyr::count(cluster.edges[, 1])$freq > 1, ]$x
  
  # Pull out the node scores for nodes connecting the leafy cluster nodes
  leafy.node.scores = node.scores[node.scores$node %in% leafy.connecting.nodes, ]
  top.node = leafy.node.scores[1, ]
  return(top.node)
}

