

##########################################################################
### Calculate marker (differentially expressed) genes for each cluster ###
##########################################################################
calculateClusterSpecificExpression<-function(geneList, seurat.object, exp.threshold=0, 
                                             threshold=0.5, cluster.set = NULL) {
  geneList = unique(geneList)
  results = data.frame()
  expression.matrix = seurat.object@data[geneList, ]
  pb <- progress_bar$new(format = "Processing [:bar] :percent eta: :eta", 
                         total = length(names(table(seurat.object@ident))), clear=FALSE) 
  if (is.null(cluster.set)) {
    cluster.set = names(table(seurat.object@ident))
  }
  for (i in cluster.set) {
    
    # Get the set of cells for the target cluster
    cluster = i
    cell.set = names(seurat.object@ident[seurat.object@ident==cluster])
    remainder.set = names(seurat.object@ident[seurat.object@ident!=cluster])
    
    # Get genes expressed in the foreground set
    gene.proportions.cluster = apply(expression.matrix[, cell.set], 1, function(x)
    {return(sum(x > exp.threshold)/length(x))})
    names(gene.proportions.cluster) = rownames(expression.matrix)
    genes.cluster = names(gene.proportions.cluster[which(gene.proportions.cluster > threshold)])
    exp.percentages = gene.proportions.cluster[genes.cluster]
    
    # Calculate average log2 fold-changes for the genes in the cluster
    log2.fc = getLog2FCList(seurat.object, genes.cluster, cluster)
    
    # Calculate averaged Log2 expression for genes in the cluster
    aveLog2Expression = getLog2AvereragedRNA(seurat.object = seurat.object, thisCluster = cluster, c.id.used = TRUE)
    aveLog2Expression = aveLog2Expression[genes.cluster]
    
    # build a data-frame for the results
    this.result = data.frame(Cluster = rep(i, length(genes.cluster)), Gene=genes.cluster, 
                             Ave_Log2_exp = aveLog2Expression, Pct_expressed = exp.percentages, Log2FC = log2.fc)
    rownames(this.result) = paste0(i, ".", genes.cluster)
    results = rbind(results, this.result)
    pb$tick()
  }
  
  return(results)
}

#########################################################################################
### For each gene pair, test if one of them is expressed above threshold in a cluster ###
#########################################################################################
getGenePairExpressionValues <- function(ligand.receptor.pairs, mapping.table, cluster.names, gene.de.results) {
  
  all.de.results = data.frame()
  pb <- progress_bar$new(format = "Processing [:bar] :percent eta: :eta", 
                         total = nrow(ligand.receptor.pairs), clear=FALSE) 
  for (n in 1:nrow(ligand.receptor.pairs)){
    ligand = as.character(mapping.table[as.character(ligand.receptor.pairs[n, 1]), 2])
    receptor = as.character(mapping.table[as.character(ligand.receptor.pairs[n, 2]), 2])
    
    results = data.frame()
    for (i in cluster.names) {
      
      cluster1 = i
      gene1.label = paste0(i,".",ligand)
      for (j in cluster.names) {
        
        gene2.label = paste0(j, ".", receptor)
        
        ### get log2 fold-change and P-value
        aveLog2Exp = gene.de.results[rownames(gene.de.results) %in% c(gene1.label, gene2.label), ]$Ave_Log2_exp
        if (length(aveLog2Exp) == 2) {
          #res = gene.de.results[c(gene1.label, gene2.label), ]$Pvalue
          
          this.result = data.frame(Cluster1=i, Gene1=ligand, Gene1.Ave_Log2_exp=gene.de.results[gene1.label, ]$Ave_Log2_exp,
                                   Gene1.Log2_fold_change=gene.de.results[gene1.label, ]$Log2FC, 
                                   Gene1.Pct_expressed = gene.de.results[gene1.label, ]$Pct_expressed, Cluster2=j,
                                   Gene2=receptor, Gene2.Ave_Log2_exp=gene.de.results[gene2.label, ]$Ave_Log2_exp,
                                   Gene2.Log2_fold_change=gene.de.results[gene2.label, ]$Log2FC, 
                                   Gene2.Pct_expressed = gene.de.results[gene2.label, ]$Pct_expressed
          )
          results = rbind(results, this.result)
          
        }
      }
    }
    
    all.de.results = rbind(all.de.results, results)
    pb$tick()
  }
  
  return(all.de.results)
}


#################################################################################
### Get STRINGdb scores: given a list of ligands and list of receptors, pulls ###
### out associations from STRINGdb and builds a ligand-receptor scoring table ###
#################################################################################
get_STRING_table <- function(ligands, receptors) {
  ## Gene names that aren't automatically mapped to STRING and need to be mapped to alternative identifier
  ## Ackr3 -> Cxcr7
  string.chromium.map = c("Ackr3")
  names(string.chromium.map) = c("Cxcr7")
  
  string_db <- STRINGdb$new( version="10", species=10090, score_threshold=0, input_directory="STRINGdb")
  
  ## For this analysis only one gene that needs to be mapped
  gene.table = data.frame(Gene = as.character(unique(ligands)), Class = "ligand")
  gene.table = rbind(gene.table, (data.frame(Gene = as.character(unique(receptors)), Class = "receptor")))
  gene.table = rbind(gene.table, (data.frame(Gene = "Cxcr7", Class = "receptor")))
  
  gene.table.mapped <- string_db$map(gene.table, "Gene", removeUnmappedRows = TRUE )
  rownames(gene.table.mapped) = gene.table.mapped$STRING_id
  
  ## Get the list of STRING identifiers
  hits <- gene.table.mapped$STRING_id
  
  ## Interaction table for identifier list
  gene.interactions.table = string_db$get_interactions(hits)
  print(paste0(nrow(gene.interactions.table), " interactions identified from STRINGdb"))
  
  ## gene interactions table are using STRING ID so need to convert back to gene name
  source.names = unlist(lapply(gene.interactions.table$from, function(x) gene.table.mapped[x, ]$Gene))
  target.names = unlist(lapply(gene.interactions.table$to, function(x) gene.table.mapped[x, ]$Gene))
  
  gene.interactions.table$from.gene = source.names
  gene.interactions.table$to.gene = target.names
  
  gene.interactions.table.subset = gene.interactions.table[, c("from.gene", "to.gene", "combined_score")]
  
  ## Replace STRING gene names with original gene names in the from column
  replace.index = gene.interactions.table.subset$from.gene %in% names(string.chromium.map)
  if (sum(replace.index) > 0) {
    gene.interactions.table.subset[replace.index, 1] = as.character(string.chromium.map[gene.interactions.table.subset[replace.index, 1]])
  }
  
  ## Replace STRING gene names with original gene names in the to column
  replace.index = gene.interactions.table.subset$to.gene %in% names(string.chromium.map)
  if (sum(replace.index) > 0) {
    gene.interactions.table.subset[replace.index, 2] = as.character(string.chromium.map[gene.interactions.table.subset[replace.index, 1]])
  }
  
  ## Now identify all ligand-receptor interactions from the input table
  lr_score_table = data.frame()
  
  ## re-order the to-from relationships to indicate ligand-receptor direction
  for (i in c(1:nrow(gene.interactions.table.subset))) {
    gene1 = gene.interactions.table.subset[i, 1]
    gene2 = gene.interactions.table.subset[i, 2]
    score = gene.interactions.table.subset[i, 3]
    if ((gene1 %in% ligands) & (gene2 %in% receptors)) {
      ## Keep order
      this.row = data.frame(Ligand = gene1, Receptor = gene2, Combined_score = score)
      rownames(this.row) = paste0(gene1, ".", gene2)
      lr_score_table = rbind(lr_score_table, this.row)
    } else if ((gene2 %in% ligands) & (gene1 %in% receptors)) {
      ## Reverse the order
      this.row = data.frame(Ligand = gene2, Receptor = gene1, Combined_score = score)
      rownames(this.row) = paste0(gene2, ".", gene1)
      lr_score_table = rbind(lr_score_table, this.row)
    } ## Ignore other combinations
  }
  nrow(lr_score_table)
  head(lr_score_table)
  
  return(lr_score_table)
}

##############################################################################
### Given a set of weighted edges, build a network and find weighted paths ###
### between query source nodes and all other nodes                         ###
##############################################################################
getWeightedPaths <- function(score.table, source.population, target.population, min.weight = 2, print.num.connections = F) {
  ## put edges in format for creating a graph
  all.weights = score.table$weight
  all.edges = score.table[, c("source", "target")]
  
  edges.list = unlist(lapply(t(all.edges), function(x) c(x[1])))
  
  ## Make directed graph and add edge weights
  lr.plot <- make_graph(edges.list, directed = TRUE)
  E(lr.plot)$weight <- all.weights
  
  ## Calculate all shortest paths (ignoring weight) between source and target node
  all.paths = all_shortest_paths(lr.plot, from = source.population, to = target.population, mode = "out", weights = NA)
  
  ## Calculate sum of weights for all shortest paths
  path.weight.table = matrix(0, ncol = 5, nrow = length(all.paths$res))
  colnames(path.weight.table) = c("Source", "Ligand", "Receptor", "Target", "Weight")
  path.weight.table = as.data.frame(path.weight.table)
  for (i in c(1:length(all.paths$res))) {
    node.names = names(unlist(all.paths$res[i]))
    this.epath = E(lr.plot, path=unlist(all.paths$res[i]))
    this.weight = sum(E(lr.plot)$weight[this.epath])
    this.elem = c(node.names[1], node.names[2], node.names[3], node.names[4], this.weight)
    path.weight.table[i, ] = this.elem
  }
  path.weight.table$Weight = as.numeric(path.weight.table$Weight)
  
  ## Filter paths by a minimum weight
  path.weight.table = path.weight.table[path.weight.table$Weight >= min.weight, ]
  if (print.num.connections) {
    print(paste0(nrow(path.weight.table), " paths for ", source.population, " to ", target.population))
  }
  
  return(path.weight.table)
}

####################################################################################
### Given a set of real ligand-receptor connections and weights, randomly        ###
### select fold-changes from the cluster-ligand and receptor-cluster connections ###
####################################################################################
getRandomFCWeights <- function(ppi.weights, cluster.ligand.table, receptor.cluster.table) {
  
  sample.size = length(ppi.weights)
  
  ## Get the set of randomly selected ligands
  rand.weights1 = sample(cluster.ligand.table$weight, size = sample.size, replace = FALSE)
  
  ## Get the set of randomly selected receptors
  rand.weights2 = sample(receptor.cluster.table$weight, size = sample.size, replace = FALSE)
  
  ## Add the weights together and return
  combined.weights = rand.weights1 + ppi.weights + rand.weights2
  
  return(combined.weights)
}

##################################
### Calculate Log2 fold-change ###
##################################
getLog2FC = function(x, y) {
  ## Assumes data has previous been processed as log(x+1)
  return(log2(mean(exp(x))/mean(exp(y))))
}

#########################################################################################
### Get a list of Log2 fold-change values for a gene list between a specified cluster ###  
### and either all remaining cells (default) or an alternative cluster                ###
#########################################################################################
getLog2FCList <- function(seurat.object, geneList, cluster1, cluster2=NULL) {
  
  ### Need to check whether the input is a cluster label or vector of cell names
  if (length(cluster1) == 1){ # cluster identity used as input
    foreground.set = names(seurat.object@ident[seurat.object@ident==cluster1])
  } else { # cell identity used as input
    foreground.set = cluster1
  }
  if (is.null(cluster2)) {
    remainder.set = setdiff(seurat.object@cell.names, foreground.set)
  } else {
    if (length(cluster2) == 1) { # cluster identity used as input
      remainder.set = names(seurat.object@ident[seurat.object@ident==cluster2])
    } else { # cell identity used as input
      remainder.set = cluster2
    }
  }
  log2.fc.values = apply(seurat.object@data[geneList, c(foreground.set, remainder.set)], 1, function(x) 
    getLog2FC(x[foreground.set], x[remainder.set]))
  return(log2.fc.values)
}

Log2ExpMean <- function (x) {
  return(log2(x = mean(x = exp(x = x) - 1) + 1))
}

####################################################################
### return an average of the log2-transformed expression data ###
####################################################################
getLog2AvereragedRNA <- function(seurat.object, thisCluster, c.id.used = FALSE){
  
  if (length(thisCluster) == 0) {
    detection.rate = rep(0, nrow(seurat.object@data))
    names(detection.rate) = rownames(seurat.object@data)
    return(detection.rate)
  }
  
  if (c.id.used == TRUE){ # cluster identity used as input
    cell.set = names(seurat.object@ident[seurat.object@ident==thisCluster])
  } else { # cell identity used as input
    cell.set = thisCluster
  }
  
  if (length(cell.set) == 1) {
    return(seurat.object@data[, cell.set])
  }
  
  average.rna = apply(seurat.object@data[, cell.set], 1, function(x) Log2ExpMean(x))
  return(average.rna)
}