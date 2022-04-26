
# SNN Script

snn <- function(seurat, k.param = 30){
  
  print(paste0("Computing for: ", paste(names(seurat), collapse = "-")))
  
  set.seed(0305)
  
  seuratsnn <- FindNeighbors(seurat, k.param = k.param)
  
  return(seuratsnn)
  
}

# cluster Script

cluster <- function(seurat, resolution = 0.8){
  
  print(paste0("Computing for: ", paste(names(seurat), collapse = "-")))
  
  set.seed(0305)
  
  seuratclust <- FindClusters(seurat, resolution = resolution)
  
  return(seuratclust)
  
}


snnCluster <- function(seurat, k.param = NULL, resolution = NULL, clusters = "seurat_clusters") {
  
  print(paste0("Computing for: ", paste(names(seurat), collapse = "-")))
  
  set.seed(0305)
  
  if(is.null(k.param) & is.null(resolution)){
    
    print("k.param: 20")
    print("resolution: 0.8")
    k.param <- 20
    resolution <- 0.8
    
  } else {
    print(paste0("k.param: ", k.param))
    print(paste0("resolution: ", resolution))
  }
  
  seuratsnnCluster <- FindNeighbors(seurat, k.param = k.param)
  seuratsnnCluster <- FindClusters(seuratsnnCluster, resolution = resolution)
  
  seuratsnnCluster[[paste0("res.", resolution)]] <- seuratsnnCluster[[clusters]]
  
  return(seuratsnnCluster)  
  
  
}

# addClusAnno <- function(seurat, prefix = "seurat_clusters", new = "res.0.6"){
#   
#   seurat[[new]] <- seurat[[prefix]]
#   
#   return(seurat)
# }

clustreeRun <- function(seurat, prefix = "res."){
  
  ctree <- clustree(seurat, prefix = prefix)
  
  return(ctree)
}

silhouette <- function(seurat, slot = "data", clustername = "res.0.4" ,
                       clusters = NULL, sample = NULL){
  
  # for information on slots https://www.biostars.org/p/9484293/
  
  mat <- t(as.matrix(GetAssayData(object = seurat, slot = slot)))
  
  if(is.null(clusters)){
    clusters <- seurat[[clustername]]
  }
  
  stopifnot(rownames(clusters) == rownames(mat))
  stopifnot(length(clusters[,1]) == nrow(mat))
  
  sil <- bluster::approxSilhouette(mat, clusters[,1])
  
  if(is.null(sample)){
    sil$Sample <- rep(unique(seurat$Sample), nrow(sil))
  } else {
    sil$Sample <- rep(sample, nrow(sil))
  }

  return(sil)
}

plotSil <- function(sil, clusters = NULL, type = "boxplot", title = "res 0.2", 
                    color = NULL){
  
  if (type == "boxplot"){
    
    bp <- ggplot(as.data.frame(sil), aes(x=cluster, y=width, fill=cluster)) + 
      geom_boxplot()+
      labs(title=paste0(unique(sil$Sample)," ",title),x="Cluster", y = "Silhouette") + theme_classic()
    
    return(bp)
    
  } 
  
  if (type == "heatmap"){
    
    if(!is.null(color)){
      color = colorRampPalette(rev(brewer.pal(n = 7, name = color)))(100)
    } else {
      color =  colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)
    }
    
    best.choice <- ifelse(sil$width > 0, sil$cluster, sil$other)
    
    pheatmap(table(Assigned=sil$cluster, Closest=best.choice), cluster_cols = F,
             cluster_rows = F, color = color)
  }

  
}


purity <- function(seurat, slot = "data", clustername = "res.0.4", clusters = NULL){
  
  # for information on slots https://www.biostars.org/p/9484293/
  
  mat <- t(as.matrix(GetAssayData(object = seurat, slot = slot)))

  if(is.null(clusters)){
    clusters <- seurat[[clustername]]
  }
  
  stopifnot(rownames(clusters) == rownames(mat))
  stopifnot(length(clusters[,1]) == nrow(mat))
  
  pur <- bluster::neighborPurity(mat, clusters[,1])
  
  return(pur)
}

plotPurity <- function(pur, clusters){

  if (type == "boxplot"){
    
    bp <- ggplot(as.data.frame(pur), aes(x=cluster, y=purity, fill=cluster)) + 
      geom_boxplot()+
      labs(title=paste0("Purity ", title),x="Cluster", y = "Purity") + theme_classic()
    
    return(bp)
    
  } 
  
  if (type == "heatmap"){
    
    if(!is.null(color)){
      color = colorRampPalette(rev(brewer.pal(n = 7, name = color)))(100)
    } else {
      color =  colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)
    }
    
    #table(Assigned=clusters, Max=pure$maximum)
    
    pheatmap(table(Assigned=pur$cluster, Max=pur$maximum), cluster_cols = F,
             cluster_rows = F, color = color)
  }
  
}

## does not work with seurat graph output - g needs to be in igraph format
graphMod <- function(seurat, clusters = NULL,clustername = "res.0.4", g = NULL,
                     graph = "SCT_nn"){
  
  if(is.null(g)){
    g <- seurat@graphs[[graph]]
  }
  
  if(is.null(clusters)){
    clusters <- seurat[[clustername]]
  }
  print(head(g))
  ratio <- bluster::pairwiseModularity(g, clusters, as.ratio=TRUE)
  
  p <- pheatmap(log10(ratio+1), cluster_cols=FALSE, cluster_rows=FALSE,
           col=rev(heat.colors(100)))
  return(p)
}

addLabel <- function(seurat, chosenCluster = "res.0.2"){
  
  seurat[["Cluster"]] <- seurat[[chosenCluster]]
  
  return(seurat)
}
