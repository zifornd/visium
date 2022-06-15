
#' Wrapper for Seurat FindNeighbors function
#' embeds cells in a graph structure
#' @param seurat Seurat object.
#' @param k.param k.param for k nearest neighbours 
#' @return Seurat object with graph
snn <- function(seurat, k.param = 30) {

  print(paste0("Computing for: ", paste(names(seurat), collapse = "-")))

  set.seed(0305)

  seuratsnn <- Seurat::FindNeighbors(seurat, k.param = k.param)

  return(seuratsnn)

}

#' Wrapper for Seurat FindClusters function
#' Partition the graph into highly connected communities
#' @param seurat Seurat object.
#' @param resolution res param that sets the 'granularity' of the downstream clustering
#' @return Seurat object with clusters
cluster <- function(seurat, resolution = 0.8) {

  print(paste0("Computing for: ", paste(names(seurat), collapse = "-")))

  set.seed(0305)

  seuratclust <- Seurat::FindClusters(seurat, resolution = resolution)

  return(seuratclust)

}

#' Full cluster for seurat object - FindNeighbors and FindClusters
#' Embeds cells in a graph structure
#' Partition the graph into highly connected communities
#' @param seurat Seurat object.
#' @param k.param k.param for k nearest neighbours.
#' @param resolution res param that sets the 'granularity' of the downstream clustering.
#' @param clusters Name of clusters saved into meta data.
#' @return Seurat object with clusters
snnCluster <- function(seurat, k.param = NULL, resolution = NULL,
                       clusters = "seurat_clusters") {

  print(paste0("Computing for: ", paste(names(seurat), collapse = "-")))

  set.seed(0305)

  if (is.null(k.param) & is.null(resolution)) {

    print("k.param: 20")
    print("resolution: 0.8")
    k.param <- 20
    resolution <- 0.8

  } else {

    print(paste0("k.param: ", k.param))

    print(paste0("resolution: ", resolution))
  }

  seuratsnnCluster <- Seurat::FindNeighbors(seurat, k.param = k.param)

  seuratsnnCluster <- Seurat::FindClusters(seuratsnnCluster,
                                           resolution = resolution)

  seuratsnnCluster[[paste0("res.", resolution)]] <- seuratsnnCluster[[clusters]]

  return(seuratsnnCluster)

}

#' Plot clusters overlayed onto tissue image
#' @param seurat Seurat object.
#' @param group Name of group/cluster to plot
#' @param label Whether to plot cluster labels.
#' @param label.size Size of labels for clusters.
#' @return Seurat::SpatialDimPlot
plotSpatialCluster <- function(seurat, group, label = TRUE, label.size = 3) {

  Seurat::Idents(seurat) <- group

  sdimplot <- Seurat::SpatialDimPlot(seurat, label = TRUE, label.size = 3)

  return(sdimplot)
}

#' Wrapper for clustree 
#' Creates a cluster dendrogram acros resolution params
#' @param seurat Seurat object.
#' @param prefix name of cluster prefix - default = res.
#' @return Clustree plot
clustreeRun <- function(seurat, prefix = "res.") {

  ctree <- clustree::clustree(seurat, prefix = prefix)

  return(ctree)
}

#' Run Silhouette scores for clusters using bluster
#' @param seurat Seurat object.
#' @param slot name of cluster prefix - default = res.
#' @param clustername name of cluster to assess.
#' @param clusters user supplied clusters instead of retrieving from seurat obj
#' @param sample sample name - if not supplied - taken from seurat object
#' @return Silhouette scores
silhouette <- function(seurat, slot = "data", clustername = "res.0.4",
                       clusters = NULL, sample = NULL) {

  # for information on slots https://www.biostars.org/p/9484293/

  mat <- t(as.matrix(Seurat::GetAssayData(object = seurat, slot = slot)))

  if (is.null(clusters)) {

    clusters <- seurat[[clustername]]

  }

  stopifnot(rownames(clusters) == rownames(mat))
  stopifnot(length(clusters[, 1]) == nrow(mat))

  sil <- bluster::approxSilhouette(mat, clusters[, 1])

  if (is.null(sample)) {

    sil$Sample <- rep(unique(seurat$Sample), nrow(sil))

  } else {

    sil$Sample <- rep(sample, nrow(sil))

  }

  return(sil)
}

#' Plot Silhouette scores
#' @param sil Silhouette scores from bluster
#' @param type type of plot.
#' @param title title of plot.
#' @param color Color paletter for heatmap
#' @return Ggplot object
plotSil <- function(sil, type = "boxplot", title = "res 0.2",
                    color = NULL) {

  if (type == "boxplot") {

    bp <- ggplot(as.data.frame(sil), aes(x = cluster, y = width,
                                         fill = cluster)) +
          geom_boxplot() +
          labs(title = paste0(unique(sil$Sample), " ", title),
               x = "Cluster", y = "Silhouette") +
          theme_classic()

    return(bp)

  }

  if (type == "heatmap") {

    if (!is.null(color)) {

      pal <- rev(RColorBrewer::brewer.pal(n = 7, name = color))

      color <- colorRampPalette(pal)(100)

    } else {

      pal <- rev(RColorBrewer::brewer.pal(n = 7, name = "YlOrRd"))

      color <-  colorRampPalette(pal)(100)

    }

    best.choice <- ifelse(sil$width > 0, sil$cluster, sil$other)

    pheatmap::pheatmap(table(Assigned = sil$cluster, Closest = best.choice),
                       cluster_cols = F, cluster_rows = F, color = color)
  }

}

#' Run Purity scores for clusters using bluster
#' @param seurat Seurat object.
#' @param slot name of cluster prefix - default = res.
#' @param clustername name of cluster to assess.
#' @param clusters user supplied clusters instead of retrieving from seurat obj
#' @param sample sample name - if not supplied - taken from seurat object
#' @return Purity scores
purity <- function(seurat, slot = "data", clustername = "res.0.4",
                   clusters = NULL, sample = NULL) {

  # for information on slots https://www.biostars.org/p/9484293/

  mat <- t(as.matrix(Seurat::GetAssayData(object = seurat, slot = slot)))

  if (is.null(clusters)) {

    clusters <- seurat[[clustername]]

  }

  stopifnot(rownames(clusters) == rownames(mat))

  stopifnot(length(clusters[, 1]) == nrow(mat))

  pur <- bluster::neighborPurity(mat, clusters[,1])

  if (is.null(sample)) {

    pur$Sample <- rep(unique(seurat$Sample), nrow(pur))

  } else {

    pur$Sample <- rep(sample, nrow(pur))

  }

  return(pur)
}

#' Plot Purity scores
#' @param pur Purity scores from bluster
#' @param type type of plot.
#' @param title title of plot.
#' @param color Color paletter for heatmap
#' @return Ggplot object
plotPurity <- function(pur, type = "boxplot", title = "res 0.2",
                       color = NULL) {

  if (type == "boxplot") {

    bp <- ggplot(as.data.frame(pur), aes(x = cluster, y = purity,
                                         fill = cluster)) +
          geom_boxplot() +
          labs(title = paste0("Purity ", title), x = "Cluster", y = "Purity") +
          theme_classic()

    return(bp)

  }

  if (type == "heatmap") {

    if (!is.null(color)) {

      pal <- rev(RColorBrewer::brewer.pal(n = 7, name = color))

      color <- colorRampPalette(pal)(100)

    } else {

      pal <- rev(RColorBrewer::brewer.pal(n = 7, name = "YlOrRd"))

      color <-  colorRampPalette(pal)(100)

    }

    pheatmap(table(Assigned = pur$cluster, Max = pur$maximum), cluster_cols = F,
             cluster_rows = F, color = color)

  }
}

#' Add Chosen Cluster label to Seurat meta data as "Cluster"
#' @param seurat Seurat object.
#' @param chosenCluster Chosen cluster.
#' @return Seurat object.
addLabel <- function(seurat, chosenCluster = "res.0.2") {

  seurat[["Cluster"]] <- seurat[[chosenCluster]]

  return(seurat)
}
