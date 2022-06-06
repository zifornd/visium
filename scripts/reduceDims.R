
# PCA Script

pca <- function(seurat, assay = "SCT") {

  print(paste0("Computing for: ", paste(names(seurat), collapse = "-")))

  seuratpca <- Seurat::RunPCA(seurat, assay = "SCT", verbose = FALSE)

  return(seuratpca)

}

# UMAP Script

umap <- function(seurat, reduction = "pca", dims = NULL,
                 n.neighbors = NULL, min.dist = NULL) {

  set.seed(0305)

  print(paste0("Computing for: ", paste(names(seurat), collapse = "-")))

  if (is.null(dims)) {

    dims <- as.character(unique(seurat[["Elbow"]]))

    print(paste0("PCs: ", dims))

    dims <- as.numeric(strsplit(dims, ":")[[1]])[1]:as.numeric(strsplit(dims, ":")[[1]])[2]

  } else {

    print(paste0("PCs: ", dims))

  }

  if (is.null(n.neighbors) & is.null(min.dist)) {

    print("n.neighbors: 30")
    print("min.dist: 0.3")
    n.neighbors <- 30
    min.dist <- 0.3

  } else {

    print(paste0("n.neighbors: ", n.neighbors))
    print(paste0("min.dist: ", min.dist))

  }

  seuratumap <- Seurat::RunUMAP(seurat, reduction = reduction, dims = dims,
                                n.neighbors = n.neighbors, min.dist = min.dist)

  return(seuratumap)

}

# TSNE Script

tsne <- function(seurat, reduction = "pca", dims = NULL,
                 perplexity = NULL, max_iter = 10000) {

  set.seed(0305)

  print(paste0("Computing for: ", paste(names(seurat), collapse = "-")))

  if (is.null(dims)) {

    dims <- as.character(unique(seurat[["Elbow"]]))

    print(paste0("PCs: ", dims))

    dims <- as.numeric(strsplit(dims, ":")[[1]])[1]:as.numeric(strsplit(dims, ":")[[1]])[2]

  } else {

    print(paste0("PCs: ", dims))

  }

  if (is.null(perplexity)) {

    print(paste0("perplexity: ", round(ncol(seurat)^ (1 / 2), digits = -1)))

    perplexity <- round(ncol(seurat)^ (1 / 2), digits = -1)

  } else {

    print(paste0("perplexity: ", perplexity))

    if (perplexity == "N^1/2"){

      print(paste0("perplexity: ", round(ncol(seurat)^ (1 / 2), digits = -1)))

      perplexity <- round(ncol(seurat)^ (1 / 2), digits = -1)

    }

  }

  seurattsne <- Seurat::RunTSNE(seurat, reduction = reduction, dims = dims, 
                                perplexity = perplexity, max_iter = max_iter)

  return(seurattsne)

}


# dimension plot for seurat

dimplot <- function(seurat, reduction = "umap", dims = c(1,2), label = TRUE, 
                    group.by = "orig.ident", cols = NULL, 
                    groups = NULL, split.by = NULL, title = NULL) {

  print(paste0("Computing for: ", paste(names(seurat), collapse = "-")))

  if (!is.null(groups)) {

    cols <- cols

    names(cols) <- groups

  }

  if (is.null(title)) {

    title <- unique(seurat$Sample)

  }

  seuratdim <- Seurat::DimPlot(seurat, reduction = reduction,
                               dims = dims, label = label,
                               group.by =  group.by, cols = cols,
                               split.by = split.by) + 
               labs(title = title)

  return(seuratdim)
}

# feature plot

featureplot <- function(seurat, reduction = "umap", dims = c(1,2), label = TRUE,
                        features = "nCount_Spatial", title = NULL) {

  print(paste0("Computing for: ", paste(names(seurat), collapse = "-")))

  if (is.null(title)) {

    title <- unique(seurat$Sample)

  }

  seuratfeature <- Seurat::FeaturePlot(seurat, reduction = reduction,
                                       dims = dims, label = label,
                                       features =  features) +
                   scale_color_viridis(option = "magma") +
                   labs(title = title, colour = features)

  return(seuratfeature)
}

# feature plot in Vln format - can access any metadata col or Gene of interest

featurevln <- function(seurat, features = "nCount_Spatial",
                       slot = NULL, log = T,
                       group.by = NULL, flip = FALSE,
                       xlabs = TRUE, title = NULL) {

  print(paste0("Computing for: ", paste(names(seurat), collapse = "-")))

  if (is.null(title)) {

    title <- paste0(unique(seurat$Sample), " ", features)

  } else {

    title <- paste0(title, " ", features)

  }

  if (is.null(slot)) {

    seuratvln <- VlnPlot(seurat, features = features, group.by = group.by) +
      labs(title = title)

  } else {

    # for raw counts slot = "counts" else normalised values will be used above
    seuratvln <- VlnPlot(seurat, features = features, slot = slot, log = log) +
      labs(title = title)
  }

  if (flip == TRUE) {

    seuratvln <- seuratvln + ggplot2::coord_flip()

  }

  if (xlabs == FALSE) {

    seuratvln <- seuratvln + ggplot2::theme(axis.text.x = ggplot2::element_blank())

  }
  return(seuratvln)

}

# Elbow plot

elbow <- function(seurat, ndims = 20, reduction = "pca", vline1 = 5,
                  vline2 = 10, title = NULL) {

  if (is.null(title)) {

    title <- unique(seurat$Sample)

  }

  elbowplot <- ElbowPlot(seurat, ndims = ndims, reduction = reduction) +
               geom_vline(xintercept = vline1, linetype = "dashed",
                          color = "red", size = 0.5) +
               geom_vline(xintercept = vline2, linetype = "dashed",
                          color = "red", size = 0.5) +
               geom_rect(aes(xmin = vline1, xmax = vline2, ymin = -Inf, ymax = Inf),
                         alpha = 0.002) +
               labs(title = title)

  return(elbowplot)
}

# Add elbow to seurat

addElbow <- function(seurat, elbow_vector) {

  dims <- elbow_vector[[unique(seurat$Sample)]]

  seurat[["Elbow"]] <- rep(dims, ncol(seurat))

  return(seurat)
}
