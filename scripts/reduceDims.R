
#' Wrapper for Seurat RunPCA function
#'
#' @param seurat Seurat object.
#' @param assay Seurat assay to run PCA on
#' @return Seurat object with PCA
pca <- function(seurat, assay = "SCT") {

  print(paste0("Computing for: ", paste(names(seurat), collapse = "-")))

  seuratpca <- Seurat::RunPCA(seurat, assay = "SCT", verbose = FALSE)

  return(seuratpca)

}

#' Wrapper for Seurat RunUMAP function
#'
#' @param seurat Seurat object.
#' @param reduction name of dimensionality reduction to run UMAP on.
#' @param dims Dimensions to run UMAP on - if NULL will take from Seurat object
#' @param n.neighbors important param for UMAP plotting
#' @param min.dist important param for UMAP plotting
#' @return Seurat object with UMAP
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

#' Wrapper for Seurat RunTSNE function
#'
#' @param seurat Seurat object.
#' @param reduction name of dimensionality reduction to run TSNE on.
#' @param dims Dimensions to run TSNE on - if NULL will take from Seurat object
#' @param n.neighbors important param for TSNE plotting
#' @param max_iter maximum number of iterations for TSNE to converge.
#' @return Seurat object with TSNE
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

    print(paste0("perplexity: ", round(ncol(seurat)^(1 / 2), digits = -1)))

    perplexity <- round(ncol(seurat)^(1 / 2), digits = -1)

  } else {

    print(paste0("perplexity: ", perplexity))

    if (perplexity == "N^1/2") {

      print(paste0("perplexity: ", round(ncol(seurat)^(1 / 2), digits = -1)))

      perplexity <- round(ncol(seurat)^(1 / 2), digits = -1)

    }

  }

  seurattsne <- Seurat::RunTSNE(seurat, reduction = reduction, dims = dims, 
                                perplexity = perplexity, max_iter = max_iter)

  return(seurattsne)

}


#' Dimension plot used for any dimensionality reduction of seurat object (e.g. PCA, UMAP, TSNE)
#'
#' @param seurat Seurat object.
#' @param reduction name of dimensionality reduction to plot (e.g. pca, tsne, umap)
#' @param dims Dimensions to plot e.g. c(1,2)
#' @param label Whether to label plots
#' @param group.by How to colour plots - should be a header in meta data
#' @param cols custom colours if groups are supplied
#' @param groups name of groups within group.by meta data column
#' @param split.by if the plot should be split by a particular meta data variable
#' @param title title of plot - default to sample name
#' @return Seurat::DimPlot
dimplot <- function(seurat, reduction = "umap", dims = c(1, 2), label = TRUE,
                    group.by = "orig.ident", cols = NULL,
                    groups = NULL, split.by = NULL, title = NULL, title_size = 12) {

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
               labs(title = title) +
               theme(plot.title = element_text(size = title_size))


  return(seuratdim)
}

#' Feature plot used for any dimensionality reduction of seurat object with Feature overlayed
#'
#' @param seurat Seurat object.
#' @param reduction name of dimensionality reduction to plot (e.g. pca, tsne, umap)
#' @param dims Dimensions to plot e.g. c(1,2)
#' @param label Whether to label plots
#' @param features names of feature to overlay onto dimensionality reduction plot - from meta data
#' @param title title of plot - default to sample name
#' @return Seurat::FeaturePlot
featureplot <- function(seurat, reduction = "umap", dims = c(1, 2),
                        label = TRUE, features = "nCount_Spatial",
                        title = NULL, title_size = 12) {

  print(paste0("Computing for: ", paste(names(seurat), collapse = "-")))

  if (is.null(title)) {

    title <- unique(seurat$Sample)

  }

  seuratfeature <- Seurat::FeaturePlot(seurat, reduction = reduction,
                                       dims = dims, label = label,
                                       features =  features) +
                   scale_color_viridis(option = "magma") +
                   labs(title = title, colour = features) +
                   theme(plot.title = element_text(size = title_size))

  return(seuratfeature)
}

#' Violin Feature plot to look across groups of interest e.g. batches or samples
#' For example to see if % Mitochondria or Cell cycle activation much higher in a subset of samples
#'
#' @param seurat Seurat object.
#' @param features names of feature to to plot - from meta data
#' @param slot Slot used for plotting
#' @param log Should the value be in log scale?
#' @param group.by How to colour plots - should be a header in meta data
#' @param flip Should the plot be flipped?
#' @param xlabs label for x axis
#' @param title title of plot - default to sample name
#' @return Seurat::VlnPlot
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

    seuratvln <- Seurat::VlnPlot(seurat, features = features, group.by = group.by) +
      labs(title = title)

  } else {

    # for raw counts slot = "counts" else normalised values will be used above
    seuratvln <- Seurat::VlnPlot(seurat, features = features, slot = slot, log = log) +
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

#' Elbow plot of explained variance per principal component
#'
#' @param seurat Seurat object.
#' @param ndims Number of dimensions to plot
#' @param reduction name of dimension reduction used for plotting (e.g. pca)
#' @param vline1 first vertical line to plot
#' @param vline2 second vertical line to plot
#' @param title Title of plot
#' @return Seurat::ElbowPlot
elbow <- function(seurat, ndims = 40, reduction = "pca", vline1 = 5,
                  vline2 = 10, title = NULL, title_size = 12) {

  if (is.null(title)) {

    title <- unique(seurat$Sample)

  }

  elbowplot <- Seurat::ElbowPlot(seurat, ndims = ndims, reduction = reduction) +
               geom_vline(xintercept = vline1, linetype = "dashed",
                          color = "red", size = 0.5) +
               geom_vline(xintercept = vline2, linetype = "dashed",
                          color = "red", size = 0.5) +
               geom_rect(aes(xmin = vline1, xmax = vline2, ymin = -Inf, ymax = Inf),
                         alpha = 0.002) +
               labs(title = title) +
               theme(plot.title = element_text(size = title_size))

  return(elbowplot)
}

#' Add chosen elbow range to a seurat object
#'
#' @param seurat Seurat object.
#' @param elbow_vector character vector of elbow ranger (in string - c("1:20")).
#' @return Seurat object with elbow range added to meta data
addElbow <- function(seurat, elbow_vector) {

  dims <- elbow_vector[[unique(seurat$Sample)]]

  seurat[["Elbow"]] <- rep(dims, ncol(seurat))

  return(seurat)
}
