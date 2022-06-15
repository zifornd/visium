
#' Wrapper script for SCTransform - to normalise seurat object.
#'
#' @param seurat Seurat object.
#' @param variable.features.n Number of highly variable gene features to use for normalisation.
#' @param vars.to.regress Variable to regress out and include in SCT normalisation - must be in seurat meta.data.
#' @param verbose
#' @param min_cells Only use genes that have been detected in at least this many cells; default is 5.
#' @return Seurat object normalised
normalise <- function(seurat, variable.features.n = 3000,
                      vars.to.regress = NULL, verbose = FALSE,
                      min_cells = 5) {

  print(paste0("Computing for: ", paste(names(seurat), collapse = "-")))

  seuratsct <- SCTransform(seurat,
                           assay = "Spatial",
                           verbose = verbose,
                           variable.features.n = variable.features.n,
                           vars.to.regress = vars.to.regress,
                           min_cells = min_cells)

  return(seuratsct)

}

#' Wrapper script for FindVariableFeatures
#'
#' @param seurat Seurat object.
#' @param nfeatures Number of variable features to select for.
#' @param selection.method Method to find variable features.
#' @return Seurat object with variable features applied
findFeatures <- function(seurat, nfeatures = 2000, selection.method = "vst") {

  seurat_varfeatures <- FindVariableFeatures(seurat,
                                             selection.method = selection.method,
                                             nfeatures = nfeatures)

  return(seurat_varfeatures)

}

#' Plot the highly variable features and label the top points
#'
#' @param seurat Seurat object.
#' @param top Number of variable features to label.
#' @param selection.method Method to find variable features.
#' @return Seurat::VariableFeaturePlot
plotVariableGenes <- function(seurat, top = 10, selection.method = "sct") {

  # Identify the 10 most variable genes
  top10 <- head(Seurat::VariableFeatures(seurat), top)

  # plot variable features with and without labels
  plot1 <- Seurat::VariableFeaturePlot(seurat,
                                       selection.method = selection.method)

  plot2 <- Seurat::LabelPoints(plot = plot1, points = top10, repel = TRUE)

  return(plot2)
}


#' Wrapper for ScaleData function to regress out variable of interest
#'
#' @param seurat Seurat object.
#' @param vars.to.regress Variable to regress - must be in seurat meta.data.
#' @param do.scale Scale data.
#' @param do.center Center data.
#' @return Seurat object with populate scale.data slot
regressOut <- function(seurat, vars.to.regress = "batch",
                       do.scale = FALSE, do.center = FALSE) {

  seuratscaled <- Seurat::ScaleData(seurat, vars.to.regress = vars.to.regress,
                                    do.scale = do.scale, do.center = do.center)

  return(seuratscaled)

}

#' Normalise with SCtransform using the vst package on assay data
#'
#' @param seurat Seurat object.
#' @param slot Slot to retrieve data from
#' @return Output from sctransform::vst
normaliseSctransform <- function(seurat, slot = "counts") {

  set.seed(0305)

  s <- Seurat::GetAssayData(object = seurat, slot = slot)

  vst_out <- sctransform::vst(s, latent_var = c("log_umi"),
                              return_gene_attr = TRUE,
                              return_cell_attr = TRUE, verbosity = 1)

  return(vst_out)

}

#' Plot an example of the SCT model fit for a given set of features.
#'
#' @param seurat Seurat object.
#' @param features Example features
#' @param vst_outlist Output from sctransform::vst
#' @param plot_residual Whether to plot residuals
#' @param slot Slot to retrieve data from
#' @return Output from sctransform::plot_model
plotModelSctransform <- function(seurat, features, vst_outlist,
                                 plot_residual = TRUE,
                                 slot = "counts") {

  set.seed(0305)

  vst_out <- vst_outlist[[unique(seurat$Sample)]]

  s <- Seurat::GetAssayData(object = seurat, slot = slot)

  p <- sctransform::plot_model(vst_out, s, features,
                               plot_residual = plot_residual)

  return(p)
}
