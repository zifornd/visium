
# Normalise Script

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

findFeatures <- function(seurat, nfeatures = 2000, selection.method = "vst") {

  seurat_varfeatures <- FindVariableFeatures(seurat,
                                             selection.method = selection.method,
                                             nfeatures = nfeatures)

  return(seurat_varfeatures)

}

# Plot variable features, selection.method can be sct or vst

plotVariableGenes <- function(seurat, top = 10, selection.method = "sct") {

  # Identify the 10 most variable genes
  top10 <- head(Seurat::VariableFeatures(seurat), top)

  # plot variable features with and without labels
  plot1 <- Seurat::VariableFeaturePlot(seurat,
                                       selection.method = selection.method)

  plot2 <- Seurat::LabelPoints(plot = plot1, points = top10, repel = TRUE)

  return(plot2)
}

# Merge per sample (pairwise merge)
# Will require downstream investigation for potential regressing out of batch

mergeSamples <- function(seurat) {

  for (i in 1:(length(seurat) - 1)) {

    print(paste0("Merging: ", names(seurat[[i]])[3],
                 " ", names(seurat[[i + 1]])[3]))

    if (i == 1) {

      seurat_merged <- merge(seurat[[i]], seurat[[i + 1]])

    }
    if (i != 1) {

      seurat_merged <- merge(seurat_merged, seurat[[i + 1]])

    }
  }
  return(seurat_merged)
}

# Regress out batch effect

regressOut <- function(seurat, vars.to.regress = "batch",
                       do.scale = FALSE, do.center = FALSE) {

  seuratscaled <- Seurat::ScaleData(seurat, vars.to.regress = vars.to.regress,
                                    do.scale = do.scale, do.center = do.center)

  return(seuratscaled)

}


normaliseSctransform <- function(seurat, slot = "counts") {

  set.seed(0305)

  s <- Seurat::GetAssayData(object = seurat, slot = slot)

  vst_out <- sctransform::vst(s, latent_var = c("log_umi"),
                              return_gene_attr = TRUE,
                              return_cell_attr = TRUE, verbosity = 1)

  return(vst_out)

}

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
