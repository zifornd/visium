
#' Load in Spatial data via Seurat helper functions
#'
#' @param samplename individual sample name to load in.
#' @param prefix folder location for samples.
#' @param metadata Experiment metadata from samplesheet which includes image and sample information
#' @param saveout Whether the files should be saved out as rds objects instead of returned
#' @param filename Feature bc matrix in h5 format intended to load in.
#' @param path path of spatial data folder which holds files such tissue_hires_image.jpg
#' @return Seurat object
loadInSpatial <- function(samplename, prefix, metadata, saveout = T,
                          filename = "filtered_feature_bc_matrix.h5",
                          path = "/outs/spatial/") {

  suppressMessages(library(Seurat))

  samplefolder <- paste0(prefix, samplename)

  print(paste0("Reading in: ", samplefolder))

  image <- Seurat::Read10X_Image(image.dir = paste0(samplefolder, path))

  slice <- metadata$image[metadata$sample == samplename]

  seurat <- Seurat::Load10X_Spatial(data.dir = paste0(samplefolder, "/outs/"),
                                    filename = filename,
                                    assay = "Spatial",
                                    slice = slice,
                                    filter.matrix = TRUE,
                                    to.upper = FALSE,
                                    image = image)

  seurat <- addMeta(seurat, metadata, samplename)

  # save out
  if (saveout) {

    saveRDS(seurat, file = paste0("../output/01-data-loading-",
                                  samplename, ".rds"))

    print(paste0("Saved out as: ", "../output/01-data-loading-",
                 samplename, ".rds"))

  } else {

    return(seurat)

  }

}

#' Wrapper for loadInSpatial to Load in Spatial data list
#'
#' @param idx sample index.
#' @param samples character vector of sample names.
#' @param prefix folder location for samples.
#' @param metadata Experiment metadata from samplesheet which includes image and sample information
#' @param saveout Whether the files should be saved out as rds objects instead of returned
#' @param filenames character vector of feature bc matrix files in h5 format intended to load in.
#' @return Seurat object
loadInSpatialList <- function(idx, samples, prefix,
                              metadata, saveout, filenames) {

  seurat <- loadInSpatial(samples[[idx]], prefix,
                          metadata, saveout, filenames[[idx]])

  return(seurat)

}

#' Add Meta data to Seurat object
#'
#' @param seurat Seurat object.
#' @param metadata Experiment metadata from samplesheet which includes image and sample information.
#' @param samplename Name of sample.
#' @return Seurat object
addMeta <- function(seurat, metadata, samplename) {

  # Add meta data to seurat objects

  seurat@meta.data$Sample <- rep(samplename, nrow(seurat@meta.data))

  seurat@meta.data$Image <- rep(metadata$image[metadata$sample == samplename],
                               nrow(seurat@meta.data))

  seurat@meta.data$Slide <- rep(metadata$slide[metadata$sample == samplename],
                               nrow(seurat@meta.data))

  seurat@meta.data$Group <- rep(metadata$group[metadata$sample == samplename],
                               nrow(seurat@meta.data))

  seurat@meta.data$Area <- rep(metadata$area[metadata$sample == samplename],
                              nrow(seurat@meta.data))

  return(seurat)
}


#' Collect metrics which are output from 10x spaceranger run
#'
#' @param samplename Seurat object.
#' @param prefix folder location for samples.
#' @param filename Name of filename suffix.
#' @return metrics dataframe of QC metrics.
collectMetrics <- function(samplename, prefix,
                           filename = "_metrics_summary.csv") {

  samplefolder <- paste0(prefix, samplename)

  datadir <- paste0(samplefolder, "/outs/", samplename, filename)

  print(paste0("Reading in: ", datadir))

  # Columns made R safe
  metrics <- read.csv(datadir)

  # Rename last column to original if needed
  if ('Number.of.Panel.Genes....10.UMIs' %in% colnames(metrics)) {

    colnames(metrics)[colnames(metrics) == 'Number.of.Panel.Genes....10.UMIs'] <- 'Number.of.Panel.Genes.Greater.Or.Equal.10.UMIs'

  }

  return(metrics)

}

#' Collect QC metrics from seurat metadata into single dataframe for plotting
#'
#' @param seuratList List of Seurat objects.
#' @return meta dataframe of combined Seurat meta data.
collectQC <- function(seuratList) {

  meta <- lapply(seuratList, function(seurat) {

    meta <- seurat@meta.data

    meta
  })

  meta <- do.call("rbind.data.frame", meta)

  return(meta)
}


#' Merge samples based on seurat supplied extension of merge
#'
#' @param seuratList List of Seurat objects.
#' @param project Name of project used for new merged seurat object.
#' @return seurat Seurat object merged.
seuratMerge <- function(seuratList, project = "SeuratProject") {

  if (length(seuratList) > 1) {

    seurat <- merge(seuratList[[1]], y = seuratList[2:length(seuratList)],
                    add.cell.ids = names(seuratList), project = project,
                    merge.data = TRUE)

    Seurat::VariableFeatures(seurat) <- unique(unlist(lapply(seuratList,
                                                             VariableFeatures)))

  } else {

    seurat <- seuratList[[1]]

  }

  return(seurat)
}

#' Integrate seurat objects which have been normalised with SCT
#'
#' @param seuratList List of Seurat objects.
#' @param nfeatures Name of project used for new merged seurat object.
#' @param normalization.method Method of normalisation for integration.
#' @return seurat Seurat object integrated.
seuratIntegrate <- function(seuratList, nfeatures = 2000,
                            normalization.method = "SCT") {

  if (length(seuratList) > 1) {

    features <- Seurat::SelectIntegrationFeatures(object.list = seuratList,
                                                  nfeatures = nfeatures)

    seuratList <- Seurat::PrepSCTIntegration(object.list = seuratList,
                                             anchor.features = features)

    anchors <- Seurat::FindIntegrationAnchors(object.list = seuratList,
                                              normalization.method = normalization.method,
                                              anchor.features = features)

    combined.sct <- Seurat::IntegrateData(anchorset = anchors,
                                          normalization.method = normalization.method)

  } else {

    combined.sct <- seuratList[[1]]

  }

  return(combined.sct)

}

#' Rename cells post integration to ensure naming correct
#'
#' @param seurat Seurat object.
#' @return seurat Seurat object.
parseIntegrate <- function(seurat) {

  char_remove <- unlist(gregexpr("_", colnames(seurat)))

  clnames <- substr(colnames(seurat), 1, (char_remove - 1))

  clnames <- paste(seurat@meta.data$Sample, clnames, sep = "_")

  seurat <- Seurat::RenameCells(seurat, new.names = clnames)

  return(seurat)
}

#' Pairwise merge command - DEPRECATED
#'
#' @param seurat Seurat object.
#' @return Seurat object merged
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

#' Prepare seurat object prior to FindMarkers, or Spatial analysis
#'
#' @param seurat Single seurat object or list of seurat objects.
#' @return List of seurat objects
prepSeuratLoad <- function(seurat) {

  if (class(seurat) == "Seurat") {

    seuratList <- list(seurat)

    # If merged or combined will have multiple samples in meta data
    # so save into list as combined seurat

    if (length(unique(seurat@meta.data$Sample)) > 1) {

      names(seuratList) <- "Combined.Seurat"

      # Run PrepSCTFindMarkers
      seuratprep <- lapply(seuratList, Seurat::PrepSCTFindMarkers)

    }

    # If single sample save into list as sample name
    if (length(unique(seurat@meta.data$Sample)) == 1) {

      names(seuratList) <- unique(seurat@meta.data$Sample)

      # If only a single un merged/combined experiment
      # - do not need to run PrepSCTFindMarkers

      seuratprep <- seuratList

    }

  } else {

      # If only list of un merged/combined experiments
      # - do not need to run PrepSCTFindMarkers

      stopifnot(class(seurat) == "list")

      seuratprep <- seurat
  }

  return(seuratprep)
}