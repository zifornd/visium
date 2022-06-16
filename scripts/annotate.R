
#' Annotate seurat obj metadata (per spot) with T/F if meets thres of expression for a given set of genes.
#' Labels column based on some user defined name
#' Can take multiple genes at once - T when all genes meet > threshold
#'
#' @param seurat Seurat object.
#' @param features Name of gene features to check.
#' @param anno Name of annotation column to apply to seurat meta data and name of annotation given to spot if it meets threshold.
#' @param assay Name of Seurat assay to retrieve feature expression values.
#' @param slot Name of data slot to retrieve feature expression values.
#' @param threshold threshold for annotation.
#' @param other Name of applied to spot if it does not meet threshold.
#' @return Seurat object with additional annotation.
annotate <- function(seurat, features, anno, assay = "Spatial",
                     slot = "data", threshold = 1, other = "Other") {

  # ready data
  seuratdata <- Seurat::FetchData(seurat, vars = features, slot = slot)

  seurat[[features]] <- as.data.frame(as.matrix(seuratdata) > threshold)

  seurat[[anno]] <- apply(seurat[[features]], 1, function(x) ifelse(T %in% x & !F %in% x, anno, other))

  # remove feature cols from metadata -
  # as this can affect downstream retrieval of features

  seurat@meta.data[, features] <- NULL

  return(seurat)
}

#' Annotate seurat obj metadata (per spot) with T/F if meets thres of expression
#' And if other gene feature is lower than some threshold
#' Labels column based on some user defined name
#' Can take multiple genes at once - T when both genes meet threshold
#'
#' @param seurat Seurat object.
#' @param features_high Name of gene features to check > threshold.
#' @param features_low Name of gene features to check < threshold.
#' @param anno Name of annotation column to apply to seurat meta data and name of annotation given to spot if it meets threshold.
#' @param slot Name of data slot to retrieve feature expression values.
#' @param threshold threshold for annotation.
#' @param other Name of applied to spot if it does not meet threshold.
#' @return Seurat object with additional annotation.
annotate_highlow <- function(seurat, features_high, features_low, anno,
                             slot = "data", threshold = 0.5, other = "Other") {

  seuratdata <- Seurat::FetchData(seurat, vars = features_high, slot = slot)

  seurat[[features_high]] <- as.data.frame(as.matrix(seuratdata) > threshold)

  seurat[[features_low]] <- as.data.frame(as.matrix(seuratdata) < threshold)

  seurat_high_low <- seurat[[c(features_high, features_low)]]

  seurat[[anno]] <- apply(seurat_high_low, 1, function(x) ifelse(T %in% x & !F %in% x, anno, other))

  # remove feature cols from metadata
  # - as this can affect downstream retrieval of features

  seurat@meta.data[, c(features_high, features_low)] <- NULL

  return(seurat)

}

#' Parses label for combining
#'
#' @param rowlabel Name of a annotation label.
#' @param other Name of annotation label which defines spots not defined by thresholds.
#' @return Parsed label.
defineLabel <- function(rowlabel, other = "Other") {

  if (length(unique(rowlabel)) == 1 & unique(rowlabel)[1] == other) {

    return(other)

  }

  if (other %in% rowlabel) {

    label <- unique(rowlabel)

    label <- paste0(label[!label %in% other], collapse = "_")

    return(label)

  } else {

    label <- unique(rowlabel)

    label <- paste0(label, collapse = "_")

    return(label)
  }
}

#' Combine annotation together
#'
#' @param seurat Seurat object.
#' @param anno Name of annotation column in seurat meta data
#' @return Seurat object with additional annotation.
combineAnno <- function(seurat, anno) {

  anno_interest <- seurat[[anno]]

  seurat[["Annotation"]] <- apply(anno_interest, 1, defineLabel)

  return(seurat)
}

#' Plots the distribution of expression for a marker of interest to help define what is "expressed" vs "not expressed"
#' Use this visual guide to feed into functions above.
#'
#' @param seurat Seurat object.
#' @param features Name of gene features to check.
#' @param assay Name of Seurat assay to retrieve feature expression values.
#' @param slot Name of data slot to retrieve feature expression values.
#' @param threshold threshold to plot on x axis of distribution.
#' @return Seurat object with additional annotation.
defineExpressed <- function(seurat, features, assay = "Spatial", slot = "data",
                            threshold = NULL) {

  # ready data

  if (features %in% rownames(GetAssayData(object = seurat, slot = slot))) {

    seuratdata <- Seurat::GetAssayData(object = seurat, slot = slot)

    df <- as.data.frame(t(as.data.frame(seuratdata[features, , drop = F])))

    xlabel <- "Normalised Counts"

  } else {

    # try set default assay back to spatial if gene has been filtered by SCT 

    Seurat::DefaultAssay(seurat) <- assay

    seuratdata <- Seurat::GetAssayData(object = seurat, slot = slot)

    df <- as.data.frame(t(as.data.frame(seuratdata[features, , drop = F])))

    # but default assay back to previous SCT

    Seurat::DefaultAssay(seurat) <- "SCT"

    xlabel <- "Counts"
  }

  df2 <- reshape2::melt(df)

  colnames(df2) <- c("Feature", "Counts")

  # plot distribution

  if (!is.null(threshold)) {

    p <- ggplot(df2, aes(x = Counts, color = Feature)) +
      geom_density() +
      geom_vline(aes(xintercept = threshold, color = Feature),
                     linetype = "dashed") +
      labs(title = paste0(unique(seurat$Sample)), x = xlabel, y = "Density") +
      theme_classic()

  } else {

    p <- ggplot(df2, aes(x = Counts, color = Feature)) +
      geom_density() +
      labs(title = paste0(unique(seurat$Sample)), x = xlabel, y = "Density") +
      theme_classic()
  }

  return(p)

  # decide on threshold
  # ? gaussian mixture?
}