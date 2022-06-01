
# Annotate seurat obj metadata (per spot) with T/F if meets thres of expression
# Labels column based on some user defined name
# Can take multiple genes at once - T when both genes meet threshold

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


#  Annotate seurat obj metadata (per spot) with T/F if meets thres of expression
# And if other gene feature is lower than some threshold
# Labels column based on some user defined name
# Can take multiple genes at once - T when both genes meet threshold

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

# combing annotation

combineAnno <- function(seurat, anno) {

  anno_interest <- seurat[[anno]]

  seurat[["Annotation"]] <- apply(anno_interest, 1, defineLabel)

  return(seurat)
}

# Defines the returned label for a set of labels for combining

defineLabel <- function(rowlabel, other = "Other") {

  if (length(unique(rowlabel)) == 1 & unique(rowlabel)[1] == other) {

    return(other)

  }

  if (other %in% rowlabel){

    label <- unique(rowlabel)

    label <- paste0(label[!label %in% other], collapse = "_")

    return(label)

  } else {

    label <- unique(rowlabel)

    label <- paste0(label, collapse = "_")

    return(label)
  }
}

# Define expressed genes and plot distribution

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