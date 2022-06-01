
## parse cloupe csv - adding in sample name to cell ids and standardising colnames

parseCloupe <- function(cloupe, sample = NULL) {

  if (is.null(sample)) {

    sample <- colnames(cloupe)[2]

  }

  rownames(cloupe) <- paste(rep(sample, nrow(cloupe)), cloupe$Barcode,
                             sep = "_")

  cloupe$Barcode <- NULL

  colnames(cloupe) <- "Loupe"

  return(cloupe)
}

## Read in loupe annotation and return as merged single col annotation

readCloupe <- function(cloupeLoc) {

  suppressMessages(library(dplyr))

  cloupes <- list.files(path = cloupeLoc, pattern = "*_cloupe.csv")

  cloupes <- lapply(paste0(cloupeLoc, cloupes), read.csv, check.names = F)

  cloupes <- lapply(cloupes, parseCloupe)

  print(head(cloupes[[1]]))

  cloupes <- dplyr::bind_rows(cloupes)

  return(cloupes)
}

## Add annotation to seurat object

addCloupe <- function(seurat, cloupe) {

  # find those spots not in annotation

  unknw <- rownames(seurat@meta.data)[!rownames(seurat@meta.data) %in% rownames(cloupe)]

  noinanno <- data.frame(Loupe = rep("Unknown", length(unknw)))

  rownames(noinanno) <- unknw

  # noinanno all new cell ids
  stopifnot(!rownames(noinanno) %in% rownames(cloupe))

  # check length ok
  stopifnot((length(rownames(seurat@meta.data)) - nrow(cloupe)) == length(unknw))

  # add not in anno spots to loupe anno
  cloupe <- rbind.data.frame(cloupe, noinanno)

  stopifnot(nrow(cloupe) == nrow(seurat@meta.data))

  cloupe <- cloupe[match(rownames(seurat@meta.data), rownames(cloupe)), ,drop = F]

  stopifnot(rownames(cloupe) == rownames(seurat@meta.data))

  seurat$Loupe <- cloupe$Loupe

  print(table(seurat$Loupe))

  return(seurat)
}


testDoHeatmap <- function(doHeatmapCell, doHeatmap, metadata) {

  doHeatmapLoupe <- unique(doHeatmap$Identity[doHeatmap$Cell == doHeatmapCell])

  if (grep("^[0-9]", doHeatmapCell)) {

    print(metadata$Loupe[rownames(metadata) == doHeatmapCell])

    print(doHeatmapLoupe)

    stopifnot(metadata$Loupe[rownames(metadata) == doHeatmapCell] == doHeatmapLoupe)

    }
}
