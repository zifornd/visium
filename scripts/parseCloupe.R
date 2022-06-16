
#' Read in loupe annotation and return as merged single column annotation
#'
#' @param cloupeLoc Location of cloupe files.
#' @return Merged single column cloupe annotation.
readCloupe <- function(cloupeLoc) {

  suppressMessages(library(dplyr))

  cloupes <- list.files(path = cloupeLoc, pattern = "*_cloupe.csv")

  cloupes <- lapply(paste0(cloupeLoc, cloupes), read.csv, check.names = F)

  cloupes <- lapply(cloupes, parseCloupe)

  print(head(cloupes[[1]]))

  cloupes <- dplyr::bind_rows(cloupes)

  return(cloupes)
}

#' Parse cloupe csv - adding in sample name to cell ids and standardising colnames
#'
#' @param cloupe Output from cloupe manual annotation of interesting cell types.
#' @param sample Name of sample - if NULL taken from cloupe output.
#' @return Parsed cloupe output dataframe.
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

#' Add cloupe annotation to seurat object
#'
#' @param seurat Seurat object.
#' @param cloupe Parsed cloupe output dataframe.
#' @return Seurat object with cloupe annotation added.
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
