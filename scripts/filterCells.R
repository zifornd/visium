#' Filter Seurat object by threshold, direction, and feature
#'
#' @param seurat Seurat object.
#' @param threshold Threshold to filter seurat object by.
#' @param by Feature to filter by e.g. nCount_Spatial
#' @param direction Direction of filter - up or down
#' @return Seurat object filtered.
filterBy <- function(seurat, threshold, by = "nCount_Spatial",
                     direction = "down") {

  # fetch variable of interest from seurat object as data frame
  df <- Seurat::FetchData(seurat, vars = by)

  # keep only cells of interest
  if (direction == "up") {

    keep <- df[, by] > threshold

  }
  if (direction == "down") {

    keep <- df[, by] < threshold

  }

  seurat_filt <- seurat[, keep]

  return(seurat_filt)

}

#' Plot filtered Seurat object spots by threshold, direction, and feature
#'
#' @param seurat Seurat object.
#' @param threshold Threshold to filter seurat object by.
#' @param by Feature to filter by e.g. nCount_Spatial
#' @param direction Direction of filter - up or down
#' @return Seurat::SpatialDimPlot
filterByPlot <- function(seurat, threshold, by = "nCount_Spatial",
                         direction = "down") {

  # fetch variable of interest from seurat object as data frame
  df <- Seurat::FetchData(seurat, vars = by)

  # keep only cells of interest
  if (direction == "up") {

    keep <- df[, by] > threshold

  }
  if (direction == "down") {

    keep <- df[, by] < threshold

  }

  seurat_filt <- seurat[, keep]

  plot <- Seurat::SpatialDimPlot(seurat,
                                 cells.highlight = colnames(seurat_filt))

  return(plot)
}

#' Filter Seurat object by x and y coordinates
#' Intersection of coordinates defines filtered area
#' @param seurat Seurat object.
#' @param x1 x axis coord 1
#' @param x2 x axis coord 2
#' @param y1 y axis coord 1
#' @param y2 y axis coord 2
#' @return Seurat object filtered.
filterByXY <- function(seurat, x1, x2, y1, y2) {

  # presuming only one image
  stopifnot(length(seurat@images) == 1)

  image_key <- seurat@images[[1]]@key

  # fetch xy axis data frame
  spatial_df <- Seurat::FetchData(seurat,
                                  vars = c(paste0(image_key, "imagerow"),
                                           paste0(image_key, "imagecol")))

  # keep only cells of interest
  keep <- spatial_df[, paste0(image_key, "imagerow")] >= y1 &
          spatial_df[, paste0(image_key, "imagerow")] <= y2 &
          spatial_df[, paste0(image_key, "imagecol")] >= x1 &
          spatial_df[, paste0(image_key, "imagecol")] <= x2

  seurat_filt <- seurat[, keep]

  return(seurat_filt)

}
