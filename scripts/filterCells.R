

# filter by value in seurat meta data table

filterBy <- function(seurat, threshold, by = "nCounts_Spatial", 
                     direction = "down") {
  
  # fetch variable of interest from seurat object as data frame
  df <- FetchData(seurat, vars = by)
  
  # keep only cells of interest
  if(dir == "up"){
    keep <- df[,by] > threshold
  }
  if(dir == "down"){
    keep <- df[,by] <= threshold
  }

  SpatialDimPlot(cortex,cells.highlight = keep)
  
  seurat_filt <- seurat[,keep]
  
  return(seurat_filt)
  
}


# filter spots spatially 

filterByXY <- function(seurat, x1, x2, y1, y2){
  
  # presuming only one image
  stopifnot(length(seurat@images) == 1)
  image_key <- seurat@images[[1]]@key
  
  # fetch xy axis data frame
  spatial_df <- FetchData(seurat, vars = c(paste0(image_key, "imagerow"), 
                                         paste0(image_key, "imagecol")))
  
  # keep only cells of interest
  keep <- spatial_df[,paste0(image_key, "imagerow")] >= y1 & 
          spatial_df[,paste0(image_key, "imagerow")] <= y2 &
          spatial_df[,paste0(image_key, "imagecol")] >= x1 &
          spatial_df[,paste0(image_key, "imagecol")] <= x2
  
  seurat_filt <- seurat[,keep]
  
  return(seurat_filt)
  
}

