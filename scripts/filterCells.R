
# filter by value in seurat meta data table

filterBy <- function(seurat, threshold, by = "nCount_Spatial", 
                     direction = "down") {
  
  # fetch variable of interest from seurat object as data frame
  df <- FetchData(seurat, vars = by)
  
  # keep only cells of interest
  if(direction == "up"){
    keep <- df[,by] > threshold
  }
  if(direction == "down"){
    keep <- df[,by] < threshold
  }

  #SpatialDimPlot(seurat,cells.highlight = keep)
  
  seurat_filt <- seurat[,keep]
  
  return(seurat_filt)
  
}

filterByPlot <- function(seurat, threshold, by = "nCount_Spatial", 
                     direction = "down") {
  
  # fetch variable of interest from seurat object as data frame
  df <- FetchData(seurat, vars = by)
  
  # keep only cells of interest
  if(direction == "up"){
    keep <- df[,by] > threshold
  }
  if(direction == "down"){
    keep <- df[,by] < threshold
  }

  seurat_filt <- seurat[,keep]

  plot <- SpatialDimPlot(seurat,cells.highlight = colnames(seurat_filt))
    
  return(plot)
  
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

