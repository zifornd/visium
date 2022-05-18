### RScript for loading in Spatial Transcriptomic data

## Load in Spatial data via Seurat helper functions

loadInSpatial <- function(samplename, prefix, metadata, saveout = T,
                          filename = "filtered_feature_bc_matrix.h5",
                          path = "/outs/spatial/") {
  
  suppressMessages(library(Seurat))

  samplefolder <- paste0(prefix, samplename)
  
  print(paste0("Reading in: ", samplefolder))
  
  image <- Read10X_Image(image.dir = paste0(samplefolder, path))
  
  slice <- metadata$image[metadata$sample == samplename]

  seurat <- Load10X_Spatial(
    data.dir = paste0(samplefolder, "/outs/"),
    filename = filename,
    assay = "Spatial",
    slice = slice,
    filter.matrix = TRUE,
    to.upper = FALSE,
    image = image
    )
  
  seurat <- addMeta(seurat, metadata, samplename)
  
  # save out
  if(saveout){
    
    saveRDS(seurat, file = paste0("../output/01-data-loading-",
                                  samplename ,".rds"))
    
    print(paste0("Saved out as: ", "../output/01-data-loading-",
                 samplename ,".rds"))
  
  } else {
    
    return(seurat)
    
  }

}

## Add on Meta data to Seurat object

addMeta <- function(seurat, metadata, samplename){
  
  # Add meta data to seurat objects
  
  seurat@meta.data$Sample = rep(samplename, nrow(seurat@meta.data))
  
  seurat@meta.data$Image = rep(metadata$image[metadata$sample == samplename], 
                               nrow(seurat@meta.data))
  
  # seurat@meta.data$Slide = as.character(sapply(seurat@meta.data$Image, 
  #                                              function(x) {unlist(strsplit(x, "_"))[[1]]}))
  seurat@meta.data$Slide = rep(metadata$slide[metadata$sample == samplename], 
                               nrow(seurat@meta.data))
  
  seurat@meta.data$Group = rep(metadata$group[metadata$sample == samplename], 
                               nrow(seurat@meta.data))
  
  seurat@meta.data$Area = rep(metadata$area[metadata$sample == samplename], 
                              nrow(seurat@meta.data))
  
  return(seurat)
}


## Collect metrics which are output from 10x spaceranger run

collectMetrics <- function(samplename, prefix, 
                           filename = "metrics_summary.csv"){
  
  samplefolder <- paste0(prefix, samplename)
  
  datadir <- paste0(samplefolder, "/outs/", filename)
  
  print(paste0("Reading in: ", datadir))
  
  # Columns made R safe
  metrics <- read.csv(datadir)
  
  # Rename last column to original 
  colnames(metrics)[colnames(metrics) == 'Number.of.Panel.Genes....10.UMIs'] <- 'Number.of.Panel.Genes.Greater.Or.Equal.10.UMIs'
  
  return(metrics)

}

# Collect QC metrics from seurat metadata into single dataframe for plotting

collectQC <- function(seuratList){
  
  meta <- lapply(seuratList, function(seurat){

    meta <- seurat@meta.data

    meta
  })
  
  meta <- do.call("rbind.data.frame", meta)

  return(meta)
}


# merge samples based on seurat supplied extension of merge

seuratMerge <- function(seuratList, project = "SeuratProject"){
  
  seurat <- merge(seuratList[[1]], y = seuratList[2:length(seuratList)], 
                  add.cell.ids = names(seuratList), project = project,
                  merge.data = TRUE)
  
  VariableFeatures(seurat) <- unique(unlist(lapply(seuratList, VariableFeatures)))
  
  
  return(seurat)
}

# run integration of seurat objects rather than straight merge

seuratIntegrate <- function(seuratList, nfeatures = 2000, normalization.method = "SCT"){
  
  
  features <- SelectIntegrationFeatures(object.list = seuratList, nfeatures = nfeatures)
  
  seuratList <- PrepSCTIntegration(object.list = seuratList, anchor.features = features)
  
  anchors <- FindIntegrationAnchors(object.list = seuratList, normalization.method = normalization.method,
                                           anchor.features = features)
  
  combined.sct <- IntegrateData(anchorset = anchors, normalization.method = normalization.method)
  
  # Reset variable features from individual ones...not sure if this is correct or not
  # But this is what we were doing for the simple merge above
  # Other option is to use the features variable as specified above ^
  # VariableFeatures(combined.sct) <- unique(unlist(lapply(seuratList, VariableFeatures)))
  
  # VariableFeatures(combined.sct) <- features
  
  return(combined.sct)
  
}

# Rename cells so they have the same names as in merge

parseIntegrate <- function(seurat){
  
  char_remove <- unlist(gregexpr("_", colnames(seurat)))
  
  clnames <- substr(colnames(seurat),1,(char_remove - 1))
  
  clnames <- paste(seurat@meta.data$Sample, clnames, sep = "_")
  
  seurat <- RenameCells(seurat, new.names = clnames)
  
  return(seurat)
}
