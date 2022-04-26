
# Plot heatmap

plotpHeatmap <- function(genes_int, seuratprep, seurat_sub, order = "Loupe", 
                         gaps_row = NULL, max.exprs = 3.5, min.exprs = -3.5, 
                         assay = "SCT", slot = "data"){
  
  # genes_int <- unique(c(head(rownames(upgenes),25), head(rownames(downgenes),25)))
  
  DefaultAssay(seuratprep) <- assay
  
  r <- as.data.frame(t(FetchData(seuratprep, vars = genes_int, slot = slot)))
  
  mycolan <- data.frame(Cluster = seuratprep@meta.data$Cluster, Loupe = seuratprep@meta.data$Loupe, Sample = seuratprep@meta.data$Sample)
  rownames(mycolan) <- rownames(seuratprep@meta.data)

  order_heatmap <- seurat_sub@meta.data[order(seurat_sub@meta.data[[order]]),]
  cell_order <- as.character(unique(rownames(order_heatmap)))
  
  gaps <- cumsum(as.numeric(table(order_heatmap[[order]])))
  
  pheatmap(r[,cell_order], color = col, cluster_rows = F, cluster_cols = F, 
           scale = "row", breaks = seq(min.exprs, max.exprs, length.out = 101), 
           show_colnames = F, annotation_col = mycolan, gaps_col = gaps, 
           gaps_row = gaps_row)
  
}


