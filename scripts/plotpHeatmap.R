
# Plot heatmap

plotpHeatmap <- function(genes_int, seuratprep, seurat_sub, order = "Cluster",
                         gaps_row = NULL, max.exprs = 3.5, min.exprs = -3.5,
                         assay = "SCT", slot = "data") {


  Seurat::DefaultAssay(seuratprep) <- assay

  r <- as.data.frame(t(Seurat::FetchData(seuratprep,
                                         vars = genes_int,
                                         slot = slot)))

  mycolan <- data.frame(Cluster = seuratprep@meta.data$Cluster,
                        Sample = seuratprep@meta.data$Sample)

  rownames(mycolan) <- rownames(seuratprep@meta.data)

  order_heatmap <- seurat_sub@meta.data[order(seurat_sub@meta.data[[order]]), ]

  cell_order <- as.character(unique(rownames(order_heatmap)))

  gaps <- cumsum(as.numeric(table(order_heatmap[[order]])))

  # Take [[4]], which is the gtable of the plot.
  # You have to index it specifically since it's not going to trigger R's print method lookup/execution for that object inside the function call.
  # https://www.biostars.org/p/128229/

  p <- pheatmap(r[, cell_order], color = col, cluster_rows = F,
                cluster_cols = F, scale = "row",
                breaks = seq(min.exprs, max.exprs, length.out = 101),
                show_colnames = F, annotation_col = mycolan, gaps_col = gaps,
                gaps_row = gaps_row)[[4]]

}