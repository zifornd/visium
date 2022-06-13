
# Scale seurat object prior to visualisation

scaleData <- function(seuratprep, params) {

    if (params$regress != FALSE) {

        seuratscale <- Seurat::ScaleData(seuratprep@assays$SCT,
                                         features = rownames(seuratprep),
                                         vars.to.regress = params$regress)

    } else {

        seuratscale <- Seurat::ScaleData(seuratprep@assays$SCT,
                                         features = rownames(seuratprep))

    }

        seuratprep[["ScaleData"]] <- seuratscale

    return(seuratprep)

}

# Run Seurat FindAllMarkers

findAllMarkers <- function(seuratprep, params,
                           sample_name = "Combined.Seurat") {

    Seurat::Idents(seuratprep) <- seuratprep$Cluster

    res <- Seurat::FindAllMarkers(seuratprep, only.pos = params$only.pos,
                                  min.pct = params$min.pct,
                                  logfc.threshold = params$logfc.threshold)

    rownames(res) <- NULL

    res <- res[res$p_val_adj < params$pval.adj.thres, ]

    # If combined will have multiple sample names - so need to check

    if (length(unique(seuratprep@meta.data$Sample)) == 1) {

      sample_name <- unique(seuratprep@meta.data$Sample)

    }

    res$Sample <- rep(sample_name, nrow(res))

    return(res)
}

# Get top genes per cluster from results table

topFeatClust <- function(clust, res, head_num = 10) {

    features <- head(res$gene[res$cluster == clust], head_num)

    return(features)
}

# Get top genes per seurat object from results table

topFeatClustList <- function(seurat, resList, head_num = 10,
                             sample_name = "Combined.Seurat") {

    # If combined will have multiple sample names - so need to check

    if (length(unique(seurat@meta.data$Sample)) == 1) {

      sample_name <- unique(seurat@meta.data$Sample)

    }

    # get results from results list
    res <- resList[[sample_name]]

    clusters <- sort(unique(seurat@meta.data$Cluster))

    features <- lapply(clusters, topFeatClust, res = res, head_num = head_num)

    features <- unlist(features)

    return(features)
}

# Wrapper for DoHeatmap and pheatmap functions

plotHeatmap <- function(seuratprep, featureList, head_num = 10,
                        type = "DoHeatmap", legend = TRUE,
                        sample_name = "Combined.Seurat") {

    if (length(unique(seuratprep@meta.data$Sample)) == 1) {

      sample_name <- unique(seuratprep@meta.data$Sample)

    }

    # Get top features
    features <- featureList[[sample_name]]

    # Plot from Seurat using DoHeatmap function
    if (type == "DoHeatmap") {

        Seurat::Idents(seuratprep) <- seuratprep$Cluster

        p <- Seurat::DoHeatmap(seuratprep, features = features,
                               assay = "ScaleData", slot = "scale.data")

        if (!legend) {

            p <- p + Seurat::NoLegend()

        }
    }

    # Plot from pheatmap
    if (type == "pheatmap") {

        Seurat::Idents(seuratprep) <- seuratprep$Cluster

        # get gaps for heatmap
        gaps <- length(unique(seuratprep@meta.data$Cluster)) - 1

        p <- plotpHeatmap(features, seuratprep, seuratprep, order = "Cluster",
                          gaps_row = cumsum(c(rep(head_num, (gaps)))),
                          min.exprs = -3, max.exprs = 3)

    }

    return(p)
}