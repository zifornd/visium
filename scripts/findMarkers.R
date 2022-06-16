
#' Scale seurat object prior to visualisation
#' Have the option to regress out variable
#' @param seuratprep Seurat object.
#' @param regress name of variable to regress out - if FALSE will not regress out just Scale.
#' @return Seurat object
scaleData <- function(seuratprep, regress) {

    if (regress != FALSE) {

        seuratscale <- Seurat::ScaleData(seuratprep@assays$SCT,
                                         features = rownames(seuratprep),
                                         vars.to.regress = regress)

    } else {

        seuratscale <- Seurat::ScaleData(seuratprep@assays$SCT,
                                         features = rownames(seuratprep))

    }

        seuratprep[["ScaleData"]] <- seuratscale

    return(seuratprep)

}

#' Run Seurat FindAllMarkers
#' @param seuratprep Seurat object.
#' @param params full params from 08-marker-detection.qmd (only.pos, min.pct, logfc.threshold, pval.adj.thres).
#' @param sample_name User defined sample name will be used where multiple samples are saved within seurat object (merged/integrated)
#' @return Result table containing markers per cluster.
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

#' Get top genes per cluster from results table
#' @param clust Which cluster of interest to select from
#' @param res Results table supplied
#' @param head_num Number of genes to retrieve
#' @return Top genes of interest
topFeatClust <- function(clust, res, head_num = 10) {

    features <- head(res$gene[res$cluster == clust], head_num)

    return(features)
}

#' Get top genes per seurat object from results table
#' @param seurat Seurat object
#' @param resList List of results tables
#' @param head_num Number of genes to retrieve
#' @param sample_name User defined sample name will be used where multiple samples are saved within single seurat object (merged/integrated)
#' @return Top genes of interest
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
#' Wrapper for DoHeatmap and pheatmap functions
#' @param seuratprep Seurat object prepared for plotting
#' @param featureList List of features of interest
#' @param head_num Number of genes to plot per group
#' @param type Type of plot to complete
#' @param legend Whether to plot legend
#' @param sample_name User defined sample name will be used where multiple samples are saved within single seurat object (merged/integrated)
#' @return Seurat::DoHeatmap or pheatmap plot
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