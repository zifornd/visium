
#' Load OrgDb package
#'
#' @param x Organism - Mm for mouse, Hs for Human
#' @return None
OrgDb <- function(x) {

    pkg <- paste0("org.", x, ".eg.db")

    obj <- getFromNamespace(pkg, pkg)

}

#' Parse results table to retrieve a filtered vector of genes names or ranked named vector.
#' Splits results table by group e.g. cluster annotation.
#' Filters results with filter_df function ready for downstream overrep or gsea.
#' Then outputs vector of gene names or ranked named vector.
#'
#' @param res Results table (with padj, lfc, gene symbol) typically from Seurat::FindAllMarkers
#' @param group Name of group column that table would be split by if split.res = TRUE
#' @param p_adj Adjusted p-value threshold for filtering results table (set as Inf if no filtering wanted)
#' @param lfc Log2fold change threshold for filtering results table (set as 0 if no filtering wanted)
#' @param type type of parsing required - for GSEA ranked named vector is outputed.
#' @param lfc_name Name of lfc col in results table
#' @param padj_name Name of padj col in results table
#' @param gene_id_name Name of gene symbol col in results table
#' @param split.res Whether to split table based on group.
#' @return parsed results
parse_res <- function(res, group = "cluster", p_adj = Inf, lfc = 0,
                      type = "overrep", lfc_name = "avg_log2FC",
                      padj_name = "p_val_adj", gene_id_name = "gene",
                      split.res = TRUE) {

    if (split.res) {

        resList <- split_res(res, group = group)

    } else {

        resList <- list(res)

    }

    filtered_res <- lapply(resList, filter_df, p_adj = p_adj, lfc = lfc,
                           type = type, lfc_name = lfc_name,
                           padj_name = padj_name, gene_id_name = gene_id_name)

    if (type == "overrep") {

        vector_res <- lapply(filtered_res, function(x, gene_id_name) {
          x[, gene_id_name] }, gene_id_name = gene_id_name)

    }

    if (type == "gsea") {

        vector_res <-  lapply(filtered_res, make_lfc_vector,
                              col_name = gene_id_name, lfc_name = lfc_name)

    }

    return(vector_res)
}

#' Splits results into list of results by given group column
#'
#' @param res Results table (with padj, lfc, gene symbol) typically from Seurat::FindAllMarkers
#' @param group Name of group column that table would be split
#' @return List of results split by group
split_res <- function(res, group = "cluster") {

    groups <- unique(res[[group]])

    resList <- lapply(groups, function(x, res, group) {

      res[res[[group]] == x, ]

    }, res = res, group = group)

    names(resList) <- groups

    return(resList)
}


#' Filter Results table prior to overrepresentation functional enrichment analysis
#' If you dont want to filter on lfc - leave as 0
#' If you dont want to filter on p value leave as Inf
#'
#' @param df Results table (with padj, lfc, gene symbol) typically from Seurat::FindAllMarkers
#' @param p_adj Adjusted p-value threshold for filtering results table (set as Inf if no filtering wanted)
#' @param lfc Log2fold change threshold for filtering results table (set as 0 if no filtering wanted)
#' @param type type of parsing required - for GSEA ranked named vector is outputed.
#' @param lfc_name Name of lfc col in results table
#' @param padj_name Name of padj col in results table
#' @param gene_id_name Name of gene symbol col in results table
#' @return parsed results
filter_df <- function(df, p_adj = Inf, lfc = 0, type = "overrep",
                      lfc_name = "avg_log2FC", padj_name = "p_val_adj",
                      gene_id_name = "gene") {

  # if we are doing overrepresentation and need to filter on p_adj
  if (type == "overrep") {

    if (lfc > 0) {

      filtered_df <- df[df[[lfc_name]] >= lfc & df[[padj_name]] <  p_adj, ]

    }

    if (lfc < 0) {

      filtered_df <- df[df[[lfc_name]] <= lfc & df[[padj_name]] <  p_adj, ]

    }

    if (lfc == 0) {

      filtered_df <- df[df[[padj_name]] <  p_adj, ]

    }

    filtered_df <- filtered_df[!is.na(filtered_df[[gene_id_name]]) &
                                !is.na(filtered_df[[gene_id_name]]), ]

    filtered_df <- filtered_df[!is.na(filtered_df[[padj_name]]), ]

    # rank df by p value
    filtered_df <- filtered_df[order(filtered_df[[padj_name]]), ]

  }

  if (type == "gsea") {

    filtered_df <- df[complete.cases(df), ]

    # rank by lfc
    filtered_df <- filtered_df[order(-filtered_df[[lfc_name]]), ]

  }

  return(filtered_df)
}

#' Make lfc vector for gsea
#'
#' @param filtered_df Results table parsed by filter_df.
#' @param col_name Name of column with gene symbols.
#' @param lfc_name Name of lfc col in results table.
#' @return Ranked named vector.
make_lfc_vector <- function(filtered_df, col_name = "gene",
                            lfc_name = "avg_log2FC") {

  named_lfc <- filtered_df[, lfc_name]

  names(named_lfc) <- filtered_df[, col_name]

  return(named_lfc)
}

#' Map Symbol to Entrez
#'
#' @param keys Vector of gene symbols.
#' @param org Org db.
#' @return Vector of Entrez IDs.
map_symbols_entrez <- function(keys, org = org.Hs.eg.db) {

  entrez_ids <- mapIds(x = org, keys = keys,
                      keytype = "SYMBOL", column = "ENTREZID",
                      multiVals = "first")


  return(entrez_ids)

}

#' Wrapper for clusterprofiler for running with custom supplied annotations.
#'
#' @param gene_ids_vector Vector of gene symbols.
#' @param term2gene User input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene.
#' @param species Species name for Org db (Mm or Hs).
#' @param org Org db.
#' @param analysis_type Analysis type - GSEA or overrep.
#' @param pval p-value threshold cutoff for tests to report.
#' @param padj p-adjusted threshold cutoff for tests to report.
#' @param pmethod p-adjust method - one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none" (Default: BH)
#' @return clusterProfiler results object
run_clusterprofiler_custom <- function(gene_ids_vector, term2gene,
                                       species = "Mm", org = org.Hs.eg.db,
                                       analysis_type = "overrep",
                                       pval = 0.05, padj = 0.25,
                                       pmethod = "BH") {

    if (analysis_type == "overrep") {

        gene_ids_vector <- map_symbols_entrez(gene_ids_vector, org = org)

        gene_ids_vector <- gene_ids_vector[complete.cases(gene_ids_vector)]

        cp_obj <- clusterProfiler::enricher(gene = as.character(gene_ids_vector),
                                            TERM2GENE = term2gene,
                                            pvalueCutoff = pval,
                                            qvalueCutoff = padj,
                                            pAdjustMethod = pmethod)

    }

    if (analysis_type == "gsea") {

        names(gene_ids_vector) <- map_symbols_entrez(names(gene_ids_vector),
                                                     org = org)

        gene_ids_vector <- gene_ids_vector[complete.cases(names(gene_ids_vector))]

        # for GSEA pvalueCutoff is adjusted p value cutoff
        cp_obj <- clusterProfiler::GSEA(gene = gene_ids_vector,
                                        TERM2GENE = term2gene,
                                        pvalueCutoff = padj,
                                        pAdjustMethod = pmethod)

    }

  return(cp_obj)

}

#' Wrap clusterprofiler calling
#'
#' @param gene_ids_vector List of gene symbol vectors.
#' @param term2gene User input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene.
#' @param species Species name for Org db (Mm or Hs).
#' @param analysis_type Analysis type - GSEA or overrep.
#' @param pval p-value threshold cutoff for tests to report.
#' @param padj p-adjusted threshold cutoff for tests to report.
#' @param pmethod p-adjust method - one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none" (Default: BH)
#' @return List of clusterProfiler results objects.
wrap_cp_custom <- function(gene_ids_vector, term2gene,
                           species = "Mm",
                           analysis_type = "overrep",
                           pval = 0.01, padj = 0.05, pmethod = "BH") {

    org <-  OrgDb(species)

    cp_obj_list <- lapply(gene_ids_vector, run_clusterprofiler_custom,
                          term2gene, species = species, org = org,
                           analysis_type = analysis_type,
                           pval = pval, padj = padj, pmethod = pmethod)

    return(cp_obj_list)

}

#' Parse clusterprofiler object
#'
#' @param cp_obj Clusterprofiler object
#' @param type Type of analysis conducted.
#' @param db Type of database used for analysis - Those run with GO database will have gene symbols output. 
#' @param species Species name for Org db (Mm or Hs).
#' @return Df of clusterprofiler results.
parse_cp_object <- function(cp_obj, type = "overrep",
                            db = "GO", species = "Mm") {

  org <-  OrgDb(species)

  cp_df <- as.data.frame(cp_obj)

  #remove NA

  cp_df <- cp_df[complete.cases(cp_df), ]

  # parse clusterprofiler object

  if (type == "overrep") {

    if (db != "GO") {

      cp_df$gene_name  <- sapply(cp_df$geneID, parse_entrez_col,
                                 map = T, org.db = org)

     }
  }
  if (type == "gsea") {

    if (db != "GO") {

      cp_df$gene_name <- sapply(cp_df$core_enrichment, parse_entrez_col,
                                map = T, org.db = org)

     }
  }

  return(cp_df)
}

#' Wrap clusterprofiler parse
#'
#' @param cp_obj_list List of clusterprofiler objects
#' @param type Type of analysis conducted.
#' @param db Type of database used for analysis - Those run with GO database will have gene symbols output.
#' @param species Species name for Org db (Mm or Hs).
#' @return List of clusterprofiler results as data frames.
wrap_parse_cp_obj <- function(cp_obj_list, type = "overrep",
                              db = "Custom", species = "Mm") {

    cp_df_list <- lapply(cp_obj_list, parse_cp_object,
                         type = type, db = db, species = species)

    return(cp_df_list)
}

#' Parse the entrez column of some of the outputs for clusterprofiler into gene symbols.
#'
#' @param entrez_id_str Entrez ids
#' @param map Does this col need to be mapped into symbols or not?
#' @param org.db Org db.
#' @return Str Gene symbols
parse_entrez_col <- function(entrez_id_str, map = F, org.db = org.Hs.eg.db) {

  symbol_str <- strsplit(entrez_id_str, "/")[[1]]

  if (map) {

    symbol_str <- suppressMessages(mapIds(x = org.db, keys = symbol_str,
                      keytype = "ENTREZID", column = "SYMBOL",
                      multiVals = "first"))

    symbol_str <- symbol_str[complete.cases(symbol_str)]

  }

  symbol_str_for_csv <- gsub(pattern = ", ", replacement = "/",
                             toString(symbol_str))

  return(symbol_str_for_csv)
}
  
#' Plot clusterprofiler object.
#'
#' @param cp_obj ClusterProfiler results object.
#' @param title Title of plot.
#' @param plot_type Type of plot (barplot, dotplot, gsea_enrichment).
#' @param x_axis What to plot on x axis (must be in ClusterProfiler results object).
#' @param row_gsea Row number of GSEA results to plot - default top.
#' @param toptable Number of terms to plot for barplot or dotplot.
#' @param size What variable to base the dotplot point size on (e.g. Count).
#' @param color What variable to base the dotplot point colour on (e.g. p.adjust).
#' @return Plot object
plot_fea_results <- function(cp_obj, title, plot_type = "barplot",
                             x_axis = "GeneRatio", row_gsea = 1,
                             toptable = 25, size = "Count",
                             color = "p.adjust") {

  if (plot_type == "barplot") {

        if (nrow(cp_obj@result) > 0) {

            # From clusterprofiler
            p <- graphics::barplot(cp_obj, showCategory = toptable) +
                 ggtitle(title)

        } else {

            # return a null ggplot if no results
            p <- ggplot() +
                 theme_void() +
                 geom_text(aes(0, 0, label = "N/A")) +
                 xlab(NULL)
        }
  }

  if (plot_type == "dotplot") {

    if (nrow(cp_obj@result) > 0) {

        # From clusterprofiler
        p <- enrichplot::dotplot(cp_obj, showCategory = toptable,
                                 x = x_axis, size = size,
                                 color = color) + ggtitle(title)

    } else {

        # return a null ggplot if no results
        p <- ggplot() +
            theme_void() +
            geom_text(aes(0, 0, label = "N/A")) +
            xlab(NULL)
    }
  }

  if (plot_type == "gsea_enrichment") {

    # From clusterprofiler
    if (nrow(cp_obj@result) > 0) {

        title <- paste0(title, "\n", trim_label(cp_obj@result[row_gsea, "ID"]))

        p <- enrichplot::gseaplot2(cp_obj, geneSetID = row_gsea, title = title)

    } else {

        # return a null ggplot if no results
        p <- ggplot() +
             theme_void() +
             geom_text(aes(0, 0, label = "N/A")) +
             xlab(NULL)
    }

  }
  return(p)
}

#' Generic barplot function
#'
#' @param cp_df ClusterProfiler result data frame.
#' @param colnames_select Name of columns to take from the dataframe.
#' @param colour Colour of plot.
#' @param title Title of plot.
#' @param x_log Should x axis be log scale?
#' @param z_log Should z axis be log scale?
#' @param x Name of x axis.
#' @param y Name of y axis.
#' @param x_thres Threshold to draw on x axis plot.
#' @param top Number of terms to plot.
#' @return Plot object
plot_barplot <- function(cp_df,
                         colnames_select = c("Description", "padj", "GeneRatio"),
                         colour = "steelblue", title = "Barplot GO Analysis",
                         x_log = T, z_log = F,
                         x = "-log10(Adjusted P.value)", y = "GO Term",
                         x_thres = 0.05, top = 10) {

  suppressMessages(library(ggplot2))

  # if given empty df
  if (nrow(cp_df) > 0) {

    df <- cp_df[, c(colnames_select)]

    df <- head(df, top)

    colnames(df) <- c("y", "x", "z")

    if (x_log) {

      df$x <- -log10(df$x)

    }

    if (x_log) {

      x_thres <- -log10(x_thres)

    }

    if (z_log) {

      df$z <- -log10(df$z)

    }

    if (class(df$z) == "numeric") {

        df$z <- round(df$z, 2)

    }

    df$y <- factor(df$y, levels = rev(df$y))

    p <- ggplot(data = df, aes(x = x, y = y)) +
         geom_bar(stat = "identity", fill = colour) +
         geom_text(aes(label = z), hjust = 1.6, color = "white", size = 3.5) +
         theme_minimal() +
         labs(title = title, x = x, y = y) +
         geom_vline(xintercept = x_thres, linetype = "dashed", color = "red") +
         theme_classic() +
         labs(fill = colour, x = trim_label(x), y = trim_label(y)) +
         scale_y_discrete(labels = function(y) trim_label(y))


  } else {

    # return a null ggplot if no results
    p <- ggplot() +
      theme_void() +
      geom_text(aes(0, 0, label = "N/A")) +
      xlab(NULL)
  }

  return(p)

}

#' Trim labels of barplot
#'
#' @param lab label
#' @param width Width threshold for trimming
#' @param delim Returns original where delim = NA
#' @return label trimmed
trim_label <- function(lab, width = 20, delim = NULL) {

  if (!is.null(delim)) {

    if (is.na(delim)) {

      return(lab)

    }
  }

  if (nchar(lab) > width) {

    lab <- paste0(substring(lab, 1, width), "...")

  }
  return(lab)
}