
#' Ggplot2 style violin plotting function for QC from 10x output
#' Alternative to built in seurat functions
#'
#' @param seurat Seurat object
#' @param ylab Y axis label intended to plot present within seurat object meta.data
#' @param colour_by variable to colour violin plot by - present within seurat meta.data
#' @param ylim Sets Y axis limits
#' @param show_median Adds median to violin plot
#' @param xlab X axis label intended to plot present within seurat object meta.data (Default: Sample)
#' @param clean Maximum width of the clean label function - if set labels will be managed by function supplied to clean_label
#' @param xdelim Choice of deliminator to wrap_label by for xaxis label - if set to NA will return original label
#' @param ydelim Choice of deliminator to wrap_label by for yaxis label - if set to NA will return original label
#' @param delim Choice of deliminator to wrap_label by for xaxis tick labels - if set to NA will return original label
#' @param custom_lab custom labels of choice if want to overrite xaxis labels.
#' @param angle Angle for x axis tick labels
#' @param vjust Vertical adjustment for x axis tick labels
#' @param clean_label function used to tidy up plot labels (trim_label shortens, wrap_label forces line breaks)
#' @return ggplot object.
violingg <- function(seurat, ylab, colour_by = NULL,
                     ylim = NULL, show_median = F, xlab = "Sample",
                     clean = NULL, xdelim = NA, ydelim = NA, delim = NULL, custom_lab = NULL,
                     angle = 45, vjust = 1, clean_label = trim_label) {

  suppressMessages(library(ggbeeswarm))

  df_to_plot <- data.frame(X = seurat[[xlab]], Y = seurat[[ylab]])

  if (!is.null(colour_by) && colour_by != "Sample") {

    df_to_plot[, colour_by] <- seurat[[colour_by]]

  }

  if (colour_by == "Sample") {

    df_to_plot[, "label"] <- seurat[[colour_by]]

  }

  if (angle == 90) {

    vjust <- 0.5
  }

  colnames(df_to_plot) <- c("X", "Y", "label")

  plot_out <- ggplot(df_to_plot, aes(x = X, y = Y)) +
    geom_violin(fill = "grey90", alpha = 0.2, scale = "width", width = 0.8) +
    xlab(xlab) + ylab(ylab) +
    geom_quasirandom(aes(color=label), fill = "grey20", alpha = 0.6, width=0.4, groupOnX=TRUE, bandwidth=1) +
    labs(col = colour_by)

  if (show_median) {

    plot_out <- plot_out + stat_summary(fun = median,
                                        fun.min = median, fun.max = median,
                                        geom = "crossbar",
                                        width = 0.3, alpha = 0.8)
  }


  if (!is.null(ylim)) {

    plot_out <- plot_out +
                theme_classic() +
                theme(axis.text.x = element_text(angle = angle,
                                                 vjust = vjust, hjust = 1)) +
                ylim(ylim)

  } else {

    plot_out <- plot_out +
                theme_classic() +
                theme(axis.text.x = element_text(angle = angle,
                                                 vjust = vjust, hjust = 1))

  }

  if (!is.null(clean)) {

    plot_out <- plot_out +
                labs(fill = colour_by,
                     x = clean_label(xlab, width = clean, delim = xdelim),
                     y = clean_label(ylab, width = clean, delim = ydelim, custom_lab = custom_lab)) +
                scale_x_discrete(labels = function(x) clean_label(x, width = as.numeric(clean), delim = delim, custom_lab = custom_lab))

  }

  return(plot_out)

}

#' Ggplot2 style Barplots for QC from 10x output
#'
#' @param seurat Seurat object
#' @param ylab Y axis label intended to plot present within seurat object meta.data
#' @param colour_by variable to colour violin plot by - present within seurat meta.data
#' @param xlab X axis label intended to plot present within seurat object meta.data (Default: Sample)
#' @param ylim Sets Y axis limits
#' @param clean Maximum width of the clean label function - if set labels will be managed by function supplied to clean_label
#' @param xdelim Choice of deliminator to wrap_label by for xaxis label - if set to NA will return original label
#' @param ydelim Choice of deliminator to wrap_label by for yaxis label - if set to NA will return original label
#' @param delim Choice of deliminator to wrap_label by for xaxis tick labels - if set to NA will return original label
#' @param custom_lab custom labels of choice if want to overrite xaxis labels.
#' @param angle Angle for x axis tick labels
#' @param vjust Vertical adjustment for x axis tick labels
#' @param clean_label function used to tidy up plot labels (trim_label shortens, wrap_label forces line breaks)
#' @return ggplot object.
barplotgg <- function(seurat, ylab, colour_by = NULL,
                      xlab = "Sample", ylim = NULL, clean = NULL,
                      xdelim = NULL, ydelim = NULL, delim = NULL, custom_lab = NULL,
                      angle = 45, vjust = 1, clean_label = trim_label) {

  df_to_plot <- data.frame(X = seurat[[xlab]], Y = seurat[[ylab]])

  if (!is.null(colour_by) && colour_by != "Sample") {

    df_to_plot[, colour_by] <- seurat[[colour_by]]

  }
  if (colour_by == "Sample") {

    df_to_plot[, "label"] <- seurat[[colour_by]]

  }
  if (angle == 90) {

    vjust <- 0.5
  }

  colnames(df_to_plot) <- c("X", "Y", "label")

  plot_out <- ggplot(df_to_plot, aes(x = X, y = Y, fill = label)) +
    geom_bar(stat = "identity") + labs(fill = colour_by, x = xlab, y = ylab)

  if (!is.null(ylim)) {

  plot_out <- plot_out +
              theme_classic() +
              theme(axis.text.x = element_text(angle = angle, vjust = vjust, hjust = 1)) + 
              ylim(ylim)

  } else {

  plot_out <- plot_out + theme_classic() +
      theme(axis.text.x = element_text(angle = angle, vjust = vjust, hjust = 1))

  }

  if (!is.null(clean)) {

  plot_out <- plot_out +
              labs(fill = colour_by,
                   x = clean_label(xlab, width = clean, delim = xdelim),
                   y = clean_label(ylab, width = clean, delim = ydelim, custom_lab = custom_lab)) +
              scale_x_discrete(labels = function(x) clean_label(x, width = clean, delim = delim, custom_lab = custom_lab))

  }

return(plot_out)

}

#' Try and deal with long labels by returning custom labels
#'
#' @param lab string label
#' @param custom_lab named vector of custom labels of choice (names are based on original labels)
#' @return custom labels of choice
custom_label <- function(lab, custom_lab = NULL, ...) {

  if (!is.null(custom_lab)) {

    stopifnot((lab %in% names(custom_lab)))

    new_labs <- as.character(custom_lab[lab])

    return(new_labs)

  } else {

    return(lab)

  }

}

#' Prepares custom label character vector prior to plotting.
#'
#' @param custom_names character vector names of plot y axis/ x axis ticks to change to.
#' @param original_names character vector original names of plot y axis/ x axis ticks. 
#' @param custom_sample_names character vector of custom sample names.
#' @param sample_ids character vector of original sample ids.
#' @return named character vector of custom names.
get_custom_names <- function(custom_names, original_names,
                             custom_sample_names = NULL, sample_ids = NULL) {

  # if custom sample names supplied - add these to the named vector
  # else supply original sample ids.
  if (!is.null(custom_sample_names)) {

    if (length(custom_sample_names) > 0) {

      custom_names_samp <- c(custom_sample_names, custom_names)

    # if custom_sample_names is empty list
    } else {

      custom_names_samp <- c(sample_ids, custom_names)

    }
    
  # if no custom sample names supplied
  } else {

    custom_names_samp <- c(sample_ids, custom_names)

  }

  # name the vector as original names (original sample ids plus original labels)
  names(custom_names_samp) <- c(sample_ids, original_names)

  return(custom_names_samp)

}


#' Try and deal with long labels by trimming them
#'
#' @param lab string label
#' @param width Maximum width for a label before being trimmed and "..." appended
#' @param delim If NA then original label will be returned
#' @return trimmed label.
trim_label <- function(lab, width = 20, ...) {

  # lab may be character vector with > 1 element
  if (T %in% (nchar(lab) > width)) {

    lab <- paste0(substring(lab, 1, width), "...")

  }

  return(lab)
}

#' Try and deal with long labels by wrapping them
#' Returns original where delim = NA
#'
#' @param lab string label
#' @param width Maximum width for a label before being wrapped
#' @param delim Deliminator used to split text - If NA then original label will be returned
#' @return wrapped label.
wrap_label <- function(lab, width = 20, delim = NULL, ...) {

  library(stringr)

  if (!is.null(delim)) {

    if (is.na(delim)) {

      return(lab)

    }
  }

  if (nchar(lab) > width) {

    if (" " %in% lab) {

      lab <- stringr::str_wrap(lab, width = width)

    }

    if (!is.null(delim)) {

      if (delim == ".") {

        delim <- "[.]"

      }

      lab <- stringr::str_wrap(paste(unlist(strsplit(lab, delim)), collapse = " "),
                               width = width)

    }

    if (is.null(delim) & !(" " %in% lab)) {

      lab <- paste0(substring(lab, 1, width), "\n", substring(lab, (width + 1), nchar(lab)))

    }

  }
  return(lab)
}