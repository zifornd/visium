
# Ggplot2 style violin plotting function for QC from 10x output
# Alternative to built in seurat functions
# See below for details on clean label function 

violingg <- function(seurat, ylab, colour_by= NULL, 
                     point_alpha = 0.6, point_size = NULL, ylim = NULL, 
                     show_median = F, xlab = "Sample",
                     clean = NULL, xdelim = NA, ydelim = NA, delim = NULL,
                     angle = 45, vjust = 1, clean_label = trim_label){
  
  suppressMessages(library(ggbeeswarm))
  
  df_to_plot <- data.frame(X = seurat[[xlab]], Y = seurat[[ylab]])
  
  if (!is.null(colour_by) && colour_by != "Sample") {

    df_to_plot[,colour_by] <- seurat[[colour_by]]

  }

  if (colour_by == "Sample") {

    df_to_plot[,"label"] <- seurat[[colour_by]]

  }

  if (angle == 90){

    vjust <- 0.5
  }

  colnames(df_to_plot) <- c("X", "Y", "label")
  
  plot_out <- ggplot(df_to_plot, aes(x=X, y=Y)) + 
    geom_violin(fill="grey90", alpha = 0.2, scale = "width", width = 0.8) + 
    xlab(xlab) + ylab(ylab) +
    geom_quasirandom(aes(color=label), fill = "grey20", alpha = 0.6, width=0.4, groupOnX=TRUE, bandwidth=1) +
    labs(col=colour_by)

  if (show_median) {

    plot_out <- plot_out + stat_summary(fun = median, fun.min = median, fun.max = median,
                                        geom = "crossbar", width = 0.3, alpha = 0.8)
  }
  
  
  if(!is.null(ylim)){

    plot_out <- plot_out + theme_classic() +
      theme(axis.text.x = element_text(angle = angle, vjust = vjust, hjust=1)) + 
      ylim(ylim)

  } else{

    plot_out <- plot_out + theme_classic() +
      theme(axis.text.x = element_text(angle = angle, vjust = vjust, hjust=1))

  }

  if(!is.null(clean)){

    plot_out <- plot_out + labs(fill=colour_by, x = clean_label(xlab, width = clean, delim = xdelim), 
                            y = clean_label(ylab, width = clean, delim = ydelim)) + 
      scale_x_discrete(labels = function(x) clean_label(x, width = as.numeric(clean), delim = delim))

  } 

  return(plot_out)
  
}

## Ggplot2 style barplot for 10x QC output

barplotgg <- function(seurat, ylab, colour_by = NULL, 
                      xlab = "Sample", ylim = NULL, clean = NULL, 
                      xdelim = NULL, ydelim = NULL, delim = NULL,
                      angle = 45, vjust = 1, clean_label = trim_label){
  
  df_to_plot <- data.frame(X = seurat[[xlab]], Y = seurat[[ylab]])
  
  if (!is.null(colour_by) && colour_by != "Sample") {

    df_to_plot[,colour_by] <- seurat[[colour_by]]

  }
  if (colour_by == "Sample") {

    df_to_plot[,"label"] <- seurat[[colour_by]]

  }
  if (angle == 90){

    vjust <- 0.5
  }

  colnames(df_to_plot) <- c("X", "Y", "label")

  plot_out <- ggplot(df_to_plot, aes(x=X, y=Y, fill=label)) + 
    geom_bar(stat="identity") + labs(fill=colour_by, x = xlab, y = ylab)
  
  if(!is.null(ylim)){

  plot_out <- plot_out + theme_classic() +
    theme(axis.text.x = element_text(angle = angle, vjust = vjust, hjust = 1)) + 
      ylim(ylim)

  } else {

  plot_out <- plot_out + theme_classic() +
      theme(axis.text.x = element_text(angle = angle, vjust = vjust, hjust = 1))

  }

  if(!is.null(clean)){

  plot_out <- plot_out + labs(fill=colour_by, x = clean_label(xlab, width = clean, delim = xdelim), 
                              y = clean_label(ylab, width = clean, delim = ydelim)) + 
            scale_x_discrete(labels = function(x) clean_label(x, width = clean, delim = delim))

  } 

return(plot_out)

}


# Try and deal with long labels by trimming them
# Returns original where delim = NA

trim_label <- function(lab, width = 20, delim = NULL){

  if(!is.null(delim)){

    if(is.na(delim)){
      
      return(lab)
    
    }
  } 
  
  if(nchar(lab) > width){

    lab <- paste0(substring(lab, 1, width), "...")

  }

  return(lab)
}

# Try and deal with long labels by wrapping them
# Returns original where delim = NA

wrap_label <- function(lab, width = 20, delim = NULL){
  
  library(stringr)

  if(!is.null(delim)){

    if(is.na(delim)){
      
      return(lab)
    
    }
  } 

  if(nchar(lab) > width){

    if(" " %in% lab){

      lab <- stringr::str_wrap(lab, width = width)

    }

    if(!is.null(delim)){
      
      if(delim == "."){

        delim <- "[.]"

      }

      lab <- stringr::str_wrap(paste(unlist(strsplit(lab, delim)), collapse = " "), width = width)

    }

    if(is.null(delim) & !(" " %in% lab)) {
      
      lab <- paste0(substring(lab, 1, width), "\n", substring(lab, (width + 1), nchar(lab)))
      
    }
    
  }
  return(lab)
}