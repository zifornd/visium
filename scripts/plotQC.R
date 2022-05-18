
# some plotting functions for qc plotting 

violingg <- function(seurat, ylab, colour_by= NULL, 
                             point_alpha = 0.6, point_size = NULL, ylim = NULL, 
                             show_median = F, xlab = "Sample"){
  
  suppressMessages(library(ggbeeswarm))
  
  df_to_plot <- data.frame(X = seurat[[xlab]], Y = seurat[[ylab]])
  
  if (!is.null(colour_by) && colour_by != "Sample") {

    df_to_plot[,colour_by] <- seurat[[colour_by]]

  }

  if (colour_by == "Sample") {

    df_to_plot[,"label"] <- seurat[[colour_by]]

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

    plot_out + theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      ylim(ylim)

  } else{

    plot_out + theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  }
  
}


barplotgg <- function(seurat, ylab, colour_by= NULL, 
                      xlab = "Sample", ylim=NULL){
  
  df_to_plot <- data.frame(X = seurat[[xlab]], Y = seurat[[ylab]])
  
  if (!is.null(colour_by) && colour_by != "Sample") {

    df_to_plot[,colour_by] <- seurat[[colour_by]]

  }
  if (colour_by == "Sample") {

    df_to_plot[,"label"] <- seurat[[colour_by]]

  }
  
  colnames(df_to_plot) <- c("X", "Y", "label")
  
  plot_out <- ggplot(df_to_plot, aes(x=X, y=Y, fill=label)) + 
    geom_bar(stat="identity") + xlab(xlab)  + ylab(ylab) +
    labs(fill=colour_by)
  
  if(!is.null(ylim)){

  plot_out + theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      ylim(ylim)

  } else{

    plot_out + theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  }
}