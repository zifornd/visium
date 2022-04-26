
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
    geom_violin(fill="grey90", alpha = 0.2, scale = "width", width = 0.8) + xlab(xlab)  + ylab(ylab) +
    geom_quasirandom(aes(color=label), fill = "grey20", alpha = 0.6, width=0.4, groupOnX=TRUE, bandwidth=1) +
    labs(col=colour_by)
  if (show_median) {
    plot_out <- plot_out + stat_summary(fun = median, fun.min = median, fun.max = median,
                                        geom = "crossbar", width = 0.3, alpha = 0.8)
  }
  
  #point_out <- .get_point_args(colour_by, shape_by, size_by, alpha = point_alpha, size = point_size)
  
  #print(point_out)
  
  # point_FUN <- function(...) geom_quasirandom(..., width=0.4, groupOnX=TRUE, bandwidth=1)
  
  
  # plot_out <- plot_out + do.call(point_FUN, point_out$args)
   
  # plot_out  <-  plot_out + geom_quasirandom(color=label, fill = "grey20", alpha = 0.6, width=0.4, groupOnX=TRUE, bandwidth=1)
  
  
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





#' @importFrom ggplot2 aes_string
.get_point_args <- function(colour_by, shape_by, size_by, alpha=0.65, size=NULL) 
  ## Note the use of colour instead of fill when shape_by is set, as not all shapes have fill.
  ## (Fill is still the default as it looks nicer.)
{
  aes_args <- list()
  fill_colour <- FALSE
  if (!is.null(shape_by)) {
    aes_args$shape <- "shape_by"
    fill_colour <- FALSE
  }
  if (!is.null(colour_by)) {
    if (fill_colour) {
      aes_args$fill <- "colour_by"
    } else {
      aes_args$colour <- "colour_by"
    }
  }
  if (!is.null(size_by)) {
    aes_args$size <- "size_by"
  }
  new_aes <- do.call(aes_string, aes_args)
  
  geom_args <- list(mapping=new_aes, alpha=alpha)
  if (is.null(colour_by) || fill_colour) {
    geom_args$colour <- "grey70"
  }
  if (is.null(colour_by) || !fill_colour) { # set fill when there is no fill colour, to distinguish between e.g., pch=16 and pch=21.
    geom_args$fill <- "grey20"
  }
  if (is.null(shape_by)) {
    geom_args$shape <- 19
  }
  if (is.null(size_by)) {
    geom_args$size <- size
  }
  return(list(args=geom_args, fill=fill_colour))
}