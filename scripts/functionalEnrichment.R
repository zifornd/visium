
# Load OrgDb package
# e.g. x = "Hs"

OrgDb <- function(x) {

    pkg <- paste0("org.", x, ".eg.db")

    obj <- getFromNamespace(pkg, pkg)

}

# Parses results table
# Splits results table by group e.g. cluster annotation
# Filters results with filter_df function ready for downstream ovverep or gsea
# Then outputs vector of gene names or ranked named vector 

parse_res <- function(res, group = "cluster", p_adj = Inf, lfc = 0, type= "overrep", 
                      lfc_name = "avg_log2FC", padj_name = "p_val_adj", gene_id_name = "gene",
                      split.res = TRUE){
    if(split.res){
        resList <- split_res(res, group = group)
    } else {
        resList <- list(res)
    }                     

    filtered_res <- lapply(resList, filter_df, p_adj = p_adj, lfc = lfc, type= type, 
                           lfc_name = lfc_name, padj_name = padj_name, gene_id_name = gene_id_name)
    
    if(type == "overrep"){
        vector_res <- lapply(filtered_res, function(x, gene_id_name){x[,gene_id_name]}, gene_id_name = gene_id_name)
    }

    if(type == "gsea"){
        vector_res <-  lapply(filtered_res, make_lfc_vector, col_name = gene_id_name, lfc_name = lfc_name)
    }

    return(vector_res)

}

# Splits results into list of results by given group column

split_res <- function(res, group = "cluster"){

    groups <- unique(res[[group]])

    resList <- lapply(groups, function(x, res, group){ res[res[[group]] == x,]}, res = res, group = group)

    names(resList) <- groups

    return(resList)
}


# Filter DESeq2 df prior to functional enrichment analysis (e.g. p_adj, LFC, no gene symbol etc.)

# Note - GSEA does not require pre-filtering based on LFC or p-value thresholding - but will need valid gene sybmols or ens_ids
# Standard options here for filtering are adjusted p value of 0.05 
# LFC is something people often do as well. 
# Statistically speaking genes "of interest" should be defined by a threshold defined a priori - based on study design and rationale
# this is classically an alpha value of 0.05 

# If you dont want to filter on lfc - leave as 0
# if you dont want to filter on p value leave as Inf

filter_df <- function(df, p_adj = Inf, lfc = 0, type= "overrep", 
                      lfc_name = "avg_log2FC", padj_name = "p_val_adj", gene_id_name = "gene"){
  
  # if we are doing overrepresentation and need to filter on p_adj 
  if(type == "overrep"){

    if(lfc > 0){
      filtered_df <- df[df[[lfc_name]] >= lfc & df[[padj_name]] <  p_adj,]
    }
    if(lfc < 0){
      filtered_df <- df[df[[lfc_name]] <= lfc & dedfseq_df[[padj_name]] <  p_adj,]
    }
    if(lfc == 0){
      filtered_df <- df[df[[padj_name]] <  p_adj,]
      
    }

    filtered_df <- filtered_df[ !is.na(filtered_df[[gene_id_name]]) & !is.na(filtered_df[[gene_id_name]]) ,]

    filtered_df <- filtered_df[ !is.na(filtered_df[[padj_name]]),]

    # rank df by p value
    filtered_df <- filtered_df[order(filtered_df[[padj_name]]),]

  }
  
  if(type == "gsea"){
    
    filtered_df <- df[complete.cases(df),]
   
    # rank by lfc
    filtered_df <- filtered_df[order(-filtered_df[[lfc_name]]),]
    
  }

  return(filtered_df)
}

# Make lfc vector for gsea 

make_lfc_vector <- function(filtered_df, col_name = "gene", lfc_name = "avg_log2FC"){
  
  named_lfc = filtered_df[,lfc_name]

  names(named_lfc) = filtered_df[,col_name]
  
  return(named_lfc)
}



# Map Symbol to Entrez

map_symbols_entrez <- function(keys, org = org.Hs.eg.db){
  
  entrez_ids = mapIds(x = org, keys = keys, 
                      keytype = "SYMBOL", column = "ENTREZID", 
                      multiVals = "first")
  
  #entrez_ids = entrez_ids[complete.cases(entrez_ids)] 
  
  return(entrez_ids)
  
}

# Wrapper for clusterprofiler for running with custom supplied annotations

run_clusterprofiler_custom <- function(gene_ids_vector, term2gene, 
                                       species = "Mm", org = org.Hs.eg.db,
                                       analysis_type = "overrep", 
                                       pval = 0.01, padj = 0.05, pmethod = "BH"){

    if(analysis_type == "overrep"){
        
        gene_ids_vector <- map_symbols_entrez(gene_ids_vector, org = org)

        gene_ids_vector <- gene_ids_vector[complete.cases(gene_ids_vector)] 

        cp_obj <- enricher(gene = as.character(gene_ids_vector), TERM2GENE = term2gene)

    }

    if(analysis_type == "gsea"){
        
        names(gene_ids_vector) <- map_symbols_entrez(names(gene_ids_vector), org = org)
        
        gene_ids_vector <- gene_ids_vector[complete.cases(names(gene_ids_vector))]

        cp_obj <- GSEA(gene = gene_ids_vector, TERM2GENE = term2gene)

    }

  return(cp_obj)

}

# Wrap clusterprofiler calling and parse output 

wrap_cp_custom <- function(gene_ids_vector, term2gene, 
                           species = "Mm",
                           analysis_type = "overrep", 
                           pval = 0.01, padj = 0.05, pmethod = "BH"){
    
    org <-  OrgDb(species)

    cp_obj_list <- lapply(gene_ids_vector, run_clusterprofiler_custom, term2gene, 
                           species = species, org = org,
                           analysis_type = analysis_type, 
                           pval = pval, padj = padj, pmethod = pmethod)
    
    return(cp_obj_list)

}


# Wrap clusterprofiler calling and parse output 

wrap_parse_cp_obj <- function(cp_obj_list, type = "overrep", db = "Custom", species = "Mm"){
    
    cp_df_list <- lapply(cp_obj_list, parse_cp_object, type = type, db = db, species = species)

    return(cp_df_list)

}

#############################################
####  Parse clusterprofiler object ##########
#############################################

# Input: clusterprofiler object
# Output: Df of clusterprofiler results

parse_cp_object <- function(cp_obj, type = "overrep", db = "GO", species = "Mm"){
  
  org <-  OrgDb(species)

  cp_df = as.data.frame(cp_obj)
  
  #remove NA

  cp_df = cp_df[complete.cases(cp_df),]
  
  # parse clusterprofiler object 

  if(type == "overrep"){

    if(db != "GO"){

      cp_df$gene_name  = sapply(cp_df$geneID, parse_entrez_col, map=T, org.db = org)

     }
  }
  if(type == "gsea"){

    if(db != "GO"){

      cp_df$gene_name = sapply(cp_df$core_enrichment, parse_entrez_col, map=T, org.db = org)
      
     }
    
  }

  return(cp_df)
}


parse_entrez_col <- function(entrez_id_str, map=F, org.db = org.Hs.eg.db){
  
  symbol_str = strsplit(entrez_id_str, "/")[[1]]
  
  if(map){

    symbol_str = suppressMessages(mapIds(x = org.db, keys = symbol_str, 
                      keytype = "ENTREZID", column = "SYMBOL", 
                      multiVals = "first"))

    symbol_str = symbol_str[complete.cases(symbol_str)] 

  }
  
  symbol_str_for_csv = gsub(pattern = ", ", replacement="/", toString(symbol_str))
  
  return(symbol_str_for_csv)
}
  
#############################################
####  Write clusterprofiler object ##########
#############################################

# Input: Df of clusterprofiler results
# Output: Text file of enriched go/pathways/terms 

write_fea_results <- function(cp_df, filename, file = "csv"){
  
  # This includes both the terms themselves which are ranked
  # But also the genes of interest
  
  if(file == "csv"){separator = ","}
  if(file == "tsv"){separator = "\t"}
  
  write.table(x=cp_df, file = filename, quote = F, sep = separator, col.names = T, row.names = F)
  
}

###########################################
####  Plot clusterprofiler object #########
###########################################

# Input: Df of clusterprofiler results following parsing
# Output: User specified plot

plot_fea_results <- function(cp_obj, title, plot_type = "barplot", top_num=15, x_axis="GeneRatio", 
                             row_gsea = 1, toptable=25, size="Count", color = "p.adjust", width=5, height=10){
  library(enrichplot)
  
  if(plot_type == "barplot"){
    
        if(nrow(cp_obj@result) > 0){

            # From clusterprofiler
            p <- barplot(cp_obj, showCategory=toptable)  + ggtitle(title)

        } else {

            # return a null ggplot if no results
            p <- ggplot() +
                theme_void() +
                geom_text(aes(0,0,label='N/A')) +
                xlab(NULL)
        }
  }

  if(plot_type == "dotplot"){
    
    if(nrow(cp_obj@result) > 0){

        # From clusterprofiler
        p <- dotplot(cp_obj, showCategory=toptable, x = x_axis, size =size, color=color ) + ggtitle(title)

    } else {

        # return a null ggplot if no results
        p <- ggplot() +
            theme_void() +
            geom_text(aes(0,0,label='N/A')) +
            xlab(NULL)
    }
  }

  if(plot_type == "gsea_enrichment"){
    
    # From clusterprofiler
    if(nrow(cp_obj@result) > 0){

        title <- paste0(title, "\n", trim_label(cp_obj@result[row_gsea,"ID"]))
        
        p <- gseaplot2(cp_obj, geneSetID = row_gsea, title = title) 

    } else {

        # return a null ggplot if no results
        p <- ggplot() +
             theme_void() +
             geom_text(aes(0,0,label='N/A')) +
             xlab(NULL)
    }

  }
  return(p)
}

######### Generic barplot function

plot_barplot <- function(cp_df, colnames_select = c("Description", "padj", "GeneRatio"), colour = "steelblue", 
                         title = "Barplot GO Analysis", x_log = T, z_log=F,
                         x = "-log10(Adjusted P.value)", y = "GO Term", 
                         x_thres=0.05, top = 10){
  
  suppressMessages(library(ggplot2))
  
  df = cp_df[,c(colnames_select)]
  
  df <- head(df, top)

  colnames(df)  = c("y", "x", "z")
  
  if(x_log){df$x = -log10(df$x)}
  if(x_log){x_thres = -log10(x_thres)}
  if(z_log){df$z = -log10(df$z)}
  
  if(class(df$z) == "numeric"){

      df$z <- round(df$z, 2)

  }

  df$y = factor(df$y, levels = rev(df$y))

  p = ggplot(data=df, aes(x=x, y= y)) +
    geom_bar(stat="identity", fill=colour)+
    geom_text(aes(label=z), hjust=1.6, color="white", size=3.5)+
    theme_minimal() + labs(title=title, x=x, y = y) + 
    geom_vline(xintercept = x_thres, linetype="dashed", color = "red") +
    theme_classic() + labs(fill=colour, x = trim_label(x), y = trim_label(y)) +
    scale_y_discrete(labels = function(y) trim_label(y))
    
    #scale_fill_continuous(low="#fee0d2", high="#de2d26")
  
  return(p)

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

##########################################################################
#### Wrapper for clusterprofiler which does all the actual analysis ######
##########################################################################

# Input: deseq_df, filtered DF
# Output: clusterprofiler object

# Other Options 
# Analysis type: gsea, overrep
# DB choice: GO,KEGG, REACTOME, Hallmark
# DB sub cat fo GO as ont:  ALL, BP, CC, MF
# P.val threshold: standard as 0.05 (set to "Inf" if no filter wanted)

# org.Hs.eg.db as default is org.Hs.eg.db and should be loaded prior 
# - can be alternative species if needed

# Background is imporant to define in function enrichment and can effect results
# by default all expressed genes is often a standard approach 

run_clusterprofiler <- function(deseq_df, filtered_df, col_name = "gene_name", 
                                analysis_type = "overrep", 
                                database = "GO", ont = "ALL", kegg_org = "hsa",reactomepa = "human",
                                org.db = org.Hs.eg.db, keytype = "SYMBOL",
                                pval = 0.01, padj = 0.05, pmethod = "BH"){
  
  # for range of overrepresentation and GSEA approaches:
  suppressMessages(library(clusterProfiler))
  
  # DB (Can be changed for other species):
  suppressMessages(library(org.Hs.eg.db))

  if(analysis_type == "overrep"){
    
    print("Running Overrepresentation test...")
    
    #### In case a run needs entrez ids:
    
    entrez_ids = mapIds(x = org.Hs.eg.db, keys = filtered_df[,col_name], 
                        keytype = "SYMBOL", column = "ENTREZID", 
                        multiVals = "first")
    
    entrez_ids_universe = mapIds(x = org.Hs.eg.db, keys = deseq_df[,col_name], 
                                 keytype = "SYMBOL", column = "ENTREZID", 
                                 multiVals = "first")
    
    if(database == "GO"){
      
      cp_obj <- enrichGO(gene       = filtered_df[,col_name],
                      universe      = deseq_df[,col_name],
                      keyType       = keytype,
                      OrgDb         = org.db,
                      ont           = ont,
                      pAdjustMethod = pmethod,
                      pvalueCutoff  = pval,
                      qvalueCutoff  = padj)
    }

    if(database == "KEGG"){
      
      # Needs entrez id
      cp_obj <- enrichKEGG(gene     = entrez_ids[complete.cases(entrez_ids)],
                           universe = entrez_ids_universe[complete.cases(entrez_ids_universe)],
                           organism = kegg_org,
                           pAdjustMethod = pmethod,
                           pvalueCutoff  = pval,
                           qvalueCutoff  = padj)
        
    }
    
    if(database == "ReactomePA"){
      
      suppressMessages(library(ReactomePA))
      
      # Needs entrez id
      cp_obj <- enrichPathway(gene = entrez_ids[complete.cases(entrez_ids)], 
                              universe = entrez_ids_universe[complete.cases(entrez_ids_universe)],
                              organism = reactomepa,
                              pAdjustMethod = pmethod,
                              pvalueCutoff  = pval,
                              qvalueCutoff  = padj)
      
    }
    
  }else{
    
    print("Running GSEA...")
    
    # must be ordered by LFC 
    # if filtered by prior function - should already be sorted by descending lfc
    named_lfc = filtered_df[,"log2FoldChange"]
    names(named_lfc) = filtered_df[,col_name]
    
    ## Incase needs as entrez ids
    named_lfc_entrez = named_lfc
    entrez_ids = mapIds(x = org.Hs.eg.db, keys = names(named_lfc), 
                        keytype = "SYMBOL", column = "ENTREZID", 
                        multiVals = "first")
    names(named_lfc_entrez) = entrez_ids
    # remove those without a entrez id
    named_lfc_entrez = named_lfc_entrez[complete.cases(names(named_lfc_entrez))]
    
    if(database == "GO"){
      
      cp_obj <- gseGO(geneList     = named_lfc,
                      OrgDb        = org.db,
                      keyType       = keytype,
                      ont          = ont,
                      minGSSize    = 50,
                      maxGSSize    = 500,
                      eps = 0,
                      nPermSimple = 20000,
                      pvalueCutoff = pval,
                      pAdjustMethod = pmethod,
                      verbose      = FALSE) 
      
    }
    
    if(database == "KEGG"){
      
      # KEYTYPE HERE HAS TO BE SOMETHING ELSE!
      # one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
      
      cp_obj <- gseKEGG(geneList     = named_lfc_entrez,
                        organism     = kegg_org,
                        minGSSize    = 50,
                        maxGSSize    = 500,
                        eps = 0,
                        nPermSimple = 20000,
                        pvalueCutoff = pval,
                        pAdjustMethod = pmethod,
                        verbose      = FALSE)
        
    }
    
    if(database == "ReactomePA"){
      
      # I think gene might need to be as entrez gene id...
      
      cp_obj <- gsePathway(named_lfc_entrez, 
                           organism = reactomepa,
                           minGSSize    = 50,
                           maxGSSize    = 500,
                           eps = 0,
                           nPermSimple = 20000,
                           pvalueCutoff = pval,
                           pAdjustMethod = pmethod,
                           verbose = FALSE)
        
    }
  }
  
  return(cp_obj)
}