volcanoPlot <- function(res, 
                        clip_nlog10p = 10, # Upper limit for negative log10-transformed p-value
                        clip_log2fc = 2.5, # Upper limit for absolute log2-fold change
                        min_log2FoldChange = 0, # Minimum absolute log2-fold change to consider DE
                        padj_threshold = 0.05,# Adjusted p-value threshold for differential expression
                        n_label_genes = NULL
) {
  
  # Add a new column for negative log10-transformed p-value
  # Clip values in nlog10p column at clip_nlog10p
  # Add new column for log2-fold change threshold
  # Clip values in log2FoldChange column at clip_log2fc
  # Create a new column for differential expression status
  # Add external gene name column if gene is differentially expressed
  plotting_data <- res %>% 
    mutate(
      nlog10p = ifelse(nlog10p > clip_nlog10p, clip_nlog10p, nlog10p),
      fc_threshold = abs(log2FoldChange) > min_log2FoldChange,
      log2FoldChange =  ifelse(abs(log2FoldChange) > abs(clip_log2fc), sign(log2FoldChange)*clip_log2fc, log2FoldChange),
      de_status = NULL,
      de_status = case_when(
        padj < padj_threshold & (log2FoldChange >= min_log2FoldChange) ~ "upregulated",
        padj < padj_threshold & (log2FoldChange <= -min_log2FoldChange) ~ "downregulated"
      ),
      sig_external_gene_name = ifelse(!is.na(de_status), external_gene_name, NA)
    )
  
  # Calculate the number of genes for each differential expression status
  # Concatenate the counts to the plotting_data table
  holder <- plotting_data %>% 
    summarise(n = n(), .by = "de_status") %>% 
    mutate(Status = paste0(de_status, " (n = ", n,")"))
  plotting_data <- plotting_data %>% 
    left_join(., holder, by= "de_status")
  
  if(!is.null(n_label_genes)){
    topbottom <- rbind(slice_max(res, n = n_label_genes, order_by = log2FoldChange),
                         slice_min(res, n = n_label_genes, order_by = log2FoldChange)) %>%
      pull(external_gene_name)
    
    plotting_data <- plotting_data %>% mutate(label_genes = ifelse(external_gene_name %in% topbottom, external_gene_name, NA))
  }
  # Generate the volcano plot using ggplot2
  # X-axis: log2-fold change
  # Y-axis: negative log10-transformed p-value
  # Points colored based on differential expression status
  # Add vertical lines at +/- min_log2FoldChange
  # Add horizontal line at -log10(padj_threshold)
  # Use a black and white theme
  pl <- ggplot(plotting_data) + 
    geom_point(aes(x = log2FoldChange, y = nlog10p, color = Status)) +
    geom_vline(xintercept = min_log2FoldChange) + 
    geom_vline(xintercept = -min_log2FoldChange) + 
    geom_hline(yintercept = -log10(padj_threshold)) + 
    scale_color_manual(values = c("blue","grey", "red")) + 
    theme_bw() 
  
  # Return the ggplot2 object
  return(pl)
}
