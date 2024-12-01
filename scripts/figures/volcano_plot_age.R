library(DESeq2)
library(tidyverse)
library(here)
source("/mnt/morbo/Data/Users/jslosberg/aged_bulkrna/scripts/accessory/volcanoPlot.R", echo=TRUE)

dds <- readRDS(here("data/bulkrna_dds.rds"))

testnames <- resultsNames(dds)[-c(1)]
gene_df <- rowData(dds)[,c("ensembl_gene_id_trimmed","external_gene_name")] %>% as.data.frame()

res_all <- map_df(testnames, function(comp){
  res <- results(dds, name = comp) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ensembl_gene_id_trimmed") %>% left_join(., gene_df) %>%
    mutate(comparison = comp,
           nlog10p = -log10(padj),
           is_sig = padj < 0.05) %>% 
    mutate(sig_external_gene_name = ifelse(is_sig, external_gene_name, NA))
  return(res)
})

res_all_shrink <- map_df(testnames, function(comp){
  res <- lfcShrink(dds, coef = comp) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ensembl_gene_id_trimmed") %>% left_join(., gene_df) %>%
    mutate(comparison = comp,
           nlog10p = -log10(padj),
           is_sig = padj < 0.05)
  
  return(res)
})

res_17mo_vs_6mo <- results(dds, contrast=c("age", "17mo", "6mo"), alpha=0.1, pAdjustMethod="BH", parallel=TRUE, tidy = T)  %>%
  rownames_to_column(var = "ensembl_gene_id") %>% left_join(., gene_df) %>%
  mutate(comparison = "age_17mo_vs_6mo",
         nlog10p = -log10(padj),
         is_sig = padj < 0.05)

colData(dds)$age <- factor(colData(dds)$age, levels = c("6mo","17mo","P30"))
dds <- DESeq(dds)
res_17mo_vs_6mo_shrink <- lfcShrink(dds, coef = "age_17mo_vs_6mo") %>%
  as.data.frame() %>%
  rownames_to_column(var = "ensembl_gene_id_trimmed") %>% left_join(., gene_df) %>%
  mutate(comparison = "age_17mo_vs_6mo",
         nlog10p = -log10(padj),
         is_sig = padj < 0.05)

#add this comparison to the others
assertthat::assert_that(all(colnames(res_all_shrink) == colnames(res_17mo_vs_6mo_shrink)))

res_6mo <- res_all %>% filter(comparison == "age_6mo_vs_P30")
res_6mo_shrink <- res_all_shrink %>% filter(comparison == "age_6mo_vs_P30") 

res_17mo <- res_all %>% filter(comparison == "age_17mo_vs_P30")
res_17mo_shrink <- res_all_shrink %>% filter(comparison == "age_17mo_vs_P30") %>% arrange(desc(log2FoldChange)) 

#res_17mo_vs_6mo <- res_all %>% filter(comparison == "age_17mo_vs_6mo")
res_17mo_vs_6mo_shrink <- res_17mo_vs_6mo_shrink %>% arrange(desc(log2FoldChange)) 

#plotting params
clip_nlog10p <- 10
clip_log2fc <- 2.5
min_log2FoldChange <- 0.1 #default 0

p1 <-volcanoPlot(res_6mo_shrink, min_log2FoldChange = 0.1,clip_nlog10p=10, clip_log2fc=2.5, n_label_genes = 15) + 
  ggrepel::geom_label_repel(aes(x = log2FoldChange,  y= nlog10p, label = label_genes),size = 2.5, force = 20)+
  ggtitle("6mo vs 1mo shrunken LFC")

p2<-volcanoPlot(res_6mo, min_log2FoldChange = 0.1,clip_nlog10p=10, clip_log2fc=2.5, n_label_genes = 15) + 
  ggrepel::geom_label_repel(aes(x = log2FoldChange,  y= nlog10p, label = label_genes),size = 2.5, force = 20)+
  ggtitle("6mo vs 1mo LFC")

p3<- volcanoPlot(res_17mo_shrink, min_log2FoldChange = 0.1,clip_nlog10p = 10, clip_log2fc = 2.5, n_label_genes = 15) + 
  ggrepel::geom_label_repel(aes(x = log2FoldChange,  y= nlog10p, label = label_genes),size = 2.5, force = 20)+
  ggtitle("17mo vs 1mo shrunken LFC")

p4<- volcanoPlot(res_17mo, min_log2FoldChange = 0.1,clip_nlog10p=10, clip_log2fc=2.5, n_label_genes = 15) + 
  ggrepel::geom_label_repel(aes(x = log2FoldChange,  y= nlog10p, label = label_genes),size = 2.5, force = 50)+
  ggtitle("17mo vs 1mo LFC") +
  xlim(c(-5,5))

p5 <-volcanoPlot(res_17mo_vs_6mo_shrink, min_log2FoldChange = 0.1,clip_nlog10p=10, clip_log2fc=2.5, n_label_genes = 15) + 
  ggrepel::geom_label_repel(aes(x = log2FoldChange,  y= nlog10p, label = label_genes),size = 2.5, force = 20)+
  ggtitle("17mo vs 6mo shrunken LFC")

p6<-volcanoPlot(res_17mo_vs_6mo, min_log2FoldChange = 0.1,clip_nlog10p=10, clip_log2fc=2.5, n_label_genes = 15) + 
  ggrepel::geom_label_repel(aes(x = log2FoldChange,  y= nlog10p, label = label_genes),size = 2.5, force = 20)+
  ggtitle("17mo vs 6mo LFC")

pl_grid <- ggpubr::ggarrange(plotlist = list(p1,p2,p3,p4,p5,p6),ncol = 2)
pdf(here("plots/figures/6mo_vs_1mo_volcano_plot_shrunken_lfc.pdf"))
p1
dev.off()



pdf(here("plots/figures/17mo_vs_1mo_volcano_plot_shrunken_lfc.pdf"))
p3
dev.off()

pdf(here("plots/figures/17mo_vs_6mo_volcano_plot_shrunken_lfc.pdf"))
p5
dev.off()


pdf(here("plots/figures/age_all_comps_volcano_plots.pdf"), width = 15)
pl_grid
dev.off()


