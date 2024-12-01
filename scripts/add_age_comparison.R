library(here)
library(DESeq2)

dds <- readRDS(here("data/bulkrna_dds.rds"))
gene_df <- rowData(dds)[,c("ensembl_gene_id_trimmed","external_gene_name")] %>% as.data.frame() %>%
  transmute(ensemblgene_id = ensembl_gene_id_trimmed, external_gene_name = external_gene_name)

res <- results(dds, contrast=c("age", "17mo", "6mo"), alpha=0.1, pAdjustMethod="BH", parallel=TRUE, tidy = F) %>% as.data.frame() %>%
  rownames_to_column(var = "ensemblgene_id") %>%
  left_join(., gene_df, by= "ensemblgene_id") %>%
  relocate(external_gene_name) %>%
  arrange(desc(log2FoldChange))

write.csv(res, here("results/diff_exp/age_17mo_vs_6mo.csv"))


#testing accuracy of de
res_all <- map(testnames, function(comp){
  print(comp)
  res <- read.csv(here(paste0("results/diff_exp/",comp,".csv"))) %>%
    mutate(ensembl_gene_id = ensemblgene_id) %>% left_join(., gene_df) 
  
  lfc_str <- paste0(comp, "_log2FoldChange")
  padj_str <- paste0(comp, "_padj")
  
  res[,lfc_str] <- res$log2FoldChange
  res[,padj_str] <- res$padj
  
  res <- res[,c("ensembl_gene_id", "baseMean",
                "external_gene_name", lfc_str, padj_str)]
  return(res)
}) %>% reduce(left_join) %>%
  mutate(sum_sig = (age_17mo_vs_P30_padj < 0.05) + (age_6mo_vs_P30_padj < 0.05) + (age_6mo_vs_P30_padj < 0.05))

ngenes <- nrow(dds)
gene_df <- rowData(dds) %>% as.data.frame() %>%
  transmute(ensembl_gene_id = ensembl_gene_id_trimmed, external_gene_name =external_gene_name)
sample_ind <- sample(ngenes, 50)
ntd <- normTransform(dds[sample_ind,]) %>% assay()


#Prep count matrix by sample and merged by age

ntd_long <- ntd %>% as.data.frame() %>% rownames_to_column(var = "ensembl_gene_id") %>%
  pivot_longer(!ensembl_gene_id, names_to = "sample_name", values_to = "counts") %>%
  separate(., "sample_name", into = c("age","sex","rep"), sep = "-", remove = F) %>%
  mutate(age = factor(age, levels = c("P30","6mo","17mo")),
         age_numeric = case_when(
           age == "P30" ~ 30,
           age == "6mo" ~ 180, 
           age == "17mo" ~ 510)) %>% 
  left_join(gene_df)


ntd_avg <- ntd_long %>% group_by(ensembl_gene_id) %>%
  summarise(avg_expression = mean(counts), n = n()) %>%
  left_join(gene_df) %>%
  left_join(res_all)



# pdf(here("plots/diff_exp/collagen_expression.pdf"), width = 16, height = 8)
pl <- ggplot(ntd_long,aes(x = age_numeric, y = counts))+ geom_point() +
  facet_wrap(~external_gene_name, scales = "free_y") + 
  geom_text(data = ntd_avg  , aes(x = 100, y = avg_expression, label = round(age_6mo_vs_P30_log2FoldChange, 2))) + 
  geom_text(data = ntd_avg  , aes(x = 350, y = avg_expression, label = round(age_17mo_vs_6mo_log2FoldChange, 2))) +
  geom_segment(data = ntd_avg, aes(y= avg_expression -.25, yend = avg_expression -.25, x = 0, xend = 170)) +
  geom_segment(data = ntd_avg, aes(y= avg_expression -.25, yend = avg_expression -.25, x = 190, xend = 510)) +
  geom_text(data = ntd_avg  , aes(x = 250, y = avg_expression-1, label = round(age_17mo_vs_P30_log2FoldChange, 2))) + 
  geom_segment(data = ntd_avg, aes(y= avg_expression -1.1, yend = avg_expression -1.1, x = 0, xend = 510)) + 
  ylab("expression") + 
  ggtitle("expression for 50 random genes", subtitle = "numbers are log2fc between conditions")
  


#plt
library(GGally)

lfc_age <- res_all %>% dplyr::select("age_6mo_vs_P30_log2FoldChange", "age_17mo_vs_P30_log2FoldChange", "age_17mo_vs_6mo_log2FoldChange")

# Customize your scatterplots as you wish here:
lowerfun <- function(data, mapping) {
  ggplot(data = data, mapping = mapping)+ 
    ggpointdensity::geom_pointdensity(alpha = .25) + 
    scale_color_viridis_c(trans = scales::pseudo_log_trans(sigma = 0.001)) +
    geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dashed")
}

lowerfun_squish <- function(data, mapping) {
  ggplot(data = data, mapping = mapping)+ 
    ggpointdensity::geom_pointdensity(alpha = .25) + 
    scale_color_viridis_c(trans = scales::pseudo_log_trans(sigma = 0.001)) +
    geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dashed") + 
    scale_x_continuous(limits = c(-4,4), oob = scales::oob_squish) + 
    scale_y_continuous(limits = c(-4,4), oob = scales::oob_squish)
  
}

# Plot the scatterplot matrix
pl2 <- ggpairs(lfc_age, lower = list(continuous = wrap(lowerfun))) + ggtitle("DESEQ2 log2fc for all genes")

lfc_age_sig <- res_all %>% filter(sum_sig >= 1) %>% dplyr::select("age_6mo_vs_P30_log2FoldChange", "age_17mo_vs_P30_log2FoldChange", "age_17mo_vs_6mo_log2FoldChange")
pl3 <- ggpairs(lfc_age_sig, lower = list(continuous = wrap(lowerfun))) + ggtitle("DESEQ2 log2fc for significant (any conditions) genes")

age_de_lfc_shrink_merged <- read_csv("results/diff_exp/age_de_lfc_shrink_merged.csv")

lfc_age_shrink <- age_de_lfc_shrink_merged %>% dplyr::select("age_6mo_vs_P30_log2FoldChange", "age_17mo_vs_P30_log2FoldChange", "age_17mo_vs_6mo_log2FoldChange")
pl4 <- ggpairs(lfc_age_shrink, lower = list(continuous = wrap(lowerfun_squish))) + ggtitle("DESEQ2 shrunken log2fc for all genes")

lfc_age_shrink_sig <- age_de_lfc_shrink_merged %>% filter(sum_sig >0) %>% dplyr::select("age_6mo_vs_P30_log2FoldChange", "age_17mo_vs_P30_log2FoldChange", "age_17mo_vs_6mo_log2FoldChange")
pl5 <- ggpairs(lfc_age_shrink_sig, lower = list(continuous = wrap(lowerfun_squish)))+ ggtitle("DESEQ2 shrunken log2fc for significant (any conditions) genes")

pdf(here("plots/exploratory/gene_expression_lfc.pdf"), width = 15, height = 15)
pl
pl2
pl3
pl4
pl5
dev.off()
