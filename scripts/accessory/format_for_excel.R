dds <- readRDS(here("data/bulkrna_dds.rds"))

testnames <- c("age_6mo_vs_P30","age_17mo_vs_P30","age_17mo_vs_6mo")
gene_df <- rowData(dds)[,c("ensembl_gene_id_trimmed","external_gene_name")] %>% as.data.frame()

res_all <- map(testnames, function(comp){
  print(comp)
  res <- read.csv(here(paste0("results/diff_exp/",comp,".csv"))) %>%
    mutate(ensembl_gene_id_trimmed = ensemblgene_id) %>% left_join(., gene_df) 
  
  lfc_str <- paste0(comp, "_log2FoldChange")
  padj_str <- paste0(comp, "_padj")
  
  res[,lfc_str] <- res$log2FoldChange
  res[,padj_str] <- res$padj
  
  res <- res[,c("ensembl_gene_id_trimmed", "baseMean",
                  "external_gene_name", lfc_str, padj_str)]
  return(res)
}) %>% reduce(left_join) %>% 
  mutate(sum_sig = (age_17mo_vs_P30_padj < 0.05) + (age_6mo_vs_P30_padj < 0.05) + (age_17mo_vs_6mo_padj < 0.05))


res_all_sig <- res_all %>% filter(sum_sig >= 1)

write.csv(res_all, here("results/diff_exp/age_de_merged.csv"), row.names =F )
write.csv(res_all_sig, here("results/diff_exp/age_de_merged_significant.csv"), row.names = F)

#with shrinkage
res_all_shrink <- map(testnames, function(comp){
  print(comp)
  
  if(comp == "age_17mo_vs_6mo"){
    dds_2 <- dds
    dds_2$age <- factor(dds_2$age, levels = c("6mo","P30","17mo"))
    dds_2 <- DESeq(dds_2)
    res <- lfcShrink(dds_2, coef = comp)
  }else{
  res <- lfcShrink(dds, coef = comp)
   
  }
  res <- res %>% as.data.frame() %>% rownames_to_column(var = "ensembl_gene_id_trimmed") %>% left_join(., gene_df) 
  lfc_str <- paste0(comp, "_log2FoldChange")
  padj_str <- paste0(comp, "_padj")
  
  res[,lfc_str] <- res$log2FoldChange
  res[,padj_str] <- res$padj
  
  res <- res[,c("ensembl_gene_id_trimmed", "baseMean",
                "external_gene_name", lfc_str, padj_str)]
  return(res)
}) %>% reduce(left_join) %>% 
  mutate(sum_sig = (age_17mo_vs_P30_padj < 0.05) + (age_6mo_vs_P30_padj < 0.05) + (age_17mo_vs_6mo_padj < 0.05))

res_all_shrink_sig <- res_all_shrink %>% filter(sum_sig >= 1)

write.csv(res_all_shrink, here("results/diff_exp/age_de_lfc_shrink_merged.csv"), row.names =F )

#with gene descriptions
library(biomaRt)
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl") #uses mouse ensembl annotations

gene.data <- getBM(attributes=c('external_gene_name', 'ensembl_gene_id',"description","entrezgene_description"),
                   filters = 'ensembl_gene_id', values = rownames(dds) , mart = ensembl) %>% 
  mutate(ensembl_gene_id_trimmed = ensembl_gene_id) %>% group_by(ensembl_gene_id) %>% sample_n(size = 1)

res_all_shrink_sig <- res_all_shrink_sig %>% left_join(gene.data)
write.csv(res_all_shrink_sig, here("results/diff_exp/age_de_lfc_shrink_merged_significant.csv"), row.names = F)

