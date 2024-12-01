library(monocle3)
library(tidyverse)
library(here)
library(DESeq2)

lmmp <- readRDS("/data/users/jared/ENS/Timecourse_ENS/TC_LMMP.rds")


pData(lmmp) %>% head

lmmp_clean <- lmmp[,!(pData(lmmp)$cell_type == "unknown" | pData(lmmp)$cell_type == "Unknown")]

pData(lmmp_clean) <- pData(lmmp_clean) %>% as.data.frame() %>%
  mutate(cell_type_coarse = case_when(
    cell_type == "smooth muscle or junk" ~ "smooth muscle",
    cell_type == "fibroblasts" ~ "fibroblast",
    cell_type == "fibroblast 1" ~ "fibroblast",
    cell_type == "fibroblast 2" ~ "fibroblast",
    cell_type == "neuroglia/neuroendocrine" ~ "NENs",
    .default = cell_type
  )) %>% DataFrame()

coarse_markers <- top_markers(lmmp_clean, group_cells_by = "cell_type_coarse", genes_to_test_per_group = 500, cores = 16)
write.csv(coarse_markers, here("data/timecourse_sc_markers/timecourse_sc_coarse_markers.csv"))

ct_categories <- data.frame("cell_group" = c("adipocytes","B cells","endothelial","enterocytes",
                                             "fibroblast","ICC?","macrophage","MENs",
                                             "mucosal epithelium","NENs","Neuroglia","putative glia",
                                             "RBC","S/G2M","smooth muscle","T cells")
                            , "meta_cat" = c("other","immune","other","other",
                                             "fibroblast/SMC/ICC","fibroblast/SMC/ICC","immune","MENs-ENS",
                                             "other","crest-ENS","crest-ENS","crest-ENS",
                                             "other","other","fibroblast/SMC/ICC","immune"))

#prepare coarse markers by finding cell type with highest specificity measure for each gene
ggplot(coarse_markers, aes(x = specificity, y = pseudo_R2)) + geom_point(alpha = 0.1) + geom_smooth() + facet_wrap(~cell_group)
ggplot(coarse_markers, aes(x = specificity, y = marker_score)) + geom_point(alpha = 0.1) + geom_smooth() 

ranked_markers <- coarse_markers %>% group_by(gene_id) %>%
  mutate(rank_specificity = dense_rank(dplyr::desc(specificity)),
         gene_ct = case_when(rank_specificity == 1 ~ cell_group))%>%
  full_join(ct_categories)

#possibly also filter marker q < 0.05
markers_filt <- ranked_markers %>%
  filter(marker_test_q_value < 0.05) %>%
  select(gene_id, gene_short_name, specificity, gene_ct, meta_cat) %>% 
  filter(!is.na(gene_ct)) 



dds <- readRDS(here("data/bulkrna_dds.rds"))
comp <- "age_17mo_vs_P30"
  
#probably also filter log2fc padj < 0.05
log2fc <- lfcShrink(dds, coef = comp, type="apeglm")

#these genes are filtered on being 
#1) in the top 500 specific genes for at least one cell type
#2) significant marker for at least one cell type by top_markers()
fc_df <- log2fc %>% as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  left_join(markers_filt) %>%
  filter(!is.na(specificity))

fc_df$gene_ct %>% table #number of markers per cell type ranges [34, 442]
fc_df$meta_cat %>% table #number of markers per meta category = [302, 1311]

q1 <- ggplot(fc_df, aes(x = specificity, y = log2FoldChange)) + geom_point(aes(color = meta_cat, alpha = specificity)) + 
  ggtitle(paste0("Max specificity in scRNA-seq vs log2FoldChange: ", comp), 
          subtitle = paste(nrow(fc_df), "significant marker genes"))+ 
  geom_hline(yintercept = 0, linetype = "dashed")

q2 <- ggplot(filter(fc_df, abs(log2FoldChange) > 0.1), aes(x = specificity, y = log2FoldChange, alpha = specificity)) +
  geom_point(aes(color = meta_cat)) +
  ggtitle(paste0("Max specificity in scRNA-seq vs log2FoldChange: ", comp), 
            subtitle = paste(nrow(filter(fc_df, abs(log2FoldChange) > 0.1)), "significant marker genes and log2fc > 0.1")) + 
  geom_hline(yintercept = 0, linetype = "dashed")

q3 <- ggplot(filter(fc_df, padj < 0.05), aes(x = specificity, y = log2FoldChange, alpha = specificity)) +
  geom_point(aes(color = meta_cat)) +
  ggtitle(paste0("Max specificity in scRNA-seq vs log2FoldChange: ", comp), 
          subtitle = paste(nrow(filter(fc_df, padj < 0.05)), "significant marker genes and DE padj < 0.05")) + 
  geom_hline(yintercept = 0, linetype = "dashed")

p1 <- ggplot(filter(fc_df, (padj < 0.05) & meta_cat == "immune"), aes(x = specificity, y = log2FoldChange)) +
  geom_point(aes(color = gene_ct)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_blank(aes(y = -log2FoldChange)) + #keep scales symmetrical about y
  ggtitle(paste0("Immune markers: ",nrow(filter(fc_df, (padj < 0.05) & meta_cat == "immune")), " DE genes"))

p2 <- ggplot(filter(fc_df, (padj < 0.05) & meta_cat == "crest-ENS"), aes(x = specificity, y = log2FoldChange)) +
  geom_point(aes(color = gene_ct)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_blank(aes(y = -log2FoldChange)) +#keep scales symmetrical about y
  ggtitle(paste0("Crest-ENS markers: ",nrow(filter(fc_df, (padj < 0.05) & meta_cat == "crest-ENS")), " DE genes"))

p3 <- ggplot(filter(fc_df, (padj < 0.05) & meta_cat == "fibroblast/SMC/ICC"), aes(x = specificity, y = log2FoldChange)) +
  geom_point(aes(color = gene_ct)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_blank(aes(y = -log2FoldChange)) + #keep scales symmetrical about y
  ggtitle(paste0("fibroblast/SMC/ICC markers: ",nrow(filter(fc_df, (padj < 0.05) & meta_cat == "fibroblast/SMC/ICC")), " DE genes"))

p <- ggpubr::ggarrange(p1, p2, p3, ncol = 3)

pdf(here("plots/single_cell_integration/cell_type_specificity_vs_logfc_metacat_split_padj.pdf"), width = 18, height = 7.5)
p
dev.off()

pdf(here("plots/single_cell_integration/cell_type_specificity_vs_logfc_metacat.pdf"), width = 18, height = 7.5)
ggpubr::ggarrange(q1, q2, q3, ncol = 3)
dev.off()

specificity_threshold <- c(0.2, 0.5, 0.8)

de_tbl <- map_df(specificity_threshold, function(sp){
  fc_df %>% filter(specificity >= sp, padj < 0.05) %>%
    group_by(gene_ct) %>%
    mutate(de_sign = sign(log2FoldChange)) %>% 
    summarise(specificity_threshold = sp, de_genes = n(), de_up = sum(de_sign == 1),
              de_down = sum(de_sign == -1),  direction = sum(de_sign)) %>%
    arrange(desc(direction))
})
write.csv(de_tbl, here("results/diff_exp/scRNA_integration/cell_type_markers_de_results.csv"))

de_tbl_meta <- map_df(specificity_threshold, function(sp){
  fc_df %>% filter(specificity >= sp, padj < 0.05) %>%
  group_by(meta_cat) %>% mutate(de_sign = sign(log2FoldChange)) %>% 
  summarise(specificity_threshold = sp, n_de_genes = n(), de_up = sum(de_sign == 1),
            de_down = sum(de_sign == -1),  direction = sum(de_sign)) %>%
  arrange(desc(direction))
})
write.csv(de_tbl_meta, here("results/diff_exp/scRNA_integration/cell_type_markers_meta_de_results.csv"))



#for ppt figures
genes_mark <- c("Fbn1","Mmp2","Col1a1","Col6a3","Nfasc","Ric3","Agrn")
fc_df <- fc_df %>% mutate(meta_category = factor(meta_cat, levels = c("crest-ENS","MENs-ENS","fibroblast/SMC/ICC","immune","other")),
                          marker_celltype = gene_ct,
                          gene_mark = ifelse(gene_short_name %in% genes_mark, gene_short_name, NA))
fc_df[fc_df[,"marker_celltype"] == "ICC?","marker_celltype"] <- "ICC"

q3 <- ggplot(filter(fc_df, padj < 0.05), aes(x = specificity, y = log2FoldChange, alpha = specificity)) +
  geom_point(aes(color = meta_category)) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_color_manual(values = c("deepskyblue3", "chartreuse3","chocolate","brown2","gray60")) + 
  theme(text = element_text(size = 20))

p1 <- ggplot(filter(fc_df, (padj < 0.05) & meta_cat == "immune"), aes(x = specificity, y = log2FoldChange)) +
  geom_point(aes(color =  marker_celltype)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_blank(aes(y = -log2FoldChange)) + #keep scales symmetrical about y
  ggtitle("Immune markers",subtitle= paste0(nrow(filter(fc_df, (padj < 0.05) & meta_cat == "immune")), " DE genes")) + 
  theme(text = element_text(size = 20), legend.position = "bottom", legend.title = element_blank())

p2 <- ggplot(filter(fc_df, (padj < 0.05) & meta_cat == "crest-ENS"), aes(x = specificity, y = log2FoldChange)) +
  geom_point(aes(color =  marker_celltype)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_blank(aes(y = -log2FoldChange)) +#keep scales symmetrical about y
  ggtitle("Crest-ENS markers", subtitle = paste0(nrow(filter(fc_df, (padj < 0.05) & meta_cat == "crest-ENS")), " DE genes")) + 
  theme(text = element_text(size = 20),legend.position = "bottom", legend.title = element_blank(), axis.title.y = element_blank()) + 
  ggrepel::geom_text_repel(aes(label = gene_mark), min.segment.length = 0, force = 50, size =6)

p3 <- ggplot(filter(fc_df, (padj < 0.05) & meta_cat == "fibroblast/SMC/ICC"), aes(x = specificity, y = log2FoldChange)) +
  geom_point(aes(color =  marker_celltype)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_blank(aes(y = -log2FoldChange)) + #keep scales symmetrical about y
  ggtitle("Fibroblast/SMC/ICC markers",paste0(nrow(filter(fc_df, (padj < 0.05) & meta_cat == "fibroblast/SMC/ICC")), " DE genes")) + 
  theme(text = element_text(size = 20),  legend.position = "bottom", legend.title = element_blank(), axis.title.y = element_blank()) +
  ggrepel::geom_text_repel(aes(label = gene_mark), size = 6, min.segment.length = 0, force = 50)

pdf(here("plots/figures/sc_integration_de_meta.pdf"), width = 6.5, height = 5)
q3
dev.off()

pdf(here("plots/figures/cell_type_specificity_vs_logfc_metacat.pdf"), width = 18, height = 7.5)
ggpubr::ggarrange(p1, p2, p3, ncol = 3)
dev.off()
