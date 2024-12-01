

age_6mo_vs_P30 <- read_csv("results/diff_exp/age_6mo_vs_P30.csv") %>% mutate(dir = sign(log2FoldChange))
de1 <- age_6mo_vs_P30 %>% filter(padj < 0.05) %>% dplyr::select(external_gene_name, ensemblgene_id, dir)
de1b <- age_6mo_vs_P30 %>% filter(padj < 0.05) %>%
  dplyr::select(external_gene_name, ensemblgene_id, dir) %>%
  filter(!str_detect(external_gene_name,"^Gm[0-9]+"))
age_17mo_vs_P30 <- read_csv("results/diff_exp/age_17mo_vs_P30.csv")%>% mutate(dir = sign(log2FoldChange))

de2 <- age_17mo_vs_P30 %>% filter(padj < 0.05) %>% dplyr::select(external_gene_name, ensemblgene_id,dir)
de2b <- age_17mo_vs_P30 %>% filter(padj < 0.05) %>%
  dplyr::select(external_gene_name, ensemblgene_id, dir) %>%
  filter(!str_detect(external_gene_name,"^Gm[0-9]+"))
age_17mo_vs_6mo <- read_csv("results/diff_exp/age_17mo_vs_6mo.csv")%>% mutate(dir = sign(log2FoldChange))
de3 <- age_17mo_vs_6mo %>% filter(padj < 0.05) %>% dplyr::select(external_gene_name, ensemblgene_id,dir)
de3b <- age_17mo_vs_6mo %>% filter(padj < 0.05) %>%
  dplyr::select(external_gene_name, ensemblgene_id, dir) %>%
  filter(!str_detect(external_gene_name,"^Gm[0-9]+"))


library(ggVennDiagram)
gene_list <- list("6mo_vs_P30" = de1$ensemblgene_id,
                  "17mo_vs_P30" = de2$ensemblgene_id,
                  "17mo_vs_6mo" = de3$ensemblgene_id)

ggVennDiagram(gene_list) + scale_x_continuous(expand = expansion(mult = .25))

gene_list2 <- list("6mo_vs_P30_up" = de1[de1$dir >0, "ensemblgene_id",drop = T],
                  "17mo_vs_P30_up" = de2[de2$dir >0, "ensemblgene_id",drop = T],
                  "17mo_vs_6mo_up" = de3[de3$dir >0, "ensemblgene_id",drop = T],
                  "6mo_vs_P30_down" = de1[de1$dir <0, "ensemblgene_id",drop = T],
                  "17mo_vs_P30_down" = de2[de2$dir <0, "ensemblgene_id",drop = T],
                  "17mo_vs_6mo_down" = de3[de3$dir <0, "ensemblgene_id",drop = T])

ggVennDiagram(gene_list2) + scale_x_continuous(expand = expansion(mult = .25))



gene_list1b <- list("6mo_vs_P30" = de1b$ensemblgene_id,
                  "17mo_vs_P30" = de2b$ensemblgene_id,
                  "17mo_vs_6mo" = de3b$ensemblgene_id)

ggVennDiagram(gene_list1b) + scale_x_continuous(expand = expansion(mult = .25))

unlist(gene_list) %>% unique %>% length()
unlist(gene_list1b) %>% unique %>% length()

gene_list2b <- list("6mo_vs_P30_up" = de1b[de1b$dir >0, "ensemblgene_id",drop = T],
                   "17mo_vs_P30_up" = de2b[de2b$dir >0, "ensemblgene_id",drop = T],
                   "17mo_vs_6mo_up" = de3b[de3b$dir >0, "ensemblgene_id",drop = T],
                   "6mo_vs_P30_down" = de1b[de1b$dir <0, "ensemblgene_id",drop = T],
                   "17mo_vs_P30_down" = de2b[de2b$dir <0, "ensemblgene_id",drop = T],
                   "17mo_vs_6mo_down" = de3b[de3b$dir <0, "ensemblgene_id",drop = T])

down <- gene_list2b[c("6mo_vs_P30_down","17mo_vs_P30_down","17mo_vs_6mo_down")] %>% unlist() %>% unique
up <- gene_list2b[c("6mo_vs_P30_up","17mo_vs_P30_up","17mo_vs_6mo_up")] %>% unlist() %>% unique

weird <- intersect(up,down)
gene_df[gene_df$ensembl_gene_id_trimmed %in% weird,]