
#R4.2.0
.libPaths("/home/jared/R/x86_64-pc-linux-gnu-library/4.2")

#library(tidyverse)
library(here)
library(IsoformSwitchAnalyzeR)
library(dplyr)

switchRes <- readRDS(here("isoform_level_analysis/isoformSwitchAnalyzer/DEXSeq/DEXSeq_res_all.rds"))

pl <- switchPlot(switchRes, gene = "Bdnf", plotTopology = T, condition1 = "age_P30", condition2 = "age_17mo")

pdf(here("isoform_level_analysis/isoformSwitchAnalyzer/satuRn/figures/test.pdf"))
pl
dev.off()


f1 <- switchRes$isoformFeatures %>% filter(gene_id == "Bdnf") %>%
  select(isoform_id, gene_id, condition_1, gene_value_1, gene_stderr_1, iso_value_1, iso_stderr_1, IF1)
f2 <- switchRes$isoformFeatures %>% filter(gene_id == "Bdnf") %>%
  select(isoform_id, gene_id, condition_2, gene_value_2, gene_stderr_2, iso_value_2, iso_stderr_2, IF2)
colnames(f1) <- c("isoform_id", "gene_id", "condition", "gene_value", "gene_stderr", "iso_value", "iso_stderr", "IF")
colnames(f2) <- c("isoform_id", "gene_id", "condition", "gene_value", "gene_stderr", "iso_value", "iso_stderr", "IF")


features <- rbind(f1,f2) %>% distinct() %>%
  mutate(condition =case_when(
    condition == "age_P30" ~ "1mo",
    condition == "age_6mo" ~ "6mo",
    condition == "age_17mo" ~ "17mo")) %>%
  mutate(age = factor(condition, levels = c("1mo","6mo","17mo")))

#gene plots
gene_data <- features %>% transmute(gene_id, age, gene_value, gene_stderr) %>% distinct()

pl1 <- ggplot(gene_data, aes(x = age, y = gene_value)) + geom_col(aes(fill = age), color = "black") +
  geom_errorbar(aes(ymin = gene_value - gene_stderr, ymax = gene_value + gene_stderr), width = 0.1)+
  theme_light()+
  xlab("") +
  ylab("Gene Expression (TPM)") + 
  scale_fill_manual(values = c("#006D77","#83C5BE","#EDF6F9")) +
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.position = "none")  +
  theme(plot.margin = margin(5.5, 12, 53, 7, unit = "pt"))


#isoform plots
meta_bdnf <- read.csv(here("isoform_level_analysis/data/bdnf_isoform_description.csv"))
isoform_data <- features %>% select(-c(gene_value, gene_stderr)) %>%
  tidyr::separate(isoform_id, into = c("transcript_id","version"), sep = "\\.", remove = F) %>% 
  left_join(meta_bdnf) %>% 
  mutate(isoform_label = paste0(isoform, " (", exon_id, ")")) %>%
  mutate(isoform_label = factor(isoform_label,
                                levels = c("Bdnf-205 (exon6)", "Bdnf-211 (exon5)",
                                           "Bdnf-209 (exon4)", "Bdnf-206 (exon3)","Bdnf-210 (exon2)")))

isoform_data[isoform_data[,"iso_value"] == 0,c("iso_stderr", "iso_value")] <- NA

pl2 <- ggplot(isoform_data, aes(x = isoform_label, y = iso_value, fill = age)) +
  geom_bar(position = position_dodge(width = .7), stat = "identity", width = .6, color ="black") +
  geom_errorbar(aes(ymin = iso_value - iso_stderr, ymax = iso_value + iso_stderr),
                position = position_dodge(width = .7), width = 0.4) +
  theme_light()+
  xlab("") +
  ylab("Isoform Expression (TPM)") + 
  scale_fill_manual(values = c("#006D77","#83C5BE","#EDF6F9")) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0.1))+ 
  theme(legend.direction = "horizontal", legend.position = c(0.6,0.9), legend.title = element_blank()) +
  theme(plot.margin = margin(5.5, 20, 5.5, 7, unit = "pt"))


pl3 <- ggplot(isoform_data, aes(x = isoform_label, y = IF, fill = age)) +
  geom_bar(position = position_dodge(width = .7), stat = "identity", width = 0.6, color = "black") +
  theme_light()+
  xlab("") +
  ylab("Isoform Fraction (IF)") + 
  scale_fill_manual(values = c("#006D77","#83C5BE","#EDF6F9")) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0.1))+ 
  theme(legend.direction = "horizontal", legend.position = c(0.6,0.9), legend.title = element_blank()) +
  theme(plot.margin = margin(5.5, 20, 5.5, 7, unit = "pt"))

pdf(here("isoform_level_analysis/isoformSwitchAnalyzer/satuRn/figures/supp_bdnf_gene_expression.pdf"), height = 4, width = 2)
pl1 
dev.off()

pdf(here("isoform_level_analysis/isoformSwitchAnalyzer/satuRn/figures/supp_bdnf_isoform_expression.pdf"), height = 4, width = 4)
pl2 
dev.off()

pdf(here("isoform_level_analysis/isoformSwitchAnalyzer/satuRn/figures/supp_bdnf_isoform_fraction.pdf"), height = 4, width = 4)
pl3 
dev.off()

pl4 <- ggpubr::ggarrange(pl1, pl2, pl3, ncol = 3, widths = c(1, 2, 2))
pdf(here("isoform_level_analysis/isoformSwitchAnalyzer/satuRn/figures/supp_merge.pdf"), height = 4, width = 10.3)
pl4 
dev.off()   




bdnf_isos <- isoform_data$isoform_id %>% unique
bdnf_exp <- switchRes$isoformRepExpression %>% filter(isoform_id %in% bdnf_isos)

gene_exp <- bdnf_exp %>% tibble::column_to_rownames(var = "isoform_id") %>% as.matrix() %>% colSums()
mean(test[1:4])
sd(test[1:4])/2

iso_exp <- bdnf_exp %>% tibble::column_to_rownames(var = "isoform_id") %>% as.matrix() 
