Dictionary for aged bulkRNA seq analysis

# General processing and differential expression

## Differential expression 
### scripts
- diff_exp.Rmd: Initial scripts to read in abundance matrix from kallisto, create DESeq object, filter samples and genes, conduct and save differential expression results by age.
  - output: bulkrna_dds.rds
- dynamic_gene_modules.Rmd: Identify gene modules of DE genes that are dynamic with age. 
