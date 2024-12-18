## README
TODO: Currently there is no special treatment (e.g. removal) of things like ribosomal genes or immunoglobin genes. Should they be removed before gene module clustering (yes), or before library normalization (this would be require rerunning all analyses)

### Preprocessing with kallisto:

### Preprocessing with STAR (for alignment, splicing analysis):

### Exploratory analysis:
- compare abundances and differential expression between pseudoaligment a) before or b) after collapsing libraries that were split for sequencing
- differences in diff exp estimates with modified models:
  - naive: expression ~ age
  - expression ~ factor(age, ordered = T) (will look at gradients in expression that are linear with age)
  - expression ~ age + RIN (two 17mo samples have low RIN and are outliers via PCA and by correlation with other samples)
  
Order of analysis scripts (after preprocessing):
1)  diffexp.Rmd

a1) dynamic_gene_modules.Rmd

b1) query_marker_genes.Rmd

c1) projections.Rmd

###Isoform splicing analysis: ./isoform_level_analysis
Parallel analysis pipelines with:
1) cuffdiff -> isoformSwitchAnalyzeR
- run cuffcompare.sh to generate .gtf annotation file with appropriate transcript identifiers (p_id, tss_id)
- run cufflink.sh to do isoform level DE, isoform splicing analysis (differential exons, TSS, etc...)
2) isoformSwitchAnalyzeR (with DEXSeq)


Gene signalling networks grabbed from https://baderlab.org/CellCellInteractions