# Microbiome analysis for O'Brien S2EBPR project

- Sequencing data processing run using a snakemake workflow 
- Workflow consisted of the following steps:
    - Import into a QIIME object
    - Read trimming and merging
    - Deblur to identify ASVs
    - Taxonomy assignment against MiDAS v5.1

- Analysis was performed in R with the following packages:
    - phyloseq
    - qiime2R
