# Microbiome analysis for O'Brien S2EBPR project

- Sequencing data processing run using a Snakemake workflow that automated sequential steps in QIIME2

- The Snakemake workflow consisted of the following steps:
    - Import into a QIIME object
    - Read trimming and merging
    - Deblur to identify ASVs
    - Taxonomy assignment against MiDAS v5.1
 
- The Snakemake workflow was run on the Quest High-Performance Computing Cluster

- Analysis was performed in R with some work performed on the Quest analytics node



This research was supported in part through the computational resources and staff contributions provided for the Quest high performance computing facility at Northwestern University which is jointly supported by the Office of the Provost, the Office for Research, and Northwestern University Information Technology.
