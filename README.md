# Microbiome analysis for O'Brien S2EBPR project



Brief overview of workflow
- Raw sequencing data was processed using [QIIME2](https://qiime2.org) in a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow.
    - Raw reads are stored at NCBI with accession [PRJNA1068094](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1068094)
- The Snakemake workflow consisted of the following steps:
    - Import raw reads into a QIIME object
    - Read trimming and merging
    - Identify ASVs with [Deblur](https://doi.org/10.1128%2FmSystems.00191-16)
    - Taxonomy assignment using [MiDAS](https://midasfieldguide.org/guide) v5.1
- The Snakemake workflow was run on the Quest High-Performance Computing Cluster
- Data analysis was performed in R with most work performed on the Quest analytics node

This research was supported in part through the computational resources and staff contributions provided for the Quest high performance computing facility at Northwestern University which is jointly supported by the Office of the Provost, the Office for Research, and Northwestern University Information Technology.
