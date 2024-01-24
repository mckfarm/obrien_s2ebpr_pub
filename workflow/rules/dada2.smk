rule dada2:
    input:
        "qiime_io/reads_trim.qza"
    output:
        seqs = "qiime_io/rep_seqs_dada2.qza",
        table = "qiime_io/table_dada2.qza",
        stats = "qiime_io/stats_dada2.qza"
    resources:
        time = "40:00:00",
        mem = "60gb"
    threads:
        12
    shell:
        """
        module purge all
        module load qiime2/2023.2

        qiime dada2 denoise-paired --verbose --p-n-threads {threads} \
        --i-demultiplexed-seqs {input} \
        --p-trunc-len-f 0 --p-trunc-len-r 238 \
        --o-representative-sequences {output.seqs} \
        --o-table {output.table} \
        --o-denoising-stats {output.stats}

        """

rule viz_dada2:
    input:
        stats = "qiime_io/stats_dada2.qza",
        seqs = "qiime_io/rep_seqs_dada2.qza"
    output:
        stats = "qiime_io/stats_dada2.qzv",
        seqs = "qiime_io/rep_seqs_dada2.qzv"
    resources:
        time = "00:25:00",
        mem = "5gb"
    shell:
        """
        module purge all
        module load qiime2/2023.2

        qiime metadata tabulate \
        --m-input-file {input.stats} \
        --o-visualization {output.stats}

        qiime feature-table tabulate-seqs \
        --i-data {input.seqs} \
        --o-visualization {output.seqs}

        """



rule taxonomy_dada2:
    input:
        "qiime_io/rep_seqs_dada2.qza"
    output:
        "qiime_io/taxonomy_dada2.qza"
    params:
        classifier_path = "/projects/b1052/mckenna/resources/qiime_v5.1/midas_5.1_classifier.qza"
    resources:
        time = "01:00:00",
        mem = "50gb"
    threads:
        12
    shell:
        """
        module purge all
        module load qiime2/2023.2

        qiime feature-classifier classify-sklearn --p-n-jobs {threads} \
        --i-classifier {params.classifier_path} \
        --i-reads {input} \
        --o-classification {output}
        """

rule tree_dada2:
    input:
        "qiime_io/rep_seqs_dada2.qza"
    output:
        alignment = "qiime_io/aligned_rep_seqs_dada2.qza",
        masked_alignment = "qiime_io/masked_aligned_rep_seqs_dada2.qza",
        tree = "qiime_io/unrooted_tree_dada2.qza",
        rooted_tree = "qiime_io/rooted_tree_dada2.qza"
    resources:
        time = "12:00:00",
        mem = "50gb"
    threads:
        6
    shell:
        """
        module purge all
        module load qiime2/2023.2

        qiime phylogeny align-to-tree-mafft-fasttree --p-n-threads {threads} \
        --i-sequences {input} \
        --o-alignment {output.alignment} \
        --o-masked-alignment {output.masked_alignment} \
        --o-tree {output.tree} --o-rooted-tree {output.rooted_tree}
        """

rule barplot_dada2:
    input:
        table = "qiime_io/table_dada2.qza",
        taxonomy = "qiime_io/taxonomy_dada2.qza"
    output:
        "qiime_io/barplot_dada2.qzv"
    resources:
        time = "00:25:00",
        mem = "5gb"
    shell:
        """
        module purge all
        module load qiime2/2023.2

        qiime taxa barplot \
        --i-table {input.table} \
        --i-taxonomy {input.taxonomy} \
        --o-visualization {output}
        """