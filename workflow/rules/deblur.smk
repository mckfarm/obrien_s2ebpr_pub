rule deblur:
    input:
        "qiime_io/reads_trim_merge.qza"
    output:
        seqs = "qiime_io/rep_seqs_deblur.qza",
        table = "qiime_io/table_deblur.qza",
        stats = "qiime_io/stats_deblur.qza"
    resources:
        time = "24:00:00",
        mem = "10gb"
    threads:
        16
    shell:
        """
        module purge all
        module load qiime2/2023.2

        qiime deblur denoise-16S --verbose --p-jobs-to-start {threads} \
        --i-demultiplexed-seqs {input} \
        --o-table {output.table} \
        --o-representative-sequences {output.seqs} \
        --o-stats {output.stats} \
        --p-trim-length 368 \
        --p-sample-stats
        
        """

rule viz_deblur:
    input:
        stats = "qiime_io/stats_deblur.qza",
        seqs = "qiime_io/rep_seqs_deblur.qza"
    output:
        stats = "qiime_io/stats_deblur.qzv",
        seqs = "qiime_io/rep_seqs_deblur.qzv"
    resources:
        time = "00:25:00",
        mem = "5gb"
    shell:
        """
        module purge all
        module load qiime2/2023.2

        qiime deblur visualize-stats \
        --i-deblur-stats {input.stats} \
        --o-visualization {output.stats}

        qiime feature-table tabulate-seqs \
        --i-data {input.seqs} \
        --o-visualization {output.seqs}

        """

rule taxonomy_deblur:
    input:
        "qiime_io/rep_seqs_deblur.qza"
    output:
        "qiime_io/taxonomy_deblur.qza"
    params:
        classifier_path = "/projects/b1052/mckenna/resources/qiime_v5.1/midas_5.1_classifier.qza"
    resources:
        time = "24:00:00",
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

rule tree_deblur:
    input:
        "qiime_io/rep_seqs_deblur.qza"
    output:
        alignment = "qiime_io/aligned_rep_seqs_deblur.qza",
        masked_alignment = "qiime_io/masked_aligned_rep_seqs_deblur.qza",
        tree = "qiime_io/unrooted_tree_deblur.qza",
        rooted_tree = "qiime_io/rooted_tree_deblur.qza"
    resources:
        time = "12:00:00",
        mem = "50gb"
    threads:
        12
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

rule barplot_deblur:
    input:
        table = "qiime_io/table_deblur.qza",
        taxonomy = "qiime_io/taxonomy_deblur.qza"
    output:
        "qiime_io/barplot_deblur.qzv"
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