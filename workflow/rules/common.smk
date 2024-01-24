def get_rules(self):

    all_rules = []

    if config["dada2"]:
        all_rules.append("qiime_io/rep_seqs_dada2.qza")
        all_rules.append("qiime_io/table_dada2.qza")
        all_rules.append("qiime_io/stats_dada2.qza")
        all_rules.append("qiime_io/taxonomy_dada2.qza")
        all_rules.append("qiime_io/unrooted_tree_dada2.qza")
        all_rules.append("qiime_io/rooted_tree_dada2.qza")

        all_rules.append("qiime_io/stats_dada2.qzv")
        all_rules.append("qiime_io/rep_seqs_dada2.qzv")
        all_rules.append("qiime_io/barplot_dada2.qzv")

    if config["deblur"]:
        all_rules.append("qiime_io/rep_seqs_deblur.qza")
        all_rules.append("qiime_io/table_deblur.qza")
        all_rules.append("qiime_io/stats_deblur.qza")
        all_rules.append("qiime_io/taxonomy_deblur.qza")
        all_rules.append("qiime_io/unrooted_tree_deblur.qza")
        all_rules.append("qiime_io/rooted_tree_deblur.qza")

        all_rules.append("qiime_io/stats_deblur.qzv")
        all_rules.append("qiime_io/rep_seqs_deblur.qzv")
        all_rules.append("qiime_io/barplot_deblur.qzv")

    return all_rules