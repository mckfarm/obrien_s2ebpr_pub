### alpha and beta diversity plots ###


library(tidyverse)
library(ggpubr)
library(ggcorrplot)
library(reshape2)
library(rstatix)
library(MetBrewer)
library(phyloseq)
library(qiime2R)
library(vegan)

source("analysis/plotting.R")


metadata <- read.delim("data/metadata.txt", sep = "\t")
metadata$date <- mdy(metadata$date)
metadata_ps <- metadata %>% column_to_rownames("sample")
metadata_ps <- sample_data(metadata_ps)

ps <- qza_to_phyloseq(features = "qiime_io/table_deblur.qza", 
                      taxonomy = "qiime_io/taxonomy_deblur.qza",
                      tree = "qiime_io/rooted_tree_deblur.qza")

ps <- merge_phyloseq(ps, metadata_ps)

ps_rel <- transform_sample_counts(ps, function(x) x*100/sum(x))

rm(ps)

### weighted unifrac - operation mode - all data ----
dist_all <- distance(ps_rel, method = "wunifrac")
dist_ord <- ordinate(ps_rel, "PCoA", "unifrac", weighted = TRUE)
dist_df <- as.data.frame(as.matrix(dist_all))
name_order <- colnames(dist_df)

metadata_ordered <- metadata_helper[match(name_order, metadata_helper$sample), ]
groups <- factor(metadata_ordered$op_mode)

disper_all <- betadisper(dist_all, groups, bias.adjust = TRUE)

plot(disper_all, hull = FALSE, ellipse = TRUE,
     xlab = "Axis 1 [27.6%]", ylab = "Axis2 [12.2%]", sub = NULL, main = "Weighted Unifrac")

dist_df_box <- as.data.frame(as.matrix(disper_all$distances)) %>%
  rownames_to_column(var = "sample") %>%
  left_join(metadata_ordered)

ggplot(dist_df_box, aes(x = op_mode, y = V1, color = op_mode)) +
  geom_boxplot() +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "Distance to centroid", x = "")





### weighted unifrac - carbon dosing - all data ------
ps_rel_sbr23 <- subset_samples(ps_rel, reactor != "SBR1")
sbr23_dist <- distance(ps_rel_sbr23, method = "wunifrac")
sbr23_ord <- ordinate(ps_rel_sbr23, "PCoA", "unifrac", weighted=TRUE)
dist_df <- as.data.frame(as.matrix(sbr23_dist))

name_order <- colnames(dist_df)

metadata_ordered <- metadata_helper[match(name_order, metadata_helper$sample), ] %>%
  filter(reactor != "SBR1")

groups <- factor(metadata_ordered$perf)

disper_all <- betadisper(sbr23_dist, groups, bias.adjust = TRUE)

plot(disper_all, hull = FALSE, ellipse = TRUE,
     xlab = "Axis 1 [31.2%]", ylab = "Axis2 [14.2%]", sub = NULL, main = "Weighted Unifrac")

dist_df_box <- as.data.frame(as.matrix(disper_all$distances)) %>%
  rownames_to_column(var = "sample") %>%
  left_join(metadata_ordered) %>%
  mutate(perf = factor(perf, levels = c("no C", "reduced C", "high C")))

ggplot(dist_df_box, aes(x = perf, y = V1, color = perf)) +
  geom_boxplot() +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "Distance to centroid", x = "")