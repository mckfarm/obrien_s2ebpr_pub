### alpha and beta diversity plots ###

### packages ----
library(tidyverse)
library(phyloseq)
library(qiime2R)
library(vegan)
library(patchwork)
library(ggpubr)

source("scripts/plotting_cols_shapes.R")

### data read in -----
metadata <- read.delim("data/metadata.txt", sep = "\t")
metadata$date <- mdy(metadata$date)
metadata_ps <- metadata %>% column_to_rownames("sample")
metadata_ps <- sample_data(metadata_ps)

ps <- qza_to_phyloseq(features = "data/table_deblur.qza", 
                      taxonomy = "data/taxonomy_deblur.qza",
                      tree = "data/rooted_tree_deblur.qza")

ps <- merge_phyloseq(ps, metadata_ps)

ps_rel <- transform_sample_counts(ps, function(x) x*100/sum(x))


metadata_helper <- metadata %>%
  mutate(perf = case_when(
    reactor == "SBR1" ~ "no C", 
    (reactor != "SBR1" & date < ymd("2022-10-26")) ~ "no C",
    (reactor != "SBR1" & between(date, ymd("2022-10-25"), ymd("2022-12-2"))) ~ "high C",
    (reactor != "SBR1" & between(date, ymd("2022-12-2"), ymd("2023-05-12"))) ~ "reduced C",
    (reactor != "SBR1" & date > ymd("2023-05-12")) ~ "no C"
  )) %>%
  mutate(carb = case_when(
    reactor == "SBR1" ~ "C off", 
    (reactor != "SBR1" & date < ymd("2022-10-26")) ~ "C off",
    (reactor != "SBR1" & between(date, ymd("2022-10-25"), ymd("2023-05-12"))) ~ "C on",
    (reactor != "SBR1" & date > ymd("2023-05-12")) ~ "C off"
  )) %>%
  mutate(op_mode = case_when(
    reactor == "SBR1" ~ "EBPR",
    reactor == "SBR2" ~ "S2EBPR",
    (reactor == "SBR3" & date < ymd("2023-03-06")) ~ "EBPR",
    (reactor == "SBR3" & date >= ymd("2023-03-06")) ~ "S2EBPR"
  ))


### weighted unifrac - operation mode - all data ----
dist_all <- distance(ps_rel, method = "wunifrac")
dist_ord <- ordinate(ps_rel, "PCoA", "unifrac", weighted = TRUE)
dist_df <- as.data.frame(as.matrix(dist_all))
name_order <- colnames(dist_df)

metadata_ordered <- metadata_helper[match(name_order, metadata_helper$sample), ]
groups <- factor(metadata_ordered$reactor)

disper_all <- betadisper(dist_all, groups, bias.adjust = TRUE)

# nice plot
dist_result_df <- disper_all$vectors %>% 
  as.data.frame() %>% 
  select(c(PCoA1, PCoA2)) %>%
  rownames_to_column(var = "sample") %>%
  mutate(point = "points") %>%
  left_join(metadata_ordered)

centroid_df <- disper_all$centroids %>% 
  as.data.frame() %>%
  select(c(PCoA1, PCoA2)) %>% 
  rownames_to_column(var = "reactor") %>%
  mutate(point = "centroid") 

dist_result_df <- dist_result_df %>%
  bind_rows(centroid_df)

plt1 <- 
ggplot(dist_result_df, aes(x = PCoA1, y = PCoA2, color = reactor)) + 
  geom_point(data = . %>% filter(point == "points"), 
             aes(shape = op_mode), size = 2) +
  geom_point(data = . %>% filter(point == "centroid"), 
             shape = 21, size = 3) + 
  scale_color_reactor +
  scale_shape_op_mode + 
  labs(x = "Axis 1 [27.6%]", y = "Axis 2 [12.2%]", title = "Weighted Unifrac distance") + 
  theme_black_box +
  guides(color = guide_legend(order = 1), 
         shape = guide_legend(order = 2))
  

# centroid boxplot
dist_df_box <- as.data.frame(as.matrix(disper_all$distances)) %>%
  rownames_to_column(var = "sample") %>%
  left_join(metadata_ordered)

plt2 <- 
ggplot(dist_df_box, aes(x = reactor, y = V1, color = reactor)) +
  geom_boxplot(outlier.color = NA) +
  geom_point(position = position_jitter(width = 0.1)) + 
  scale_color_reactor + 
  ylim(0, 0.04) + 
  labs(y = "Distance to centroid", x = "Reactor") +
  theme_black_box + theme(legend.position = "none")



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

dist_result_df <- disper_all$vectors %>% 
  as.data.frame() %>% 
  select(c(PCoA1, PCoA2)) %>%
  rownames_to_column(var = "sample") %>%
  mutate(point = "points") %>%
  left_join(metadata_ordered)

centroid_df <- disper_all$centroids %>% 
  as.data.frame() %>%
  select(c(PCoA1, PCoA2)) %>% 
  rownames_to_column(var = "perf") %>%
  mutate(point = "centroid") 

dist_result_df <- dist_result_df %>%
  bind_rows(centroid_df) %>%
  mutate(perf = factor(perf, levels = c("no C", "reduced C", "high C")))

plt3 <- 
ggplot(dist_result_df, aes(x = PCoA1, y = PCoA2, color = perf)) + 
  geom_point(data = . %>% filter(point == "points"), 
             aes(shape = op_mode), size = 2) +
  geom_point(data = . %>% filter(point == "centroid"), 
            shape = 21, size = 3) + 
  scale_shape_op_mode + 
  scale_color_perf + 
  labs(x = "Axis 1 [31.2%]", y = "Axis 2 [14.2%]", title = "Weighted Unifrac distance - SBR2 and SBR3 only") + 
  theme_black_box + 
  guides(color = guide_legend(order = 1), 
         shape = guide_legend(order = 2))

dist_df_box <- as.data.frame(as.matrix(disper_all$distances)) %>%
  rownames_to_column(var = "sample") %>%
  left_join(metadata_ordered) %>%
  mutate(perf = factor(perf, levels = c("no C", "reduced C", "high C")))

dist_df_box %>%
  group_by(perf) %>%
  summarise(median(V1))

plt4 <- 
ggplot(dist_df_box, aes(x = perf, y = V1, color = perf)) +
  geom_boxplot(outlier.color = NA) +
  geom_point(position = position_jitter(width = 0.1)) + 
  scale_color_perf + 
  ylim(0, 0.04) + 
  labs(y = "Distance to centroid", x = "Carbon") +
  theme_black_box + theme(legend.position = "none")

# combo
plt1 + plt3 + plt2 + plt4 + plot_layout(ncol = 2)
ggsave("results/beta_diversity_bigplot.png", units = "in", width = 12, height = 7, dpi = 300)


dist_result_df <- dist_result_df %>%
  mutate_at(vars(date), funs(year, month)) %>%
  mutate(year_month = paste0(year, "-", month))

# by date 
ggplot(dist_result_df, aes(x = PCoA1, y = PCoA2, color = year_month)) + 
  geom_point(data = . %>% filter(point == "points"), 
             aes(shape = reactor), size = 2) +
  geom_point(data = . %>% filter(point == "centroid"), 
             shape = 21, size = 3) + 
  labs(x = "Axis 1 [27.6%]", y = "Axis 2 [12.2%]", title = "Weighted Unifrac distance") + 
  theme_black_box



### alpha diversity ----
source("scripts/plotting_dates.R")

tot_counts <- as.data.frame(as.matrix(sample_sums(ps))) %>%
  rename(tot_count = V1) %>%
  rownames_to_column(var = "sample")

alpha_div <- estimate_richness(ps) %>% 
  rownames_to_column(var = "sample") %>%
  mutate(sample = str_replace(sample, "\\.", "-")) %>%
  left_join(tot_counts) %>%
  pivot_longer(-sample) %>%
  left_join(metadata_helper) 

alpha_div %>% 
  filter(name %in% c("InvSimpson", "Shannon")) %>% 
  ggplot(data = ., aes(x = date, y = value, color = reactor)) +
  facet_wrap(~name, scales = "free") + 
  geom_point(size = 2, aes(shape = carb)) +
  geom_line(linewidth = 0.3) +
  scale_color_reactor +
  scale_shape_carb + 
  x_axis_date + 
  theme_black_box_facet + 
  labs(x = "Date", y = "Alpha diversity measure")
ggsave("results/alpha_diversity.png", units = "in", width = 6, height = 3, dpi = 300)

