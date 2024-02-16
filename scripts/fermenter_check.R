# check fermenter diversity # 

# packages
library(tidyverse)
library(ggpubr)
library(vegan)
library(patchwork)
library(extrafont)
library(pheatmap)

source("scripts/plotting_cols_shapes.R")
source("scripts/plotting_dates.R")


# custom stat compare means
fixed_call <- quote(ggsignif::geom_signif(comparisons = comparisons, 
                                          y_position = label.y, 
                                          test = method, test.args = method.args, 
                                          step_increase = step.increase, 
                                          size = bracket.size, textsize = 4, color = color, 
                                          map_signif_level = map_signif_level, 
                                          tip_length = tip.length, data = data,
                                          vjust = vjust, ...))

scm <- as.list(ggpubr::stat_compare_means)
scm[[26]][[2]][[3]][[13]] <- fixed_call
stat_compare_means <- function(...) {
  ggp <- getNamespace('ggpubr')
  .method_info <- ggp$.method_info
  .add_item <- ggp$.add_item
  .is_p.signif_in_mapping <- ggp$.is_p.signif_in_mapping
  .is_empty <- ggp$.is_empty
  do.call(as.function(scm), list(...), quote = FALSE)
}

# data read in and prep
reactors <- c("SBR1", "SBR2", "SBR3")

my_comparisons <- list(c("SBR1", "SBR2"),
                       c("SBR2", "SBR3"),
                       c("SBR1", "SBR3"))

rel_genus <- readRDS("data/rel_genus.RDS") %>%
  mutate(genus_format = str_replace(Genus, "Ca_", "Ca "))

metadata <- read_delim("data/metadata.txt") %>%
  mutate(date = mdy(date)) %>%
  mutate(perf = case_when(
    reactor == "SBR1" ~ "C off", 
    (reactor != "SBR1" & date < ymd("2022-10-26")) ~ "C off",
    (reactor != "SBR1" & between(date, ymd("2022-10-25"), ymd("2023-05-12"))) ~ "C on",
    (reactor != "SBR1" & date > ymd("2023-05-12")) ~ "C off"
  ))

fermenters <- read_csv("/projects/b1052/mckenna/resources/qiime_v5.1/fermenters_5.1.csv")

rel_ferm <- rel_genus %>%
  filter(genus_format %in% fermenters$`Canonical Name`) %>%
  filter(!genus_format %in% c("Ca Competibacter", "Ca Obscuribacter")) %>%
  mutate(genus_format = str_replace(genus_format, "Ca ", "Ca. "))
  
rel_ferm_wide <- rel_ferm %>%
  ungroup() %>%
  select(sample, Genus, sum) %>%
  pivot_wider(names_from = sample, values_from = sum) %>%
  column_to_rownames("Genus")

rel_ferm_wide <- as.data.frame(apply(rel_ferm_wide, 2, function(x){x / sum(x)})) # normalize by sample


rm(fermenters)

# plotting 

# evenness qualitative ----
# boxplot ----
even_dfs <- list()
even_plts <- list()
for(i in reactors){
  
  tmp <- rel_ferm %>% 
    filter(reactor == i) %>%
    group_by(reactor, genus_format) %>%
    summarise(median = median(sum)) %>%
    ungroup() %>%
    slice_max(order_by = median, n = 10) %>%
    select(-reactor)
  
  even_df <- rel_ferm %>%
    filter(reactor == i) %>%
    filter(genus_format %in% tmp$genus_format) %>%
    left_join(tmp, by = "genus_format")
    
  even_dfs[[i]] <- even_df
  
  even_plt <- ggplot(data = even_df, aes(x = reorder(genus_format, -median), y = sum)) +
    geom_boxplot(width = 0.5) + 
    labs(x = "", y = "Median relative abundance [%]", x = "Genus") +
    theme_black_box + theme(axis.text.x = element_text(face = "italic")) +
    coord_flip()
  
  even_plts[[i]] <- even_plt
}

# barplot -----
even_dfs <- list()
even_plts <- list()
reactor_cols_mapped <- list("SBR1" = "#06d6a0", 
                            "SBR2" = "#118ab2",
                            "SBR3" = "#073b4c")
for(i in reactors){
  
  tmp <- rel_ferm %>% 
    filter(reactor == i) %>%
    group_by(reactor, genus_format) %>%
    summarise(median = median(sum)) %>%
    ungroup() %>%
    slice_max(order_by = median, n = 10) %>%
    select(-reactor)
  
  even_plt <- ggplot(data = tmp, aes(x = reorder(genus_format, -median), y = median)) +
    geom_bar(stat = "identity", fill = reactor_cols_mapped[[i]]) + 
    labs(y = "Median relative abundance [%]", x = "Genus", title = paste(i)) + 
    coord_flip() + 
    ylim(0, 1.7) + 
    theme_black_box + theme(axis.text.y = element_text(face = "italic"))
  
  even_plts[[i]] <- even_plt
}

even_plts[["SBR1"]] + even_plts[["SBR2"]] + even_plts[["SBR3"]] + plot_layout(axis_titles = "collect")
ggsave("results/fermenter_even_bar.png", width = 7, height = 2.5, units = "in", dpi = 500)

# diversity ----
div <- as.data.frame(diversity(rel_ferm_wide, index = "shannon", MARGIN = 2))
colnames(div) <- "shannon"

div <- div %>% 
  rownames_to_column(var = "sample") %>%
  left_join(metadata, by = "sample")
  
ggplot(div, aes(x = reactor, y= shannon)) +
  geom_boxplot(outlier.color = NULL) +
  geom_point(position = position_jitter(width = 0.1)) +
  stat_compare_means(method = "wilcox", comparisons = my_comparisons) + 
  theme_black_box 


# evenness ----
S <- apply( rel_ferm_wide > 0, 2 , sum )
even <- as.data.frame(diversity(rel_ferm_wide, index = "shannon", MARGIN = 2) / log(S))
colnames(even) <- "shannon"

even <- even %>% 
  rownames_to_column(var = "sample") %>%
  left_join(metadata, by = "sample")

ggplot(even, aes(x = reactor, y= shannon)) +
  geom_boxplot(outlier.color = NULL) +
  geom_point(position = position_jitter(width = 0.1)) +
  stat_compare_means(method = "wilcox", comparisons = my_comparisons) + 
  theme_black_box 

# richness ----
richness <- as.data.frame(apply( rel_ferm_wide > 0, 2 , sum ))
colnames(richness) <- "richness"

richness <- richness %>% 
  rownames_to_column(var = "sample") %>%
  left_join(metadata, by = "sample")

richness %>% 
  group_by(reactor) %>%
  summarise(median(richness))

plt_rich <- ggplot(richness, aes(x = reactor, y = richness, color = reactor)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.7) +
  stat_compare_means(method = "wilcox", comparisons = my_comparisons, 
                     family = "FreeSerif") + 
  scale_color_reactor + 
  labs(x = "Reactor", y = "Number of fermentative genera") + 
  scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.1, 0.1))) +
  theme_black_box +
  theme(legend.position = "none")

# total relative abundance
plt_relab <- rel_ferm %>% 
  group_by(reactor, date) %>%
  summarise(tot = sum(sum)) %>%
  ggplot(data = ., aes(x = reactor, y = tot, color = reactor)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.7) +
  scale_color_reactor + 
  labs(x = "Reactor", y = "Relative abundance of fermentative genera [%]") +
  stat_compare_means(method = "wilcox", comparisons = my_comparisons, 
                     family = "FreeSerif") + 
  scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.1, 0.1))) +
  theme_black_box +
  theme(legend.position = "none")

plt_relab + plt_rich + plot_layout(axis_titles = "collect_x")
ggsave("results/ferm_relab_rich.png", width = 7, height = 3.5, units = "in", dpi = 500)

# presence absence heatmap
sample_order <- metadata %>%
  select(sample) %>%
  as.list()

top_n_ferm <- rel_ferm %>%
  ungroup() %>%
  group_by(Genus) %>%
  summarise(median = median(sum)) %>%
  filter(median >= 0.001)

rel_ferm_wide_pa <- rel_ferm_wide %>%
  filter(row.names(rel_ferm_wide) %in% top_n_ferm$Genus) %>%
  mutate_all(~ifelse(. == 0, 0, 1)) %>%
  as.matrix()

rel_ferm_wide_pa <- rel_ferm_wide_pa[, sample_order[["sample"]]]

my_group <- as.numeric(as.factor(metadata$reactor))
colSide <- reactor_cols[my_group]

metadata_pheat <- metadata %>%
  column_to_rownames("sample") %>%
  select(reactor) %>%
  mutate(reactor = as.factor(reactor))

ferm_pheat <- pheatmap(rel_ferm_wide_pa, treeheight_row = 0, treeheight_col = 0, 
                       legend_breaks = c(0,1), 
                       legend_labels = c("Absent", "Present"), 
                       color = c("#e5e5e5", "#fca311"),
                       annotation_col = metadata_pheat, 
                       annotation_colors = list(reactor = c(SBR1 = reactor_cols[1], 
                                                            SBR2 = reactor_cols[2],
                                                            SBR3 = reactor_cols[3])))

png("./results/ferm_heatmap.png", width=6, height=12, units="in", res=1200)
print(ferm_pheat)  ## or just `p`
dev.off()





