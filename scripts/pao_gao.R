### pao and gao ###

### packages -----
library(tidyverse)
library(ggpubr)
library(rstatix)
library(patchwork)

source("scripts/plotting_cols_shapes.R")
source("scripts/plotting_dates.R")

### data read in ----
pao_gao_list_format <- c("Ca. Accumulibacter", "Ca. Obscuribacter", "Ca. Phosphoribacter",
                         "Dechloromonas", "Tetrasphaera", 
                         "Ca. Competibacter", "Ca. Contendobacter", "Defluviicoccus")

pao_list <- c("Ca. Accumulibacter", "Ca. Obscuribacter", "Ca. Phosphoribacter",
                         "Dechloromonas", "Tetrasphaera")
gao_list <- c("Ca. Competibacter", "Ca. Contendobacter", "Defluviicoccus")

pao_gao_list <- c("Ca_Accumulibacter", "Ca_Obscuribacter", "Ca_Phosphoribacter",
                  "Dechloromonas", "Tetrasphaera", 
                  "Ca_Competibacter", "Ca_Contendobacter", "Defluviicoccus")

rel_pao_gao <- readRDS("data/rel_pao_gao.RDS") %>%
  mutate(genus_format = str_replace(Genus, "Ca_", "Ca. ")) %>%
  filter(genus_format %in% pao_gao_list_format) %>%
  mutate(genus_format = factor(genus_format, levels = pao_gao_list_format)) %>%
  mutate(perf = case_when(
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

### kruskal stats
kruskal_list <- list()
for(i in pao_gao_list){
  diff_test <- rel_pao_gao %>%
    filter(Genus == paste(i)) %>%
    ungroup() %>%
    select(-c(type, sample, Genus, perf, op_mode, genus_format)) %>%
    pivot_wider(names_from = reactor, values_from = sum) %>% # have to do this weird transformation to get complete cases
    drop_na() %>%
    pivot_longer(!date, values_to = "sum", names_to = "reactor") %>%
    kruskal_test(sum ~ reactor)
  kruskal_list[[i]] <- diff_test
  
  # list idea from https://stackoverflow.com/questions/29402528/append-data-frames-together-in-a-for-loop
  
}

kruskal_df <- do.call(rbind, kruskal_list) %>%
  rownames_to_column(var = "Genus") %>%
  mutate(Genus = str_replace(Genus, "Ca_", "Ca. "))

kruskal_for_plot <- kruskal_df %>%
  select(Genus, statistic, p) %>%
  rename(genus_format = Genus) %>%
  mutate(p = round(p, digits = 3))

df_labs <- rel_pao_gao %>%
  filter(date <= ymd("2023-07-31")) %>%
  ungroup() %>%
  select(reactor, genus_format, sum, date) %>%
  group_by(genus_format) %>%
  slice_max(order_by = sum, n = 1) %>%
  left_join(kruskal_for_plot) %>%
  rename(lab_y = sum) %>%
  mutate(lab_y = lab_y * 1.2) %>%
  ungroup() %>% 
  group_by(genus_format) %>% 
  mutate(lab_y = max(lab_y)) %>% 
  ungroup() %>%
  mutate(label_full = paste("Kruskal Wallis p =", p))

df_pao_gao_plt <- rel_pao_gao %>%
  left_join(df_labs, by = join_by(date, genus_format, reactor))

### plots ----
# main text rel ab over time
plt_pao <-
df_pao_gao_plt %>% 
  filter(date <= ymd("2023-07-31")) %>%
  filter(type == "PAO") %>%
  ggplot(data = ., aes(x = date, y = sum, color = reactor)) +
  facet_wrap(~genus_format, scales = "free") +  
  geom_line(linewidth = 0.2) +
  geom_point(size = 2, aes(shape = perf), fill = "white", alpha = 0.8) +
  geom_text(aes(label = label_full, y = lab_y, x = ymd("2022-10-20"), family = def_font), 
            size = 3, color = "black") + 
  scale_color_reactor +
  scale_shape_carb + 
  x_axis_date + 
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.2))) +
  labs(x = "", y = "Relative abundance [%]") +
  theme_black_box_taxa + 
  theme(legend.position = c(0.85, 0.2), 
        axis.title.y = element_text(hjust = 0),
        legend.spacing.y = unit(0.01, 'cm')) + 
  guides(color=guide_legend(ncol = 3), shape = guide_legend(ncol = 2))

plt_gao <- 
  df_pao_gao_plt %>% 
  filter(date <= ymd("2023-07-31")) %>%
  filter(type == "GAO") %>%
  ggplot(data = ., aes(x = date, y = sum, color = reactor)) +
  facet_wrap(~genus_format, scales = "free") + 
  geom_line(linewidth = 0.2) +
  geom_point(size = 2, aes(shape = perf), fill = "white", alpha = 0.8) +
  geom_text(aes(label = label_full, y = lab_y, x = ymd("2022-10-20"), family = def_font), 
            size = 3, color = "black") + 
  scale_color_reactor +
  scale_shape_carb + 
  x_axis_date + 
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.2))) +
  labs(x = "", y = "") +
  theme_black_box_taxa + theme(legend.position = "none", plot.margin = unit(c(0, 0.1, 0.1, 0.1), "cm"))

plt_pao / plt_gao + plot_layout(heights = c(1, 0.5)) 

ggsave("results/pao_gao_main_text_v2.png", width = 8, height = 6, units = "in", dpi = 400) 

### sbr3 main text figure in nitrospira_corr.R


# comparisons of abundance by reactor configuration?

my_comparisons <- list(c("SBR1", "SBR2"),
                       c("SBR2", "SBR3"),
                       c("SBR1", "SBR3"))

# custom stat compare means
fixed_call <- quote(ggsignif::geom_signif(comparisons = comparisons, 
                                          y_position = label.y, 
                                          test = method, test.args = method.args, 
                                          step_increase = step.increase, 
                                          size = bracket.size, textsize = 3, color = color, 
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


compare_pao <-
rel_pao_gao %>% 
  filter(genus_format %in% pao_list) %>%
  filter(date <= ymd("2023-07-31")) %>%
  ggplot(data = ., aes(x = reactor, y = sum, color = reactor)) +
  facet_wrap(~genus_format, scales = "free") + 
  geom_boxplot(outlier.color = NA) + 
  geom_point(position = position_jitterdodge()) + 
  scale_color_reactor + 
  stat_compare_means(method = "wilcox", comparisons = my_comparisons, family = def_font) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
  labs(x = "", y = "") +
  theme_black_box_taxa + theme(legend.position = "none") 

compare_gao <- 
rel_pao_gao %>% 
  filter(genus_format %in% gao_list) %>%
  filter(date <= ymd("2023-07-31")) %>%
  ggplot(data = ., aes(x = reactor, y = sum, color = reactor)) +
  facet_wrap(~genus_format, scales = "free") + 
  geom_boxplot(outlier.color = NA) + 
  geom_point(position = position_jitterdodge()) + 
  scale_color_reactor + 
  stat_compare_means(method = "wilcox", comparisons = my_comparisons, family = def_font) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
  labs(x = "", y = "") +
  theme_black_box_taxa + theme(legend.position = "none") 

compare_pao / compare_gao + plot_layout(heights = c(2, 1))
ggsave("results/pao_gao_compare_sbr.png", width = 6, height = 6, units = "in", dpi = 500)



compare_pao <- 
  rel_pao_gao %>% 
  filter(genus_format %in% pao_list) %>%
  ggplot(data = ., aes(x = op_mode, y = sum, color = op_mode)) +
  facet_wrap(~genus_format, scales = "free") + 
  geom_boxplot(outlier.color = NA) + 
  geom_point(position = position_jitterdodge()) + 
  scale_color_manual(values = mode_cols) + 
  ggpubr::stat_compare_means(method = "wilcox", size = 3, label = "p.format", vjust = -1, family = def_font) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.2))) +
  labs(x = "", y = "") +
  theme_black_box_taxa + theme(legend.position = "none") 

compare_gao <- 
  rel_pao_gao %>% 
  filter(genus_format %in% gao_list) %>%
  filter(date <= ymd("2023-07-31")) %>%
  ggplot(data = ., aes(x = op_mode, y = sum, color = op_mode)) +
  facet_wrap(~genus_format, scales = "free") + 
  geom_boxplot(outlier.color = NA) + 
  geom_point(position = position_jitterdodge()) + 
  scale_color_manual(values = mode_cols) + 
  ggpubr::stat_compare_means(method = "wilcox", size = 3, label = "p.format", vjust = -1, family = def_font) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.2))) +
  labs(x = "", y = "") +
  theme_black_box_taxa + theme(legend.position = "none") 

compare_pao / compare_gao + plot_layout(heights = c(2, 1))

ggsave("results/pao_gao_compare_mode.png", width = 6, height = 6, units = "in", dpi = 500)


### omfg why didnt I google this sooner 
### https://stackoverflow.com/questions/67813734/adjusting-y-axis-limits-in-ggplot2-with-facet-and-free-scales

## extra
# first round main text figure 
rel_pao_gao %>% 
  filter(reactor != "SBR3") %>%
  filter(date <= ymd("2023-07-31")) %>%
  filter(genus_format %in% c("Ca. Accumulibacter", "Ca. Phosphoribacter", "Tetrasphaera", "Ca. Competibacter")) %>%
  ggplot(data = ., aes(x = date, y = sum, color = reactor)) +
  facet_wrap(~genus_format) + 
  # geom_rect(aes(xmin = phases$x0, xmax = phases$red_C, ymin = -Inf, ymax = Inf), fill = "#E8D9FC", color = NA) + 
  # geom_rect(aes(xmin = phases$red_C, xmax = phases$no_C, ymin = -Inf, ymax = Inf), fill = "#FFD6E8", color = NA) + 
  geom_line(linewidth = 0.2) +
  geom_point(size = 2.5, aes(shape = perf), fill = "white") +
  scale_color_reactor +
  scale_shape_carb + 
  x_axis_date + 
  ylim(0, 2.2) + 
  labs(x = "", y = "Relative abundance [%]") +
  theme_black_box_taxa

ggsave("results/pao_gao_maintext.png", width = 6, height = 3.5, units = "in", dpi = 300)
