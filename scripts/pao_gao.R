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

### plots ----
# main text rel ab over time
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

# supplemental all pao rel ab over time
rel_pao_gao %>% 
  filter(date <= ymd("2023-07-31")) %>%
  ggplot(data = ., aes(x = date, y = sum, color = reactor)) +
  facet_wrap(~genus_format, scales = "free") + 
  # geom_rect(aes(xmin = phases$x0, xmax = phases$red_C, ymin = -Inf, ymax = Inf), fill = "#E8D9FC", color = NA) + 
  # geom_rect(aes(xmin = phases$red_C, xmax = phases$no_C, ymin = -Inf, ymax = Inf), fill = "#FFD6E8", color = NA) + 
  geom_line(linewidth = 0.2) +
  geom_point(size = 2, aes(shape = perf), fill = "white", alpha = 0.8) +
  scale_color_reactor +
  scale_shape_carb + 
  x_axis_date + 
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.2))) +
  labs(x = "", y = "Relative abundance [%]") +
  theme_black_box_taxa
ggsave("results/pao_gao_supplemental.png", width = 10, height = 6, units = "in", dpi = 300)

### sbr3 main text figure in nitrospira_corr.R


# comparisons of abundance by reactor configuration?

my_comparisons <- list(c("SBR1", "SBR2"),
                       c("SBR2", "SBR3"),
                       c("SBR1", "SBR3"))
compare_pao <- 
rel_pao_gao %>% 
  filter(genus_format %in% pao_list) %>%
  filter(date <= ymd("2023-07-31")) %>%
  ggplot(data = ., aes(x = reactor, y = sum, color = reactor)) +
  facet_wrap(~genus_format, scales = "free") + 
  geom_boxplot(outlier.color = NA) + 
  geom_point(position = position_jitterdodge()) + 
  scale_color_reactor + 
  stat_compare_means(method = "wilcox", comparisons = my_comparisons, size = 2.8, family = "serif") +
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
  stat_compare_means(method = "wilcox", comparisons = my_comparisons, size = 2.8, family = "serif") +
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
  stat_compare_means(method = "wilcox", size = 2.8, label = "p.format", vjust = -1, family = "serif") +
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
  stat_compare_means(method = "wilcox", size = 2.8, label = "p.format", vjust = -1, family = "serif") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.2))) +
  labs(x = "", y = "") +
  theme_black_box_taxa + theme(legend.position = "none") 


compare_pao / compare_gao + plot_layout(heights = c(2, 1))
ggsave("results/pao_gao_compare_mode.png", width = 6, height = 6, units = "in", dpi = 500)


### omfg why didnt I google this sooner 
### https://stackoverflow.com/questions/67813734/adjusting-y-axis-limits-in-ggplot2-with-facet-and-free-scales

