### pao and gao ###

### packages -----
library(tidyverse)
library(ggpubr)
library(reshape2)
library(rstatix)

source("scripts/plotting_cols_shapes.R")
source("scripts/plotting_dates.R")

### data read in ----
rel_genus <- readRDS("data/rel_genus.RDS")

pao_gao_list_format <- c("Ca. Accumulibacter", "Ca. Obscuribacter", "Ca. Phosphoribacter",
                         "Dechloromonas", "Tetrasphaera", 
                         "Ca. Competibacter", "Ca. Contendobacter", "Defluviicoccus")

rel_pao_gao <- rel_genus %>%
  mutate(genus_format = str_replace(Genus, "Ca_", "Ca. ")) %>%
  filter(genus_format %in% pao_gao_list_format) %>%
  mutate(genus_format = factor(genus_format, levels = pao_gao_list_format))


rel_pao_gao %>% 
  filter(reactor != "SBR3") %>%
  filter(!genus_format %in% c("Ca. Obscuribacter", "Ca. Contendobacter")) %>%
  ggplot(data = ., aes(x = date, y = sum, color = reactor, shape = reactor)) +
  facet_wrap(~genus_format) + 
  geom_rect(aes(xmin = phases$x0, xmax = phases$red_C, ymin = -Inf, ymax = Inf), fill = "#E8D9FC", color = NA) + 
  geom_rect(aes(xmin = phases$red_C, xmax = phases$no_C, ymin = -Inf, ymax = Inf), fill = "#FFD6E8", color = NA) + 
  geom_point(size = 3) +
  geom_line(linewidth = 0.3) +
  scale_color_reactor +
  scale_shape_reactor + 
  theme_black_box 

my_comparisons <- list(c("SBR1", "SBR2"),
                       c("SBR2", "SBR3"),
                       c("SBR1", "SBR3"))

rel_pao_gao %>% 
  filter(!genus_format %in% c("Ca. Obscuribacter", "Ca. Contendobacter")) %>%
  ggplot(data = ., aes(x = reactor, y = sum, color = reactor)) +
  facet_wrap(~genus_format, scales = "free") + 
  geom_boxplot(outlier.color = NA) + 
  geom_point(position = position_jitterdodge()) + 
  scale_color_reactor + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", 
                     label = "p.format", tip.length = 0.02, size = 3) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  labs(x = "", y = "Relative abundance [%]") +
  theme_black_box + theme(legend.position = "none", 
                          strip.text = element_text(face = "italic")) 

### omfg why didnt I google this sooner 
### https://stackoverflow.com/questions/67813734/adjusting-y-axis-limits-in-ggplot2-with-facet-and-free-scales
