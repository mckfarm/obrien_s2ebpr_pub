### nitrospira corr ###

library(tidyverse)
library(rstatix)
library(MetBrewer)

source("scripts/plotting_dates.R")
source("scripts/plotting_cols_shapes.R")

rel_genus <- readRDS("data/rel_genus.RDS") %>%
  mutate(genus_format = str_replace(Genus, "Ca_", "Ca. "))


rel_genus %>% 
  filter(Genus %in% c("Nitrospira","Ca_Competibacter", "Hassallia", "midas_g_179")) %>%
  ggplot(., aes(x = date, y = sum, color = genus_format, shape = genus_format)) +
  facet_wrap(~reactor) +
  geom_line(linewidth = 0.3) +
  geom_point(size = 2.5, fill = "white", alpha = 0.8) + 
  scale_shape_manual(values = c(3, 25, 24, 22), name = "Genus") + 
  scale_color_manual(values = met.brewer("Johnson", 4), name = "Genus") + 
  theme_black_box + theme(legend.text = element_text(face = "italic")) + 
  x_axis_date +
  labs(x = "Date", y = "Relative abundance [%]")
ggsave("results/nitrospira_cor.png", width = 10, height = 3.5, units = "in", dpi = 300)


df_cor <- rel_genus %>% 
  ungroup() %>%
  filter(Genus %in% taxa_filt$Genus) %>%
  filter(date <= ymd("2023-07-31")) %>%
  filter(reactor == "SBR3") %>%
  pivot_wider(names_from = Genus, values_from = sum) %>%
  select(-c(date, reactor, sample))

test_pear <- df_cor %>%
  cor_test(vars = Ca_Competibacter, method = "pearson") 

test_pear_filt <- test_pear %>%
  filter(p <= 0.05) 

