### nitrospira corr ###

library(tidyverse)
library(rstatix)
library(MetBrewer)

source("scripts/plotting_dates.R")
source("scripts/plotting_cols_shapes.R")

rel_genus <- readRDS("data/rel_genus.RDS") %>%
  mutate(genus_format = str_replace(Genus, "Ca_", "Ca. "))


# correlation -----
reactor_list <- c("SBR1", "SBR2", "SBR3")

pearson_list <- list()
for(i in reactor_list){
  df_cor <- rel_genus %>% 
    ungroup() %>%
    select(-genus_format) %>%
    filter(reactor == paste(i)) %>%
    pivot_wider(names_from = Genus, values_from = sum) %>%
    select(-c(date, reactor, sample)) %>%
    select_if(~median(.) > 0.05)
  
  test_pear <- df_cor %>%
    cor_test(vars = Ca_Competibacter, method = "pearson") 
  
  test_pear_filt <- test_pear %>%
    mutate(reactor = paste(i))
  
  pearson_list[[i]] <- test_pear_filt
}

pearson_df <- do.call(rbind, pearson_list)

pearson_df %>%
  filter(var2 %in% c("Ca_Phosphoribacter", "Nitrospira"))

rel_genus %>%
  filter(Genus == "Nitrospira") %>%
  ungroup() %>%
  select(-c(sample, Genus, genus_format)) %>%
  pivot_wider(names_from = reactor, values_from = sum) %>% # have to do this weird transformation to get complete cases
  drop_na() %>%
  pivot_longer(!date, values_to = "sum", names_to = "reactor") %>%
  kruskal_test(sum ~ reactor)


# cor plot -----
rel_genus %>% 
  filter(Genus %in% c("Nitrospira","Ca_Competibacter", "Hassallia", "midas_g_179")) %>%
  ggplot(., aes(x = date, y = sum, color = genus_format, shape = genus_format)) +
  facet_wrap(~reactor) +
  geom_line(linewidth = 0.3) +
  geom_point(size = 2.5, fill = "white", alpha = 0.8) + 
  scale_shape_manual(values = c(3, 25, 24, 22), name = "Genus") + 
  scale_color_manual(values = met.brewer("Johnson", 4), name = "Genus") + 
  theme_black_box + theme(legend.text = element_text(face = "italic"), legend.position = "left") + 
  x_axis_date +
  labs(x = "Date", y = "Relative abundance [%]")
ggsave("results/nitrospira_cor.png", width = 10, height = 3.5, units = "in", dpi = 300)


rel_genus %>% 
  mutate(perf = case_when(
    reactor == "SBR1" ~ "C off", 
    (reactor != "SBR1" & date < ymd("2022-10-26")) ~ "C off",
    (reactor != "SBR1" & between(date, ymd("2022-10-25"), ymd("2023-05-12"))) ~ "C on",
    (reactor != "SBR1" & date > ymd("2023-05-12")) ~ "C off"
  )) %>%
  filter(genus_format %in% c("Ca. Competibacter", "Ca. Phosphoribacter", "Nitrospira")) %>%
  # filter(between(date, (phases$x0 - 10), phases$x4)) %>%
  filter(date <= ymd("2023-07-31")) %>%
  ggplot(data = ., aes(x = date, y = sum, color = genus_format)) +
  facet_wrap(~reactor) + 
  geom_rect(aes(xmin = phases$x2, xmax = phases$x4, ymin = -Inf, ymax = Inf), fill = "#FFEBF2", color = NA) +
  geom_line(linewidth = 0.2) +
  geom_point(size = 2.2, aes(shape = perf), fill = "white", alpha = 0.7) +
  scale_shape_carb + 
  scale_color_manual(values = c("#90be6d", "#be95c4", "#f9844a"), name = "Genus") +
  x_axis_date + 
  ylim(0, 2.9) +
  labs(x = "", y = "Relative abundance [%]") +
  theme_black_box + theme(axis.text = element_text(size = 8)) +
  guides(color = guide_legend(override.aes = list(size = 1.5)),
         shape = guide_legend(override.aes = list(size = 1.5)))

ggsave("results/pao_gao_sbr3_maintext.png", width = 9, height = 2.5, units = "in", dpi = 500)
