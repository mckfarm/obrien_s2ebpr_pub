### performance corr check ### 

# packages
library(tidyverse)
library(readxl)
library(rstatix)
library(MetBrewer)
library(patchwork)

source("./scripts/plotting_cols_shapes.R")

# data read in and parse
ml <- read_excel("data/obrien_data.xlsx", sheet = "mixedliquor") %>% 
  select(-source)

rel_pao_gao <- readRDS("data/rel_pao_gao.RDS")

rel_pao_gao_wide <- rel_pao_gao %>%
  ungroup() %>%
  select(-c(type, sample)) %>%
  pivot_wider(names_from = Genus, values_from = sum)


pao_gao_list_format <- c("Ca. Accumulibacter", "Ca. Obscuribacter", "Ca. Phosphoribacter",
                         "Dechloromonas", "Tetrasphaera")

reactors <- c("SBR1", "SBR2", "SBR3")

# correlation matrices ----

p_perc_list <- list()
for (i in reactors){
  ml_temp <- ml %>% 
    filter(reactor == paste(i)) %>%
    select(c(date, P_perc))
  
  df_all <- rel_pao_gao_wide %>% 
    filter(reactor == paste(i)) %>%
    select(-reactor) %>%
    left_join(ml_temp, by = join_by(date)) %>%
    select(-date)
  
  df_corr <- df_all %>%
    cor_test(vars = "P_perc", method = "pearson") %>%
    mutate(reactor = paste(i))
  
  p_perc_list[[i]] <- df_corr
  
}

p_perc_df <- do.call(rbind, p_perc_list) 


p_rel_list <- list()
for (i in reactors){
  ml_temp <- ml %>% 
    filter(reactor == paste(i)) %>%
    select(c(date, P_rel))
  
  df_all <- rel_pao_gao_wide %>% 
    filter(reactor == paste(i)) %>%
    select(-reactor) %>%
    left_join(ml_temp, by = join_by(date)) %>%
    select(-date)
  
  df_corr <- df_all %>%
    cor_test(vars = "P_rel", method = "pearson") %>%
    mutate(reactor = paste(i))
  
  p_rel_list[[i]] <- df_corr
  
}

p_rel_df <- do.call(rbind, p_rel_list) 



# plots ----
p_perc_df_plt <- p_perc_df %>%
  mutate(genus_format = str_replace(var2, "Ca_", "Ca. ")) %>%
  filter(genus_format %in% pao_gao_list_format) %>%
  mutate(genus_format = factor(genus_format, levels = pao_gao_list_format)) %>%
  mutate(cor_filt = ifelse(p <= 0.10, round(cor, 2), NA)) %>%
  mutate(p_lab = ifelse(p <= 0.10, signif(p, digits = 1), NA))

p_rel_df_plt <- p_rel_df %>%
  mutate(genus_format = str_replace(var2, "Ca_", "Ca. ")) %>%
  filter(genus_format %in% pao_gao_list_format) %>%
  mutate(genus_format = factor(genus_format, levels = pao_gao_list_format)) %>%
  mutate(cor_filt = ifelse(p <= 0.10, round(cor, 2), NA)) %>%
  mutate(p_lab = ifelse(p <= 0.10, signif(p, digits = 1), NA))

heatmap_theme <- list(
  scale_y_discrete(limits = rev),
  scale_fill_gradientn(colors = met.brewer("Cassatt1"), limits = c(-1, 1), na.value = "lightgrey"),
  theme_bw(),
  theme(legend.position = "none")
)

plt_perc <- 
ggplot(p_perc_df_plt, aes(x = reactor, y = genus_format, fill = cor_filt)) +
  geom_tile() +
  geom_text(aes(label = p_lab), nudge_y = -0.2, size = 3) + 
  geom_text(aes(label = cor_filt), nudge_y = 0.2) + 
  heatmap_theme +
  labs(title = "P content [% P in VSS]", y = "Genus", x = "Reactor") +
  theme(axis.text.y = element_text(face = "italic"))

plt_rel <- 
ggplot(p_rel_df_plt, aes(x = reactor, y = genus_format, fill = cor_filt)) +
  geom_tile() +
  geom_text(aes(label = p_lab), nudge_y = -0.2, size = 3) + 
  geom_text(aes(label = cor_filt), nudge_y = 0.2) + 
  heatmap_theme +
  labs(title = "P release [mgP/", y = "", x = "Reactor") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

plt_perc + plt_rel 
ggsave("results/biomass_corr.png", width = 6, height = 3, units = "in", dpi = 300)


# Pearson correlation with p perc by taxa
  
df_labs <- rel_pao_gao %>%
  ungroup() %>%
  select(reactor, Genus, sum, date) %>%
  group_by(reactor, Genus) %>%
  slice_max(order_by = sum, n = 1) %>%
  left_join(select(p_perc_df, var2, cor, p, reactor), by = c("Genus" = "var2",
                                                              "reactor" = "reactor")) %>%
  rename(lab_y = sum) %>%
  mutate(lab_y = lab_y * 1.15) %>%
  ungroup() %>% 
  group_by(Genus) %>% 
  mutate(lab_y = max(lab_y)) %>% 
  ungroup() %>%
  mutate(label_full = paste("r =", cor, "p =", p))
# im sure there's a better way to do this but this assigns the correlation coefs to a maximum y value
# the max y value is based on the max abundance of each genus because of the facet grid with fixed y axes


df_all <- rel_pao_gao %>%
  left_join(ml, by = join_by(date, reactor)) %>%
  left_join(df_labs, by = join_by(date, reactor, Genus)) %>%
  mutate(genus_format = str_replace(Genus, "Ca_", "Ca. "))
  

df_all %>%
  filter(type == "PAO") %>%
  ggplot(., aes(x = P_perc, y = sum)) +
  facet_grid(genus_format~reactor, scales = "free_y") +
  geom_point(size = 1) +
  geom_text(aes(label = label_full, y = lab_y, x = 0.04, family = def_font), size = 3) + 
  theme_black_box +
  scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.1, 0.2))) +
  labs(y = "Relative abundance [%]", x = "Biomass P content [mgP/mgVSS]")
ggsave("results/p_perc_rel_ab.png", width = 5, height = 6, units = "in", dpi = 300)

