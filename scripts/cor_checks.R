### cor checks ###
# packages
library(tidyverse)
library(reshape2)
library(rstatix)
library(MetBrewer)
library(reshape2)


# read in
rel_pao_gao <- readRDS("data/rel_pao_gao.RDS")

pao_gao_list <- c("Ca_Accumulibacter", "Ca_Obscuribacter", "Ca_Phosphoribacter",
                  "Dechloromonas", "Tetrasphaera", 
                  "Ca_Competibacter", "Ca_Contendobacter", "Defluviicoccus")


# kruskal test
# this measures if ALL are different from each other, similar to anova but doesnt assume normality
kruskal_list <- list()
for(i in pao_gao_list){
  diff_test <- rel_pao_gao %>%
    filter(Genus == paste(i)) %>%
    ungroup() %>%
    select(-c(type, sample, Genus)) %>%
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

kruskal_df


# wilcox for sbr1 and sbr2

wilcox_list <- list()
for(i in pao_gao_list){
  diff_test <- rel_pao_gao %>%
    filter(reactor != "SBR3") %>%
    filter(date <= ymd("2023-07-31")) %>%
    filter(Genus == paste(i)) %>%
    ungroup() %>%
    select(-c(type, sample, Genus)) %>%
    pivot_wider(names_from = reactor, values_from = sum) %>% # have to do this weird transformation to get complete cases
    drop_na() %>%
    pivot_longer(!date, values_to = "sum", names_to = "reactor") %>%
    wilcox_test(sum ~ reactor)
  wilcox_list[[i]] <- diff_test
  
}

wilcox_df <- do.call(rbind, wilcox_list) %>%
  rownames_to_column(var = "Genus") %>%
  mutate(Genus = str_replace(Genus, "Ca_", "Ca. "))

wilcox_df

# pearson for sbr1 and sbr2

pearson_list <- list()
for(i in pao_gao_list){
  diff_test <- rel_pao_gao %>%
    filter(reactor != "SBR3") %>%
    filter(date <= ymd("2023-07-31")) %>%
    filter(Genus == paste(i)) %>%
    ungroup() %>%
    select(-c(type, sample, Genus)) %>%
    pivot_wider(names_from = reactor, values_from = sum) %>% # have to do this weird transformation to get complete cases
    drop_na() %>%
    cor_test(SBR1, SBR2, method = "pearson")
  pearson_list[[i]] <- diff_test
  
}

pearson_df <- do.call(rbind, pearson_list) %>%
  rownames_to_column(var = "Genus") %>%
  mutate(Genus = str_replace(Genus, "Ca_", "Ca. "))

pearson_df


# spearman for taxa within each reactor

pao_gao_list_format <- c("Ca. Accumulibacter", "Ca. Obscuribacter", "Ca. Phosphoribacter",
                         "Dechloromonas", "Tetrasphaera", 
                         "Ca. Competibacter", "Ca. Contendobacter", "Defluviicoccus")

reactors <- c("SBR1", "SBR2", "SBR3")

rel_pao_gao_wide <- rel_pao_gao %>%
  ungroup() %>%
  mutate(Genus = str_replace(Genus, "Ca_", "Ca. ")) %>%
  select(-c(type, sample)) %>%
  pivot_wider(names_from = Genus, values_from = sum)

# heatmap_theme <- list(
#   scale_y_discrete(limits = rev),
#   scale_fill_gradientn(colors = met.brewer("Cassatt1"), limits = c(-1, 1), na.value = "lightgrey"),
#   theme_bw(),
#   theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
# )

# correlation matrices ----

pg_cor_list <- list()
for (i in reactors){
  df_corr <- rel_pao_gao_wide %>% 
    filter(reactor == paste(i)) %>%
    select(-c(reactor, date)) %>%
    cor_mat(vars = NULL, method = "pearson")
  
  pg_cor_list[[i]] <- df_corr
  
}

# pg_cor_df <- do.call(rbind, pg_cor_list) %>%
#   rename(all_of(c(var1 = "rowname", var2 = "variable"))) %>%
#   mutate(var1 = factor(var1, levels = pao_gao_list_format),
#          var2 = factor(var2, levels = pao_gao_list_format)) %>%
#   mutate(cor_filt = ifelse(p <= 0.10, round(cor, 2), NA)) %>%
#   mutate(p_lab = ifelse(p <= 0.10, signif(p, digits = 1), NA))


pg_cor_list[["SBR1"]] %>% 
  pull_lower_triangle() %>%
  cor_plot(
    method = "number",
    label = TRUE,
    palette = met.brewer("Cassatt1", 200)) 

pg_cor_list[["SBR2"]] %>% 
  pull_lower_triangle() %>%
  cor_plot(
    method = "number",
    label = TRUE,
    palette = met.brewer("Cassatt1", 200)) 

pg_cor_list[["SBR3"]] %>% 
  pull_lower_triangle() %>%
  cor_plot(
    method = "number",
    label = TRUE,
    palette = met.brewer("Cassatt1", 200)) 
