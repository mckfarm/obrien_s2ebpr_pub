### cor checks ###

library(tidyverse)
library(reshape2)
library(rstatix)

source("scripts/plotting_cols_shapes.R")

rel_pao_gao <- readRDS("data/rel_pao_gao.RDS")

pao_gao_list <- c("Ca_Accumulibacter", "Ca_Obscuribacter", "Ca_Phosphoribacter",
                  "Dechloromonas", "Tetrasphaera", 
                  "Ca_Competibacter", "Ca_Contendobacter", "Defluviicoccus")

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



wilcox_list <- list()
for(i in pao_gao_list){
  diff_test <- rel_pao_gao %>%
    filter(reactor != "SBR3") %>%
    filter(Genus == paste(i)) %>%
    ungroup() %>%
    select(-c(type, sample, Genus)) %>%
    pivot_wider(names_from = reactor, values_from = sum) %>% # have to do this weird transformation to get complete cases
    drop_na() %>%
    pivot_longer(!date, values_to = "sum", names_to = "reactor") %>%
    wilcox_test(sum ~ reactor)
  wilcox_list[[i]] <- diff_test
  
  # list idea from https://stackoverflow.com/questions/29402528/append-data-frames-together-in-a-for-loop
  
}

wilcox_df <- do.call(rbind, wilcox_list) %>%
  rownames_to_column(var = "Genus") %>%
  mutate(Genus = str_replace(Genus, "Ca_", "Ca. "))

wilcox_df


pearson_list <- list()
for(i in pao_gao_list){
  diff_test <- rel_pao_gao %>%
    filter(reactor != "SBR3") %>%
    filter(Genus == paste(i)) %>%
    ungroup() %>%
    select(-c(type, sample, Genus)) %>%
    pivot_wider(names_from = reactor, values_from = sum) %>% # have to do this weird transformation to get complete cases
    drop_na() %>%
    cor_test(SBR1, SBR2, method = "pearson")
  pearson_list[[i]] <- diff_test
  
  # list idea from https://stackoverflow.com/questions/29402528/append-data-frames-together-in-a-for-loop
  
}

pearson_df <- do.call(rbind, pearson_list) %>%
  rownames_to_column(var = "Genus") %>%
  mutate(Genus = str_replace(Genus, "Ca_", "Ca. "))

pearson_df
