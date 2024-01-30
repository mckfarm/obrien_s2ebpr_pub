### performance corr check ### 

# packages
library(tidyverse)
library(readxl)
library(reshape2)
library(rstatix)

# data read in and parse
ml <- read_excel("data/obrien_data.xlsx", sheet = "mixedliquor") %>% 
  select(-source)

rel_pao_gao <- readRDS("data/rel_pao_gao.RDS")

rel_pao_gao_wide <- rel_pao_gao %>%
  ungroup() %>%
  select(-c(type, sample)) %>%
  pivot_wider(names_from = Genus, values_from = sum)

# correalation

reactors <- c("SBR1", "SBR2", "SBR3")

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
    cor_test(vars = "P_perc", method = "spearman") %>%
    mutate(reactor = paste(i))
  
  p_perc_list[[i]] <- df_corr
  
}

p_perc_df <- do.call(rbind, p_perc_list) 
p_perc_df



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
    cor_test(vars = "P_rel", method = "spearman") %>%
    mutate(reactor = paste(i))
  
  p_rel_list[[i]] <- df_corr
  
}

p_rel_df <- do.call(rbind, p_rel_list) 
p_rel_df


