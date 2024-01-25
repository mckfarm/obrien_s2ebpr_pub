### seasonal impacts ###

library(Maaslin2)
library(tidyverse)

# prep ----
rel_genus <- readRDS("data/rel_genus.RDS")
maaslin_genus <- rel_genus %>%
  ungroup() %>%
  pivot_wider(names_from = Genus, values_from = sum) %>%
  select(-c(date, reactor)) %>%
  column_to_rownames("sample")

saveRDS(maaslin_genus, "data/maaslin_genus.RDS")


metadata <- read_delim("data/metadata.txt", delim = "\t")
metadata$date <- mdy(metadata$date)

inf <- readRDS("data/influent.RDS")
inf_wide <- inf %>% select(-source) %>%
  pivot_wider(names_from = var, values_from = value)

maaslin_metadata <- metadata %>%
  left_join(inf_wide, by = join_by(date)) %>%
  mutate(perf = case_when(
    reactor == "SBR1" ~ "no C", 
    (reactor != "SBR1" & date < ymd("2022-10-26")) ~ "no C",
    (reactor != "SBR1" & between(date, ymd("2022-10-25"), ymd("2022-12-2"))) ~ "high C",
    (reactor != "SBR1" & between(date, ymd("2022-12-2"), ymd("2023-05-12"))) ~ "reduced C",
    (reactor != "SBR1" & date > ymd("2023-05-12")) ~ "no C"
  )) %>%
  mutate(op_mode = case_when(
    reactor == "SBR1" ~ "EBPR",
    reactor == "SBR2" ~ "S2EBPR",
    (reactor == "SBR3" & date < ymd("2023-03-06")) ~ "EBPR",
    (reactor == "SBR3" & date >= ymd("2023-03-06")) ~ "S2EBPR"
  )) %>%
  column_to_rownames(var = "sample")

saveRDS(maaslin_metadata, "data/maaslin_metadata.RDS")

fit_genus_sep <- Maaslin2(input_data = maaslin_genus,
                     input_metadata = maaslin_metadata,
                     min_abundance = 0.1,
                     min_prevalence = 0.5,
                     normalization = "TSS",
                     fixed_effects = c("reactor", "temp", "flow_2weekave", "perf"),
                     reference = c(("reactor,SBR1"),("perf,no C")),
                     output = "data/maaslin_genus")

sig_table <- fit_genus_sep$results %>% 
  mutate(map_val = -1*log(qval) * sign(coef)) %>% # how they make the heatmap in tutorial
  mutate(map_val_abs = abs(map_val)) %>%
  mutate(map_val2 = pmax(-20, pmin(20, map_val))) %>% # this is in the package R code
  mutate(map_sign = ifelse(coef > 0.0, "+", ifelse(coef < 0.0, "-", ""))) %>%
  mutate(Genus = str_replace(feature, "Ca_", "Ca. ")) %>%
  filter(pval <= 0.01)

sig_table %>%
  ggplot(., aes(x=value, y = Genus, fill=coef)) +
  geom_tile() +
  geom_text(aes(label=map_sign)) +
  scale_fill_distiller(palette = "Spectral") +
  scale_y_discrete(limits = rev) +
  theme_bw() +
  labs(x = "Variable", y = "Genus")


