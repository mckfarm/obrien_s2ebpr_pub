library(tidyverse)

table_rel_long <- readRDS("./data/table_rel_long.RDS")
midas_taxonomy <- read_delim("~/b1052/mckenna/resources/qiime_v5.1/qiime_5.1.txt", 
                             col_names = c("kingdom", "phylum", "class", 
                                           "order", "family", "genus", "species")) %>%
  separate(col = kingdom, into = c("id", "kingdom"), sep = "\\s") %>%
  mutate_at(vars(kingdom:species), ~ str_replace(., ";", "")) %>%
  mutate_at(vars(kingdom:species), ~gsub(x = ., "[A-z]__", ""))




blast <- read_delim("./data/competibacter_blast.txt", 
                            col_names = c("query", "subject", "pident", "length", 
                                          "mismatch", "gapopen", "qstart", "qend",
                                          "sstart", "send", "evalue", "bitscore")) %>%
  filter(pident >= 97) %>%
  left_join(midas_taxonomy, by = c("subject" = "id"))
