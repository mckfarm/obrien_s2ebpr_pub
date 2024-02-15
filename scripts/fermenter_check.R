

library(tidyverse)
library(ggpubr)
library(vegan)
library(MetBrewer)

source("scripts/plotting_cols_shapes.R")
source("scripts/plotting_dates.R")

reactors <- c("SBR1", "SBR2", "SBR3")

rel_genus <- readRDS("data/rel_genus.RDS") %>%
  mutate(genus_format = str_replace(Genus, "Ca_", "Ca "))

metadata <- read_delim("data/metadata.txt") %>%
  mutate(date = mdy(date)) %>%
  mutate(perf = case_when(
    reactor == "SBR1" ~ "C off", 
    (reactor != "SBR1" & date < ymd("2022-10-26")) ~ "C off",
    (reactor != "SBR1" & between(date, ymd("2022-10-25"), ymd("2023-05-12"))) ~ "C on",
    (reactor != "SBR1" & date > ymd("2023-05-12")) ~ "C off"
  ))

fermenters <- read_csv("/projects/b1052/mckenna/resources/qiime_v5.1/fermenters_5.1.csv")

rel_ferm <- rel_genus %>%
  filter(genus_format %in% fermenters$`Canonical Name`) %>%
  filter(!genus_format %in% c("Ca Competibacter", "Ca Obscuribacter"))

my_comparisons <- list(c("SBR1", "SBR2"),
                       c("SBR2", "SBR3"),
                       c("SBR1", "SBR3"))

rel_ferm %>% 
  group_by(reactor, date) %>%
  summarise(tot = sum(sum)) %>%
  ggplot(data = ., aes(x = reactor, y = tot, color = reactor)) +
  geom_boxplot(outlier.color = NULL) +
  geom_point(position = position_jitter(width = 0.2)) + 
  scale_color_reactor + 
  labs(x = "", y = "Relative abundance [%]") +
  stat_compare_means(method = "wilcox", comparisons = my_comparisons) +
  theme_black_box +
  theme(legend.position = "none")


for(i in reactors){
  rel_ferm_plt <- 
    rel_ferm %>% 
    filter(reactor == i) %>%
    group_by(genus_format) %>%
    summarise(median = median(sum)) %>%
    filter(median > 0.01) %>%
    ggplot(data = ., aes(x = reorder(genus_format, -median), y = median)) +
    geom_bar(stat = "identity") + 
    labs(x = "", y = "Relative abundance [%]", title = paste(i)) +
    theme_black_box +
    theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
  
  print(rel_ferm_plt)
}

# prep for diversity and evenness calcs
rel_ferm_wide <- rel_ferm %>%
  ungroup() %>%
  select(sample, Genus, sum) %>%
  pivot_wider(names_from = sample, values_from = sum) %>%
  column_to_rownames("Genus")

rel_ferm_wide <- as.data.frame(apply(rel_ferm_wide, 2, function(x){x / sum(x)})) # normalize by sample

# diversity
div <- as.data.frame(diversity(rel_ferm_wide, index = "shannon", MARGIN = 2))
colnames(div) <- "shannon"

div <- div %>% 
  rownames_to_column(var = "sample") %>%
  left_join(metadata, by = "sample")
  
ggplot(div, aes(x = reactor, y= shannon)) +
  geom_boxplot(outlier.color = NULL) +
  geom_point(position = position_jitter(width = 0.1)) +
  stat_compare_means(method = "wilcox", comparisons = my_comparisons) + 
  theme_black_box 


# evenness
S <- apply( rel_ferm_wide > 0, 2 , sum )
even <- as.data.frame(diversity(rel_ferm_wide, index = "shannon", MARGIN = 2) / log(S))
colnames(even) <- "shannon"

even <- even %>% 
  rownames_to_column(var = "sample") %>%
  left_join(metadata, by = "sample")

ggplot(even, aes(x = reactor, y= shannon)) +
  geom_boxplot(outlier.color = NULL) +
  geom_point(position = position_jitter(width = 0.1)) +
  stat_compare_means(method = "wilcox", comparisons = my_comparisons) + 
  theme_black_box 

# richness
richness <- as.data.frame(apply( rel_ferm_wide > 0, 2 , sum ))
colnames(richness) <- "richness"

richness <- richness %>% 
  rownames_to_column(var = "sample") %>%
  left_join(metadata, by = "sample")

richness %>% 
  group_by(reactor) %>%
  summarise(median(richness))

ggplot(richness, aes(x = reactor, y = richness, color = reactor)) +
  geom_boxplot(outlier.color = NULL) +
  geom_point(position = position_jitter(width = 0.1)) +
  stat_compare_means(method = "wilcox", comparisons = my_comparisons) + 
  scale_color_reactor + 
  labs(x = "Reactor", y = "Number of genera with fermentative metabolisms") + 
  theme_black_box +
  theme(legend.position = "none")
