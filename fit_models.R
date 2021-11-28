library(tidyverse)

source('latent_class_em.R')
source('preproc.R')

haplo <- read_csv('haplotype.csv')

# debugging code
raw_data <- haplo
n_latent_classes <- 3
tol <- 1e-1


temp2 <- haplo %>% 
  select(-id) %>% 
  drop_na() %>% 
  group_by(race) %>% 
  summarise_all(mean) %>%
  select(-race) %>% 
  t()


temp <- marginal_counts %>% 
  mutate(obs_id = rep(1:nrow(id), each = n_latent_classes)) %>% 
  select(obs_id, group, n_s) %>% 
  pivot_wider(names_from = group, values_from = n_s) %>% 
  mutate(race = rep(drop_na(haplo)$race, each  = 2))
ggtern(temp, aes(group1, group2, group3)) + 
  geom_point(aes(color = race))
