# Fits latent class Bayesian model via NIMBLE on haplotype data

library(tidyverse)
library(ggtern) # for Barycentric coordinate plots

haplo <- read_csv('haplotype.csv')

load(mcmc0_out.Rdata)
all_chains <- mcmc0_out$summary$chain1

which(stringr::str_detect(rownames(all_chains), 'Q'))-> Qs
Qdraws <- all_chains[Qs,]
Qdraws_means <- Qdraws[,'Mean']
Qdraws2 <- array(Qdraws_means, c(168, 3))
#cbind(Qdraws[1:168,1], Qdraws[169:336,1], Qdraws[337:504,1])
colnames(Qdraws2) <- c('x', 'y', 'z')
Qdraws2 <- as_tibble(Qdraws2) %>% 
  mutate(race = drop_na(haplo)$race)


ggtern(Qdraws2, aes(x, y, z)) + 
  geom_point(aes(color = race))
ggsave('haplo_tern.png')


# P
which(stringr::str_detect(rownames(all_chains), 'P')) -> Ps
Pdraws <- all_chains[Ps,]
Pdraws_means <- Pdraws[,'Mean']
temp <- array(Pdraws_means, dim = c(3, 97, 2))

# Z
which(stringr::str_detect(rownames(all_chains), 'Z')) -> Zs
Zdraws <- all_chains[Zs,]
Zdraws_means <- Zdraws[,'Mean']
temp <- array(Zdraws_means, c(168, 97, 2, 3))

