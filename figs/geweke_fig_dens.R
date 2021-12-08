# creates figure for Geweke diagnostics

library(tidyverse)
library(cowplot)

load('../data/ge_diagnostics.Rdata')


# By hand

mod1 <- as.data.frame(out$ge_1)
colnames(mod1) <- 'geweke'
mod1 <- ggplot(mod1) +
  geom_histogram(aes(geweke, y=..density..), bins = 50)  +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = 'red') +
  labs(x='', y='')

mod2 <- as.data.frame(out$ge_2)
colnames(mod2) <- 'geweke'
mod2 <- ggplot(mod2) +
  geom_histogram(aes(geweke, y=..density..), bins = 50)  +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = 'red') +
  labs(x='', y='')

mod2b <- as.data.frame(out$ge_2b)
colnames(mod2b) <- 'geweke'
mod2b <- ggplot(mod2b) +
  geom_histogram(aes(geweke, y=..density..), bins = 50)  +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = 'red') +
  labs(x='', y='')

# NIMBLE

mod3 <- as.data.frame(out$ge_3)
colnames(mod3) <- 'geweke'
mod3 <- ggplot(mod3) +
  geom_histogram(aes(geweke, y=..density..), bins = 50)  +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = 'red') +
  labs(x='', y='')

mod4 <- as.data.frame(out$ge_4)
colnames(mod4) <- 'geweke'
mod4 <- ggplot(mod4) +
  geom_histogram(aes(geweke, y=..density..), bins = 50)  +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = 'red') +
  labs(x='Geweke Diagnostics', y='')

mod4b <- as.data.frame(out$ge_4b)
colnames(mod4b) <- 'geweke'
mod4b <- ggplot(mod4b) +
  geom_histogram(aes(geweke, y=..density..), bins = 50)  +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = 'red') +
  labs(x='', y='')

# make plot
plot_grid(mod1, mod2, mod2b, mod3, mod4, mod4b, nrow = 2,
          label_size = 10, label_y = 1.01, label_x = 0,
          labels = c('Model (1) by hand', 'Model (2) by hand', 'Model (3) by hand',
                     'Model (1) NIMBLE', 'Model (2) NIMBLE', 'Model (3) NIMBLE'))
ggsave('geweke_diags.png', width = 11, height = 5)
