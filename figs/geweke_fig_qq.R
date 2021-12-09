# creates figure for Geweke diagnostics

library(tidyverse)
library(cowplot)

load('../data/ge_diagnostics.Rdata')


# By hand

mod1 <- as.data.frame(out$ge_1)
colnames(mod1) <- 'geweke'
mod1 <- ggplot(mod1, aes(sample = geweke)) + 
  stat_qq() +
  stat_qq_line() +
  labs(y = '', 
       x = '')

mod2 <- as.data.frame(out$ge_2)
colnames(mod2) <- 'geweke'
mod2 <- ggplot(mod2, aes(sample = geweke)) + 
  stat_qq() +
  stat_qq_line()  +
  labs(y = '', 
       x = '')

mod2b <- as.data.frame(out$ge_2b)
colnames(mod2b) <- 'geweke'
mod2b <- ggplot(mod2b, aes(sample = geweke)) + 
  stat_qq() +
  stat_qq_line()   +
  labs(y = '', 
       x = '')

# NIMBLE

mod3 <- as.data.frame(out$ge_3)
colnames(mod3) <- 'geweke'
mod3 <- ggplot(mod3, aes(sample = geweke)) + 
  stat_qq() +
  stat_qq_line()  +
  labs(y = 'Quantiles of Geweke Diagnostics', 
       x = 'Theoretical Quantiles')

mod4 <- as.data.frame(out$ge_4)
colnames(mod4) <- 'geweke'
mod4 <- ggplot(mod4, aes(sample = geweke)) + 
  stat_qq() +
  stat_qq_line()  +
  labs(y = '', 
       x = '')

mod4b <- as.data.frame(out$ge_4b)
colnames(mod4b) <- 'geweke'
mod4b <- ggplot(mod4b, aes(sample = geweke)) + 
  stat_qq() +
  stat_qq_line() +
  labs(y = '', 
       x = '')

# make plot
plot_grid(mod1, mod2, mod2b, mod3, mod4, mod4b, nrow = 2)
