# Using Raftery-Lewis diagnostic to understand convergence criteria of the MCMC

library(tidyverse)
library(cowplot)
library(coda)

# by hand with Binary X (CASI's model); object is 'mcmc_out'
load('mcmc_BHBinX.Rdata')

draws_P <- mcmc_out$draws_P
draws_Q <- mcmc_out$draws_Q

# by hand with Multinomial X (my model); object is 'mcmc_out'
load('mcmc_BHMultiX.Rdata')

draws_P <- mcmc_out$draws_P
draws_Q <- mcmc_out$draws_Q


# NIMBLE with Binary X (CASI's model); object is 'mcmc_out'
load('mcmc_PPLBinX.Rdata')
P_index <- which(str_detect(colnames(mcmc_out$samples), 
                            'P\\[[:graph:]{2}[:blank:][:graph:]{1,4}, 1\\]'))
Q_index <- which(str_detect(colnames(mcmc_out$samples),
                 'Q\\[[:graph:]{1,3}, [1,2]\\]'))
draws <- mcmc_out$samples[,c(P_index, Q_index)]
ge_3 <- geweke.diag(draws)
ge_3_plot <- ggplot() + 
  geom_histogram(aes(ge_3$z, y = ..density..)) + 
  stat_function(fun = dnorm, args = list(mean = mean(ge_3$z), sd = sd(ge_3$z))) +
  labs(x = 'Geweke Diagnostics', y = 'Density', 
       main = 'Geweke Diagnostics of MCMC chains') + 
  theme_minimal()


# NIMBLE with Multinomial X (my model); object is 'mcmc_out'
load('mcmc_PPLMultiX.Rdata')
P_index <- which(str_detect(colnames(mcmc_out$samples), 
                            'P\\[[:graph:]{2}[:blank:][:graph:]{1,4}, [1,2]\\]'))
Q_index <- which(str_detect(colnames(mcmc_out$samples),
                 'Q\\[[:graph:]{1,3}, [1,2]\\]'))
draws <- mcmc_out$samples[,c(P_index, Q_index)]
ge_4 <- geweke.diag(draws)
ge_4_plot <- ggplot() + 
  geom_histogram(aes(ge_4$z, y = ..density..)) + 
  stat_function(fun = dnorm, args = list(mean = mean(ge_4$z), sd = sd(ge_4$z))) +
  labs(x = 'Geweke Diagnostics', y = 'Density', 
       main = 'Geweke Diagnostics of MCMC chains') + 
  theme_minimal()

plot_grid(ge_3_plot, ge_4_plot)
ggsave('geweke_plot.png')
