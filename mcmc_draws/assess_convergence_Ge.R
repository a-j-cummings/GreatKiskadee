# Using Raftery-Lewis diagnostic to understand convergence criteria of the MCMC

library(tidyverse)
library(cowplot)
library(coda)

# by hand with Binary X (CASI's model); object is 'mcmc_out'
load('mcmc_BHBinX.Rdata')
draws_P <- mcmc_out$draws_P[,,1,1,]
ge_p <- apply(draws_P, c(1,2), \(draws) geweke.diag(draws)$z)
draws_Q <- mcmc_out$draws_Q[,1:2,,]
ge_q <- apply(draws_Q, c(1,2), \(draws) geweke.diag(draws)$z)
ge_1 <- c(ge_p, ge_q)


# by hand with Multinomial X (my model); object is 'mcmc_out'
load('mcmc_BHMultiX.Rdata')
draws_P <- mcmc_out$draws_P[,,1:2,1,]
ge_p <- apply(draws_P, c(1,2,3), \(draws) geweke.diag(draws)$z)
draws_Q <- mcmc_out$draws_Q[,1:2,,]
ge_q <- apply(draws_Q, c(1,2), \(draws) geweke.diag(draws)$z)
ge_2 <- c(ge_p, ge_q)


# by hand with Multinomial X and hyperpriors
load('mcmc_BHMultiX2.Rdata')
draws_P <- mcmc_out$draws_P[,,1:2,1,]
ge_p <- apply(draws_P, c(1,2,3), \(draws) geweke.diag(draws)$z)
draws_Q <- mcmc_out$draws_Q[,1:2,,]
ge_q <- apply(draws_Q, c(1,2), \(draws) geweke.diag(draws)$z)
draws_alpha <- mcmc_out$draws_alpha[,1,]
ge_alpha <- apply(draws_alpha, 1, \(draws) geweke.diag(draws)$z)
ge_2b <- c(ge_p, ge_q, ge_alpha)


# NIMBLE with Binary X (CASI's model); object is 'mcmc_out'
load('mcmc_PPLBinX.Rdata')
P_index <- which(str_detect(colnames(mcmc_out$samples), 
                            'P\\[[:graph:]{2}[:blank:][:graph:]{1,4}, 1\\]'))
Q_index <- which(str_detect(colnames(mcmc_out$samples),
                 'Q\\[[:graph:]{1,3}, [1,2]\\]'))
draws <- mcmc_out$samples[,c(P_index, Q_index)]
ge_3 <- geweke.diag(draws)$z


# NIMBLE with Multinomial X (my model); object is 'mcmc_out'
load('mcmc_PPLMultiX.Rdata')
P_index <- which(str_detect(colnames(mcmc_out$samples), 
                            'P\\[[:graph:]{2}[:blank:][:graph:]{1,4}, [1,2]\\]'))
Q_index <- which(str_detect(colnames(mcmc_out$samples),
                 'Q\\[[:graph:]{1,3}, [1,2]\\]'))
draws <- mcmc_out$samples[,c(P_index, Q_index)]
ge_4 <- geweke.diag(draws)$z


# NIMBLE with Multinomial X and hyperpriors
load('mcmc_PPLMultiX2.Rdata')
P_index <- which(str_detect(colnames(mcmc_out$samples), 
                            'P\\[[:graph:]{2}[:blank:][:graph:]{1,4}, [1,2]\\]'))
Q_index <- which(str_detect(colnames(mcmc_out$samples),
                 'Q\\[[:graph:]{1,3}, [1,2]\\]'))
alpha_index <- which(str_detect(colnames(mcmc_out$samples), 'alpha'))
draws <- mcmc_out$samples[,c(P_index, Q_index, alpha_index)]
ge_4b <- geweke.diag(draws)$z


# save to output
out <- list(ge_1 = ge_1, ge_2 = ge_2, ge_2b = ge_2b,
            ge_3 = ge_3, ge_4 = ge_4, ge_4b = ge_4b)
save(out, file = 'ge_diagnostics.Rdata')

