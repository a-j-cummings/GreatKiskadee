# Using Raftery-Lewis diagnostic to understand convergence criteria of the MCMC

library(tidyverse)
library(coda)
library(narray)

# by hand with Binary X (CASI's model); object is 'mcmc_out'
load('mcmc_BHBinX.Rdata')
draws_P <- mcmc_out$draws_P[,,1,1,]
rl_p <- apply(draws_P, c(1,2), \(draws) raftery.diag(draws)$resmatrix)
draws_Q <- mcmc_out$draws_Q[,1:2,,]
rl_q <- apply(draws_Q, c(1,2), \(draws) raftery.diag(draws)$resmatrix)
max_p <- apply(rl_p, 1, max)
max_q <- apply(rl_q, 1, max)
rl_1_res <- apply(rbind(max_p, max_q), 2, max)
names(rl_1_res) <- c('nburn', 'ntotal', 'nmin', 'k')



# by hand with Multinomial X (my model); object is 'mcmc_out'
load('mcmc_BHMultiX.Rdata')
draws_P <- mcmc_out$draws_P[,,1:2,1,]
rl_p <- apply(draws_P, c(1,2,3), \(draws) raftery.diag(draws)$resmatrix)
draws_Q <- mcmc_out$draws_Q[,1:2,,]
rl_q <- apply(draws_Q, c(1,2), \(draws) raftery.diag(draws)$resmatrix)
max_p <- apply(rl_p, 1, max)
max_q <- apply(rl_q, 1, max)
rl_2_res <- apply(rbind(max_p, max_q), 2, max)
names(rl_2_res) <- c('nburn', 'ntotal', 'nmin', 'k')



# by hand with Multinomial X and hyperprior on concentration parameters
load('mcmc_BHMultiX2.Rdata')
draws_P <- mcmc_out$draws_P[,,1:2,1,]
rl_p <- apply(draws_P, c(1,2,3), \(draws) raftery.diag(draws)$resmatrix)
draws_Q <- mcmc_out$draws_Q[,1:2,,]
rl_q <- apply(draws_Q, c(1,2), \(draws) raftery.diag(draws)$resmatrix)
draws_alpha <- mcmc_out$draws_alpha[,1,]
rl_alpha <- apply(draws_alpha, 1, \(draws) raftery.diag(draws)$resmatrix)
max_p <- apply(rl_p, 1, max)
max_q <- apply(rl_q, 1, max)
max_alpha <- apply(rl_alpha, 1, max)
rl_2b_res <- apply(rbind(max_p, max_q, max_alpha), 2, max)
names(rl_2b_res) <- c('nburn', 'ntotal', 'nmin', 'k')


# NIMBLE with Binary X (CASI's model); object is 'mcmc_out'
load('mcmc_PPLBinX.Rdata')
P_index <- which(str_detect(colnames(mcmc_out$samples), 
                            'P\\[[:graph:]{2}[:blank:][:graph:]{1,4}, 1\\]'))
Q_index <- which(str_detect(colnames(mcmc_out$samples),
                 'Q\\[[:graph:]{1,3}, [1,2]\\]'))
draws <- mcmc_out$samples[,c(P_index, Q_index)]
rl_3 <- t(apply(draws, 2, \(col) raftery.diag(col)$resmatrix))
colnames(rl_3) <- c('nburn', 'ntotal', 'nmin', 'k')
rl_3_res <- apply(rl_3, 2, max)



# NIMBLE with Multinomial X (my model); object is 'mcmc_out'
load('mcmc_PPLMultiX.Rdata')
P_index <- which(str_detect(colnames(mcmc_out$samples), 
                            'P\\[[:graph:]{2}[:blank:][:graph:]{1,4}, [1,2]\\]'))
Q_index <- which(str_detect(colnames(mcmc_out$samples),
                 'Q\\[[:graph:]{1,3}, [1,2]\\]'))
draws <- mcmc_out$samples[,c(P_index, Q_index)]
rl_4 <- t(apply(draws, 2, \(col) raftery.diag(col)$resmatrix))
colnames(rl_4) <- c('nburn', 'ntotal', 'nmin', 'k')
rl_4_res <- apply(rl_4, 2, max)



# NIMBLE with Multinomial X and hyperprior on concentration parameters
load('mcmc_PPLMultiX2.Rdata')
P_index <- which(str_detect(colnames(mcmc_out$samples), 
                            'P\\[[:graph:]{2}[:blank:][:graph:]{1,4}, [1,2]\\]'))
Q_index <- which(str_detect(colnames(mcmc_out$samples),
                 'Q\\[[:graph:]{1,3}, [1,2]\\]'))
alpha_index <- which(str_detect(colnames(mcmc_out$samples), 'alpha'))
draws <- mcmc_out$samples[,c(P_index, Q_index, alpha_index)]
rl_4b <- t(apply(draws, 2, \(col) raftery.diag(col)$resmatrix))
colnames(rl_4b) <- c('nburn', 'ntotal', 'nmin', 'k')
rl_4b_res <- apply(rl_4b, 2, max)



# summarize it all
rl_res_all <- rbind(rl_1_res, rl_2_res, rl_2b_res,
                    rl_3_res, rl_4_res, rl_4b_res)
rownames(rl_res_all) <- c('BHBinX', 'BHMultiX', 'BHMultiX2', 
                          'PPLBinX', 'PPLMultiX', 'PPLMultiX2')
rl_res_df <- as_tibble(rl_res_all) %>%
  add_column(model = rownames(rl_res_all), .before=1) %>%
  mutate(Nthin = ceiling(k),
         Nburn = nburn,
         Ntotal = Nthin * nmin)
         
rl_res_df  %>%
   write_csv('rl_res.csv')

