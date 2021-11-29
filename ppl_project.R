# Fits latent class Bayesian model via Stan on haplotype data

library(tidyverse)
library(rstan)
library(R2jags)

source('preproc.R')

# set up parallel processing
nCores <- parallel::detectCores()
options(mc.cores = nCores) # use all available cores
rstan_options(auto_write = TRUE) # cache compiled code

haplo <- read_csv('haplotype.csv')

# 0. The model as found in CASI
haplo_preproc0 <- preproc(haplo, dist = 'bernoulli')
# prepare data objects
n_classes <- 3
snp <- haplo_preproc0$preproc_data
n_copies <- 2
n_ind <- nrow(snp)/n_copies
n_loc <- ncol(snp)
copy1 <- snp[rep(c(TRUE, FALSE), n_ind),] %>% 
  mutate_all(\(x) as.numeric(x) - 1) %>% 
  as.matrix()
copy2 <- snp[rep(c(FALSE, TRUE), n_ind),] %>% 
  mutate_all(\(x) as.numeric(x) - 1) %>% 
  as.matrix()
X <- array(dim = c(n_ind, n_loc, n_copies))
X[,,1] <- copy1
X[,,2] <- copy2
ind <- unique(haplo_preproc0$id$ind)
copy <- haplo_preproc0$id$copy
alpha_q <- c(1,1,1)
alpha_p <- c(1,1)
data0 <- c('n_classes', 'n_copies', 'n_ind', 'n_loc', 'X', 'alpha_q', 'alpha_p')
model0 <- "
  model {
    for (i in 1:n_ind){
      Q[i,] ~ ddirich(alpha_q)
      for (m in 1:n_loc){
        for (c in 1:n_copies){
          Z[i,m,c,] ~ dmulti(Q[i,], 1)
          for (j in 1:n_classes){
            if (Z[i,m,c,j] == 1){
              X[i,m,c,] ~ dmulti(P[j,m,], 1)
            }
          }
        }
      }
    }
    for (j in 1:n_classes){
      for (m in 1:n_loc){
        P[j,m,] ~ ddirich(alpha_p)
      }
    }
  }
"
params0 <- c('Q', 'Z', 'P')
writeLines(model0, 'model0.txt')

sim0 <- jags(data=data0, inits=NULL,
             parameters.to.save=params0,
             model.file='model0.txt',
             n.iter=10, n.burnin=0,
             n.chains=1, n.thin=1)







data0 <- list(n_classes = n_classes,
              X = X, 
              n_copies = n_copies,
              n_ind = n_ind,
              n_loc = n_loc,
              ind = ind,
              copy = copy,
              alpha_q = alpha_q,
              alpha_p = alpha_p)
fit0 <- stan(model_code = readLines('model0.stan'),
             data = data0,
             iter = 10,
             warmup = 0,
             thin = 1,
             chains = nCores)





# 1. The model with genotype as a multinomial

haplo_preproc1 <- preproc(haplo, dist = 'multinomial')

