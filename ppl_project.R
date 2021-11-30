# Fits latent class Bayesian model via NIMBLE on haplotype data

library(tidyverse)
library(nimble)

source('preproc.R')

# set up parallel processing
nCores <- parallel::detectCores()
options(mc.cores = nCores) # use all available cores

haplo <- read_csv('haplotype.csv')

# 0. The model as found in CASI
haplo_preproc0 <- preproc(haplo, dist = 'bernoulli')
# prepare data objects and constants
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
alpha_q <- c(1,1,1)
alpha_p <- c(1,1)
# prepare initial values
set.seed(1)
Q0 <- t(replicate(n_ind, rdirch(1, alpha_q)))
Z0 <- array(dim = c(n_ind, n_loc, n_copies, n_classes))
for (i in 1:n_ind){
  for (m in 1:n_loc){
    for (c in 1:n_copies){
      Z0[i, m, c,] <- rmulti(1, 1, Q0[i,]) 
    }
  }
}
P0 <- aperm(replicate(n_classes, replicate(n_loc, rdirch(1, alpha_p))), c(3,2,1))

model0_code <- nimbleCode({
  for (i in 1:n_ind){
    Q[i,1:3] ~ ddirch(alpha_q[1:3])
    for (m in 1:n_loc){
      for (c in 1:n_copies){
        Z[i,m,c,1:3] ~ dmultinom(Q[i,1:3], 1)
        j[i,m,c] <- 1*(Z[i,m,c,1] == 1) + 2*(Z[i,m,c,2] == 1) + 3*(Z[i,m,c,3] == 1)
        X[i,m,c] ~ dbinom(P[j[i,m,c],m,1], 1)
      }
    }
  }
  for (k in 1:n_classes){
    for (m in 1:n_loc){
      P[k,m,1:2] ~ ddirch(alpha_p[1:2])
    }
  }
})

model0_consts <- list(n_classes = n_classes,
                      n_copies = n_copies,
                      n_ind = n_ind,
                      n_loc = n_loc,
                      alpha_q = alpha_q,
                      alpha_p = alpha_p)

model0_data <- list(X = X)

model0_inits <- list(Z=Z0,
                     Q=Q0,
                     P=P0)

model0 <- nimbleModel(code = model0_code, name = 'model0', 
                     constants = model0_consts, data = model0_data,
                     inits = model0_inits)

mcmc0_out <- nimbleMCMC(code = model0_code, constants = model0_consts, 
                        data = model0_data, inits = model0_inits, 
                        nchains = 1, niter = 10, summary = TRUE, 
                        WAIC=TRUE, monitors = c('Q', 'P', 'Z'))

