# Fits latent class Bayesian model via NIMBLE on haplotype data

library(tidyverse)
library(nimble)
library(ggtern) # for Barycentric coordinate plots

set.seed(1)

# read in the data
haplo <- read_csv('../data/haplotype.csv') %>%
  drop_na()

iters <- 146094
burn <- 154
thin <- 39
# prepare data objects and constants
I <- nrow(haplo)
M <- ncol(haplo)-2
J <- 3
C <- 2
X <- array(NA, dim = c(I, M)) 
for (i in 1:I){
  for (m in 1:M){
    if (haplo[i,(m+2)] == 0){
      X[i,m] <- 1 # 1 indicates two recesive genes
    } else if (haplo[i,(m+2)] == 1){
      X[i,m] <- 2 # 2 indicates one recesive and one dominant gene 
    } else {
      X[i,m] <- 3 # 3 indicates two dominant genes 
    }
  }
}
alpha_p <- c(1,1,1)

# prepare initial values
Q0 <- matrix(c(1/3, 1/3, 1/3), nrow = I, ncol = J, byrow = TRUE)
alpha_q0 <- c(1,1,1)
P0 <- array(1/3, dim = c(J, M, 3))
Z0 <- matrix(apply(Q0, 1, \(prob) rcat(1, prob)), nrow = I, ncol = M) 

# The NIMBLE part of this function
model0_code <- nimbleCode({
 for (i in 1:I){
    Q[i,1:3] ~ ddirch(alpha_q[1:3])
    for (m in 1:M){
      Z[i,m] ~ dcat(Q[i,1:3])
      X[i,m] ~ dcat(P[Z[i,m],m,1:3])
    }
 }
 for (j in 1:J){
   alpha_q[j] ~ dgamma(2,2)
   for (m in 1:M){
      P[j,m,1:3] ~ ddirch(alpha_p[1:3])
   }
 }
})

model0_consts <- list(J = J,
                      C = C,
                      I = I,
                      M = M,
                      alpha_p = alpha_p)

model0_data <- list(X = X)

model0_inits <- list(Z=Z0,
                     Q=Q0,
                     P=P0,
                     alpha_q = alpha_q0) 

mcmc_out <- nimbleMCMC(code = model0_code, constants = model0_consts, 
                       data = model0_data, inits = model0_inits, 
                       nchains = 1, niter = iters + burn, nburnin = burn,
                       thin = thin, summary = TRUE, 
                        WAIC=TRUE, monitors = c('Q', 'P', 'Z', 'alpha_q'))

save(mcmc_out, file = '../mcmc_draws/mcmc_PPLMultiX2.Rdata')

Qs <- which(stringr::str_detect(rownames(mcmc_out$summary), 'Q'))
Qdraws <- mcmc_out$summary[Qs,]
Qdraws_means <- Qdraws[,'Mean']
Qdraws2 <- array(Qdraws_means, c(168, 3))
colnames(Qdraws2) <- c('x', 'y', 'z')
Qdraws2 <- as_tibble(Qdraws2) %>% 
  mutate(race = haplo$race)

ggtern(Qdraws2, aes(x, y, z)) + 
  geom_point(aes(color = race))
ggsave('../figs/haplo_tern_PPLMultiX2.png')  

