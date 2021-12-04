# Fits latent class Bayesian model via NIMBLE on haplotype data

library(tidyverse)
library(nimble)
library(ggtern) # for Barycentric coordinate plots

set.seed(1)

# read in the data
haplo <- read_csv('../data/haplotype.csv') %>%
  drop_na()

# 0. The model as found in CASI
iters <- 100
burn <- 10
thin <- 2
# prepare data objects and constants
I <- nrow(haplo)
M <- ncol(haplo)-2
J <- 3
C <- 2
X1 <- array(NA, dim = c(I, M)) 
X2 <- array(NA, dim = c(I, M))
for (i in 1:I){
  for (m in 1:M){
    if (haplo[i,(m+2)] == 0){
      X1[i,m] <- X2[i,m] <- 1 # 1 indicates a recesive gene
    } else if (haplo[i,(m+2)] == 1){
      u <- runif(1) # a touch of random assignment when only 1 dominant gene
      if (u >= 0.5){
        X1[i,m] <- 1 
        X2[i,m] <- 2 # 2 indicates a dominant gene
      } else {
        X1[i,m] <- 2
        X2[i,m] <- 1
      }
    } else {
        X1[i,m] <- X2[i,m] <- 2
    }
  }
}
alpha_q <- c(1,1,1)
alpha_p <- c(1,1)

# prepare initial values
set.seed(1)
Q0 <- matrix(c(1/3, 1/3, 1/3), nrow = I, ncol = J, byrow = TRUE)
P0 <- array(1/2, dim = c(J, M, 2))
Z10 <- Z20 <- matrix(apply(Q0, 1, \(prob) rcat(1, prob)), nrow = I, ncol = M) 
#Z1 for copy 1 Z2 for copy 2

# The NIMBLE part of this function
model0_code <- nimbleCode({
 for (i in 1:I){
    Q[i,1:3] ~ ddirch(alpha_q[1:3])
    for (m in 1:M){
      Z1[i,m] ~ dcat(Q[i,1:3])
      Z2[i,m] ~ dcat(Q[i,1:3])
      X1[i,m] ~ dcat(P[Z1[i,m],m,1:2])
      X2[i,m] ~ dcat(P[Z2[i,m],m,1:2])
    }
 }
 for (j in 1:J){
   for (m in 1:M){
      P[j,m,1:2] ~ ddirch(alpha_p[1:2])
   }
 }
})

model0_consts <- list(J = J,
                      C = C,
                      I = I,
                      M = M,
                      alpha_q = alpha_q,
                      alpha_p = alpha_p)

model0_data <- list(X1 = X1,
                    X2 = X2)

model0_inits <- list(Z1=Z10,
                     Z2=Z20,
                     Q=Q0,
                     P=P0) 

mcmc_out <- nimbleMCMC(code = model0_code, constants = model0_consts, 
                        data = model0_data, inits = model0_inits, 
                        nchains = 1, niter = iters + burn, nburnin = burn,
                        thin = thin, summary = TRUE, 
                        WAIC=TRUE, monitors = c('Q', 'P', 'Z1', 'Z2'))

save(mcmc_out, file = '../mcmc_draws/mcmc_PPLBinX.Rdata')

Qs <- which(stringr::str_detect(rownames(mcmc_out$summary), 'Q'))
Qdraws <- mcmc_out$summary[Qs,]
Qdraws_means <- Qdraws[,'Mean']
Qdraws2 <- array(Qdraws_means, c(168, 3))
colnames(Qdraws2) <- c('x', 'y', 'z')
Qdraws2 <- as_tibble(Qdraws2) %>% 
  mutate(race = haplo$race)

ggtern(Qdraws2, aes(x, y, z)) + 
  geom_point(aes(color = race))
ggsave('../figs/haplo_tern_PPLBinX.png')  

