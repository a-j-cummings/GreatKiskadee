library(tidyverse)
library(gtools) # for rdirichlet
library(ggtern) # for Barycentric coordinate plots

set.seed(1)

haplo <- read_csv('../data/haplotype.csv') %>% 
  drop_na()

# data preprocessing
I <- nrow(haplo)
M <- ncol(haplo) - 2
J = 3
C = 2
k_x <- J
k_z <- J

X <- array(NA, dim=c(I, M, k_x))
for (i in 1:I){
  for (m in 1:M){
    if (haplo[i,(m+2)] == 0){
      X[i,m,1:3] <- c(1,0,0)
    } else if (haplo[i,(m+2)] == 1){
      X[i,m,1:3] <- c(0,1,0)
    } else {
      X[i,m,1:3] <- c(0,0,1)
    }
  }
}

# construct updating functions

update_Z <- function(){
  Z <- array(NA, dim = c(I, M, k_z))
  for (i in 1:I){
    for (m in 1:M){
      pz <- c(Q[i,1]*dmultinom(X[i,m,], prob = P[1,m,]),
              Q[i,2]*dmultinom(X[i,m,], prob = P[2,m,]),
              Q[i,3]*dmultinom(X[i,m,], prob = P[3,m,]))
      Z[i,m,] <- rmultinom(1, 1, pz)
    }
  }
  return(Z)
}

update_P <- function(lambda = c(1, 1, 1)){
  P <- array(NA, dim = c(J, M, k_x))
  for (j in 1:J){
    for (m in 1:M){
      njm0 <- sum((X[,m,1] == 1)*(Z[,m,j] == 1))
      njm1 <- sum((X[,m,2] == 1)*(Z[,m,j] == 1))
      njm2 <- sum((X[,m,3] == 1)*(Z[,m,j] == 1))
      P[j,m,] <- rdirichlet(1, lambda + c(njm0, njm1, njm2))
    }
  }
  return(P)
}

update_Q <- function(gamma = c(1, 1, 1)){
  Q <- array(NA, dim = c(I, k_z))
  for (i in 1:I){
   mi1 <- sum(Z[i,,1] == 1)
   mi2 <- sum(Z[i,,2] == 1)
   mi3 <- sum(Z[i,,3] == 1)
   Q[i,] <- rdirichlet(1, gamma + c(mi1, mi2, mi3))
  }
  return(Q)
}

# Initial values
Q <- matrix(c(1/3, 1/3, 1/3), nrow = I, ncol = J, byrow = TRUE)
P <- array(1/2, dim = c(J, M, k_x))
# init for Z not needed (will be sampled in first step)

# run the sampler
nchains <- 1
nburn <- 40
niters <- 67428
nthin <- 18
draws_Z <- array(NA, dim = c(I, M, C, k_z, nchains, niters/nthin))
draws_P <- array(NA, dim = c(J, M, k_x, nchains, niters/nthin))
draws_Q <- array(NA, dim = c(I, J, nchains, niters/nthin))
for (chain in 1:nchains){
  print(paste0('Running chain ', chain))
  for (burn in 1:nburn){
    Z <- update_Z()
    P <- update_P()
    Q <- update_Q()  
  }
  pb = txtProgressBar(min = 0, max = niters, initial = 0)
  for (iter in 1:niters){
    setTxtProgressBar(pb,iter)
    Z <- update_Z()
    P <- update_P()
    Q <- update_Q()
    if (iter %% nthin == 0){
      loc <- iter/nthin
      draws_Z[,,,,chain,loc] <- Z
      draws_P[,,,chain,loc] <- P
      draws_Q[,,chain,loc] <- Q 
    }
  }
  close(pb)
}

mcmc_out <- list(draws_Z=draws_Z, draws_P=draws_P, draws_Q=draws_Q)
save(mcmc_out, file = '../mcmc_draws/mcmc_BHMultiX.Rdata')


post_means <- apply(draws_Q[,,,niters/nthin], c(1, 2), mean) %>%
  as_tibble() %>%
  rename(x = V1, y = V2, z = V3) %>%
  mutate(race = haplo$race)

ggtern(post_means, aes(x, y, z)) + 
  geom_point(aes(color = race))
ggsave('../figs/haplo_tern_BHMultiX.png')
