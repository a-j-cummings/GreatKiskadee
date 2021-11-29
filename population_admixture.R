library(tidyverse)
library(gtools) # for rdirichlet
library(ggtern) # for Barycentric coordinate plots

set.seed(1)

haplo <- read_csv('haplotype.csv') %>% 
  rename(id = X1) %>% 
  drop_na()

# data preprocessing
I <- nrow(haplo)
M <- ncol(haplo) - 2
J = 3
C = 2
k_x <- 2
k_z <- J

X <- array(NA, dim = c(I, M, C, k_x)) 

for (i in 1:I){
  for (m in 1:M){
    if (haplo[i,(m+2)] == 0){
      X[i,m,1,] <- c(1, 0)
      X[i,m,2,] <- c(1, 0)
    } else if (haplo[i,(m+2)] == 1){
      u <- runif(1) # a touch of random assignment when only 1 dominant gene
      if (u >= 0.5){
        X[i,m,1,] <- c(0, 1)
        X[i,m,2,] <- c(1, 0)
      } else {
        X[i,m,1,] <- c(1, 0)
        X[i,m,2,] <- c(0, 1)
      }
    } else {
        X[i,m,1,] <- c(0, 1)
        X[i,m,2,] <- c(0, 1)
    }
  }
}


# construct updating functions

update_Z <- function(){
  Z <- array(NA, dim = c(I, M, C, k_z))
  for (i in 1:I){
    for (m in 1:M){
      pz1 <- c(Q[i,1]*dmultinom(X[i,m,1,], prob = P[1,m,]),
               Q[i,2]*dmultinom(X[i,m,1,], prob = P[2,m,]),
               Q[i,3]*dmultinom(X[i,m,1,], prob = P[3,m,]))
      Z[i,m,1,] <- rmultinom(1, 1, pz1)
      pz2 <- c(Q[i,1]*dmultinom(X[i,m,2,], prob = P[1,m,]),
               Q[i,2]*dmultinom(X[i,m,2,], prob = P[2,m,]),
               Q[i,3]*dmultinom(X[i,m,2,], prob = P[3,m,]))
      Z[i,m,2,] <- rmultinom(1, 1, pz2)
    }
  }
  return(Z)
}

update_P <- function(lambda = c(1, 1)){
  P <- array(NA, dim = c(J, M, k_x))
  for (j in 1:J){
    for (m in 1:M){
      njm0 <- sum((X[,m,,1] == 1)*(Z[,m,,j] == 1))
      njm1 <- sum((X[,m,,2] == 1)*(Z[,m,,j] == 1))
      P[j,m,] <- rdirichlet(1, lambda + c(njm0, njm1))
    }
  }
  return(P)
}

update_Q <- function(gamma = c(1, 1, 1)){
  Q <- array(NA, dim = c(I, k_z))
  for (i in 1:I){
   mi1 <- sum(Z[i,,,1] == 1)
   mi2 <- sum(Z[i,,,2] == 1)
   mi3 <- sum(Z[i,,,3] == 1)
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
nburn <- 1000
nkeep <- 2000
niters <- nkeep + nburn
draws_Q <- array(NA, dim = c(I, J, nchains, niters))
for (chain in 1:nchains){
  for (iter in 1:niters){
    Z <- update_Z()
    P <- update_P()
    Q <- update_Q()
    draws_Q[,,chain,iter] <- Q
  }
}


post_means <- apply(draws_Q[,,,(nburn+1):niters], c(1, 2), mean) %>%
  as_tibble() %>%
  rename(x = V1, y = V2, z = V3) %>%
  mutate(race = haplo$race)

ggtern(post_means, aes(x, y, z)) + 
  geom_point(aes(color = race))
ggsave('haplo_tern.png')



