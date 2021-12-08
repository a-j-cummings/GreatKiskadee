# 

library(tidyverse)

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


