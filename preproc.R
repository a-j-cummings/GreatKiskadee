# functions for performing preprocessing on haplotype data

preproc <- function(raw_data, 
                    dist = 'bernoulli',
                    rand_seed = 258, 
                    na_action = 'drop'){
  "Performs preprocessing on haplotype data for latent class model fitting
  
  Keyword Arguments:
  ------------------
  raw_data -- the raw 'Human ancestry' data found at 
    https://web.stanford.edu/~hastie/CASI/data.html
  dist -- string, one of 'bernoulli' or 'multinomial'; indicates how to treat 
    encode the allele levels. CASI endorses a 'bernoulli' fit in the discussion
    on page 258; default is 'bernoulli'
  rand_seed -- numeric, the random seed that is used to do preprocessing when
    dist == 'bernoulli', default is 258
  na_action -- one of 'drop' or 'handle'; indicates how to treat missing values
    SNP measurements in the data
  "
  require('tidyverse')
  if (na_action == 'drop') {
    raw_data <- raw_data %>% 
      drop_na()
  }
  snp <- raw_data %>% 
    select(-id, -race)
  n_ind <- nrow(snp) # the number of individuals represented in the data
  if (dist == 'bernoulli') {
    n_loc <- ncol(snp) # the number of genome locations represented in the data
    copies <- array(0, dim = c(n_ind, n_loc, 2))
    # a helper function for spreading the allele levels
    helper1 <- function(x){
      x <- as.character(x)
      switch(x,
             '0' = c(0,0),
             '1' = sample(0:1, 2),
             '2' = c(1,1))
    }
    # double for loop for readability, not efficiency
    for (ind in 1:n_ind) {
      for (loc in 1:n_loc) {
        x <- snp[ind, loc]
        copies[ind, loc, ] <- helper1(x)
      }
    }
    cols <- colnames(snp)
    snp <- rbind(copies[,,1], copies[,,2])
    colnames(snp) <- cols
    snp <- snp %>% 
      as_tibble %>% 
      add_column(ind = rep(1:n_ind, 2), .before = cols[1]) %>% 
      arrange(ind) %>% 
      add_column(copy = rep(1:2, n_ind), .after = 'ind')
    id <- snp %>% 
      select(ind, copy)
    snp <- snp %>% 
      select(-ind, -copy) %>% 
      mutate_all(factor)
    snp <- snp %>% 
      .[,apply(snp, 2, \(col) length(unique(col))) > 1]
  } else if (dist == 'multinomial') {
    snp <- snp %>% 
      mutate_all(factor)
    snp <- snp %>% 
      .[,apply(snp, 2, \(col) length(unique(col))) > 1]
    id <- tibble(ind = 1:n_ind)
  }
  ret <- list(preproc_data = snp, id = id)
}



