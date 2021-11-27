# Function for fitting a latent class model via the EM algorithm

LatentClassEM <- function(raw_data, 
                          preproc,
                          n_latent_classes,
                          rand_seed = 637,
                          tol = 1e-1){
  "Fits a Latent Class model via the EM algorithm.
  
  Keyword Arguments:
  ------------------
  raw_data -- data frame of raw data
  preproc -- function that performs any necessary preprocessing; takes only the
    raw data as an input, returns (1) 'preproc_data', a tibble/data.frame, the 
    preprocessed manifest variables left unsummarised (2) 'id', a 
    tibble/data.frame containing column(s) that will be used to generate latent 
    class estimates, the rows should correspond with the rows of 'preproc_data'
  n_latent_classes -- numeric integer, the number of latent classes to be fit in
    the model
  rand_seed -- numeric integer, the random seed to fit model with; default = 637
  tol -- numeric, tolerance threshhold for determining EM algorithm convergence
  "
  require('tidyverse')
  set.seed(rand_seed)
  
  preproc_out <- preproc(raw_data)
  preproc_data <- preproc_out$preproc_data
  id <- preproc_out$id
  manifest_variables <- colnames(preproc_data)
  # get joint counts
  joint_counts <- preproc_data %>% 
    group_by_all() %>% 
    summarise(n = n(), .groups = 'drop')
  # initialize marginal counts (randomly)
  helper1 <- function(n, n_latent_classes){
    n_in_class <- numeric(n_latent_classes)
    for (k in 1:(n_latent_classes-1)){
      n_in_class[k] <- sample(0:n, 1)
      n <- n - n_in_class[k]
    }
    n_in_class[k+1] <- n
    return(n_in_class)
  }
  n_in_class <- lapply(joint_counts$n, \(n) helper1(n, n_latent_classes))
  marginal_counts <- joint_counts %>% 
    mutate(n_s = n_in_class) %>% 
    unnest(n_s) %>% 
    mutate(group = factor(rep(1:n_latent_classes, nrow(preproc_data))))
  # fit log linear model
  formula_for_log_linear_model <- formula(paste0('n_s ~ (',
     paste0(manifest_variables, collapse = " + "), ')*group' ))
  mod <- glm(formula_for_log_linear_model, 
             data = marginal_counts, 
             family = poisson)
  # obtain marginal expectations and update marginal counts
  mu_s <- predict(mod, type = 'response')
  marginal_counts <- marginal_counts %>% 
    mutate(mu_s = mu_s)
  marginal_counts <- marginal_counts %>% 
    left_join(marginal_counts %>% 
                group_by_at(manifest_variables) %>% 
                summarise(mu = sum(mu_s), .groups = 'drop'), 
              by = manifest_variables) %>% 
    mutate(n_s = n*mu_s/mu)
  # conduct EM algorithm
  iter <- 0
  diff <- 1
  options(warn = -1)
  while(diff > tol){
    iter <- iter + 1
    prev_mu_s <- mu_s
    # M
    mod <- glm(formula_for_log_linear_model, 
               data = marginal_counts, 
               family = poisson)
    mu_s <- predict(mod, type = 'response')
    marginal_counts <- marginal_counts %>%
      select(-n_s, -mu_s, -mu) %>% 
      mutate(mu_s = mu_s)    
    # E
    marginal_counts <- marginal_counts %>% 
      left_join(marginal_counts %>% 
                  group_by_at(manifest_variables) %>% 
                  summarise(mu = sum(mu_s), .groups = 'drop'), 
                by = manifest_variables) %>% 
      mutate(n_s = n*mu_s/mu)
    diff <- sqrt(sum((prev_mu_s - mu_s)^2))
    print(paste('iteration', iter, ':', round(diff, 4)))
  }
  options(warn = 0)  
}


LatentClassEM(haplo, preproc, 3)

