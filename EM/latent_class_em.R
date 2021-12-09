# Function for fitting a latent class model via the EM algorithm

library(tidyverse)

source('preproc.R')

haplo <- read_csv('../data/haplotype.csv') %>%
  drop_na()

variances <- haplo %>%
  select(-id, -race, -Snp3, -Snp8, -Snp14) %>% # these SNP loci have no variance
  as.matrix() %>%
  apply(2, var) %>%
  sort(decreasing = TRUE)
snp <- variances %>%
  names() %>%
  str_extract('[:digit:]{1,3}') %>%
  as.numeric()



LatentClassEM <- function(raw_data, 
                          preproc,
                          n_vars,
                          n_latent_classes = 3,
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
    and the first column of 'id' should uniquely identify distinct combinations
    of the manifest variables
  n_latent_classes -- numeric integer, the number of latent classes to be fit in
    the model
  rand_seed -- numeric integer, the random seed to fit model with; default = 637
  tol -- numeric, tolerance threshhold for determining EM algorithm convergence
  "
  require('tidyverse')
  require('janitor')
  require('ggtern')
  set.seed(rand_seed)
  
  #rand_choice <- c(10, 12, 18, 20, 25, 26, 35, 38, 39, 40, 56, 60, 65, 67, 86, 94, 97)
  #more_careful <- c(94, 67, # A and AA
  #                  60, 25, # E and A
  #                  58, 99, # J and A
  #                  20, 86, # E and AA
  #                  65, 40, # J and AA
  #                  10, 56) # E and J 
  preproc_out <- preproc(raw_data)
  cols <- paste0('Snp', snp[1:n_vars])
  preproc_data <- preproc_out$preproc_data[,cols]#1:20]#sample(1:97, 20)]
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
    mutate(n_s = n_in_class) %>% #rep(n_in_class, each = 2)) %>% 
    unnest(n_s) %>% 
    mutate(group = factor(paste0('group', 
                                 rep(1:n_latent_classes, nrow(joint_counts)))))
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
    #print(paste('iteration', iter, ':', round(diff, 4)))
  }
  options(warn = 0)
  # get joint probability of manifest and latent variables
  #id_cols <- colnames(id)
  fitted <- marginal_counts %>% 
    mutate(pr_joint = n_s/sum(n)*n_latent_classes) %>% 
    select(-n, -mu_s, -mu, -n_s)
  # get the marginal probability of manifest variables
  marg_probs_manifest <- fitted %>% 
    pivot_longer(manifest_variables, names_to = 'manifest_variable', 
                 values_to = 'manifest_value') %>% 
    group_by(manifest_variable, manifest_value) %>% 
    summarise(pr_marginal = sum(pr_joint), .groups = 'drop') 
  marg_probs_manifest_vec <- fitted %>% 
    select_at(manifest_variables) %>% 
    mutate(obs_id = rep(1:nrow(joint_counts), each = n_latent_classes))  %>% 
    pivot_longer(manifest_variables, names_to = 'manifest_variable',
                 values_to = 'manifest_value') %>% 
    left_join(marg_probs_manifest, 
              by = c('manifest_variable', 'manifest_value')) %>% 
    group_by(obs_id) %>% 
    summarise(log_pr_marginal = sum(log(pr_marginal)))
  # get the conditional probabilities of the latent variables
  cond_prob_latent <- fitted %>% 
    mutate(obs_id = rep(1:nrow(joint_counts), each = n_latent_classes)) %>%
    #select(obs_id, group, pr_joint) %>%
    pivot_wider(names_from = 'group', values_from = 'pr_joint') %>% 
    left_join(marg_probs_manifest_vec, by = 'obs_id') %>% 
    mutate(group1 = group1/exp(log_pr_marginal),
           group2 = group2/exp(log_pr_marginal),
           group3 = group3/exp(log_pr_marginal)) %>%
    select(-log_pr_marginal) %>% 
    adorn_percentages('row', na.rm = TRUE, c('group1', 'group2', 'group3')) %>% 
    select(-obs_id) %>%
    left_join(bind_cols(id, preproc_data), by = manifest_variables) %>% 
    group_by(ind) %>% 
    summarise(group1 = prod(group1),
              group2 = prod(group2),
              group3 = prod(group3)) %>% 
    adorn_percentages('row', na.rm = TRUE, c('group1', 'group2', 'group3')) %>% 
    bind_cols(drop_na(haplo))

  ret <- list(plot = ggtern(data = cond_prob_latent, 
                            aes(x = group1, y=group2, z=group3)) +
                            geom_point(aes(color = race)),
              model = mod,
              coords = cond_prob_latent,
              n_iters = iter)
  return(ret)
}

#init <- LatentClassEM(haplo, preproc, 3, tol = 1)
#step1 <- LatentClassEM(haplo, preproc, 3, tol = 0.15)
#step2 <- LatentClassEM(haplo, preproc, 3, tol = 1e-1)
#step3 <- LatentClassEM(haplo, preproc, 3, tol = 0.05)
#figs_list <- list(init = init, step1 = step1, step2 = step2, step3 = step3)
#save(figs_list, file = 'EM_figs_list.Rdata')

