library(tidyverse)
library(cowplot)
library(GGally)

source('latent_class_em.R')

# Types of fit
fit1 <- LatentClassEM(haplo, preproc, 1, 3, tol = 3e-1)
fit1$plot + labs(color = 'Ethnicity')
ggsave('clustering1.png')
fit20 <- LatentClassEM(haplo, preproc, 20, 3, tol = 3e-1)
fit20$plot + labs(color = 'Ethnicity')
ggsave('clustering20.png')


fit7 <- LatentClassEM(haplo, preproc, 7, 3, tol = 3e-1)
fit7$plot + labs(color = 'Ethnicity')
ggsave('clustering7.png')
fit8 <- LatentClassEM(haplo, preproc, 8, 3, tol = 3e-1)
fit8$plot + labs(color = 'Ethnicity')
ggsave('clustering8.png')
fit9 <- LatentClassEM(haplo, preproc, 9, 3, tol = 3e-1)
fit9$plot + labs(color = 'Ethnicity')
ggsave('clustering9.png')

# to understand variable selection
n_pops <- 3
preproc_out <- preproc(haplo)
out_tibble <- tibble()
for (n_vars in 1:20){
  cols <- paste0('Snp', snp[1:n_vars])
  preproc_data <- preproc_out$preproc_data[,cols]
  manifest_variables <- colnames(preproc_data)
  # get joint counts
  joint_counts <- preproc_data %>% 
    group_by_all() %>% 
    summarise(n = n(), .groups = 'drop') %>%
    mutate(n_vars = n_vars) %>%
    ungroup() %>%
    group_by(n_vars) %>%
    summarise(n_filled = n(),
              mean_count_no0 = mean(n)) %>%
    mutate(n_cells = n_pops * 2^n_vars,
           n_empty = n_cells - n_filled,
           prop_empty = n_empty / n_cells,
           mean_count = 1/n_cells * mean_count_no0*n_filled,
           n_params = (n_pops - 1) + n_vars + (n_pops - 1)*n_vars)
  out_tibble <- bind_rows(out_tibble, joint_counts)
}

vars_plot <- ggplot(out_tibble) + 
  geom_point(aes(n_vars, mean_count)) + 
  labs(y = 'Average cell count (including 0s)') + 
  scale_x_continuous(name="Number of SNP loci included in model",
                     limits=c(1, 12), breaks = 1:12)
cells_plot <- ggplot(out_tibble) +
  geom_point(aes(n_vars, n_cells)) + 
  labs(y = 'Number of cells in joint table') + 
  scale_x_continuous(name="Number of SNP loci included in model",
                     limits=c(1, 12), breaks = 1:12)
plot_grid(vars_plot, cells_plot)
ggsave('variable_selection.png')

deviances <- c()
zeros <- tibble(group = rep(c('1', '2', '3'), 4), 
                ethnicity = rep(c('African', 'African American', 'European', 'Japanese'), each = 3),
                n = 0)
preds_tables <- list()
for (n_vars in 1:20){
  fit <- LatentClassEM(haplo, preproc, n_vars, 3, tol = 3e-1)
  deviances <- c(deviances, deviance(fit$model))
  preds_tables[[n_vars]] <- fit$coords %>%
    select(group1, group2, group3) %>%
    as.matrix %>%
    apply(1, \(row_i) which(row_i == max(row_i))) %>%
    cbind(fit$coords$race) %>%
    as_tibble() %>%
    rename(group = '.', ethnicity = 'V2') %>%
    group_by(group, ethnicity) %>%
    summarise(n = n(), .groups = 'drop') %>%
    bind_rows(zeros) %>%
    group_by(group, ethnicity) %>%
    summarise(n = sum(n), .groups = 'drop') %>%
    pivot_wider(names_from = 'ethnicity', values_from = 'n') %>%
    mutate(n_vars = n_vars)
}
preds_df <- do.call(bind_rows, preds_tables) %>%
  pivot_longer(cols = c(-group, -n_vars), names_to = 'ethnicity', 
               values_to = 'count')

groupings <- ggplot(preds_df) +
  geom_point(aes(n_vars, count, color = group)) +
  geom_line(aes(n_vars, count, color = group)) + 
  facet_wrap(~ ethnicity) +
  labs(x = ' Number of SNP loci included in model',
       y = 'Number of individuals grouped in latent class',
       color = 'Latent class')

deviances_plot <- ggplot() + 
  geom_point(aes(1:20, deviances)) +
  geom_line(aes(1:20, deviances)) +
  labs(x = 'Number of SNP loci included in model',
       y = 'Model residual deviance')


get_dist <- function(matr){
  dist1 <- sqrt(sum((matr[,1] - matr[,2])^2))
  dist2 <- sqrt(sum((matr[,2] - matr[,3])^2))
  dist3 <- sqrt(sum((matr[,1] - matr[,3])^2))
  return(c(dist1, dist2, dist3))
}

distances_plot <- do.call(rbind, 
        preds_tables %>%
           lapply(\(tab) tab %>% 
              select_at(c(2, 4, 5)) %>%
              as.matrix() %>%
              get_dist())) %>%
  as_tibble() %>%
  mutate(n_vars = 1:20) %>%
  rename('A-E' = 'V1',
         'E-J' = 'V2',
         'A-J' = 'V3') %>%
  pivot_longer(cols = -n_vars, names_to = 'Comparison',
               values_to = 'Euc. Dist.') %>%
  ggplot(aes(n_vars, `Euc. Dist.`, color = Comparison)) + 
    geom_point() +
    geom_line() +
    labs(x = 'Eucledian Distance', 'Number of SNP loci included in model')

ggsave('groupings_plot.png', groupings)
ggsave('deviances.png', deviances_plot)
ggsave('distances.png', distances_plot)

# a fun result
deviances[20] - deviances[1]


plot_grid(distances_plot, deviances_plot)


