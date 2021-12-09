# assessing the similarity between "by hand" code and NIMBLE code

library(tidyverse)
library(coda)

haplo <- read_csv('../data/haplotype.csv') %>%
  drop_na()
ethnicity <- haplo$race
afam <- which(ethnicity == 'African American')
euro <- which(ethnicity == 'European')
japa <- which(ethnicity == 'Japanese')
afri <- which(ethnicity == 'African')
eth_id <- list(afam = afam, euro = euro, japa = japa, afri = afri)

dist_from_center <- function(draw){
  center <- apply(draw, 2, mean)
  mean_dist <- apply(draw, 1, \(coord) sqrt(sum((coord - center)^2))) %>%
    mean()
  return(mean_dist)
}

mc_ci <- function(index){
  bhq_eth <- apply(bhq[,index,], 1, dist_from_center)
  pplq_eth <- apply(pplq[,index,], 1, dist_from_center)
  nb <- effectiveSize(bhq_eth)
  np <- effectiveSize(pplq_eth)
  pointb <- mean(bhq_eth)
  pointp <- mean(pplq_eth)
  mccib <- pointb + c(-1,1)*qnorm(0.975)*sd(bhq_eth)/sqrt(nb)
  mccip <- pointp + c(-1,1)*qnorm(0.975)*sd(pplq_eth)/sqrt(np)
  ret <- matrix(c(pointb, mccib,
                  pointp, mccip), nrow = 2, ncol = 3, byrow = TRUE)
  return(ret)
}

make_df <- function(eth_id, model_number){
  mccis <- lapply(eth_id, mc_ci)
  mccis <- do.call(rbind, mccis)
  colnames(mccis) <- c('E(g(Q))', '95% MC CI LB', '95% MC CI UB')
  mccis_df <- as_tibble(mccis) %>%
    add_column(Ethnicity = rep(c('African American', 'European', 'Japanese', 
                                 'African'), each = 2), .before = 1) %>%
    add_column(Model = model_number, .before = 1) %>%
    add_column(Implementation = rep(c('By hand', 'NIMBLE'), 4), .after = 2)
  return(mccis_df)
}



# Model 1
load('Run2/mcmc_BHBinX.Rdata')
bh <- mcmc_out
bhq <- aperm(bh$draws_Q[,,1,], c(3, 1, 2))
load('Run2/mcmc_PPLBinX.Rdata')
ppl <- mcmc_out$samples
pplq <- array(ppl[,which(str_detect(colnames(ppl), 'Q'))], c(10000, 168, 3))
mod1 <- make_df(eth_id, 1)

# Model 2
load('Run2/mcmc_BHMultiX2.Rdata')
bh <- mcmc_out
bhq <- aperm(bh$draws_Q[,,1,], c(3, 1, 2))
load('Run2/mcmc_PPLMultiX2.Rdata')
ppl <- mcmc_out$samples
pplq <- array(ppl[,which(str_detect(colnames(ppl), 'Q'))], c(10000, 168, 3))
mod2 <- make_df(eth_id, 2)

# Model 3
load('Run2/mcmc_BHMultiX.Rdata')
bh <- mcmc_out
bhq <- aperm(bh$draws_Q[,,1,], c(3, 1, 2))
load('Run2/mcmc_PPLMultiX.Rdata')
ppl <- mcmc_out$samples
pplq <- array(ppl[,which(str_detect(colnames(ppl), 'Q'))], c(10000, 168, 3))
mod3 <- make_df(eth_id, 3)


bind_rows(mod1, mod2, mod3) %>%
  write_csv('models123summary.csv')



# Model 2 -- further investigation
load('Run2/mcmc_BHMultiX2.Rdata')
bh <- mcmc_out
alpha_bh <- t(bh$draws_alpha[,1,])
load('Run2/mcmc_PPLMultiX2.Rdata')
ppl <- mcmc_out$samples
alpha_ppl <- ppl[,which(str_detect(colnames(ppl), 'alpha'))]
conc_freq_bh <- janitor::adorn_percentages(as.data.frame(alpha_bh),
  denominator = 'row', na.rm = TRUE, c(1, 2, 3))
conc_freq_ppl <- janitor::adorn_percentages(as.data.frame(alpha_ppl),
  denominator = 'row', na.rm = TRUE, c(1, 2, 3))
# expectations
crint_lbbh <- apply(conc_freq_bh, 2, quantile, 0.025)
crint_ubbh <- apply(conc_freq_bh, 2, quantile, 0.975) 
crint_lbppl <- apply(conc_freq_ppl, 2, quantile, 0.025)
crint_ubppl <- apply(conc_freq_ppl, 2, quantile, 0.975)
n_cfbh <- effectiveSize(conc_freq_bh)
n_cfppl <- effectiveSize(conc_freq_ppl)
point_cfbh <- apply(conc_freq_bh, 2, mean)
point_cfppl <- apply(conc_freq_ppl, 2, mean)
sd_cfbh <- apply(conc_freq_bh, 2, sd)
sd_cfppl <- apply(conc_freq_ppl, 2, sd)
mccis2 <- tibble(point = c(point_cfbh, point_cfppl),
                 cri_lb = c(crint_lbbh, crint_lbppl),
                 cri_ub = c(crint_ubbh, crint_ubppl),
                 sd = c(sd_cfbh, sd_cfppl),
                 n_eff = c(n_cfbh, n_cfppl)) %>%
          mutate(mc_lb = point - qnorm(0.975)*sd/n_eff,
                 mc_ub = point + qnorm(0.975)*sd/n_eff) %>%
          round(3)
# concentration parameters
crint_bhlb <- apply(alpha_bh, 2, quantile, 0.025)
crint_bhub <- apply(alpha_bh, 2, quantile, 0.975)
crint_ppllb <- apply(alpha_ppl, 2, quantile, 0.025)
crint_pplub <- apply(alpha_ppl, 2, quantile, 0.975)
n_bh <- effectiveSize(alpha_bh)
n_ppl <- effectiveSize(alpha_ppl)
point_bh <- apply(alpha_bh, 2, mean)
point_ppl <- apply(alpha_ppl, 2, mean)
sd_bh <- apply(alpha_bh, 2, sd)
sd_ppl <- apply(alpha_ppl, 2, sd)
mccis3 <- tibble(point = c(point_bh, point_ppl),
                 cri_lb = c(crint_bhlb, crint_ppllb),
                 cri_ub = c(crint_bhub, crint_pplub),
                 sd = c(sd_bh, sd_ppl),
                 n_eff = c(n_bh, n_ppl)) %>%
          mutate(mc_lb = point - qnorm(0.975)*sd/n_eff,
                 mc_ub = point + qnorm(0.975)*sd/n_eff) %>%
          round(3)
write_csv(mccis2, 'model2_expectation.csv')
write_csv(mccis3, 'model2_concentration.csv')

write_csv(bind_cols(mccis2, mccis3), 'model2_all.csv')

