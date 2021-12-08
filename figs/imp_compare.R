# comparing implementations with a summary value

library(tidyverse)
library(cowplot)
library(latex2exp)

sumtab <- read_csv('../data/models123summary.csv') %>%
  rename(expectation = 'E(g(Q))', lb = '95% MC CI LB', ub = '95% MC CI UB') %>%
  mutate(Model = paste('Model', Model))

ggplot(sumtab) + 
  geom_errorbar(aes(Ethnicity, ymin=lb, ymax=ub, color = Implementation,
                    width = 0.3, linetype = Implementation)) +
  geom_point(aes(Ethnicity, expectation, color = Implementation)) +
  facet_grid(cols = vars(Model)) +
  labs(y = TeX("$\\phi_t$"))
ggsave('model_comp.png', width = 12, height = 5)

