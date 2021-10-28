## antioxi tara plotting stoichiometry
library(ggplot2)
library(magrittr)
library(dplyr)
library(ggridges)
library(ggpubr)

stoich_t <- read.csv(file = "data/tara_stoichiometry_antiox.csv")

levels(stoich_t$prot_name) <- c('APX', 'CAT', 'CCP', 'Ferritin', 'GPX', 'PRX', 'SOD')

stoich_t_ratio <- stoich_t %>% 
  mutate(C_N = C/N,
         H_N = H/N,
         O_N = O/N,
         S_N = S/N)

hydro_plot <- stoich_t_ratio %>% 
  filter(prot_name != "Ferritin") %>% 
  ggplot(aes(x = H_N,
             y = prot_name,
             fill = prot_name)) +
  stat_density_ridges(quantile_lines = TRUE,
                      alpha = 0.8) +
  theme_bw() +
  xlab("H:N") +
  ylab('Protein Name');hydro_plot
  
carbon_plot <- stoich_t_ratio %>% 
  filter(prot_name != "Ferritin") %>% 
  ggplot(aes(x = C_N,
             y = prot_name,
             fill = prot_name)) +
  stat_density_ridges(quantile_lines = TRUE,
                      alpha = 0.8) +
  theme_bw() +
  xlab("C:N") +
  ylab('Protein Name');carbon_plot

oxy_plot <- stoich_t_ratio %>% 
  filter(prot_name != "Ferritin") %>% 
  ggplot(aes(x = O_N,
             y = prot_name,
             fill = prot_name)) +
  stat_density_ridges(quantile_lines = TRUE,
                      alpha = 0.8) +
  theme_bw() +
  xlab("O:N") +
  ylab('Protein Name');oxy_plot

sul_plot <- stoich_t_ratio %>% 
  filter(prot_name != "Ferritin") %>% 
  ggplot(aes(x = S_N,
             y = prot_name,
             fill = prot_name)) +
  stat_density_ridges(quantile_lines = TRUE,
                      alpha = 0.8) +
  theme_bw() +
  xlab("S:N") +
  ylab('Protein Name');sul_plot


compiled_macro_stoich <- ggarrange(hydro_plot, carbon_plot, oxy_plot, sul_plot, legend = "none", 
          labels = c("a", "b", "c", "d"))

ggsave(compiled_macro_stoich, file = 'figures/compiled_macro_stoich.png', width = 7.24, height = 5.88)

# getting proetein lengths and summary statistics
summary_table_stoich <- stoich_t_ratio %>% 
  group_by(prot_name) %>% 
  summarize(median_length = median(Sequence_Length),
            sd_length = sd(Sequence_Length),
            median_c_n = median(C_N) %>% round(2),
            median_h_n = median(H_N) %>% round(2),
            median_o_n = median(O_N) %>% round(2),
            median_s_n = median(S_N, na.rm = TRUE) %>% round(2))

write.csv(summary_table_stoich, 
          'data/summary_table_stoich.csv', row.names = FALSE)

# number atoms per protein
n_atoms_per_prot <- stoich_t_ratio %>% 
  group_by(prot_name) %>% 
  summarize(median_n = median(N))

write.csv(n_atoms_per_prot, 'data/n_atoms_per_prot.csv',
          row.names = FALSE)







