### library(ggplot2)
library(magrittr)
library(ggpubr)
library(dplyr)

li2014_emp <- read.csv("data/li2014_empirical_dist.csv")
nunn2013_emp <- read.csv("data/nunn2013_empirical_dist.csv")
schmidt0216_emp <- read.csv("data/schmidt2016_empirical_dist.csv")

outcome_var <- function(number_to_gen, metal_per_antioxidant = 2,
                        anti_per_aa = 1/310,
                        antioxidant_per_prot_max = 0.04){
  
  ## Set up sampling so that equal number of samples are expected from each dataset
  lens <- c(length(schmidt0216_emp$fold_change),
            length(nunn2013_emp$fold_change),
            length(li2014_emp$fold_change))
  prob <- rep(1/(lens*3), lens)
  
  ## sample number_to_gen items from these empirical values
  expression_foldchange <- sample(x = c(schmidt0216_emp$fold_change,
                                        nunn2013_emp$fold_change,
                                        li2014_emp$fold_change),
                                  size=number_to_gen, replace=TRUE, prob=prob)
  antioxidant_per_prot <- antioxidant_per_prot_max*expression_foldchange
  
  protein_aa_per_N <- 0.699
  
  N_total_protein_per_N_total <- runif(n = number_to_gen, 
                                       min = 0.5, 
                                       max = 0.85)
  
  N_only <- truncnorm::rtruncnorm(n = number_to_gen, a = 0, b = Inf, 
                                  mean = 16, sd = 5)
  
  N_to_C <- N_only/106
  
  output <- metal_per_antioxidant*anti_per_aa*antioxidant_per_prot*protein_aa_per_N*N_total_protein_per_N_total*N_to_C
  
  transformed_output <- output*1e6
  
  return(transformed_output)
  
}


ranges_of_pars <- function(number_to_gen, metal_per_antioxidant = c(2),
                           anti_per_aa = c(1/400),
                           antioxidant_per_prot_max = c(0.04)){
  ###
  # wrapper function for above MC simulation
  ###

  c_to_fe_vec <- c()
  metal_per_antioxidant_vec <- c()
  anti_per_aa_vec <- c()
  antioxidant_per_prot_max_vec <- c()

  for(i in 1:length(metal_per_antioxidant)){
    for(j in 1:length(anti_per_aa)){
      for(k in 1:length(antioxidant_per_prot_max)){
        c_to_fe_vec_sub <- outcome_var(number_to_gen = number_to_gen,
                                       metal_per_antioxidant = metal_per_antioxidant[i],
                                       anti_per_aa = anti_per_aa[j],
                                       antioxidant_per_prot_max = antioxidant_per_prot_max[k])

        c_to_fe_vec <- c(c_to_fe_vec,
                         c_to_fe_vec_sub)

        metal_per_antioxidant_vec <- c(metal_per_antioxidant_vec,
                                       rep(metal_per_antioxidant[i],
                                           number_to_gen))
        anti_per_aa_vec <- c(anti_per_aa_vec,
                             rep(anti_per_aa[j],
                                 number_to_gen))
        antioxidant_per_prot_max_vec <- c(antioxidant_per_prot_max_vec,
                                          rep(antioxidant_per_prot_max[k],
                                              number_to_gen))

      }
    }
  }

  return(data.frame(c_to_fe = c_to_fe_vec,
                    metal_per_antiox = paste0("Metal per Antioxidant: ",
                                              metal_per_antioxidant_vec),
                    amino_per_anti = paste0("AA per Antioxidant: ",
                                            1/anti_per_aa_vec),
                    anti_per_prot = paste0("Antioxidant Proteomic Fraction: ",
                                           antioxidant_per_prot_max_vec)))
}

# reading in the protein length data for Frag
prot_length_frag <- read.csv("data/protein_expression_data/frag_anti_prot_lengths.csv")
express_frag <- read.csv("data/protein_expression_data/frag_max_antiox_expression.csv")

# reading in the protein length data for Phaeo
prot_length_phae <- read.csv("data/protein_expression_data/phae_anti_prot_lengths.csv")
express_phae <- read.csv("data/protein_expression_data/phae_max_antiox_expression.csv")

prot_length_overall <- rbind(prot_length_frag,
                             prot_length_phae)

# getting average protein length
prot_length_overall_sum <- prot_length_overall %>% 
  group_by(anti) %>% 
  summarize(mean_len = mean(sequence_lengths))


express_overall <- rbind(express_frag,
                         express_phae) %>% 
  group_by(antioxi_string) %>% 
  summarize(max_exp = max(total_val),
            mean_exp = mean(total_val))

number_to_gen_i <- 1000000

# running monte carlo -----------------------------------------------------------

mnfesod_both <- ranges_of_pars(number_to_gen_i, 
                               metal_per_antioxidant = c(1), 
                               anti_per_aa = c(1/prot_length_overall_sum[prot_length_overall_sum$anti == "MnFeSOD", ]$mean_len),
                               antioxidant_per_prot_max = express_overall[express_overall$antioxi_string == "MnFeSOD", ]$mean_exp)

cat_both <- ranges_of_pars(number_to_gen_i, 
                           metal_per_antioxidant = c(1), 
                           anti_per_aa = c(1/prot_length_overall_sum[prot_length_overall_sum$anti == "CAT", ]$mean_len),
                           antioxidant_per_prot_max = express_overall[express_overall$antioxi_string == "CAT", ]$mean_exp)


apx_val_both <- ranges_of_pars(number_to_gen_i, 
                               metal_per_antioxidant = c(1), 
                               anti_per_aa = c(1/prot_length_overall_sum[prot_length_overall_sum$anti == "APX", ]$mean_len),
                               antioxidant_per_prot_max = express_overall[express_overall$antioxi_string == "APX", ]$mean_exp)

ccp_val_both <- ranges_of_pars(number_to_gen_i, 
                               metal_per_antioxidant = c(2), 
                               anti_per_aa = c(1/prot_length_overall_sum[prot_length_overall_sum$anti == "CCP", ]$mean_len),
                               antioxidant_per_prot_max = express_overall[express_overall$antioxi_string == "CCP", ]$mean_exp)

cuznsod_val_both <- ranges_of_pars(number_to_gen_i, 
                                   metal_per_antioxidant = c(1), 
                                   anti_per_aa = c(1/prot_length_overall_sum[prot_length_overall_sum$anti == "CuZnSOD", ]$mean_len),
                                   antioxidant_per_prot_max = express_overall[express_overall$antioxi_string == "CuZnSOD", ]$mean_exp)

nisod_val_both <- ranges_of_pars(number_to_gen_i, 
                                 metal_per_antioxidant = c(1), 
                                 anti_per_aa = c(1/prot_length_overall_sum[prot_length_overall_sum$anti == "NiSOD", ]$mean_len),
                                 antioxidant_per_prot_max = express_overall[express_overall$antioxi_string == "NiSOD", ]$mean_exp)



# aggregating the results
overall_c_to_fe_variation_both <- data.frame(c_to_fe_overall = mnfesod_both$c_to_fe + 
                                               apx_val_both$c_to_fe + 
                                               ccp_val_both$c_to_fe + 
                                               cat_both$c_to_fe)
overall_c_to_fe_variation_both$mnfesod <- rep('MnFeSOD as FeSOD', nrow(overall_c_to_fe_variation_both)) 

overall_c_to_fe_variation_both_no_mnfesod <- data.frame(c_to_fe_overall = apx_val_both$c_to_fe + 
                                               ccp_val_both$c_to_fe + 
                                               cat_both$c_to_fe)
overall_c_to_fe_variation_both_no_mnfesod$mnfesod <- rep('MnFeSOD as MnSOD', nrow(overall_c_to_fe_variation_both_no_mnfesod)) 

overall_c_to_fe <- rbind(overall_c_to_fe_variation_both, 
                         overall_c_to_fe_variation_both_no_mnfesod)

# plotting  -----------------------------------------------------------

contribution_quantiles_both <- overall_c_to_fe %>%
  group_by(mnfesod) %>%
  # sample_n(1000) %>%
  summarize(lower = quantile(c_to_fe_overall, probs = .025),
            upper = quantile(c_to_fe_overall, probs = .975),
            median = quantile(c_to_fe_overall, probs = 0.5))

print(contribution_quantiles_both)

both_contribution_overall_p <- overall_c_to_fe %>% 
  ggplot(aes(x = c_to_fe_overall)) +
  geom_density(aes(fill = mnfesod),
               alpha = 0.4) +
  theme_bw() +
  xlim(0, 9) +
  xlab('Overall Antioxidant Contribution to Fe:C (umol/mol)') +
  ylab('Kernel Density') +
  theme(legend.position = c(0.8, 0.8), 
        legend.title = element_blank());both_contribution_overall_p

ccp_both_p1 <- ccp_val_both %>% 
  ggplot(aes(x = c_to_fe)) +
  geom_density(fill = 'firebrick4', alpha = 0.4) +
  theme_bw() +
  xlim(0, 9) +
  ylab('Kernel Density') +
  xlab('Cytochrome C Peroxidase Contribution to Fe:C (umol/mol)');ccp_both_p1

cat_both_p1 <- cat_both %>% 
  ggplot(aes(x = c_to_fe)) +
  geom_density(fill = 'firebrick4', alpha = 0.4) +
  theme_bw() + 
  xlim(0, 9) +
  ylab('Kernel Density') +
  xlab('Catalase Contribution to Fe:C (umol/mol)');cat_both_p1

mnfesod_both_p1 <- mnfesod_both %>% 
  ggplot(aes(x = c_to_fe)) +
  geom_density(fill = 'firebrick4', alpha = 0.4) +
  theme_bw() + 
  xlim(0, 9) +
  ylab('Kernel Density') +
  xlab('MnFeSOD Contribution to Mn, Fe:C (umol/mol)');mnfesod_both_p1

apx_both_p1 <- apx_val_both %>% 
  ggplot(aes(x = c_to_fe)) +
  geom_density(fill = 'firebrick4', alpha = 0.4) +
  theme_bw() + 
  xlim(0, 9) +
  ylab('Kernel Density') +
  xlab('Ascorbate Peroxidase Contribution to Fe:C (umol/mol)');apx_both_p1


fe_anti_plot <- ggarrange(ggarrange(ccp_both_p1, 
                    cat_both_p1, 
                    mnfesod_both_p1, 
                    apx_both_p1, align = 'hv', ncol = 2, nrow = 2,
                    labels = c('a', 'b', 'c', 'd')),
          both_contribution_overall_p, ncol = 1, labels = c('', 'e'))

ggsave(fe_anti_plot, filename = 'figures/fe_anti_plot.png', width = 9.17, height = 7.27)

# plotting other micronutrients -------------------------------------------

nisod_both_p1 <- nisod_val_both %>% 
  ggplot(aes(x = c_to_fe)) +
  geom_density(fill = 'firebrick4', alpha = 0.4) +
  # geom_histogram() +
  theme_bw() + 
  ylab('Kernel Density') +
  # xlim(0,) +
  # xlim(0, 50) +
  xlim(0, 2) +
  xlab('NiSOD Contribution to Ni:C (umol/mol)');nisod_both_p1


cuznsod_both_p1 <- cuznsod_val_both %>% 
  ggplot(aes(x = c_to_fe)) +
  geom_density(fill = 'firebrick4', alpha = 0.4) +
  # geom_histogram() +
  theme_bw() + 
  # xlim(0,) +
  xlim(0, 1) +
  ylab('Kernel Density') +
  xlab('CuZnSOD Contribution to Cu, Zn:C (umol/mol)');cuznsod_both_p1

print(quantile(mnfesod_both$c_to_fe, probs = c(0.025, 0.5, 0.975)))
print(quantile(nisod_val_both$c_to_fe, probs = c(0.025, 0.5, 0.975)))
print(quantile(cuznsod_val_both$c_to_fe, probs = c(0.025, 0.5, 0.975)))

other_micro <- ggarrange(nisod_both_p1, cuznsod_both_p1, labels = c('a', 'b'), ncol = 2)

ggsave(other_micro, filename = 'figures/ni_zn_cu_sod.png', width = 9.17, height = 3.35)



