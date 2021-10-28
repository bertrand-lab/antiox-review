# getting the empirical distribution of changes in protein expression, where the protein expression is expressed as a percentage of total protein

library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(fitdistrplus)

nunn_data <- read.csv("data/protein_expression_data/nunn_data/Nunn2013_Table_S1_annotated_formatted.csv")

# mean protein expression values for each protein for the + Fe treatments
plus_fe_cv <- nunn_data %>%
  filter(Primary.vs.Secondary.IDs.homologous.protein.sequences. == 'primary') %>% 
  rowwise() %>% 
  mutate(tp_fe_mean = mean(c(tp_fe1, tp_fe2, tp_fe3, tp_fe4), na.rm = TRUE),
         tp_fe_sd = sd(c(tp_fe1, tp_fe2, tp_fe3, tp_fe4), na.rm = TRUE),
         tp_fe_cv = tp_fe_sd/tp_fe_mean)

# mean protein expression values for each protein for the - Fe treatments
minus_fe_cv <- nunn_data %>%
  filter(Primary.vs.Secondary.IDs.homologous.protein.sequences. == 'primary') %>% 
  rowwise() %>% 
  mutate(tp_nofe_mean = mean(c(tp_nofe1, tp_nofe2, tp_nofe3, tp_nofe4), na.rm = TRUE),
         tp_nofe_sd = sd(c(tp_nofe1, tp_nofe2, tp_nofe3, tp_nofe4), na.rm = TRUE),
         tp_nofe_cv = tp_nofe_sd/tp_nofe_mean)

minus_fe_cv_percent <- minus_fe_cv %>% 
  filter(Primary.vs.Secondary.IDs.homologous.protein.sequences. == 'primary') %>% 
  rowwise() %>% 
  mutate(tp_nofe1_percent = tp_nofe1 %>% as.character() %>% as.numeric()/sum(nunn_data$tp_nofe1 %>% as.character() %>% as.numeric(), na.rm = TRUE) %>% as.numeric(),
         tp_nofe2_percent = tp_nofe2 %>% as.character() %>% as.numeric()/sum(nunn_data$tp_nofe2 %>% as.character() %>% as.numeric(), na.rm = TRUE) %>% as.numeric(),
         tp_nofe3_percent = tp_nofe3 %>% as.character() %>% as.numeric()/sum(nunn_data$tp_nofe3 %>% as.character() %>% as.numeric(), na.rm = TRUE) %>% as.numeric(),
         tp_nofe4_percent = tp_nofe4 %>% as.character() %>% as.numeric()/sum(nunn_data$tp_nofe4 %>% as.character() %>% as.numeric(), na.rm = TRUE) %>% as.numeric(),
         tp_nofe_mean_percent = mean(c(tp_nofe1_percent, tp_nofe2_percent, tp_nofe3_percent, tp_nofe4_percent), na.rm = TRUE),
         tp_nofe_sd_percent = sd(c(tp_nofe1_percent, tp_nofe2_percent, tp_nofe3_percent, tp_nofe4_percent), na.rm = TRUE),
         tp_nofe_cv_percent = tp_nofe_sd_percent/tp_nofe_mean_percent)

minus_fe_cv_percent <- minus_fe_cv %>% 
  filter(Primary.vs.Secondary.IDs.homologous.protein.sequences. == 'primary') %>% 
  rowwise() %>% 
  mutate(tp_nofe1_percent = tp_nofe1 %>% as.character() %>% as.numeric()/sum(nunn_data$tp_nofe1 %>% as.character() %>% as.numeric(), na.rm = TRUE) %>% as.numeric(),
         tp_nofe2_percent = tp_nofe2 %>% as.character() %>% as.numeric()/sum(nunn_data$tp_nofe2 %>% as.character() %>% as.numeric(), na.rm = TRUE) %>% as.numeric(),
         tp_nofe3_percent = tp_nofe3 %>% as.character() %>% as.numeric()/sum(nunn_data$tp_nofe3 %>% as.character() %>% as.numeric(), na.rm = TRUE) %>% as.numeric(),
         tp_nofe4_percent = tp_nofe4 %>% as.character() %>% as.numeric()/sum(nunn_data$tp_nofe4 %>% as.character() %>% as.numeric(), na.rm = TRUE) %>% as.numeric(),
         tp_fe1_percent = tp_fe1 %>% as.character() %>% as.numeric()/sum(nunn_data$tp_fe1 %>% as.character() %>% as.numeric(), na.rm = TRUE) %>% as.numeric(),
         tp_fe2_percent = tp_fe2 %>% as.character() %>% as.numeric()/sum(nunn_data$tp_fe2 %>% as.character() %>% as.numeric(), na.rm = TRUE) %>% as.numeric(),
         tp_fe3_percent = tp_fe3 %>% as.character() %>% as.numeric()/sum(nunn_data$tp_fe3 %>% as.character() %>% as.numeric(), na.rm = TRUE) %>% as.numeric(),
         tp_fe4_percent = tp_fe4 %>% as.character() %>% as.numeric()/sum(nunn_data$tp_fe4 %>% as.character() %>% as.numeric(), na.rm = TRUE) %>% as.numeric(),
         tp_nofe_mean_percent = mean(c(tp_nofe1_percent, tp_nofe2_percent, tp_nofe3_percent, tp_nofe4_percent), na.rm = TRUE),
         tp_fe_mean_percent = mean(c(tp_fe1_percent, tp_fe2_percent, tp_fe3_percent, tp_fe4_percent), na.rm = TRUE),
         tp_percent_foldchange = tp_nofe_mean_percent/tp_fe_mean_percent,
         tp_percent_foldchange_reverse = tp_fe_mean_percent/tp_nofe_mean_percent,
         tp_nofe_sd_percent = sd(c(tp_nofe1_percent, tp_nofe2_percent, tp_nofe3_percent, tp_nofe4_percent), na.rm = TRUE),
         tp_nofe_cv_percent = tp_nofe_sd_percent/tp_nofe_mean_percent,
         spectral_counts_mean = mean(c(tp_nofe1, tp_nofe2, tp_nofe3, tp_nofe4), na.rm = TRUE) > 10)


minus_fe_cv_percent_non_na_or_inf <- minus_fe_cv_percent %>% 
  filter(!is.na(tp_percent_foldchange),
         !is.infinite(tp_percent_foldchange))

# looking at the distribution of values
nunn2013 <- fitdistr(minus_fe_cv_percent_non_na_or_inf$tp_percent_foldchange, 'log-normal')
write.csv(data.frame(meanlog = nunn2013$estimate[1], 
                     meansd = nunn2013$estimate[2]),
          file = 'data/nunn2013_lgnormal_dist.csv')

# writing out a dataframe of empirical values
write.csv(data.frame(fold_change = minus_fe_cv_percent_non_na_or_inf$tp_percent_foldchange, 
                     data_set = rep('nunn2013', length(minus_fe_cv_percent_non_na_or_inf$tp_percent_foldchange))),
          file = 'data/nunn2013_empirical_dist.csv')

