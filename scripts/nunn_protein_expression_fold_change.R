

library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(fitdistrplus)

nunn_data <- read.csv("data/protein_expression_data/nunn_data/Nunn2013_Table_S1_annotated_formatted.csv")

# getting the sum of spectral counts per coarse grained protein group for + Fe
plus_fe_cv <- nunn_data %>%
  filter(Primary.vs.Secondary.IDs.homologous.protein.sequences. == 'primary') %>% 
  rowwise() %>% 
  mutate(tp_fe_mean = mean(c(tp_fe1, tp_fe2, tp_fe3, tp_fe4), na.rm = TRUE),
         tp_fe_sd = sd(c(tp_fe1, tp_fe2, tp_fe3, tp_fe4), na.rm = TRUE),
         tp_fe_cv = tp_fe_sd/tp_fe_mean)

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
         tp_nofe_sd_percent = sd(c(tp_nofe1_percent, tp_nofe2_percent, tp_nofe3_percent, tp_nofe4_percent), na.rm = TRUE),
         tp_nofe_cv_percent = tp_nofe_sd_percent/tp_nofe_mean_percent,
         spectral_counts_mean = mean(c(tp_nofe1, tp_nofe2, tp_nofe3, tp_nofe4), na.rm = TRUE) > 10)


minus_fe_cv_percent_non_na_or_inf <- minus_fe_cv_percent %>% 
  filter(!is.na(tp_percent_foldchange),
         !is.infinite(tp_percent_foldchange),
         tp_percent_foldchange > 0)

minus_fe_cv_percent %>% 
  ggplot(aes(x = tp_nofe_cv_percent, y = tp_nofe_mean_percent)) + 
  geom_point(aes(shape = spectral_counts_mean)) +
  scale_y_log10() +
  facet_grid(~spectral_counts_mean)
  # scale_x_log10()

nunn2013 <- fitdistr(minus_fe_cv_percent_non_na_or_inf$tp_percent_foldchange, 'log-normal')

write.csv(data.frame(meanlog = nunn2013$estimate[1], 
                     meansd = nunn2013$estimate[2]),
          file = 'data/nunn2013_lgnormal_dist.csv')

# 
# par(mfrow = c(2, 1))
# 
# blah <- minus_fe_cv_percent %>% filter(ec_no1 == "1.15.1.1" | ec_no1 == "1.11.1.7")
# hist(blah$tp_percent_foldchange, xlim = c(0, 30), breaks = 100)
# 
# minus_fe_cv_percent_non_na_or_inf$tp_percent_foldchange %>% 
#   hist(breaks = 1000, 
#        main = 'Fold Change Distribution from low to high Fe (Nunn et al 2013)', 
#        xlab = 'Fold Change',  xlim = c(0, 30))
# minus_fe_cv_percent_non_na_or_inf$tp_percent_foldchange %>% quantile()
# rgamma(10000, shape = 1.15, rate = 0.66) %>% hist(breaks = 1000)
# rlnorm(1000, meanlog = 0.0634, sdlog = 0.91827) %>% 
#   hist(breaks = 1000, main = 'Fold Change Distribution (Log Normal Distribution)', 
#        xlab = 'Fold Change', xlim = c(0, 30))
# rlnorm(10000, meanlog = 0.0634, sdlog = 0.91827) %>% quantile()
# 
# minus_fe_cv_percent %>% 
#   ggplot(aes(x = tp_nofe_sd_percent, tp_nofe_mean_percent)) +
#   geom_point() +
#   scale_y_log10() +
#   scale_x_log10()
# 
# minus_fe_cv_percent %>% 
#   ggplot(aes(x = tp_percent_foldchange)) +
#   geom_histogram() +
#   scale_x_log10()
# 
# minus_fe_cv_percent %>% 
#   filter(tp_nofe_cv_percent < 0.02) %>% 
#   ggplot(aes(x = tp_nofe_mean_percent, y = tp_nofe_cv_percent)) +
#   geom_point() +
#   scale_y_log10() +
#   scale_x_log10() +
#   geom_smooth(method = "lm")
# 
# minus_fe_cv_percent %>% filter(ec_no1 == "1.11.1.7") %>% as.data.frame()
# minus_fe_cv_percent %>% filter(ec_no1 == "1.15.1.1" | ec_no1 == "1.11.1.7") %>% 
#   ggplot(aes(x = tp_percent_foldchange)) +
#   geom_histogram()
# 
# 
# lm_out <- lm(log(minus_fe_cv_percent$tp_nofe_cv_percent) ~ log(minus_fe_cv_percent$tp_nofe_mean_percent))
# 
# summary(lm_out)
# 
# 
# # histogram of all cv for both minus and plus cv
# histogram_vals_cv <- data.frame(cv_proteins = c(plus_fe_cv$tp_fe_cv,
#                                                 minus_fe_cv$tp_nofe_cv),
#                                 mean_proteins = c(plus_fe_cv$tp_fe_mean,
#                                                   minus_fe_cv$tp_nofe_mean),
#                                 sd_proteins = c(plus_fe_cv$tp_fe_sd,
#                                                 minus_fe_cv$tp_nofe_sd),
#                                 cv_condition = c(rep('Fe', nrow(plus_fe_cv)),
#                                                  rep('noFe', nrow(minus_fe_cv))))
# 
# 
# histogram_vals_cv2 <- histogram_vals_cv %>% filter(cv_proteins > 0)
# 
# lm_out <- lm(histogram_vals_cv2$cv_proteins ~ histogram_vals_cv2$mean_proteins)
# 
# summary(lm_out)
# 
# par(mfrow = c(2, 2))
# plot(lm_out)
# 
# minus_fe_cv_percent %>% 
#   ggplot(aes(y = tp_nofe_cv_percent, x = tp_nofe_mean_percent)) +
#   geom_point() +
#   scale_x_log10() +
#   scale_y_log10() +
#   geom_smooth()
# 
# histogram_vals_cv2 %>% 
#   group_by(cv_condition) %>% 
#   summarize(mean_proteins_sum = sum(mean_proteins))
# 
# histogram_vals_cv2 %>% 
#   filter(cv_condition == 'Fe') %>%
#   ggplot(aes(x = mean_proteins, y = cv_proteins)) +
#   geom_jitter() +
#     geom_smooth(method = "gam")
#   # scale_x_log10() +
#   # scale_y_log10()
# 
# nunn_protein_histogram_cv <- histogram_vals_cv2 %>% 
#   ggplot(aes(x = cv_proteins)) +
#   geom_histogram() +
#   theme_bw() +
#   ylab('Count') +
#   xlab('Protein-specific Coefficient of Variation\nacross constant conditions');nunn_protein_histogram_cv
# # facet_wrap(~cv_condition, nrow = 2)
# 
# mean_cv_cor_nunn <- histogram_vals_cv2 %>% 
#   ggplot(aes(x = mean_proteins, y = cv_proteins)) +
#   geom_point(alpha = 0.8) +
#   theme_bw() +
#   ylab('Protein-specific Coefficient of Variation') +
#   xlab('Mean Protein Abundance Value')
# 
# 
# nunn_supp_cv_plot <- ggarrange(nunn_protein_histogram_cv, 
#                                mean_cv_cor_nunn, nrow = 1,
#                                labels = c('a', 'b'))
# 
# ggsave(nunn_supp_cv_plot, filename = 'figures/nunn_supp_cv_plot.png',
#        width = 9.02, height = 5.75)
# 
