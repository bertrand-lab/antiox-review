library(ggplot2)
library(readxl)

# read in data
gene_data <- read_xlsx("data/NIHMS570024-supplement-02 (1).xlsx", sheet = 1)

# cleaning data
gene_data$mops_complete_no_br <- gsub(pattern = '\\[|\\]', 
                                            replacement = '',
                                            x = gene_data$`MOPS complete`) %>% as.numeric()
gene_data$mops_minimal_no_br <- gsub(pattern = '\\[|\\]', 
                                            replacement = '',
                                            x = gene_data$`MOPS minimal`) %>% as.numeric()
gene_data$mops_complete_wo_met_no_br <- gsub(pattern = '\\[|\\]', 
                                            replacement = '',
                                            x = gene_data$`MOPS complete without methionine`) %>% as.numeric()

gene_data_above_zero <- gene_data %>%
  filter(mops_complete_no_br > 0,
         mops_minimal_no_br > 0,
         mops_complete_wo_met_no_br > 0)

# getting distribution parameters
gene_data_pars <- fitdistr(x = c(gene_data_above_zero$mops_complete_no_br/gene_data_above_zero$mops_minimal_no_br,
               gene_data_above_zero$mops_complete_no_br/gene_data_above_zero$mops_complete_wo_met_no_br), 'lognormal')

write.csv(data.frame(meanlog = gene_data_pars$estimate[1], meansd = gene_data_pars$estimate[2]),
          file = 'data/li2014_lgnormal_dist.csv')

# writing file of empirical distribution
write.csv(data.frame(fold_change = c(gene_data_above_zero$mops_complete_no_br/gene_data_above_zero$mops_minimal_no_br,
                                     gene_data_above_zero$mops_complete_no_br/gene_data_above_zero$mops_complete_wo_met_no_br), data_set = rep('li2014', length(c(gene_data_above_zero$mops_complete_no_br/gene_data_above_zero$mops_minimal_no_br,
                                                                                                                                                                  gene_data_above_zero$mops_complete_no_br/gene_data_above_zero$mops_complete_wo_met_no_br)))),
          file = 'data/li2014_empirical_dist.csv')


