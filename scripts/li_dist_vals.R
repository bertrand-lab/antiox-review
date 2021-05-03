library(ggplot2)
library(readxl)


gene_data <- read_xlsx("data/NIHMS570024-supplement-02 (1).xlsx", sheet = 1)

head(gene_data)

gene_data$mops_complete_no_br <- gsub(pattern = '\\[|\\]', 
                                            replacement = '',
                                            x = gene_data$`MOPS complete`) %>% as.numeric()
gene_data$mops_minimal_no_br <- gsub(pattern = '\\[|\\]', 
                                            replacement = '',
                                            x = gene_data$`MOPS minimal`) %>% as.numeric()
gene_data$mops_complete_wo_met_no_br <- gsub(pattern = '\\[|\\]', 
                                            replacement = '',
                                            x = gene_data$`MOPS complete without methionine`) %>% as.numeric()

gene_data %>% 
  mutate(mops_complete_min = mops_complete_no_br/mops_minimal_no_br) %>% 
  ggplot(aes(x = mops_complete_min)) +
  geom_histogram(binwidth = 0.1) +
  xlim(0, 30)

gene_data %>% 
  mutate(mops_complete_min = mops_complete_no_br/mops_minimal_no_br) %>% 
  ggplot(aes(x = mops_complete_min)) +
  geom_histogram(binwidth = 0.1) +
  xlim(0, 30) +
  geom_vline(xintercept = 1)


# gene_data_above_zero <- gene_data %>% 
#   filter(mops_complete_no_br > 0,
#          mops_minimal_no_br > 0,
#          mops_complete_wo_met_no_br > 0)
# 
# 
#   ggplot(aes(x = mops_complete_no_br, y = mops_minimal_no_br)) +
#   geom_point()
# 
#   
# gene_data %>% 
#   filter(mops_complete_no_br == 0) %>% dim()


gene_data_pars <- fitdistr(x = c(gene_data_above_zero$mops_complete_no_br/gene_data_above_zero$mops_minimal_no_br,
               gene_data_above_zero$mops_complete_no_br/gene_data_above_zero$mops_complete_wo_met_no_br), 'lognormal')

write.csv(data.frame(meanlog = gene_data_pars$estimate[1], meansd = gene_data_pars$estimate[2]),
          file = 'data/li2014_lgnormal_dist.csv')

# hist(gene_data_above_zero$mops_complete_no_br/gene_data_above_zero$mops_minimal_no_br,
#      xlim = c(0, 30), breaks = 3000)
# hist(gene_data_above_zero$mops_complete_no_br/gene_data_above_zero$mops_complete_wo_met_no_br,
#      xlim = c(0, 30), breaks = 3000)


gene_data %>% 
  filter(grepl('sod', Gene)) %>% 
  ggplot(aes(x = mops_complete_no_br, y = mops_minimal_no_br)) +
  geom_text(aes(label = Gene))






