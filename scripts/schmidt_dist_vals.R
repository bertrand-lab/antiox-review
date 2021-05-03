#### schmidt et al data

library(ggplot2)
library(readxl)

schmidt <- read_xlsx('data/41587_2016_BFnbt3418_MOESM18_ESM.xlsx', sheet = 23, skip = 2)
all_schmidt <- read_xlsx('data/41587_2016_BFnbt3418_MOESM18_ESM.xlsx', sheet = 5, skip = 2)

schmidt_foldchange <- schmidt %>% dplyr::select(contains('medianRatio'))

# schmidt %>% 
#   ggplot(aes(x = medianRatio_Glycerin.AA_vs_Glucose)) +
#   geom_histogram(binwidth = 0.1) +
#   xlim(0, 30)

all_schmidt_vals <- all_schmidt[, c(7:25)]

empty_matrix <- matrix(nrow = 2039)

for(column_i in 1:ncol(all_schmidt_vals)){
  
  # column_i <- 1
  print(column_i)
  column_i_sub <- all_schmidt_vals[ ,column_i]
  for(div_i in 1:ncol(all_schmidt_vals)){
    # div_i <- 2
    print(div_i)
    vec_subset <- as.vector(column_i_sub/all_schmidt_vals[, div_i])
    empty_matrix <- cbind(empty_matrix, vec_subset[,1])
  }
}

empty_matrix_g <- empty_matrix[,-1]
values <- empty_matrix_g[upper.tri(empty_matrix_g, diag = FALSE)]

# fort_m <- fortify(empty_matrix_g)

# autot(log(empty_matrix_g))


vector_val <- as.vector(values)
# length(vector_val)

hist(vector_val, xlim = c(0, 30), breaks = 1000, main = 'E. coli Distribution of Fold Changes')


schmidt2016_dis <- fitdistr(vector_val, 'log-normal')

write.csv(data.frame(meanlog = schmidt2016_dis$estimate[1], 
                     meansd = schmidt2016_dis$estimate[2]),
          file = 'data/schmidt2016_lgnormal_dist.csv')



