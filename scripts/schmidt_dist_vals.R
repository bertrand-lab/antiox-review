#### schmidt et al data

library(ggplot2)
library(readxl)

schmidt <- read_xlsx('data/41587_2016_BFnbt3418_MOESM18_ESM.xlsx', sheet = 23, skip = 2)
all_schmidt <- read_xlsx('data/41587_2016_BFnbt3418_MOESM18_ESM.xlsx', sheet = 5, skip = 2)

# all conditions subsetted
all_schmidt_vals <- all_schmidt[, c(7:25)]

# getting an all-vs-all comparison
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

vector_val <- as.vector(values)

schmidt2016_dis <- fitdistr(vector_val, 'log-normal')
write.csv(data.frame(meanlog = schmidt2016_dis$estimate[1], 
                     meansd = schmidt2016_dis$estimate[2]),
          file = 'data/schmidt2016_lgnormal_dist.csv')

# writing empirical distribution for Schmidt et al 2016
write.csv(data.frame(fold_change = vector_val, data_set = rep('schmidt2016', length(vector_val))),
          file = 'data/schmidt2016_empirical_dist.csv')


