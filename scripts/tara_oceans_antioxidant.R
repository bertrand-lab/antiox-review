library(ggplot2)

eek <- read.delim(file = "../../antiox_prelim/data/tara_ocean_smags/SMAGs_v1_EggNog.tsv", 
                  col.names = c("Gene_ID", "unsure_1", "unsure_2", "unsure_3", "Taxonomic_affiliation_unsure", "unsure_4", "GO_id", "EC_number", "KO_number", "KO_number2", "unsure_5", "KEGG_Reaction_Number", "KEGG_Reaction_Class", "KO_number3", "unsure_8", "unsure_9", "unsure_10", "unsure_11", "COG", "unsure_12", "unsure_13", "description"))

## gettin antioxidant EC numbers
superoxide_dismutase_ec <- "1.15.1.1" #contains different SODs, not just FeSOD e.g.
ascorbate_peroxidase_ec <- "1.11.1.11"
cytochrome_c_peroxidase_ec <- "1.11.1.5"
glutathione_peroxidase_ec <- "1.11.1.9"

# peroxiredoxins
# in BRENDA, EC 1.11.1.15 entry says: "deleted. Now described by EC 1.11.1.24, thioredoxin-dependent peroxiredoxin; EC 1.11.1.25, glutaredoxin-dependent peroxiredoxin; EC 1.11.1.26, NADH-dependent peroxiredoxin; EC 1.11.1.27, glutathione-dependent peroxiredoxin; EC 1.11.1.28, lipoyl-dependent peroxiredoxin; and EC 1.11.1.29, mycoredoxin-dependent peroxiredoxin"
# but it also says " references in articles please use BRENDA:EC1.11.1.15 "
# also, the tara data only has the 1.11.1.15 number, so going forward with that one.
# peroxiredoxins <- c("1.11.1.15", "1.11.1.24", "1.11.1.25", "1.11.1.26", "1.11.1.27", "1.11.1.28", "1.11.1.29")
peroxiredoxins_ec <- c("1.11.1.15")
catalase_ec <- "1.11.1.6"

antiox_ec_all <- c(superoxide_dismutase_ec,
                   ascorbate_peroxidase_ec,
                   cytochrome_c_peroxidase_ec,
                   glutathione_peroxidase_ec,
                   peroxiredoxins_ec,
                   catalase_ec)

eek[grep(pattern = "glutathione", 
         x = eek$description, 
         ignore.case = TRUE), ]

eek[grep(pattern = "1.11.1.6", 
         x = eek$EC_number, 
         ignore.case = TRUE), ]$Gene_ID

antioxi_subset <- eek %>% 
  dplyr::filter(EC_number %in% antiox_ec_all)

write.csv(antioxi_subset, file = '../../antiox_prelim/data/antiox_subset_tara.csv',
          row.names = FALSE)


############

antioxi_subset <- read.csv("data/antiox_subset_tara.csv")

antioxi_subset$Gene_ID


