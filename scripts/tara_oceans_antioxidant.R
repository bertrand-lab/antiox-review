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

library(readxl)
library(magrittr)

meta_data_genomes <- read_excel('data/tara_ocean_smags/Table_S03_statistics_nr_SMAGs_METdb.xlsx',
                                sheet = 1,
                                skip = 2)

# subsetting only the photosynthetic phytoplankton
genomes_photosynthetic_w_na <- meta_data_genomes[meta_data_genomes$Phytoplankton == "Phytoplankton", ]$`Genome_Id final names`
genomes_photosynthetic <- genomes_photosynthetic_w_na[!is.na(genomes_photosynthetic_w_na)]

matches <- unique (grep(paste(toMatch,collapse="|"), 
                        myfile$Letter, value=TRUE))

# getting the genes that are antixodiants (above)
antioxi_subset <- read.csv("data/antiox_subset_tara.csv")

# subsetting the proteins that belong to photosytnehtic organisms
grep(pattern = "TARA_AON_82_MAG_00189", 
     x = antioxi_subset$Gene_ID)

antioxi_subset_photo <- dplyr::filter(antioxi_subset, 
                              grepl(paste(genomes_photosynthetic, 
                                     collapse="|"), Gene_ID))

## making txt files for each of the lists of gene names
sod_gene_id <- antioxi_subset_photo[antioxi_subset_photo$EC_number == superoxide_dismutase_ec, ]$Gene_ID %>% gsub(pattern = "mRNA.", replacement = "")
apx_gene_id <- antioxi_subset_photo[antioxi_subset_photo$EC_number == ascorbate_peroxidase_ec, ]$Gene_ID %>% gsub(pattern = "mRNA.", replacement = "")
cyto_c_gene_id <- antioxi_subset_photo[antioxi_subset_photo$EC_number == cytochrome_c_peroxidase_ec, ]$Gene_ID %>% gsub(pattern = "mRNA.", replacement = "")
gpx_gene_id <- antioxi_subset_photo[antioxi_subset_photo$EC_number == glutathione_peroxidase_ec, ]$Gene_ID %>% gsub(pattern = "mRNA.", replacement = "")
prx_gene_id <- antioxi_subset_photo[antioxi_subset_photo$EC_number == peroxiredoxins_ec, ]$Gene_ID %>% gsub(pattern = "mRNA.", replacement = "")
cat_gene_id <- antioxi_subset_photo[antioxi_subset_photo$EC_number == catalase_ec, ]$Gene_ID %>% gsub(pattern = "mRNA.", replacement = "")

# getting gene lists for getting their actual protein sequences
write.table(sod_gene_id, file = "data/antiox_gene_name_lists/sod_gene_id.txt", 
            row.names = FALSE, 
            quote = FALSE)
write.table(apx_gene_id, file = "data/antiox_gene_name_lists/APX_gene_id.txt", 
            row.names = FALSE, 
            quote = FALSE)
write.table(cyto_c_gene_id, file = "data/antiox_gene_name_lists/cyto_c_gene_id.txt", 
            row.names = FALSE, 
            quote = FALSE)
write.table(gpx_gene_id, file = "data/antiox_gene_name_lists/gpx_gene_id.txt", 
            row.names = FALSE, 
            quote = FALSE)
write.table(cat_gene_id, file = "data/antiox_gene_name_lists/cat_gene_id.txt", 
            row.names = FALSE, 
            quote = FALSE)
write.table(prx_gene_id, file = "data/antiox_gene_name_lists/prx_gene_id.txt", 
            row.names = FALSE, 
            quote = FALSE)

