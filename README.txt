`scripts` includes all the scripts for analyzing the mass spectrometry data, and doing the Monte Carlo estimates, and making all the figures.

These scripts assume the following file structure

```
/scripts/
/data/
     /tara_ocean_smags/
     /antiox_gene_name_lists/
     /protein_expression_data/
                             /phytoplankton_genomes/
                                                   /frag_genome/
                                                   /phaeo_genome/
                             /mass_spec_data/
                                            /frag_data/
                                                      /mzml_converted/
                                                      /fdr_idxml/
                                            /phaeo_data_3/
                                                      /mzml_converted/
                                                      /fdr_idxml/

```
### Data

For the `frag_data`, I got data from PRIDE project with the ID: PXD007098. I specifically used these data files:

```
D05_1.raw
D120_1.raw
D1_2.raw
L120_3.raw
T0_3.raw
```

For the `phaeo_data_3`, I got data from the PRIDE project with the ID: PXD014877. I used the script `download_phaeo_proteome_3.sh` for this.

Data for the phytoplankton genomes for protein stoichiometry were downloaded from this awesome paper by Delmont et al (2021)[https://www.biorxiv.org/content/10.1101/2020.10.15.341214v2] specifically with:

```
# annotation files
wget https://www.genoscope.cns.fr/tara/localdata/data/SMAGs-v1/SMAGs_v1_EggNog.tar.gz

# protein sequences
wget https://www.genoscope.cns.fr/tara/localdata/data/SMAGs-v1/SMAGs_v1_concat.faa.tar.gz

```

Other files ('Table_S03_statistics_nr_SMAGs_METdb.xlsx') were manually downloaded from https://www.genoscope.cns.fr/tara/, and the ferritin sequence used is from https://www.uniprot.org/uniprot/B6DMH6.fasta
.

### Genomic Data Processing

`converting_tara_to_single_line_fasta.sh`

Converts the large file of SMAGs to a single line fasta file.

`tara_oceans_antioxidant.R`

Script that gets all the antioxidant protein sequence IDs and makes txt files of lists, for eventually subsetting the large fasta file.

`antioxidant_stoich_from_seqs.ipynb`

Jupyter notebook that calculates stoichiometric composition of various phytoplankton proteins.

### Proteomic Data Processing

Once the raw MS data have been downloaded, they need to be converted to mzML. They were converted with `convert_raw_to_mzml.sh`, `frag_converting_file.sh`, and `phaeo_converting_file.sh`.

Then they are searched using the corresponding genomes with `database-searching-openms.sh`, `database-searching-frag.sh`, and `database-searching-phaeo.sh`. Before they are searched the genomes are appended with a contaminant database (CRAP), with `adding_crap_to_genomes.sh`.

Peptides are quantified using the FeatureFinderIdentification method (Weisser et al 2017) with the script `feature-finder-general.sh`, `feature-finder-phyto-proteomes.sh`, and `feature-finder-phyto-proteomes-phaeo.sh`.

The output of this is then converted to csv with `protein-quant-openms.sh`, `protein-quant-frag.sh`, and `protein-quant-phaeo.sh`.

### Monte Carlo Estimation



