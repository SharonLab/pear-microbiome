# 01.load_data_to_ps_object_and_rarefy
# this file contains the code that ;
# 1. takes the .qza files created by qiime2 
# 1a. imports the data into phyloseq
# 1b. manipulates the metadata to fit the analyses that we will perform
# 1c. takes all the Erwiniaceae, retreive their sequences, performs blast and 
#     addition of Erwinia Amylovora according to the blast results
# 2. subset the data into phyloseq object that fit the analyses:
# 2a. All data
# 2b. Coscia 
# 2c. Coscia flowers only
# 2d. Leaves only
# 2e. Spadona and Coscia early stages ( for comparisons of the cultivars)
# 2f. Spadona alone
# 2. Creates rarefaction curves for each phylose qbject
# 3. rarefies the phyloseq objects
# 4. saves all of them as RDS



# load_packages -----------------------------------------------------------
options(getClass.msg=FALSE)

library(phyloseq)
library(dplyr)
# library(tidyr)
library(stringr)
library(tidyverse)
library(vegan)
library(readxl)
library(ggpubr)
library(qiime2R)
library(ggplot2)
library(RColorBrewer)
library(microViz)
detach("package:tidytree", unload = TRUE)

# load_data_onto_phyloseq -------------------------------------------------

ps <- qza_to_phyloseq(
  features = "Data/input_data/06.table-no-mitochondria-no-chloroplast-silva.qza",
  tree = "Data/input_data/03.rooted-tree.qza",
  "Data/input_data/04.taxonomy-silva.qza",
  metadata = "Data/input_data/metadata-filtered-copy.tsv"
)
ps
## this changes the name of the sample column to a name that does not conflict
## with phyloseq functions downstream
metadata <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
names(metadata)[names(metadata) == "Sample"] <- "Sample.id"
## correct the spelling error
metadata[metadata == "Costia"] <- "Coscia"
metadata <- data.frame(sample_data(ps)) %>%
  mutate(Phenological.stage = case_when(
    Phenological.stage.old == "E2" ~ "White buds",
    Phenological.stage.old == "F" ~ "Initial bloom",
    Phenological.stage.old == "F2-G" ~ "Full bloom",
    Phenological.stage.old == "G" ~ "Initial petal fall",
    Phenological.stage.old == "H" ~ "Complete petal fall",
    Phenological.stage.old == "Leaf" ~ "Leaf"
  )) %>%
  mutate(ps.group = case_when(
    Phenological.stage.old == "E2" | Phenological.stage.old == "F" ~ "Early",
    Phenological.stage.old == "F2-G" |
      Phenological.stage.old == "G" |
      Phenological.stage.old == "H" ~ "Late",
    Phenological.stage.old == "Leaf" ~ "Leaf"
  )) %>% 
  mutate(Group = case_when(Tree.number==7 | Tree.number==8 | Tree.number == 9 ~ "Group 1",
                           Tree.number==10 | Tree.number==11 | Tree.number == 12 ~ "Group 2",
                           Tree.number==1 | Tree.number==2 | Tree.number == 3 ~ "Group 3",
                           Tree.number==5 | Tree.number==6 | Tree.number == 4 ~ "Group 4"
                           )) %>%
  mutate(Height = case_when(Location.on.tree == "U" ~ "High",
                            Location.on.tree == "D" ~ "Low"
  ))
metadata$cd.ps <- paste0(metadata$Collection.date, " ", metadata$Phenological.stage.old)
## create a new column where we only have the day and mobth without the year
metadata$Collection.date <- as.character(metadata$Collection.date)
metadata$Short.collection.date <- gsub("\\.\\d{4}$", "", metadata$Collection.date)
metadata$CPS <- paste(metadata$Cultivar , metadata$Phenological.stage)
sample_data(ps) <- metadata

# Save_ps_into_file -------------------------------------------------------

saveRDS(ps, "Data/phyloseq_all_data")
ps <- readRDS("Data/phyloseq_all_data")

ps


# Replace the taxonomy of verified Erwinia amylovora  ---------------------

# erwinia amylovora phyloseq ----------------------------------------------

ps.e <- subset_taxa(ps, Family == "Erwiniaceae") %>%
  subset_taxa(is.na(Genus) | Genus == "Erwinia")
erwinia_tax <- data.frame(tax_table(ps.e), stringsAsFactors = FALSE)
erwinia_tax$hash <- rownames(erwinia_tax)

erwinia_otu <- data.frame(otu_table(ps.e), stringsAsFactors = FALSE)
erwinia_otu$hash <- rownames(erwinia_otu)
hash.list <- rownames(erwinia_otu)
write_delim(erwinia_tax, "Data/Erwinia/Erwinia_taxonomy_table.tsv", delim = "\t")
write_delim(erwinia_otu, "Data/Erwinia/Erwinia_otu_table.tsv", delim = "\t")

## read all rep-seq from the fasta file
library(phylotools)
rep_seqs <- read.fasta("Data/input_data/all_sequences.fasta")

## keep only the ones in the hash.list
erwinia_seqs <- subset(rep_seqs, seq.name %in% hash.list)
erwinia_seqs <- erwinia_seqs[order(match(erwinia_seqs$seq.name, erwinia_tax$hash)), ]
writeFasta <- function(data, filename) {
  fastaLines <- c()
  for (rowNum in 1:nrow(data)) {
    fastaLines <- c(fastaLines, as.character(paste(">", data[rowNum, "seq.name"], sep = "")))
    fastaLines <- c(fastaLines, as.character(data[rowNum, "seq.text"]))
  }
  fileConn <- file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}
writeFasta(erwinia_seqs, "Data/Erwinia/erwiniacea-rep-seqs.fasta")

### I manually blasted the sequences against nr and added the taxonomy to the taxonomy table.
# read the new taxonomy table
true_e <- read.delim("Data/Erwinia/Erwinia_taxonomy_table_manual.tsv")
# change it to the format needed by phyloseq
rownames(true_e) <- true_e$hash
true_e$X <- NULL
true_e$hash <- NULL
taxa_names(ps.e)
ps.e.1 <- ps.e
# change the taxonomy in the phyloseq to the new manual taxonomy
tax_table(ps.e.1) <- as.matrix(true_e)

# subset and leave only the erwinia amylovora
ps.e.2 <- subset_taxa(ps.e.1, Species == "amylovora")
otu_table(ps.e.2)
tax_table(ps.e.2)
## save the phyloseq object
saveRDS(ps.e.2, "Data/phyloseq_erwinia_amylovora")


# load the phyloseq with the manually verified EA.
# the sequences were blasted against nr and has 8 sequences
ps.e.2 <- readRDS("Data/phyloseq_erwinia_amylovora")
ps.e.2
tax_table(ps.e.2)
EA_hashes <- rownames(data.frame(tax_table(ps.e.2), stringsAsFactors = FALSE))
# extract the taxonomy table of ps
TAX <- data.frame(tax_table(ps), stringsAsFactors = FALSE)
TAX[rownames(TAX) %in% EA_hashes, ]
TAX$Genus[rownames(TAX) %in% EA_hashes] <- "Erwinia"
TAX$Species[rownames(TAX) %in% EA_hashes] <- "amylovora"
# test to see that this is corrected
TAX[rownames(TAX) %in% EA_hashes, ]

TAX2 <- TAX %>%
  # mutate_all(as.character) %>%
  as.matrix()
tax_table(ps) <- TAX2

# fix the taxonomy table to remove NAs
ps <- microViz::tax_fix(ps)

# check
tax.ps <- data.frame(tax_table(ps))
View(tax.ps)

saveRDS(ps, "Data/phyloseq_all_data")
ps <- readRDS("Data/phyloseq_all_data")
# otu_table()   OTU Table:          [ 18507 taxa and 221 samples ]:
#   sample_data() Sample Data:        [ 221 samples by 30 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 18507 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 18507 tips and 18400 internal nodes ]:
sum(sample_sums(ps)) # 4190267

# rarefaction of the entire dataset-------------------------------------------------------------

# the sample size is based on rarefaction curves from qiime, done on the filtered tables.
# this rarefied phyloseq contains all of the dataincluding Buchnera which will 
# later be removed from some further analyses. 

ps.rarefied <- rarefy_even_depth(ps, rngseed = 1, sample.size = 5400, replace = F)
# 16 samples were removed and 5015OTUs
saveRDS(ps.rarefied, "Data/rarefied.phyloseq")

ps.rarefied <- readRDS("Data/rarefied.phyloseq")


#  Remove Buchnera --------------------------------------------------------

ps <- readRDS(file = "Data/phyloseq_all_data")

ps.no.buchnera <- ps %>%
  subset_taxa(Genus != "Buchnera")

# check
ps.no.buchnera
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 18451 taxa and 221 samples ]:
#   sample_data() Sample Data:        [ 221 samples by 30 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 18451 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 18451 tips and 18344 internal nodes ]:
#   
sum(sample_sums(ps.no.buchnera)) # 4146694

hist(sample_sums(ps.no.buchnera), breaks = 50)
sorted_sample_sums <- data.frame(sample_size = sort(sample_sums(ps.no.buchnera)))

saveRDS(ps.no.buchnera, "Data/phyloseq_no_buchnera")
ps.no.buchnera <- readRDS("Data/phyloseq_no_buchnera")
ps.no.buchnera.rarefied <- rarefy_even_depth(ps.no.buchnera, sample.size = 2432, replace = FALSE)
# 8 samples and 8237 ASVs removed
ps.no.buchnera.rarefied
# otu_table()   OTU Table:          [ 10214 taxa and 213 samples ]:
#   sample_data() Sample Data:        [ 213 samples by 30 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 10214 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 10214 tips and 10130 internal nodes ]

sum(sample_sums(ps.no.buchnera.rarefied)) # 518016

saveRDS(ps.no.buchnera.rarefied, file = "Data/phyloseq_all_data_rarefied_noB")




rarefaction_curve_data.obs <- calculate_rarefaction_curves(ps.no.buchnera, c("Observed"), rep(c(50, 100, 200, 500, 1000, 2000, 3000, 4000, 7000, 10000), each = 50))
# remove the automatic X generated by r
rarefaction_curve_data.obs <- rarefaction_curve_data.obs %>% mutate(Sample = as.numeric(gsub("X", "", Sample)))

rarefaction_curve_data_summary.obs <- ddply(rarefaction_curve_data.obs, c("Depth", "Sample", "Measure"), summarise,
  Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity)
)
# Add sample data
rarefaction_curve_data_summary_verbose.obs <- merge(rarefaction_curve_data_summary.obs, data.frame(sample_data(ps.no.buchnera)),
  by.x = "Sample", by.y = "row.names"
)
# Plot
flowers <- ggplot(data = rarefaction_curve_data_summary_verbose.obs, aes(
  x = Depth, y = Alpha_diversity_mean, # axis label
  ymin = Alpha_diversity_mean - Alpha_diversity_sd,
  ymax = Alpha_diversity_mean + Alpha_diversity_sd,
  group = Sample, colour = CPS
)) +
  geom_line() +
  labs(colour = "Cultivar and \nPhenological stage") +
  ylab("Observed Features") +
  xlab(" ") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16)) +
  theme_minimal() +
  ggtitle("Rarefaction curve no Buchnera")
print(flowers)

# Re-do taxonomy ----------------------------------------------------------
# Given the subsequent use of "-ota" phyla names, Proteobacteria should be "Pseudomonadota" to be consistent. Also appears elsewhere in the manuscript.
# Aia: check if all Proteobacteria go to Pseudomonadota. Change also Firmicutes
## above are the reviewer's comments. 
## steps:
# 1. make a phyloeseq with only proteobacteria
# 2. check if they are all or nearly all pseudomonadota
# 3. change the taxonomy accordingly
# 4. apply the changes to the taxonomy of the main phyloseq
# 5. change all firmicutes to Bacillota
# 6. re-do the graphs. if it looks well, save te changes to "Data/phyloseq_all_data and re-run this script. 

ps.all <- readRDS("Data/phyloseq_no_buchnera")
# 1. make a phyloeseq with only proteobacteria
ps.proteobacteria <- ps.all %>%
  subset_taxa(Phylum == "Proteobacteria")
unique(data.frame(tax_table(ps.proteobacteria))$Class)
pseudo.tax.table <- data.frame(tax_table(ps.proteobacteria))
pseudo.tax.table[pseudo.tax.table$Class == "Proteobacteria Phylum",]
ps.proteobacteria_phylum <- subset_taxa(ps.proteobacteria, Class == "Proteobacteria Phylum")
sum(taxa_sums(ps.proteobacteria_phylum))
## this is equal to 133 reads, with 2-5 reads for each ASV. so we can also classify them as Pseudomonadota. 
tax_all <- data.frame(tax_table(ps.all))
tax_all$Phylum[tax_all$Phylum == "Proteobacteria"] <- "Pseudomonadota"
tax_all$Phylum[tax_all$Phylum == "Firmicutes"] <- "Bacillota"
tax_all$Phylum[tax_all$Phylum == "Cyanobacteria"] <- "Cyanobacteriota"
tail(taxa_names(ps.all))
tail(rownames(tax_all))
tax_table(ps.all) <- as.matrix(tax_all)
saveRDS(ps.all, "Data/phyloseq_no_buchnera")

ps.no.buchnera.rarefied <- rarefy_even_depth(ps.all, sample.size = 2432, replace = FALSE)

saveRDS(ps.no.buchnera.rarefied, file = "Data/phyloseq_all_data_rarefied_noB")


# Coscia only ----------------------------------------------------------
# remove spadona
ps.no.buchnera <- readRDS("Data/phyloseq_no_buchnera")
ps.coscia <- subset_samples(ps.no.buchnera, Cultivar == "Coscia") ## 170 samples left
# check
ps.coscia
# otu_table()   OTU Table:          [ 18451 taxa and 170 samples ]:
#   sample_data() Sample Data:        [ 170 samples by 30 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 18451 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 18451 tips and 18344 internal nodes ]:
View(sample_data(ps.coscia))
sum(sample_sums(ps.coscia))


saveRDS(ps.coscia, "Data/phyloseq_costia")

sorted_sample_sums <- data.frame(sample_size = sort(sample_sums(ps.coscia)))
ps.coscia.rarefied <- rarefy_even_depth(ps.coscia, sample.size = 6110, rngseed = 300, replace = FALSE)
ps.coscia.rarefied
# 1 samples removed, 7096 ASVs
# otu_table()   OTU Table:          [ 11355 taxa and 169 samples ]:
#   sample_data() Sample Data:        [ 169 samples by 30 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 11355 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 11355 tips and 11275 internal nodes ]:
saveRDS(ps.coscia.rarefied, "Data/phyloseq_coscia_rarefied")
sum(sample_sums(ps.coscia.rarefied)) # 1032590

## Coscia flowers only
ps.coscia.flowers <- subset_samples(ps.coscia, Flower.or.leaf == "Flower") ## 141 samples left
# otu_table()   OTU Table:          [ 18451 taxa and 141 samples ]:
#   sample_data() Sample Data:        [ 141 samples by 30 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 18451 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 18451 tips and 18344 internal nodes ]:
sum(sample_sums(ps.coscia.flowers)) # 3165831


sample_data(ps.coscia.flowers)
hist(sample_sums(ps29.coscia.flowers), breaks = 50)
sorted_sample_sums <- data.frame(sample_size = sort(sample_sums(ps.coscia.flowers)))
View(sorted_sample_sums)
# rarefying to 5900 will only remove 2 samples. so:
ps.coscia.rarefied <- rarefy_even_depth(ps.coscia.flowers, rngseed = 300, sample.size = 6110, replace = F)
# 1 samples removed and 9114 ASVs
# otu_table()   OTU Table:          [ 9337 taxa and 140 samples ]:
#   sample_data() Sample Data:        [ 140 samples by 30 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 9337 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 9337 tips and 9267 internal nodes ]:
sum(sample_sums(ps.coscia.rarefied)) # 855400
saveRDS(ps.coscia.rarefied, "Data/phyloseq_costia_flowers_rarefied")


# Coscia Leaves only
ps.coscia.leaves <- subset_samples(ps.coscia, Flower.or.leaf == "Leaf") ## 29 samples left
ps.coscia.leaves
# otu_table()   OTU Table:          [ 18451 taxa and 29 samples ]:
#   sample_data() Sample Data:        [ 29 samples by 30 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 18451 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 18451 tips and 18344 internal nodes ]:
sum(sample_sums(ps.coscia.leaves))
View(sample_data(ps.coscia.leaves))
hist(sample_sums(ps.coscia.leaves), breaks = 50)
sorted_sample_sums <- data.frame(sample_size = sort(sample_sums(ps.coscia.leaves)))
sorted_sample_sums
# rarefying to 5900 will only remove 2 samples. so:
ps.coscia.rarefied <- rarefy_even_depth(ps.coscia.leaves, rngseed = 300, sample.size = 11070, replace = F)
sample_data(ps.coscia.rarefied)
# 3 samples removed and 13668 ASVs

# otu_table()   OTU Table:          [ 4783 taxa and 26 samples ]:
#   sample_data() Sample Data:        [ 26 samples by 30 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 4783 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 4783 tips and 4742 internal nodes ]:
sum(sample_sums(ps.coscia.rarefied)) # 287820
saveRDS(ps.coscia.rarefied, "Data/phyloseq_costia_leaves_rarefied")



# Early stages ------------------------------------------------------------
# this data set is used for comparing the sapdona and coscia cultivars
ps.initial <- ps.no.buchnera %>% 
  subset_samples(Phenological.stage.old == "E2" | Phenological.stage.old == "F") %>%
  subset_samples(Collection.date != "15.04.2019")
ps.initial
# otu_table()   OTU Table:          [ 18451 taxa and 81 samples ]:
#   sample_data() Sample Data:        [ 81 samples by 30 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 18451 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 18451 tips and 18344 internal nodes ]:
sum(sample_sums(ps.initial)) # 958252
hist(sample_sums(ps.initial), breaks = 50)
sorted_sample_sums <- data.frame(sample_size = sort(sample_sums(ps.initial)))
sorted_sample_sums
ps.initial.rarefied <- rarefy_even_depth(ps.initial, rngseed = 300, sample.size = 2530, replace = F)
# 6 samples removed and 12322OTUs, 768502 reads removed leaving 189750 reads
ps.initial.rarefied
# otu_table()   OTU Table:          [ 6129 taxa and 75 samples ]:
#   sample_data() Sample Data:        [ 75 samples by 30 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 6129 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 6129 tips and 6077 internal nodes ]:
sum(sample_sums(ps.initial.rarefied)) # 189750
sum(sample_sums(ps.initial))
saveRDS(ps.initial.rarefied, "Data/phyloseq_initial_stages_rarefied")



#  Leaves from both cultivars-----------------------------------------------------------------

ps.leaves <- ps.no.buchnera %>% 
  subset_samples(Flower.or.leaf == "Leaf")
ps.leaves  ## 46 samples
sum(sample_sums(ps.leaves))
hist(sample_sums(ps.leaves), breaks = 30)
sort(sample_sums(ps.leaves))

ps.leaves.rarefied <- rarefy_even_depth(ps.leaves, rngseed = 300, sample.size = 5250, replace = F)
# 6 samples removed and 12669OTUs
ps.leaves.rarefied
# otu_table()   OTU Table:          [ 5782 taxa and 40 samples ]:
#   sample_data() Sample Data:        [ 40 samples by 30 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 5782 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 5782 tips and 5733 internal nodes ]:
sum(sample_sums(ps.leaves.rarefied)) # 210000
sum(sample_sums(ps.leaves))   # 725847

saveRDS(ps.leaves, file = "Data/phyloseq_leaves")
saveRDS(ps.leaves.rarefied, file = "Data/phyloseq_leaves_rarefied")


#  Spadona ----------------------------------------------------------------

# remove Coscia
ps.spadona <- subset_samples(ps.no.buchnera, Cultivar == "Spadona") ## 51 samples left
ps.spadona
# otu_table()   OTU Table:          [ 18451 taxa and 51 samples ]:
#   sample_data() Sample Data:        [ 51 samples by 30 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 18451 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 18451 tips and 18344 internal nodes ]:
sum(sample_sums(ps.spadona)) # 418619

saveRDS(ps.spadona, "Data/phyloseq_spadona")
sorted_sample_sums <- data.frame(sample_size = sort(sample_sums(ps.spadona)))
sorted_sample_sums
hist(sample_sums(ps.spadona), breaks = 20)
ps.spadona.rarefied <- rarefy_even_depth(ps.spadona, sample.size = 2400, rngseed = 300, replace = FALSE)
# 8 samples removed, 13873 ASVs
ps.spadona.rarefied
# otu_table()   OTU Table:          [ 4578 taxa and 43 samples ]:
#   sample_data() Sample Data:        [ 43 samples by 30 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 4578 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 4578 tips and 4539 internal nodes ]:
sum(sample_sums(ps.spadona.rarefied)) # 103200


saveRDS(ps.spadona.rarefied, "Data/phyloseq_spadona_rarefied")

## spadona flowers only
ps.spadona.flowers <- subset_samples(ps.spadona, Flower.or.leaf == "Flower") ## 34 samples left
ps.spadona.flowers
# otu_table()   OTU Table:          [ 18451 taxa and 34 samples ]:
#   sample_data() Sample Data:        [ 34 samples by 30 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 18451 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 18451 tips and 18344 internal nodes ]:
hist(sample_sums(ps.spadona.flowers), breaks = 30)
View(sample_data(ps.spadona.flowers))
sorted_sample_sums <- data.frame(sample_size = sort(sample_sums(ps.spadona.flowers)))
sorted_sample_sums
# rarefying to  will only remove 2 samples. so:
ps.spadona.rarefied <- rarefy_even_depth(ps.spadona.flowers, rngseed = 300, sample.size = 2400, replace = F)
# 4 samples removed and 14906 OTUs
ps.spadona.rarefied
# otu_table()   OTU Table:          [ 3545 taxa and 30 samples ]:
#   sample_data() Sample Data:        [ 30 samples by 30 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 3545 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 3545 tips and 3514 internal nodes ]:
sum(sample_sums(ps.spadona.flowers)) # 255016
sum(sample_sums(ps.spadona.rarefied)) # 72000
saveRDS(ps.spadona.rarefied, "Data/phyloseq_spadona_flowers_rarefied")


# spadona Leaves only
ps.spadona.leaves <- subset_samples(ps.spadona, Flower.or.leaf == "Leaf") ## 17 samples left
ps.spadona.leaves
# otu_table()   OTU Table:          [ 18451 taxa and 17 samples ]:
#   sample_data() Sample Data:        [ 17 samples by 30 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 18451 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 18451 tips and 18344 internal nodes ]:
hist(sample_sums(ps.spadona.leaves), breaks = 30)
sorted_sample_sums <- data.frame(sample_size = sort(sample_sums(ps.spadona.leaves)))
sorted_sample_sums
sum(sample_sums(ps.spadona.leaves)) # 163603

ps.spadona.rarefied <- rarefy_even_depth(ps.spadona.leaves, rngseed = 300, sample.size = 3400, replace = F)
ps.spadona.rarefied
# 4 samples removed and 15679 OTUs
# otu_table()   OTU Table:          [ 2772 taxa and 13 samples ]:
#   sample_data() Sample Data:        [ 13 samples by 30 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 2772 taxa by 7 taxonomic ranks ]:
#   phy_tree()    Phylogenetic Tree:  [ 2772 tips and 2759 internal nodes ]:
sum(sample_sums(ps.spadona.rarefied)) # 44200
saveRDS(ps.spadona.rarefied, "Data/phyloseq_spadona_leaves_rarefied")


#  Look at the Buchnera samples -------------------------------------------
# start with the non-rarefied data
ps <- readRDS("Data/phyloseq_all_data")
ps.buchnera.all <- ps %>% 
  subset_taxa(Genus == "Buchnera")
ps.buchnera.all # 56 ASVs
# check that we really have the right sequences
tax_table(ps.buchnera.all)
taxa_sums(ps.buchnera.all)
# we want to look at the taxonomy of the high abundance Buchnera
# so we subset the phyloseq to remove the ones with the low counts

ps.buchnera <- ps.buchnera.all %>% 
  subset_taxa(taxa_sums(ps.buchnera.all) > 100)
# this leaves us with the 7 most abundant taxa
# check:
tax_table(ps.buchnera)
taxa_sums(ps.buchnera)

## export the taxonomy tables, the otu tables and the sequences in fasta format. 
buchnera_tax <- data.frame(tax_table(ps.buchnera), stringsAsFactors = FALSE)
buchnera_tax$hash <- rownames(buchnera_tax)

buchnera_otu <- data.frame(otu_table(ps.buchnera), stringsAsFactors = FALSE)
buchnera_otu$hash <- rownames(buchnera_otu)
hash.list <- rownames(buchnera_otu)
write_delim(buchnera_tax, "Data/Buchnera/buchnera_taxonomy_table.tsv", delim = "\t")
write_delim(buchnera_otu, "Data/Buchnera/buchnera_otu_table.tsv", delim = "\t")

## read all rep-seq from the fasta file
library(phylotools)
rep_seqs <- read.fasta("Data/input_data/all_sequences.fasta")

## keep only the ones in the hash.list
buchnera_seqs <- subset(rep_seqs, seq.name %in% hash.list)
buchnera_seqs <- buchnera_seqs[order(match(buchnera_seqs$seq.name, buchnera_tax$hash)), ]
writeFasta <- function(data, filename) {
  fastaLines <- c()
  for (rowNum in 1:nrow(data)) {
    fastaLines <- c(fastaLines, as.character(paste(">", data[rowNum, "seq.name"], sep = "")))
    fastaLines <- c(fastaLines, as.character(data[rowNum, "seq.text"]))
  }
  fileConn <- file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}
writeFasta(buchnera_seqs, "Data/Buchnera/buchnera-rep-seqs.fasta")
## results:
## All seven sequences have 100% identity to B. aphidicola species in NCBI
## except for the sequence starting with ac305 who had 270 reads and has 99.7% 
## identity to B. aphidicola

# look at the samples where we have b. aphidicola
zero_samples <- data.frame(samples = sample_sums(ps.buchnera)[sample_sums(ps.buchnera) == 0])
# check that the numbers make sense
length(zero_samples$samples)
zero_samples$name <- rownames(zero_samples)

ps.buchnera <- microbiome::remove_samples(zero_samples$name, ps.buchnera)
ps.buchnera
# this leaves 7ASVs in 51 samples
my_palette <- colorRampPalette(brewer.pal(9, "Set1"))(7)
plot_ordination_pears(ps.buchnera, OrdMeth = "bray", metadataColumn = "CPS", titleAddition = "CPS")
p1 <- plot_bar(ps.buchnera, x = "CPS", fill = "OTU", facet_grid = ~ Collection.date) +
  scale_fill_manual(values =  my_palette)  
png("images/buchnera_abundances.png", height = 3000, width = 3000, res = 0)
print(p1)
dev.off()


p2 <- plot_bar(ps.buchnera, x = "Collection.date", fill = "OTU") +
  scale_fill_manual(values =  my_palette) 
  



# count samples by collection dates and phenological stages ---------------
metadata.rarefied <- data.frame(sample_data(ps.rarefied), stringsAsFactors = FALSE)
sample_counts <- metadata.rarefied %>%
  group_by(Collection.date, Phenological.stage) %>%
  summarise(n_samples = n()) %>%
  tidyr::pivot_wider(names_from = Collection.date, values_from = n_samples)
data_bars <- make_data_bars1(tax_glom(ps.rarefied, taxrank = "Phylum"), 0.000001, "CPS")
taxa_bar_plot(data_bars = data_bars)

# rarefaction function ----------------------------------------------------


# Code from https://github.com/joey711/phyloseq/issues/143
set.seed(42)

calculate_rarefaction_curves <- function(psdata, measures, depths) {
  library("plyr") # ldply
  library("reshape2") # melt

  estimate_rarified_richness <- function(psdata, measures, depth) {
    if (max(sample_sums(psdata)) < depth) {
      return()
    }
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)

    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)

    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)

    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c("Sample", "Measure"), value.name = "Alpha_diversity")

    molten_alpha_diversity
  }

  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = "Depth", .progress = ifelse(interactive(), "text", "none"))

  # convert Depth from factor to numeric

  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]

  rarefaction_curve_data
}

# check buchnera ----------------------------------------------------------
library(ggtree)
ps.e <- subset_taxa(ps, Genus == "Buchnera")
buchnera.tree <- phyloseq::phy_tree(ps.e)
ggtree::ggtree(buchnera.tree) + theme_tree2() +
  geom_tiplab(align = TRUE, linesize = .5)
erwinia_tax <- data.frame(tax_table(ps.e))
erwinia_tax$hash <- rownames(erwinia_tax)

erwinia_otu <- data.frame(otu_table(ps.e, taxa_are_rows = TRUE))
erwinia_otu$hash <- rownames(erwinia_otu)
hash.list <- rownames(erwinia_otu)
## read all rep-seq from the fasta file usint the phytools package
rep_seqs <- phylotools::read.fasta("all_sequences.fasta")

## keep only the ones in the hash.list
erwinia_seqs <- subset(rep_seqs, seq.name %in% hash.list)
erwinia_seqs <- erwinia_seqs[order(match(erwinia_seqs$seq.name, erwinia_tax$hash)), ]

writeFasta(erwinia_seqs, "buchnera-rep-seqs.fasta")


