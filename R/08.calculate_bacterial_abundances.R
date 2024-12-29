# calculate bacteria abundances, as well as the number of reads in each sample, the percentages of mitochondria and chloroplast  the number of microbiome reads
# in each sample and other statistics.
library(phyloseq)
library(ggplot2)
library(dplyr)
library(qiime2R)
library(purrr)
library(tidyverse)
source("R/02.general_definitions.R")
source("R/02.helper_functions.R")

metadata <- read.csv( "Data/metadata_with_chloroplast_reads.csv")


# number of bacterial flowers in the literature: --------------------------

# pusey 2009 10^2 - 10^6
# the graph goes up then down. first you have little, it peaks after 3 days then drops in 1994. apples. 
# in 1995 the pattern is similar to what we see. went up slower or remained the same. 1995 was colder and rainier than 1994
# 

# mitochondria ------------------------------------------------------------


n_mito <- 400 
# arabidopsis has 200-300 mito per cell
# tobacco has 500-600 mito per cell

n_mito_genomes <- 50

size_mito_genome <- 459000

# Fewer genes than organelles: Extremely low and variable gene copy numbers in mitochondria of somatic plant cells
# they show that plants have 50-500 copies per cell but that flowers have 50 in both tobacco and arabidopsis.
# they claim that some mitochondria do not have genomes!

# https://www.nature.com/articles/s41438-019-0177-8
# Copy numbers of mitochondrial genes change during melon leaf development and are lower than the numbers of mitochondria
# The gene copy numbers varied from approximately 36.8 ± 4.5 (atp6 copies in the 15th leaf) 
# to approximately 82.9 ± 5.7 (nad9 copies in the 9th leaf), 
# while the mean number of mitochondria was approximately 416.6 ± 182.7 in the 15th leaf and
# 459.1 ± 228.2 in the 9th leaf.

# chloroplasts ------------------------------------------------------------


n_chloroplasts <- mean(3:200)
n_chloro_genome <- mean(10:400) * n_chloroplasts
size_of_chloroplast_genome <- 160000
# the size is for pear

# Chloroplast genome copy numbers per cell are highest in young photosynthetically active leaves. 
# Chloroplast genome copy number varies widely between tissues, 
# ranging from 3 to 275 copies per plastid in leaf cells of different developmental stages (Zoschke et al., 2007; Liere and Borner, 2013). 
# For other species there are 10–400 copies of the chloroplast genome per plastid,
# translating to 1000–50,000 genome copies per plant cell (also see Boffey and Leech, 1982; Tymms et al., 1983).

## as for the numbers of chloroplasts: these range from 3 to 200 in a cell. 


# pear genome -------------------------------------------------------------
# pyrus communis genome size. in bp. 
size_pear_genome <- 577000000

# 16S copies --------------------------------------------------------------

n_16S <- 6 # assuming 6 16S per bacterial genome

bacterial_reads <-  mean(! is.na(metadata$reads_after_filtration) - ! is.na(metadata$number_of_chloro_mito))





# genome_sizes ------------------------------------------------------------
bacterial_genome_size <- 3
pear_genome_size <- 527
expected_ratio <- bacterial_genome_size / pear_genome_size
# model -------------------------------------------------------------------

one_plant_genome_expected <- n_mito_genomes + n_chloro_genome 
plant_genomes <- mean(metadata$number_of_chloro_mito) / one_plant_genome
bacterial_genomes <- bacterial_reads / n_16S

bact_to_pear_genome_ratio <- bacterial_genomes / plant_genomes

expected_to_seen_ratio <- bact_to_pear_genome_ratio / expected_ratio

mean_DNA <- mean(metadata$concentration)



# for each CPS separately -------------------------------------------------
# change the format to a list of tibbles by CPS
metadata_nested <- metadata %>%
  group_by(CPS) %>%
  nest()

# I want to find the expected_to_seen_ratio for each CPS separately

# see if the literature agrees with our data ------------------------------

# does the ratio of mito/chloro stay the same? seems so
# the mean ratio of mito/chloro is 0.14
# which means that for every mitochondria genome there are 7 chloroplast genomes
# 

metadata$mito_over_chloro <- metadata$mito_reads/metadata$chloro_reads
ggplot(metadata, aes(y = mito_over_chloro, x = CPS)) +
  geom_boxplot() + 
  scale_x_discrete(limits = positions_cps)


display_stats( metadata$mito_over_chloro)
# Mean: 0.1423537 
# Median: 0.1447555 
# Standard Deviation: 0.05847977 

plot(density(metadata$mito_over_chloro), main = "Density Plot of Numbers", xlab = "Value", ylab = "Density", col = "blue")



# helper functions --------------------------------------------------------

display_stats <- function(x) {
  mean_value <- mean(x)
  median_value <- median(x)
  sd_value <- sd(x)
  
  cat("Mean:", mean_value, "\n")
  cat("Median:", median_value, "\n")
  cat("Standard Deviation:", sd_value, "\n")
}
