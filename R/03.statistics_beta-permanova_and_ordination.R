## PERMANOVA using vegan::adonis function, on bray curtis beta matrix
#

# load_packages -----------------------------------------------------------
library(phyloseq)
library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)
library(vegan)
library(ggplot2)
library(writexl)
library(ggpubr)

setwd("~/qiime/fire-blight/pear.flower/20200125/qiime2.aia.filter/R_pears")
## tested on 13-3-2023 till line 258
## re-tested on 8-8-23
## the different tests are performed on phyloseqs with different rarefaction depths
## as described in file 01.load_data_to_ps_object_rarefy
## in order to preserve as much data as possible

## within coscia checks
ps.coscia.rarefied <- readRDS("Data/phyloseq_coscia_rarefied")
## within coscia flowers
ps.coscia.flowers <- readRDS("Data/phyloseq_costia_flowers_rarefied")
## comparison of Coscia and Spadona initial stages
ps.init <- readRDS("Data/phyloseq_initial_stages_rarefied")
ps.coscia.leaves <- readRDS("Data/phyloseq_costia_leaves_rarefied")


# statistical testing -----------------------------------------------------------------

# transform phyloseq to vegan format
ps_otu2veg <- function(physeq) {
  library(phyloseq)
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  OTU <- as.data.frame(OTU)
  return(OTU)
}


# function that arranges results ------------------------------------------
# this function takes the adonis results and inserts them into a
# data frame that can be saved at the end of the script.

insert_result_to_df <- function(r1) {
  # take the adonis results and parses it into a data frame.
  res <- NULL
  res$call <- as.data.frame(paste(as.character(result$call), collapse = " , "))
  x <- result$aov.tab
  res$test <- rownames(as.data.frame(x))
  res$DF <- as.data.frame(result$aov.tab$Df)
  res$SumsOfSqs <- as.data.frame(result$aov.tab$SumsOfSqs)
  res$MeanSqs <- as.data.frame(result$aov.tab$MeanSqs)
  res$FModel <- as.data.frame(result$aov.tab$F.Model)
  res$R2 <- as.data.frame(result$aov.tab$R2)
  res$Pvalue <- as.data.frame(result$aov.tab$`Pr(>F)`)
  res <- as.data.frame(res)
}


# Coscia amongst itself ---------------------------------------------------
# each test need a vegan type abundance table and the corresponding metadata
# note that the ps_otu2veg function is a little slow
# we need the general Coscia ps to compare the flowers to the leaves on the same collection date

metadata.coscia <- data.frame(sample_data(ps.coscia.rarefied), stringsAsFactors = FALSE)
vegan.coscia <- ps_otu2veg(ps.coscia.rarefied) ## this function takes some time to run


vegan.ps.flowers <- ps_otu2veg(ps.coscia.flowers)
metadata.flowers <- data.frame(sample_data(ps.coscia.flowers), stringsAsFactors = FALSE)


vegan.ps.coscia.leaves <- ps_otu2veg(ps.coscia.leaves)
metadata.coscia.leaves <- data.frame(sample_data(ps.coscia.leaves), stringsAsFactors = FALSE)


# choose beta matrix ------------------------------------------------------
# we need a beta matrix for each test.
method <- "bray"
beta.otu.coscia <- vegdist(vegan.coscia, method = "bray")
beta.otu.flowers <- vegdist(vegan.ps.flowers, method = "bray")
beta.otu.leaves <- vegdist(vegan.ps.coscia.leaves, method = "bray")


# adonis tests ------------------------------------------------------------
## flowers vs. leaves - same collection dates
rm(adonis_results)
metadata.coscia$Flower.vs.leaf.vs.collection.date <- paste0(metadata.coscia$Flower.or.leaf, "_", metadata.coscia$Collection.date)
result <- adonis(beta.otu.coscia ~ Flower.vs.leaf.vs.collection.date, metadata.coscia, perm = 1000)
temp_res <- insert_result_to_df(result)
adonis_results <- temp_res

result <- adonis(beta.otu.coscia ~ Flower.or.leaf + Collection.date, metadata.coscia, perm = 1000)
temp_res <- insert_result_to_df(result)
adonis_results <- temp_res
# leaves only
# leaves vs. collection dates
result <- adonis(beta.otu.leaves ~ Collection.date, metadata.coscia.leaves, perm = 1000)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)

# leaves: location in the orchard
metadata.coscia.leaves
result <- adonis(beta.otu.leaves ~ Group, metadata.coscia.leaves, perm = 1000, strata = metadata.coscia.leaves$Collection.date)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)

## flowers coscia
## collection date, phenological stage as strata
result <- adonis(beta.otu.flowers ~ Collection.date, metadata.flowers, perm = 1000, strata = metadata.flowers$Phenological.stage)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)


result <- adonis(beta.otu.flowers ~ Phenological.stage, metadata.flowers, perm = 1000, strata = metadata.flowers$ollection.date)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)

metadata.flowers$CD.PS <- with(metadata.flowers, interaction(Collection.date, Phenological.stage))
result <- adonis(beta.otu.flowers ~ Group, metadata.flowers, perm = 1000, strata = metadata.flowers$CD.PS)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)

result <- adonis(beta.otu.flowers ~ Location.on.tree, metadata.flowers, perm = 1000, strata = metadata.flowers$CD.PS)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)

result <- adonis(formula = beta.otu.flowers ~ Collection.date + Phenological.stage, data = metadata.flowers, permutations = 1000)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)

result <- adonis(formula = beta.otu.flowers ~ Collection.date * Phenological.stage, data = metadata.flowers, permutations = 1000)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)

result <- adonis(formula = beta.otu.flowers ~ Phenological.stage * Collection.date, data = metadata.flowers, permutations = 1000)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)

# the effect of cultivar --------------------------------------------------
# flowers of both cultivars, first phenological stages only
metadata.ps.init <- data.frame(sample_data(ps.init), stringsAsFactors = FALSE)
vegan.ps.init <- ps_otu2veg(ps.init)
beta.ps.init <- vegdist(vegan.ps.init, method = "bray")


result <- adonis(beta.ps.init ~ Cultivar, metadata.ps.init, perm = 1000, strata = metadata.ps.init$Phenological.stage)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)


result <- adonis(formula = beta.ps.init ~ Phenological.stage + Cultivar, data = metadata.ps.init, permutations = 1000)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)

# add a combined collection date and phenological stage block
metadata.ps.init$CD.CU <- with(metadata.ps.init, interaction(Collection.date, Cultivar))
adonis(beta.ps.init ~ CD.CU, data = metadata.ps.init, permutations = 1000, strata = metadata.ps.init$Phenological.stage)

homo <- with(metadata.ps.init, betadisper(beta.ps.init, CD.CU))
homo
plot(homo)
boxplot(homo)
betadisper(beta.ps.init, metadata.ps.init$CD.CU)
permutest(homo)
# test the effect of cultivar on leaves -----------------------------------
# I think that the results will not be valid since the collection dates were different.


# spadona flowers ---------------------------------------------------------

ps.spadona.flowers <- readRDS("Data/phyloseq_spadona_flowers_rarefied")
metadata.spadona.flowers <- data.frame(sample_data(ps.spadona.flowers), stringsAsFactors = FALSE)
vegan.spadona.flowers <- ps_otu2veg(ps.spadona.flowers)
beta.sf <- vegdist(vegan.spadona.flowers, method = "bray")
metadata.spadona.flowers$CD.PS <- with(metadata.spadona.flowers, interaction(Collection.date, Phenological.stage))
metadata.spadona.flowers
# location in the orchard
result <- adonis(beta.sf ~ Group, metadata.spadona.flowers, perm = 1000, strata = metadata.spadona.flowers$CD.PS)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)

# effect of phenological stage, same dates
result <- adonis(beta.sf ~ Collection.date, metadata.spadona.flowers, perm = 1000, strata = metadata.spadona.flowers$Phenological.stage)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)

# effect of collection date, same phenological stage. this is for one phenological stage (F?)
result <- adonis(beta.sf ~ Phenological.stage, metadata.spadona.flowers, perm = 1000)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)

result <- adonis(beta.sf ~ CD.PS, metadata.spadona.flowers, perm = 1000)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)

result <- adonis(beta.sf ~ Collection.date * Phenological.stage, metadata.spadona.flowers, perm = 1000)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)

result <- adonis(beta.sf ~ Phenological.stage * Collection.date, metadata.spadona.flowers, perm = 1000)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)

# spadona leaves ----------------------------------------------------------
ps.spadona <- readRDS(file = "Data/phyloseq_spadona_rarefied")
ps.spadona.leaves <- ps.spadona %>%
  subset_samples(Flower.or.leaf == "Leaf")
metadata.spadona.leaves <- data.frame(sample_data(ps.spadona.leaves))
vegan.spadona.leaves <- ps_otu2veg(ps.spadona.leaves)
beta.sl <- vegdist(vegan.spadona.leaves, method = "bray")

metadata.spadona.leaves$CD.PS <- with(metadata.spadona.leaves, interaction(Collection.date, Phenological.stage))
metadata.spadona.leaves
# spadona leaves by date
result <- adonis(beta.sl ~ Collection.date, metadata.spadona.leaves, perm = 1000)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)

# spadona leaves by location in the orchard
result <- adonis(beta.sl ~ Group, metadata.spadona.leaves, perm = 1000, strata = metadata.spadona.leaves$Collection.date)
result
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)


# collection date on leaves from both cultivars ---------------------------
ps.leaves <- readRDS("Data/phyloseq_leaves_rarefied")
metadata.leaves <- data.frame(sample_data(ps.leaves))
vegan.leaves <- ps_otu2veg(ps.leaves)
beta.l <- vegdist(vegan.leaves, method = "bray")
metadata.leaves
result <- adonis(beta.l ~ Collection.date, metadata.leaves, perm = 1000)
temp_res <- insert_result_to_df(result)
adonis_results <- rbind(adonis_results, temp_res)


# write to txt ------------------------------------------------------
colnames(adonis_results) <- c("model", "test", "DF", "SumOfSquares", "MeanSquare", "F", "RSquare", "PValue")
write.table(adonis_results, file = "PERMANOVA_pears_6th_edition_no_bucnera.txt", sep = "\t")

# test group homogeneities using betadisper -------------------------------


homo <- with(metadata.flowers, betadisper(beta.otu.flowers, Collection.date))
homo
plot(homo)
boxplot(homo)
# the anova and permutation tests check if multivariate dispersions are different
anova(homo)
permutest(homo)
TukeyHSD(homo) # tests the differences in the variances amongst these groups

homo <- with(metadata.flowers, betadisper(beta.otu.flowers, Phenological.stage))
homo
plot(homo)
boxplot(homo)
# the anova and permutation tests check if multivariate dispersions are different
anova(homo)
permutest(homo)
TukeyHSD(homo) # tests the differences in the variances amongst these groups


homo <- with(metadata.coscia.leaves, betadisper(beta.otu.leaves, Collection.date))
homo
plot(homo)
boxplot(homo)
# the anova and permutation tests check if multivariate dispersions are different
anova(homo)


homo <- with(metadata.spadona.leaves, betadisper(beta.sl, Collection.date))
homo
plot(homo)
boxplot(homo)
# the anova and permutation tests check if multivariate dispersions are different
anova(homo)


homo <- with(metadata.ps.init, betadisper(beta.ps.init, Collection.date))
homo
plot(homo)
boxplot(homo)
# the anova and permutation tests check if multivariate dispersions are different
anova(homo)
