library(phyloseq)
library(dplyr)
library(tidyr)
library(reshape2)

ps.rarefied <- readRDS(file = "Data/phyloseq_all_data_rarefied_noB")

metadata <- as.data.frame(sample_data(ps.rarefied))
metadata$cd.ps <- paste0(metadata$Collection.date, "_", metadata$Phenological.stage)
sample_data(ps.rarefied) <- metadata

data_bars <- ps.rarefied %>%
  tax_glom(taxrank = "Phylum") %>%
  speedyseq::merge_samples2(group = "cd.ps", fun_otu = mean) %>%
  transform_sample_counts(function(x) {
    x / sum(x)
  }) %>% # Transform to rel. abundance
  psmelt()
data_bars$cd.ps.short <- gsub("\\.\\d{4}", "", data_bars$Sample)
phyla_table <- reshape2::dcast(data_bars, cd.ps.short ~ Phylum, value.var = "Abundance", fun.aggregate = mean)

phyla_table <- phyla_table %>%
  separate(cd.ps.short, into = c("Date", "Stage"), sep = "_", remove = FALSE, extra = "merge")

date_order <- c("24.03", "27.03", "31.03", "11.04", "13.04", "15.04", "18.04", "21.04")
stage_order <- c("White buds", "Initial bloom", "Full bloom", "Initial petal fall", "Complete petal fall", "Leaf")

phyla_table$Date <- factor(phyla_table$Date, levels = date_order)
phyla_table$Stage <- factor(phyla_table$Stage, levels = stage_order)

phyla_table <- phyla_table %>%
  arrange(Date, Stage)
rownames(phyla_table) <- phyla_table$cd.ps.short
phyla_table <- subset(phyla_table, select = -c(1:3))
# arrange the phyla by abundances
phyla_table <- phyla_table[, order(-colSums(phyla_table, na.rm = TRUE))]
phyla_table[1:ncol(phyla_table)] <- round(phyla_table[1:ncol(phyla_table)], 5)

write.table(phyla_table, file = "images/phyla_table.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
vars.to.remove <- c("data_bars", "metadata", "phyla_table")
rm(list = vars.to.remove)
rm("vars.to.remove")

