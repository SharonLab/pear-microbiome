# 
## this contains the scripts for the combined bar plot and box plot figures
## these are:
# 1. figure 3 comparison of the early stages of spadona and coscia. it also has a pcoa
# 2. figure S2: the most abundant genera with their alpha diversities by time and phenological stages
# figure S3: relative abundances of phyla
# tables S4 and S5

setwd("~/qiime/fire-blight/pear.flower/20200125/qiime2.aia.filter/R_pears")
library(phyloseq)
library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(reshape2)
# library(speedyseq)
source("R/02.helper_functions.R")
source("R/02.general_definitions.R")
ps <- readRDS("Data/phyloseq_no_buchnera")
ps
metadata <- data.frame(sample_data(ps))



# agglomerate the phyloseq object by taxonomic ranks ---------------------------------

ps.genera <- tax_glom(ps, taxrank = "Genus")
ps.phyla <- tax_glom(ps, taxrank = "Phylum")
ps.family <- tax_glom(ps, taxrank = "Family")
ps.order <- tax_glom(ps, taxrank = "Order")


# Figure S2 bar plot ----------------------------------------------------

ps.current <- ps.genera # for Fig S2
ps.current <- ps.phyla  # for figure S3
metadata.current <- sample_data(ps.current)
metadata.current$cd.ps <- paste0(metadata.current$Collection.date, " ", metadata.current$Phenological.stage)
sample_data(ps.current) <- metadata.current

db.current1 <- make_data_bars1(ps.current, 0.005, "cd.ps") # for genus figure S2
# for phyla use the following:
db.current1 <- make_data_bars1(ps.current, 0.00005, "cd.ps") # for phyla figure S3
# split the sample into phenological stage and collection date
db.current <- separate(db.current1, Sample, into = c("Collection.date", "Phenological.stage"), sep = "^\\S*\\K\\s+")


order.ps <- c("White buds", "Initial bloom", "Full bloom", "Initial petal fall", "Complete petal fall", "Leaf")
positions <- c("24.03.2019", "27.03.2019", "31.03.2019", "11.04.2019", "13.04.2019", "15.04.2019", "18.04.2019", "21.04.2019")

tbpf <-
  taxa_bar_plot(db.current, "Collection.date")
tbpf
# tbpf + geom_point(data = data_bars.ra1.e, aes(x = CD, y = total_ra))
p1 <-
  tbpf + # theme(legend.position = "None") +
  scale_x_discrete(limits = positions) +
  facet_wrap(~ factor(Phenological.stage, levels = order.ps), nrow = 6, ncol = 1) +
  labs(x = "\nCollection date", title = "Taxonomy bar plots") +
  geom_vline(xintercept = 3.5, color = "mediumorchid", linewidth = 1.5)

png(filename = "images/alpha and barplot phyla S3.png", height = 4500, width = 2500, res = 400)
print(p1)
dev.off()



# Figure S2 alpha diversity  ---------------------------------------------------------

ps.rarefied <- readRDS("Data/phyloseq_all_data_rarefied_noB")
metadata.rarefied <- data.frame(sample_data(ps.rarefied), stringsAsFactors = FALSE)
metadata.rarefied$Shannon <- estimate_richness(ps.rarefied, measures = "Shannon")$Shannon
sample_data(ps.rarefied) <- metadata.rarefied


rects <- data.frame(xstart = c("24.03.2019", "11.04.2019"), xend = c("24.03.2019", "11.04.2019"), colss = c("Spadona", "Coscia"))
p2 <- ggplot(metadata.rarefied, aes(x = Collection.date, y = Shannon)) +
  geom_boxplot() +
  facet_wrap(~ factor(Phenological.stage, levels = order.ps), nrow = 6, ncol = 1) +
  # scale_color_gradient(palette = "Spectral") +
  scale_x_discrete(limits = positions) +
  geom_vline(xintercept = 3.5, color = "mediumorchid", linewidth = 1.5) +
  theme_minimal() +
  labs(y = "Shannon Diversity Index\n", x = "\nCollection date", title = "Shannon diversity index") +
  theme(legend.position = "None", axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = -0.1))

p2
pdf(file = "images/bar_plots_with_alpha.pdf", height = 10, width = 8)
## or
png(filename = "images/alpha and barplot family.png", height = 4500, width = 4000, res = 500)
ggarrange(p1, NULL, p2, NULL,
  labels = c("A", " ", "B", " "),
  ncol = 4, nrow = 1,
  widths = c(10, 1, 5, 1)
) + ggtitle("0.01 per sample")


dev.off()




# alpha stats -------------------------------------------------------------
## kruskal walis and post hoc using the agricolae package on shannon diversity

ps.rarefied <- readRDS("Data/phyloseq_coscia_rarefied")
metadata.rarefied <- data.frame(sample_data(ps.rarefied), stringsAsFactors = FALSE)
metadata.rarefied$Shannon <- estimate_richness(ps.rarefied, measures = "Shannon")$Shannon

# arrange the dataframe for the analysis

alpha_for_k <- cbind(metadata.rarefied$cd.ps, as.numeric(metadata.rarefied$Shannon)) %>%
  data.frame()
alpha_for_k$X1 <- gsub(" ", "_", x = alpha_for_k$X1)
alpha_for_k <- alpha_for_k[order(alpha_for_k$X1), ]
levels(alpha_for_k$X1)

kw <- kruskal.test(X2 ~ X1, data = alpha_for_k)
kw
# Kruskal-Wallis rank sum test
#
# data:  X2 by X1
# Kruskal-Wallis chi-squared = 135.49, df = 25, p-value < 2.2e-16
# These vaues are for Coscia only!!

# manipulate the data frame so that we can test each phenological stage and collection date separately.
#
alpha_for_k1 <- alpha_for_k %>%
  separate(X1, c("CD", "PS"), "_")
dates <- as.vector(unique(alpha_for_k1$CD))
Phenological.stages <- as.vector(unique(alpha_for_k1$PS))

library(FSA)
date_kw <- data.frame()

for (i in dates) {
  print(i)
  temp_df <- alpha_for_k1[alpha_for_k1$CD == i, ]
  kw <- kruskal.test(X2 ~ PS, data = temp_df)
  line1 <- cbind(i, kw$statistic[1], round(kw$p.value[1], digits = 4))
  date_kw <- rbind(date_kw, line1)

}
date_kw$padj <- p.adjust(date_kw$V3, method = "fdr") %>%
  round(digits = 3)
write.csv(date_kw, file = "Results/kruskal walis test results collection dates Table S4.csv")
ps_kw <- data.frame()
for (i in Phenological.stages) {
  print(i)
  temp_df <- alpha_for_k1[alpha_for_k1$PS == i, ]
  kw <- kruskal.test(X2 ~ CD, data = temp_df)
  line1 <- cbind(i, round(kw$p.value[1], digits = 4))
  ps_kw <- rbind(ps_kw, line1)
}

ps_kw$padj <- p.adjust(ps_kw$V2, method = "fdr") %>%
  round(digits = 3)
write.csv(ps_kw, file = "Results/kruskal walis test results Phenological stages table S5.csv")

rm(line1, temp_df)



# Figure 3 Spadona vs Coscia white bud and initial bloom --------------------------
# compare the same phenological stages in both cultivars
# -

## We are comparing tthe white buds and Initial bloom stages from coscia and spadona
## In Coscia we analyze only the first collection dates before Erwinia Amylovora took over the population

# load the rarefied initial stages phyloseq object
ps.initial <- readRDS("Data/phyloseq_initial_stages_rarefied")
ps.initial
sample_sums(ps.initial)

p.init.pca <- plot_ordination_pears(ps.initial, metadataColumn = "Short.collection.date") 
p.init.pca  
p.init.pca.pe <- plot_ordination_pears(ps.initial, metadataColumn = "Phenological.stage") 

ps.initial.g <- tax_glom(ps.initial, taxrank = "Genus")
data_bars <- make_data_bars1(ps.initial.g,  grouping_var = "CPS")
data_bars<- data_bars %>%
  mutate(Genus = ifelse(Genus == "uncultured", paste("Uncultured", Family), Genus))
data_bars <- add_threshold_to_databars(data_bars, 0.01)

positions <- c("Spadona White buds", "Spadona Initial bloom", "", "Coscia White buds", "Coscia Initial bloom")
p <- PlotBarsGenus(data_bars)
print(p)
p1 <- p + scale_x_discrete(
  limits = positions,
  labels = c("White buds", "Initial bloom", " ", "White buds", "Initial bloom")) +
  labs(y = "Relative Abundance\n", x = "\n") +
  theme(legend.key.size = unit(0.37, "cm"))
print(p1)

positions <- c("White buds", "Initial bloom", "Full bloom", "Initial petal fall", "Complete petal fall", "Leaf")
positions <- c("24.03.2019", "27.03.2019", "31.03.2019", "11.04.2019", "13.04.2019")
p2.0 <- plot_ordination_pears(ps.initial, OrdMeth = "bray", metadataColumn = "CPS", titleAddition = "cultivar and phenological stage")
# the manual colours are the first four elements of "Set1"
colorc("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
p2 <- p2.0 +
  scale_colour_manual(values = c("#E41A1C", "#FF7F00", "#377EB8", "#4DAF4A"), name = NULL) +
  xlab(paste0("12%\n\n\n\n"))
print(p2)


alpha_diversity_melt_s <- make_alpha_melted_df(ps.initial, "CPS", index = "Shannon")
alpha_diversity_melt_c <- make_alpha_melted_df(ps.initial, "CPS", index = "Chao1")

## consider refactorin:
# 1. use the helper funtion plot_alpha_boxplot
# 2. use the color schemes from the general definitions file
positions <- c("Spadona White buds", "Spadona Initial bloom", " ", "Coscia White buds", "Coscia Initial bloom")
positions_for_legend <- c("Spadona White buds", "Spadona Initial bloom", "Coscia White buds", "Coscia Initial bloom")
p3 <-
  ggplot(alpha_diversity_melt_c, aes(x = Group, y = value, fill = Group)) +
  geom_boxplot() +
  labs(y = "Chao1 index", x = "") +
  theme_minimal() +
  scale_x_discrete(limits = positions, labels = c("White buds", "Initial bloom", " ", "White buds", "Initial bloom")) +
  scale_fill_manual(values = color_values_cps_init) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), legend.position = "off") +
  ggtitle(" ")
print(p3)

p4 <-
  ggplot(alpha_diversity_melt_s, aes(x = factor(Group), y = value, fill = Group)) +
  geom_boxplot() +
  labs(y = "Shannon", x = "") +
  theme_minimal() +
  scale_x_discrete(limits = positions, labels = c("White buds", "Initial bloom", " ", "White buds", "Initial bloom")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), legend.title = element_blank()) +
  scale_fill_manual(
    values = color_values_cps_init,
    name = NULL, limits = positions_for_legend
  ) +
  guides(color = guide_legend(reverse = TRUE)) +
  ggtitle(" ")
print(p4)

# Test differences of Shannon 
kruskal.test(value ~ Group, data = alpha_diversity_melt_s)
# Kruskal-Wallis rank sum test
# 
# data:  value by Group
#Kruskal-Wallis chi-squared = 4.2444, df = 3, p-value = 0.2363
## 

#Test Chao1 
kruskal.test(value ~ Group, data = alpha_diversity_melt_c)
# data:  value by Group
# Kruskal-Wallis chi-squared = 14.3, df = 3, p-value = 0.002524

shannon_wilcoxon <- pairwise.wilcox.test(alpha_diversity_melt_s$value, alpha_diversity_melt_s$Group, p.adjust.method = "fdr")
# Pairwise comparisons using Wilcoxon rank sum exact test 
# 
# data:  alpha_diversity_melt_s$value and alpha_diversity_melt_s$Group 
# 
#                 Coscia Initial bloom Coscia White buds Spadona Initial bloom
# Coscia White buds     0.35                 -                 -                    
# Spadona Initial bloom 0.68                 0.35              -                    
# Spadona White buds    0.35                 0.90              0.42                 

# P value adjustment method: fdr 

Chao_wilcoxon <- pairwise.wilcox.test(alpha_diversity_melt_c$value, alpha_diversity_melt_s$Group, p.adjust.method = "none")
Chao_wilcoxon
# Pairwise comparisons using Wilcoxon rank sum exact test 
# 
# data:  alpha_diversity_melt_c$value and alpha_diversity_melt_s$Group 
# 
#               Coscia Initial bloom Coscia White buds Spadona Initial bloom
# Coscia White buds     0.00119              -                 -                    
# Spadona Initial bloom 0.77869              0.26018           -                    
# Spadona White buds    0.00027              0.01801           0.04476              

# P value adjustment method: none 

top_row <- ggpubr::ggarrange(p1, NULL,
  labels = c("D", ""),
  nrow = 1, ncol = 2,
  widths = c(20, 6)
)
bottom_row <- ggpubr::ggarrange(p2, NULL, p3, NULL, p4,
  labels = c("A", " ", "B", " ", "C"),
  nrow = 1, ncol = 5,
  widths = c(20,0.5, 6, 0.5, 12)
)
p5 <- ggpubr::ggarrange( bottom_row, top_row,
  nrow = 2
)


p5
png(file = "images/Fig 3 Cultivar comparison.png", width = 5000, height = 4000, res = 500)
print(p5)
dev.off() 

# alpha plot coscia no leaves ---------------------------------------------
# a comparison of alpha diversities of Coscia flowers 
# by early and late phenological stages.
# this was not included in the manuscript. 

order.ps1 <- c("Early", "Late", "Leaf")
metadata.rarefied.c <- metadata.rarefied[metadata.rarefied$Cultivar == "Coscia", ]
metadata.rarefied.c.f <- metadata.rarefied.c[metadata.rarefied.c$Flower.or.leaf == "Flower", ]

p3 <- ggplot(metadata.rarefied.c.f, aes(x = Collection.date, y = Shannon, color = Phenological.stage)) +
  geom_boxplot() +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(y = "Shannon Diversity Index\n", x = "\nCollection date", title = "Shannon diversity index") +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = -0.1))
p3

