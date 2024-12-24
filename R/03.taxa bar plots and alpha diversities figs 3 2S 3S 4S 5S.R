# 
## this contains the scripts for the combined bar plot and box plot figures
## these are:
# 1. figure 3 comparison of the early stages of spadona and coscia. it also has a pcoa
# 2. figure S2: the most abundant genera with their alpha diversities by time and phenological stages
# figure S3: relative abundances of phyla
# tables S4 and S5

library(phyloseq)
library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(patchwork)
source("R/02.helper_functions.R")
source("R/02.general_definitions.R")
source("R/BarPlots.R")
ps <- readRDS("Data/phyloseq_no_buchnera")
ps
metadata <- data.frame(sample_data(ps))



# agglomerate the phyloseq object by taxonomic ranks ---------------------------------

ps.genera <- tax_glom(ps, taxrank = "Genus")
ps.phyla <- tax_glom(ps, taxrank = "Phylum")
ps.family <- tax_glom(ps, taxrank = "Family")
ps.order <- tax_glom(ps, taxrank = "Order")


# Figure S2 ----------------------------------------------------


ps.current <- ps.phyla  # for figure S2
metadata.current <- sample_data(ps.current)
metadata.current$cd.ps <- paste0(metadata.current$Collection.date, " ", metadata.current$Phenological.stage)
sample_data(ps.current) <- metadata.current

# for phyla use the following:
db.current1 <- make_data_bars1(ps.current, 0.00005, "cd.ps") # for phyla figure S2
# split the sample into phenological stage and collection date
db.current <- separate(db.current1, Sample, into = c("Collection.date", "Phenological.stage"), sep = "^\\S*\\K\\s+", remove = FALSE)
#db.current <- db.current1
db.current <- db.current %>% filter(!is.na(Abundance))
order.ps <- c("White buds", "Initial bloom", "Full bloom", "Initial petal fall", "Complete petal fall", "Leaf")
positions <- c("24.03", "27.03", "31.03", "11.04", "13.04", "15.04", "18.04", "21.04")

# Adjusted plot
number_of_phyla <- length(unique(db.current$Phylum))

tbpf <- ggplot(db.current, aes(x = Short.collection.date, y = Abundance, fill = forcats::fct_reorder(Phylum, Abundance))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = taxrank, values = rev(getPalette(number_of_phyla))) +
  theme_minimal() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) + 
  facet_wrap(~ factor(Phenological.stage, levels = order.ps), nrow = 6, ncol = 1) 
tbpf
# Additional adjustments
p1 <- tbpf +
  scale_x_discrete(limits = positions) +
  labs(x = "\n\nCollection date", title = "Phylum bar plots", y = "Relative abundance\n") +
  geom_vline(xintercept = 3.5, color = "mediumorchid", linewidth = 1.5)

p1

png(filename = "images/Fig S2 barplot phyla.png", height = 4000, width = 3000, res = 450)
print(p1)
dev.off()


# figure S4 bar plot ------------------------------------------------------

ps.current <- ps.genera # for Fig S4
metadata.current <- sample_data(ps.current)
metadata.current$cd.ps <- paste0(metadata.current$Collection.date, " ", metadata.current$Phenological.stage)
sample_data(ps.current) <- metadata.current

db.current1 <- make_data_bars1(ps.current, 0.005, "cd.ps") # for genus figure S4
db.current <- separate(db.current1, Sample, into = c("Collection.date", "Phenological.stage"), sep = "^\\S*\\K\\s+", remove = FALSE)

db.current <- db.current %>% filter(!is.na(Abundance))
order.ps <- c("White buds", "Initial bloom", "Full bloom", "Initial petal fall", "Complete petal fall", "Leaf")
positions <- c("24.03", "27.03", "31.03", "11.04", "13.04", "15.04", "18.04", "21.04")

# Adjusted plot
number_of_phyla <- length(unique(db.current$Genus))
missing_name <- "Genus < 0.005"
all_taxa <- unique(db.current$Genus)
all_taxa_except_Erwinia <- all_taxa[all_taxa != "Erwinia" & all_taxa != missing_name]

getPalette <- colorRampPalette(brewer.pal(9, "Set1"))


my_palette <- setNames(getPalette(length(all_taxa_except_Erwinia)), all_taxa_except_Erwinia)
my_palette["Pseudomonas"] <- "cornflowerblue"
my_palette["Erwinia"] <- "red"
my_palette[missing_name] <- "grey70"

db.current$Genus <- factor(db.current$Genus, levels = c(missing_name, sort(all_taxa_except_Erwinia), "Erwinia"))

tbpf <- ggplot(db.current, aes(x = Short.collection.date, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Genus", values = rev(getPalette(number_of_phyla))) +
  theme_minimal() +
  theme(legend.position = "right", axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = -0.1))+ 
  facet_wrap(~ factor(Phenological.stage, levels = order.ps), nrow = 6, ncol = 1) 
tbpf
# Additional adjustments
p1 <- tbpf +
  scale_x_discrete(limits = positions) +
  labs(x = "\n\nCollection date", title = "Genera bar plots", y = "Relative abundance\n", fill = "Genus") +
  geom_vline(xintercept = 3.5, color = "mediumorchid", linewidth = 1.5) +
  guides(fill = guide_legend(ncol = 1))

p1


# Figure S4 alpha diversity  ---------------------------------------------------------

ps.rarefied <- readRDS("Data/phyloseq_all_data_rarefied_noB")
metadata.rarefied <- data.frame(sample_data(ps.rarefied), stringsAsFactors = FALSE)
metadata.rarefied$Shannon <- estimate_richness(ps.rarefied, measures = "Shannon")$Shannon
sample_data(ps.rarefied) <- metadata.rarefied


rects <- data.frame(xstart = c("24.03", "11.04"), xend = c("24.03", "11.04"), colss = c("Spadona", "Coscia"))
p2 <- ggplot(metadata.rarefied, aes(x = Short.collection.date, y = Shannon)) +
  geom_boxplot() +
  facet_wrap(~ factor(Phenological.stage, levels = order.ps), nrow = 6, ncol = 1) +
  # scale_color_gradient(palette = "Spectral") +
  scale_x_discrete(limits = positions) +
  geom_vline(xintercept = 3.5, color = "mediumorchid", linewidth = 1.5) +
  theme_minimal() +
  labs(y = "Shannon Diversity Index\n", x = "\nCollection date", title = "Shannon diversity index") +
  theme(legend.position = "None", axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = -0.1))

p2


combined_plot <- wrap_plots(
  p1 + labs(tag = "A"),
  p2 + labs(tag = "B"),
  ncol = 2
) + 
  plot_annotation(
    title = " "
  ) & 
  theme(plot.margin = margin(25, 25, 5, 5))  # Left padding only for all plots

print(combined_plot)
# pdf(file = "images/bar_plots_with_alpha.pdf", height = 10, width = 8)
## or
png(filename = "images/Fig S4 alpha and barplot genera alphabetically.png", height = 4500, width = 4500, res = 500)
print(combined_plot)
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



# Figure 2E Spadona vs Coscia white bud and initial bloom, height, location in the orchard --------------------------

# A. hight alpha diversities
# B. location in the orchard alpha diversities
# C. Initial stages PCoA
# D Initial stages alpha diversities
# E. Initial stages bar plot


# compare the same phenological stages in both cultivars
# -

## We are comparing the white buds and Initial bloom stages from coscia and spadona
## In Coscia we analyze only the first collection dates before Erwinia Amylovora took over the population

# load the rarefied initial stages phyloseq object
ps.initial <- readRDS("Data/phyloseq_initial_stages_rarefied")
ps.initial
sample_sums(ps.initial)
ps.initial.g <- tax_glom(ps.initial, taxrank = "Genus")

# make initial ordination plots
p.init.pca <- plot_ordination_pears(ps.initial, metadataColumn = "Short.collection.date") 
p.init.pca  
p.init.pca.pe <- plot_ordination_pears(ps.initial, metadataColumn = "Phenological.stage") 

## make the bar plot for the initial stages
data_bars <- make_data_bars1(ps.initial.g,  grouping_var = "CPS")
data_bars<- data_bars %>%
  mutate(Genus = ifelse(Genus == "uncultured", paste("Uncultured", Family), Genus))
data_bars <- add_threshold_to_databars(data_bars, 0.01)

positions <- c("Spadona White buds", "Spadona Initial bloom", "", "Coscia White buds", "Coscia Initial bloom")
p_bar <- PlotBarsGenus(data_bars)
print(p_bar)


p_bar1 <- p_bar + scale_x_discrete(
  limits = positions,
  labels = c("WB\n\nSpadona", "IB\n", " ", "WB\n\nCoscia", "IB\n")
) +
  labs(y = "Relative Abundance\n", x = NULL) +
  theme(
    legend.key.size = unit(0.37, "cm"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, lineheight = 0.8, size = 10)
  )

# Add extra margin at the bottom if needed
# p_bar1 <- p_bar1 + theme(plot.margin = unit(c(5.5, 5.5, 10, 5.5), "points"))

print(p_bar1)

png("images/figure_2E_initial_stages_barplot_with alphabetised_order.png", height = 2500, width = 3000, res = 400)
print(p_bar1)
dev.off()

## make the final ordination plot for the initial stages
positions <- c("White buds", "Initial bloom", "Full bloom", "Initial petal fall", "Complete petal fall", "Leaf")
positions <- c("24.03.2019", "27.03.2019", "31.03.2019", "11.04.2019", "13.04.2019")
p2.0 <- plot_ordination_pears(ps.initial, OrdMeth = "bray", metadataColumn = "CPS", titleAddition = "cultivar and phenological stage")
# the manual colours are the first four elements of "Set1"
colorc("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
p_ordination_init <- p2.0 +
  scale_colour_manual(values = c("#E41A1C", "#FF7F00", "#377EB8", "#4DAF4A"), name = NULL) 
print(p_ordination_init)

# make the alpha diversity plots
alpha_diversity_melt_s <- make_alpha_melted_df(ps.initial, "CPS", index = "Shannon")
alpha_diversity_melt_c <- make_alpha_melted_df(ps.initial, "CPS", index = "Chao1")
alpha_diversity_init_df <- rbind(alpha_diversity_melt_c, alpha_diversity_melt_s)

p_alpha_init <- ggplot(alpha_diversity_init_df, aes(x = Group, y = value, fill = Group)) +
  geom_boxplot() +
  theme_bw() +
  geom_jitter(alpha = 0.5, size = 1) +
  facet_wrap(~measure, scales = "free_y") +
  scale_x_discrete(limits = c("Spadona White buds", "Spadona Initial bloom", "Coscia White buds", "Coscia Initial bloom"), labels = c("WB\n\nSpadona", "IB", " ", "WB\n\nCoscia", "IB\n")) +
  scale_fill_manual(values = color_values_cps_init) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10), legend.position = "off", 
        axis.title.x = element_blank(),
        plot.title = element_text(size = 10)) +
  stat_compare_means(comparisons = list(
    c("Spadona White buds", "Spadona Initial bloom"),
    c("Coscia White buds", "Coscia Initial bloom")
  ),
  method = "wilcox.test", label = "p.signif", size = 4) +  # Adjust label.y to position the labels
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))



## we now add the following (created in 03.ordination_plots):
# Chao_Group
# Shannon_Group
# Chao_Height
# Shannon_Height

## arrange with patchworks package. 
## 
#Extract legend from p_bar1
p_bar1_legend <- get_legend(p_bar1 + theme(legend.position = "right"))

# Remove legend from p_bar1
p_bar1_plot <- p_bar1 + theme(legend.position = "none")
top_row <- p_alpha_location | p_alpha_height | p_ordination_init
bottom_row <- p_alpha_init |  p_bar1 | plot_spacer() 
# Create the layout
figure_2 <- 
  top_row / bottom_row +
    plot_annotation(tag_levels = "A")

# Print the figure
print(figure_2)

png("images/figure_2.png", height = 3000, width = 5500, res = 400)
print(figure_2)
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


# bar plot of the last 3 pjhenological stages -----------------------------

ps.end <- ps.genera %>%
  subset_samples(CPS == "Coscia Full bloom" | CPS == "Coscia Initial petal fall" | CPS == "Coscia Complete petal fall")
meta.end <- as.data.frame(sample_data(ps.end)) 
unique(meta.end$CPS)
 
data_bars <- psmelt(ps.end)
taxa_bar_plot(ps.end, "sample_Sample")
  

number_of_phyla <- length(unique(data_bars$Genus))
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))(number_of_phyla)
ggplot(data_bars, aes(x = sample_Sample, y = Abundance, fill = forcats::fct_reorder(Genus, Abundance)) )+
  # geom_bar(position = "stack") +
  geom_col() +
  scale_fill_manual(values = rev(getPalette)) +
  theme_minimal() +
  theme(legend.position = "right", axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = -0.1)) +
  labs(y = "Relative Abundance\n", x = "\n", fill = taxrank) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(1, "lines")
  )


# bar plots of most abundant genera in leaves and flowers -----------------

ps.genera <- ps %>%
  tax_glom(taxrank = "Genus")
ps.genera_flowers <- ps.genera %>%
  subset_samples(Flower.or.leaf == "Flower")
ps.genera_leaves <- ps.genera %>%
  subset_samples(Flower.or.leaf == "Leaf")



db.all <- rbind(db.leaves, db.flowers)

db.all$Sample <- factor(db.all$Sample, levels = positions_cd_s)


tbpf <- basic_bar_plot(db.all) +
  facet_wrap(~ factor(Flower.or.leaf)) +
  scale_x_discrete(limits = positions_cd_s) +
  ggtitle("The most abundant Genera in pear leaves and flowers\n") +
  labs(x = "\n", y = "Relative Abundance\n") +
  theme(
    axis.text = element_text(size = 16),          # Axis tick labels
    axis.title = element_text(size = 18),         # Axis titles
    strip.text = element_text(size = 18),         # Facet titles
    legend.text = element_text(size = 16),        # Legend text
    legend.title = element_text(size = 16),       # Legend title
    plot.title = element_text(size = 20, hjust = 0.5) # Plot title, centered
  )

print(tbpf)

png("images/Genera in flowers and leaves.png", height = 4500, width = 10000, res = 750)
print(tbpf)
dev.off()


db.leatbpfdb.leaves <- make_data_bars1(ps.genera_leaves, 0.005, "Short.collection.date", RA = TRUE) 

tbpf <- basic_bar_plot(db.leaves) +
  scale_x_discrete(limits = positions_cd_s) +
  ggtitle("The most abundant Genera in leaves") +
  labs(x = "\n", y = "Relative Abundance\n")  
 
  # guides(fill = guide_legend(ncol = 1))
print(tbpf)

db.flowers <- make_data_bars1(ps.genera_flowers, 0.005, "Short.collection.date", RA = TRUE) 

tbpf2 <- basic_bar_plot(db.flowers) +
  scale_x_discrete(limits = positions_cd_s) +
  ggtitle("The most abundant Genera in flowers") +
  labs(x = "\n", y = "Relative Abundance\n")
print(tbpf2)

combined_plot <-  wrap_plots(tbpf2, tbpf, ncol = 2)
print(combined_plot)
