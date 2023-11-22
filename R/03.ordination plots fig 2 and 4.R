# figures 2 and the bottom half of 4 in the manuscript
# load libraries
library(phyloseq)
library(ggplot2)
library(dplyr)

library(stringr)
library(tidyverse)
library(vegan)
library(plyr)
library(plotly)
library(RColorBrewer)
library(ggpubr)
library(ggsignif)
source("R/02.general_definitions.R")
source("R/02.helper_functions.R")

setwd("~/qiime/fire-blight/pear.flower/20200125/qiime2.aia.filter/R_pears")

## load data
ps <- readRDS("Data/phyloseq_no_buchnera")
ps.rarefied <- readRDS("Data/phyloseq_all_data_rarefied_noB")

ps.flowers.rarefied <- readRDS("Data/phyloseq_costia_flowers_rarefied")
ps.leaves.rarefied <- readRDS("Data/phyloseq_costia_leaves_rarefied")



metadata <- data.frame(sample_data(ps.rarefied), stringsAsFactors = FALSE)
## why is this?
metadata$Date <- as.Date(metadata$Collection.date, format = "%d.%m.%Y")
sample_data(ps.rarefied) <- metadata



# Coscia flowers height and location: figure 2 in the manuscript, and pcoas for supplementary ---------------------------------------------
# this figure contains four plots:
# A. PCoA colored by the flower's height - moved to supp
# B. PcoA colored by the location in the orchard - moved to supp
# C. Alpha diversity comparison of the height.
# D alpha diversity of the location in the orchard. 

ps.coscia.flowers <- readRDS("Data/phyloseq_costia_flowers_rarefied")

meta_cosc <- data.frame(sample_data(ps.coscia.flowers), stringsAsFactors = FALSE)
# change the group and height parameters to factors
meta_cosc$Group <- as.factor(meta_cosc$Group)
meta_cosc$Height <- as.factor(meta_cosc$Height)
sample_data(ps.coscia.flowers) <- meta_cosc
purrr::map(meta_cosc, class)
# rm(meta_cosc)
## plot
# these are the colors of Set1. 
# c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
# plot_ordination_pears is in 02.helper_functions

# plot the ordination of the location on the tree and the height
p.L1 <- plot_ordination_pears(ps.coscia.flowers, OrdMeth = "bray", metadataColumn = "Height")
p.L <- p.L1 +
  scale_color_manual(values = c("#FF7F00", "#4DAF4A"), labels = c("Low", "High")) +
  theme(legend.title = element_blank()) 
print(p.L)
p1 <- plot_ordination_pears(ps.coscia.flowers, OrdMeth = "bray", metadataColumn = "Group")
p <- p1 + scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  theme(legend.title = element_blank())
print(p)
# use ggpubr to combine both plots
Location_pcoa <- ggpubr::ggarrange(p, NULL, p.L,
  labels = c("A", "", "B"),
  ncol = 3,
  widths = c(8, 1, 8)
)
print(Location_pcoa)
png(file = "images/Coscia ordination location and height.png", height = 4500, width = 5000, res = 550)
print(Location_pcoa)
dev.off()

# fig 3A CCA------------------------------------------------------------------


## make a CCA with the location in the orchard, location on the tree, ps and cd
meta_cosc$ps <- meta_cosc$Phenological.stage
meta_cosc$cd <- meta_cosc$Short.collection.date


sample_data(ps.coscia.flowers) <- meta_cosc
color_values_ps_mod <- c(
  ": White buds" = "#E41A1C",
  ": Initial bloom" = "#377EB8",
  ": Full bloom" = "gold",
  ": Initial petal fall" = "#984EA3",
  ": Complete petal fall" = "#FF7F00",
  ": Leaf" = "#4DAF4A"
)
ord_rda <- ordinate(ps.coscia.flowers, method = "CAP", distance = "bray", formula = ~ ps + cd + Height + Group)
rda_plot <- plot_ordination(ps.coscia.flowers, ordination = ord_rda, color = "ps") +
  scale_color_manual(values=color_values_ps) + # the color values are defined in 02.general_definitions.R
  #scale_fill_manual(values=color_values_ps_mod, labels = c("White buds", "Initial bloom", "Full bloom", "Initial petal fall", "Complete petal fall")) +
  theme_minimal() + geom_point(size = 3, alpha = 0.5 )
arrowmat <- vegan::scores(ord_rda, display = "bp")
  arrow_df <- data.frame(labels = rownames(arrowmat), arrowmat)
arrow_map <- aes(xend = CAP1*1.5, 
                 yend = CAP2*1.5, 
                 x = 0, y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels) 
label_map <- aes(x = 2 * CAP1, y = 1.9 * CAP2, 
                 shape = NULL, 
                 color = NULL,
                 label = labels)  

arrowhead = arrow(length = unit(0.02, "npc")) 

rda_plot1 <- rda_plot + geom_segment(mapping = arrow_map, size = .7, 
                                     data = arrow_df, color = "black", 
                                     arrow = arrowhead) + 
 # geom_text(mapping = label_map, size = 4,  fontface = "bold", 
   #         data = arrow_df, 
  #          show.legend = TRUE) + 
  ggtitle("CAP_Plot")  + theme_bw() 
plot(rda_plot1)

png("images/Constrained ordination with arrows high vs low.png", height = 3000, width = 3500, res = 400)
print(rda_plot1)
dev.off()


# make a dataframe with the alpha diversity measures and melt it by the clustering/questioning variable
# and keep just one relevant alpha diversity measurement
variables_for_boxplots <- c("Group", "Height")
colorList <- list(
  "Group" = c("#E41A1C", "#377EB8"),
  "Height" = c("#FF7F00", "#4DAF4A")
)
label_list <- list(
  "Group" = c("C", "D"),
  Height = c("E", "F")
)
# RDA by height and location in the orchard

for (Grouping in variables_for_boxplots) {
  alpha_diversity_melt_coscia_flowers_s <- make_alpha_melted_df(ps.coscia.flowers, Grouping, index = "Shannon")
  alpha_diversity_melt_coscia_flowers_c <- make_alpha_melted_df(ps.coscia.flowers, Grouping, index = "Chao1")
  # combine the two plots. it is done not with facet_wrap so that we cna use two ylabs and remove the title
  # plot_alpha_boxplot is in the helper function files
  p1 <- plot_alpha_boxplot(alpha_diversity_melt_coscia_flowers_s, colorList[[Grouping]]) 
  p2 <- plot_alpha_boxplot(alpha_diversity_melt_coscia_flowers_c, colorList[[Grouping]])
  assign(paste0("p_", Grouping), ggpubr::ggarrange(p2, p1, ncol = 2, labels = label_list[[Grouping]]), envir = .GlobalEnv)
}
print(p_Group)
print(p_Height)

alpha_location_height <- ggpubr::ggarrange(p_Group, NULL, p_Height,
  nrow = 1,
  widths = c(10, 1.5, 10)
)
print(alpha_location_height)
png(file = "images/Location on tree and in orchard.png", height = 3000, width = 3500, res = 350)
ggpubr::ggarrange(Location_pcoa, NULL, alpha_location_height,
  nrow = 3,
  heights = c(10, 1, 10)
)
dev.off()

### test the alpha diversities 
diversity_data <- estimate_richness(ps.coscia.flowers, measures = c("Shannon", "Chao1"))
sample_data(ps.coscia.flowers)$Shannon <- diversity_data$Shannon
sample_data(ps.coscia.flowers)$Chao1 <- diversity_data$Chao1
meta.cosc <- data.frame(sample_data(ps.coscia.flowers))
meta.cosc$cd.ps.g <- paste0(meta.cosc$cd.ps, "_", meta.cosc$Group)
#  a function to perform Wilcoxon test for each combination of phenological stage and date
# Create a function to perform Wilcoxon test for each combination of phenological stage and date
compare_diversity_wilcoxon <- function(data, diversity_measure = "Shannon") {
  unique_combinations <- unique(data$cd.ps)
  results <- list()
  
  for (combination in unique_combinations) {
    subset_data <- data[data$cd.ps == combination,]
    
    group1_data <- subset_data[subset_data$Group == "Group 1", diversity_measure]
    group2_data <- subset_data[subset_data$Group == "Group 2", diversity_measure]
    
    if(length(group1_data) > 0 & length(group2_data) > 0) {
      test <- wilcox.test(group1_data, group2_data)
      results[[combination]] <- test$p.value
    } else {
      results[[combination]] <- NA
    }
  }
  
  return(results)
}

#  a function to perform Wilcoxon test for each combination of phenological stage and date
# Create a function to perform Wilcoxon test for each combination of phenological stage and date
compare_diversity_wilcoxon_UD <- function(data, diversity_measure = "Shannon") {
  unique_combinations <- unique(data$cd.ps)
  results <- list()
  
  for (combination in unique_combinations) {
    subset_data <- data[data$cd.ps == combination,]
    
    group1_data <- subset_data[subset_data$Location.on.tree == "U", diversity_measure]
    group2_data <- subset_data[subset_data$Location.on.tree == "D", diversity_measure]
    
    if(length(group1_data) > 0 & length(group2_data) > 0) {
      test <- wilcox.test(group1_data, group2_data)
      results[[combination]] <- test$p.value
    } else {
      results[[combination]] <- NA
    }
  }
  
  return(results)
}

# Use the function for Shannon diversity
wilcoxon_results_shannon <- compare_diversity_wilcoxon(meta.cosc, "Shannon") %>%
  as.data.frame
print(wilcoxon_results_shannon)

# Use the function for Chao1 diversity
wilcoxon_results_chao1 <- compare_diversity_wilcoxon(meta.cosc, "Chao1")
print(wilcoxon_results_chao1)

wilcoxon_results_shannon_UD <- compare_diversity_wilcoxon_UD(meta.cosc, "Shannon")
print(wilcoxon_results_shannon_UD)
# one pair came out significant (p=0.017) so we will do FDR
wilcoxon_results_shannon_UD <- wilcoxon_results_shannon_UD %>%
  as.data.frame()
wilcoxon_results_shannon_UD <- data.frame(t(wilcoxon_results_shannon_UD))
wilcoxon_results_shannon_UD$padj <- p.adjust(wilcoxon_results_shannon_UD$V1, method = "fdr")
## and now everything is insignificant. 
# Use the function for Chao1 diversity
wilcoxon_results_chao1_UD <- compare_diversity_wilcoxon_UD(meta.cosc, "Chao1")
print(wilcoxon_results_chao1_UD)
# again, here 13.04 F2-G is significant (p = 0.04)
# # let's adjust using fdr

wilcoxon_results_chao1_UD <- data.frame(t(wilcoxon_results_chao1_UD))
wilcoxon_results_chao1_UD$padj <- p.adjust(wilcoxon_results_chao1_UD$t.wilcoxon_results_chao1_UD., method = "fdr")


# arrange the four pcoas for coscia flowers (location, height, ps, cd) --------


A <- ggpubr::ggarrange(p.ps, p.cd,
  widths = c(1.1, 1), nrow = 1,
  labels = c("A", "B")
)
B <- ggpubr::ggarrange(p.L, NULL, p, NULL,
  widths = c(0.95, 0.05, 0.80, 0.1), nrow = 1,
  labels = c("C", " ", "D", " ")
)
png(file = "images/Coscia ordination four coloring.png", height = 4500, width = 5000, res = 550)
ggpubr::ggarrange(A, B, nrow = 2)
dev.off()



# Figure 4  bottom half Coscia PCoA ----------------------------------------------------

# this is in 04.erwiniaceae abundances heatmap with the rest of the figure


# supplementary figure - five pcoas, one for each collection date  --------

  ps.rarefied <- readRDS("Data/phyloseq_coscia_rarefied")
metadata.rarefied <- data.frame(sample_data(ps.rarefied))
Collection.dates <- metadata.rarefied$Collection.date %>%
  unique()
CDPlots <- list()
for(i in 1:length(Collection.dates)){
  CD <- Collection.dates[i]
  print(i)
  ps.1 <- subset_samples(ps.rarefied, Collection.date == CD)
  print(ps.1)
  CDPlot <- plot_ordination_pears(ps.1, OrdMeth = "bray", metadataColumn = "Phenological.stage", "all data") +
    stat_ellipse() + ggtitle(CD)
  print(CDPlot)
  CDPlots[[i]] <- CDPlot
}
png(
  file = "images/PCoAs by collection dates.png",
  height = 2500, width = 3500, res = 350
)
gridExtra::grid.arrange(CDPlots[[1]], CDPlots[[2]],CDPlots[[4]], CDPlots[[5]], CDPlots[[3]])
dev.off()

# now do the same for phenological stages
PS <- metadata.rarefied$Phenological.stage %>%
  unique()
CDPlots <- list()
for(i in 1:length(PS)){
  PSCurr <- PS[i]
  print(i)
  ps.1 <- subset_samples(ps.rarefied, Phenological.stage == PSCurr)
  print(ps.1)
  CDPlot <- plot_ordination_pears(ps.1, OrdMeth = "bray", metadataColumn = "Collection.date", "all data") +
    stat_ellipse() + ggtitle(PSCurr)
  print(CDPlot)
  CDPlots[[i]] <- CDPlot
}
png(
  file = "images/PCoAs by Phenological stage.png",
  height = 2500, width = 3500, res = 350
)
gridExtra::grid.arrange(CDPlots[[1]], CDPlots[[2]], CDPlots[[3]],CDPlots[[6]] ,CDPlots[[5]], CDPlots[[4]])
dev.off()


# all data with things ----------------------------------------------------
## color by collection date
## add different shapes to the phenological stages
## add circles by the dates alone

plot_ordination_pears(ps.rarefied, OrdMeth = "bray", metadataColumn = "CPS", "all data")
p1 <- plot_ordination_pears(ps.rarefied, OrdMeth = "bray", metadataColumn = "Collection.date", "all data collection date")
plot_ordination_pears(ps.flowers.rarefied, OrdMeth = "bray", metadataColumn = "Collection.date", "flowers data collection date")
plot_ordination_pears(ps.leaves.rarefied, OrdMeth = "bray", metadataColumn = "Collection.date", "leaves data collection date")


# combine the phyla bar plot from 03. file with the pcoa and save --------------
png(
  file = "images/PCoA and phyla bar plot.png",
  height = 2500, width = 5500, res = 350
)
ggpubr::ggarrange(p1, p_phyla_bar_plot)


dev.off()

