## figure 4 and supplementary figures of specific taxa.
## NOTE that the corncob isn't well implemented and was used for exploratory data analysis only.
## It was not included in the final manuscript.
## It should be changed to differential abundance rather than what it is now (repeat bbdml without FDR)

# load required functions
source("R/02.helper_functions.R")
source("R/02.general_definitions.R")
source("R/BarPlots.R")
source("R/DataBars.R")

# load packages
library(phyloseq)
library(dplyr)
library(corncob)
library(forcats)
library(RColorBrewer)
library(ggplot2)

# load the phyloseq objects that were created in 01. load_data_to_ps_object_rarefy.R
ps.rarefied <- readRDS("Data/phyloseq_all_data_rarefied_noB") %>%
  subset_samples(Flower.or.leaf == "Flower")
ps.rarefied
taxa_names(ps.rarefied) <- paste0("ASV.", 1:ntaxa(ps.rarefied))
## positions specifies the order of the data
positions1 <- c(
  "Spadona White buds", "Spadona Initial bloom", "Coscia White buds",
  "Coscia Initial bloom", "Coscia Full bloom", "Coscia Initial petal fall",
  "Coscia Complete petal fall"
)

## this is a spaghetti code: at the moment this cannot be inside a function unless we subset the
# phyloseq outside the function because phyloseq::subset_taxa has some issues
# with the global environment vs. the local one and cannot be used
# inside functions.

# More functions used in this script
bar_plot_of_taxa <- function(data_bars, positions = NULL, colorPalette = NULL, abundance_type = "RA") {
  taxrank <- "OTU"
  p <- PlotBarDefault(data_bars = data_bars, x_group = "Sample", taxrank = taxrank, colorPalette = colorPalette, abundance_type = abundance_type) +
    theme(legend.position = "none")
  if (!is.null(positions)) {
    p <- p + scale_x_discrete(limits = positions)
  }
  return(p)
}


plot_alpha_diversity <- function(ps) {
  positions1 <- c("Spadona White buds", "Spadona Initial bloom", "Coscia White buds", "Coscia Initial bloom", "Coscia Full bloom", "Coscia Initial petal fall", "Coscia Complete petal fall")

  alpha_diversity_melt_ps <- make_alpha_melted_df(ps, "CPS", index = "Shannon")
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

  make_plot <- function(df, ylabel, legend_position = "off") {
    kruskal_result <- try(kruskal.test(value ~ Group, df), silent = TRUE)
    if (inherits(kruskal_result, "try-error")) {
      annotation_text <- "Kruskal-Wallis test failed"
    } else if (kruskal_result$p.value >= 0.05) {
      annotation_text <- paste("p-value =", format(kruskal_result$p.value, digits = 3), "\nNon significant")
    } else {
      annotation_text <- paste("p-value =", format(kruskal_result$p.value, digits = 3), "\nSignificant")
    }

    p <- ggplot(df, aes(x = Group, y = value, fill = Group)) +
      geom_boxplot() +
      labs(y = ylabel, x = "") +
      theme_minimal() +
      scale_x_discrete(limits = positions1) +
      scale_fill_manual(values = getPalette(length(unique(df$Group))), breaks = positions1) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), legend.position = "none") +
      annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, label = annotation_text, size = 3)

    return(p)
  }
  p4 <- make_plot(alpha_diversity_melt_ps, ylabel = "Shannon index", "right")


  return(p4)
}


change_metadata <- function(ps) {
  # this function addresses issues created by one of the phyloseq functions
  ps.pseudo.meta <- data.frame(sample_data(ps)) %>%
    mutate(corncob_groups = case_when(
      CPS == "Spadona White buds" | CPS == "Spadona Initial bloom" ~ "Spadona_early",
      CPS == "Coscia White buds" | CPS == "Coscia Initial bloom" ~ "Coscia_early",
      CPS == "Coscia Full bloom" | CPS == "Coscia Initial petal fall" | CPS == "Coscia Complete petal fall" ~ "Coscia_late"
    ))
  # ps.pseudo.meta <- ps.pseudo.meta[ps.pseudo.meta$Flower.or.leaf == "Flower", ]
  ps.pseudo.meta$corncob_groups <- as.character(ps.pseudo.meta$corncob_groups)
  ps.pseudo.meta$CPS <- factor(ps.pseudo.meta$CPS, levels = positions1)
  sample_data(ps) <- ps.pseudo.meta
  purrr::map(data.frame(sample_data(ps)), class)
  return(ps)
}

# perform_bbdml(taxa_list, ps.pseudo.reduced, "Pseudomonas")
# function to analyze each significant taxon separately
perform_bbdml <- function(df, ps, taxa_name = " ") {
  # ps is the phyloseq object
  # df is a dataframe with the significant taxa where the rownames are the
  # ASVs
  library(gridExtra)
  library(corncob)

  # List to store plots
  plots <- list()
  # soil_type <- unique(sample_data(ps)$Soil_type)
  # Loop over each taxa
  for (taxa in df$ASV) {
    print(paste0("Performing bbdml ", taxa))
    # Perform bbdml analysis
    corncob_da <- bbdml(
      formula = as.formula(paste(taxa, " ~ corncob_groups")),
      phi.formula = ~corncob_groups,
      data = ps
    )
    # extract the intercept pvalue
    summary_da <- summary(corncob_da)
    pvalue_intercept <- summary_da$coefficients[1, 4]


    # taxa_name <- get_taxa_name(taxa, ps)
    p <- plot(corncob_da, color = "corncob_groups", total = TRUE, B = 50) +
      theme(axis.text.x = element_blank(), legend.position = "None") +
      ggtitle(paste(taxa, "  ", format(pvalue_intercept, digits = 3)))

    plots[[taxa]] <- p
  }

  print(do.call(gridExtra::grid.arrange, c(plots, ncol = 2, top = paste0("", ""))))
  # Arrange plots
  png(filename = paste0("images/bbdml", taxa_name, ".png"), width = 2750, height = 1500 * length(plots) / 6, res = 300)
  print(do.call(gridExtra::grid.arrange, c(plots, ncol = 2, top = paste0(taxa_name))))
  dev.off()
}



# Pseudomonas -------------------------------------------------------------


ps.pseudo <- subset_taxa(ps.rarefied, Genus == "Pseudomonas")
data_bars.pseudo <- MakeDataBarsOfTaxon(ps.pseudo, grouping = "CPS")
number_of_phyla <- length(unique(data_bars.pseudo$OTU))
colorPalette1 <- colorRampPalette(brewer.pal(9, "Set1"))(number_of_phyla)
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

my_palette <- setNames(getPalette(number_of_phyla), unique(data_bars.pseudo$OTU))

p.pseudo <- bar_plot_of_taxa(data_bars.pseudo, positions = positions1, colorPalette = colorPalette1) +
  ggtitle("Pseudomonas")
print(p.pseudo)
png(filename = "images/pseudomonas bar plot RA.png", height = 3000, width = 4500, res = 450)
print(p.pseudo)
dev.off()

p.pseudo.counts <- bar_plot_of_taxa(data_bars.pseudo, positions = positions1, colorPalette = colorPalette1, abundance_type = "Counts") +
  ggtitle("Pseudomonas ") +
  ylab("Number of reads")
print(p.pseudo.counts)

pp.pseudo <- ggpubr::ggarrange(p.pseudo, NULL, p.pseudo.counts,
  ncol = 3, widths = c(30, 2, 30)
)
# widths = c(17, 2, 30))
# save pseudomonas bar plot
png(filename = "images/Pseudomonas ASV barplot.png", height = 3000, width = 4500, res = 450)
print(p.pseudo.counts)
dev.off()
p.a.p <- plot_alpha_diversity(ps = ps.pseudo) +
  ggtitle("Pseudomonas")
print(p.a.p)
png(file = "images/Pseudomonas alpha diversity.png", height = 3000, width = 3500, res = 450)
print(p.a.p)
dev.off()
hist(sample_sums(ps.pseudo), main = "Pseudomonas", breaks = 30)
hist(taxa_sums(ps.pseudo), main = "Pseudomonas", breaks = 100)


ps.pseudo <- change_metadata(ps.pseudo)
ps.pseudo.reduced <- subset_samples(ps.pseudo, Flower.or.leaf == "Flower") %>%
  subset_samples(Cultivar == "Coscia")
ps.pseudo.reduced
ntaxa <- length(unique(data_bars.pseudo$OTU))
taxa_list <- data.frame(ASV = unique(data_bars.pseudo$OTU), stringsAsFactors = FALSE)
# remove the  row in which the abundance < 0.01
taxa_list <- taxa_list[!grepl("< 0.01 abund.", taxa_list$ASV), ]
corncob::differentialTest(formula = )

perform_bbdml(taxa_list, ps.pseudo.reduced, "Pseudomonas")


# Lactobacillus -----------------------------------------------------------

ps.lacto <- subset_taxa(ps.rarefied, Genus == "Lactobacillus")
data_bars.lacto <- MakeDataBarsOfTaxon(ps.lacto, grouping = "CPS")
number_of_phyla <- length(unique(data_bars.lacto$OTU))
colorPalette1 <- colorRampPalette(brewer.pal(9, "Set1"))(number_of_phyla)
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

my_palette <- setNames(getPalette(number_of_phyla), unique(data_bars.lacto$OTU))

p.lacto <- bar_plot_of_taxa(data_bars.lacto, positions = positions1, colorPalette = colorPalette1) +
  ggtitle("Lactobacillus")
#+
# ggplot2::theme(legend.position = "none")
print(p.lacto)
p.lacto.counts <- bar_plot_of_taxa(data_bars.lacto, positions = positions1, colorPalette = colorPalette1, abundance_type = "Counts") +
  ggtitle("Lactobacillus ") +
  ylab("Number of reads")
print(p.lacto.counts)

pp.lacto <- ggpubr::ggarrange(p.lacto, NULL, p.lacto.counts,
  ncol = 3, widths = c(30, 2, 30)
)
# widths = c(17, 2, 30))
# save lactobacillus bar plot
png(filename = "images/Lactobacillus ASV barplot.png", height = 3000, width = 4500, res = 450)
print(p.lacto.counts)
dev.off()

p.a.l <- plot_alpha_diversity(ps = ps.lacto) +
  ggtitle("Lactobacillus")
print(p.a.l)
png(file = "images/lactobacillus alpha diversity.png", height = 3000, width = 3500, res = 450)
print(p.a.l)
dev.off()
hist(sample_sums(ps.lacto), main = "Lactobacillus", breaks = 30)
hist(taxa_sums(ps.lacto), main = "Lactobacillus", breaks = 100)


ps.lacto <- change_metadata(ps.lacto)
ps.lacto.reduced <- subset_samples(ps.lacto, Flower.or.leaf == "Flower") %>%
  subset_samples(Cultivar == "Coscia")
ps.lacto.reduced
ntaxa <- length(unique(data_bars.lacto$OTU))
taxa_list <- data.frame(ASV = unique(data_bars.lacto$OTU), stringsAsFactors = FALSE)
# remove the bottom row which is the abundance < 0.01
taxa_list <- data.frame(ASV = taxa_list[taxa_list$ASV != "< 0.01 abund.", ])

perform_bbdml(taxa_list, ps.lacto.reduced, "Lactobacillus")




# Pantoea -----------------------------------------------------------------

ps.Pantoea <- subset_taxa(ps.rarefied, Genus == "Pantoea")
data_bars.Pantoea <- MakeDataBarsOfTaxon(ps.Pantoea, grouping = "CPS")
number_of_phyla <- length(unique(data_bars.Pantoea$OTU))
colorPalette1 <- colorRampPalette(brewer.pal(9, "Set1"))(number_of_phyla)
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

my_palette <- setNames(getPalette(number_of_phyla), unique(data_bars.Pantoea$OTU))

p.Pantoea <- bar_plot_of_taxa(data_bars.Pantoea, positions = positions1, colorPalette = colorPalette1) +
  ggtitle("Pantoea")

print(p.Pantoea)
p.Pantoea.counts <- bar_plot_of_taxa(data_bars.Pantoea, positions = positions1, colorPalette = colorPalette1, abundance_type = "Counts") +
  ggtitle("Pantoea ") +
  ylab("Number of reads")
print(p.Pantoea.counts)

pp.Pantoea <- ggpubr::ggarrange(p.Pantoea, NULL, p.Pantoea.counts,
  ncol = 3, widths = c(30, 2, 30)
)

# save Pantoea bar plot
png(filename = "images/Pantoea ASV barplot.png", height = 3000, width = 4500, res = 450)
print(p.Pantoea.counts)
dev.off()

p.a.pa <- plot_alpha_diversity(ps = ps.Pantoea) +
  ggtitle("Pantoea")

png(file = "images/Pantoea alpha diversity.png", height = 3000, width = 3500, res = 450)
print(p.a.pa)


dev.off()
hist(sample_sums(ps.Pantoea), main = "Pantoea", breaks = 30)
hist(taxa_sums(ps.Pantoea), main = "Pantoea", breaks = 100)


ps.Pantoea <- change_metadata(ps.Pantoea)
ps.Pantoea.reduced <- subset_samples(ps.Pantoea, Flower.or.leaf == "Flower") %>%
  subset_samples(Cultivar == "Coscia")
ps.Pantoea.reduced
ntaxa <- length(unique(data_bars.Pantoea$OTU))
taxa_list <- data.frame(ASV = unique(data_bars.Pantoea$OTU), stringsAsFactors = FALSE)
# remove the bottom row which is the abundance < 0.01
taxa_list <- data.frame(ASV = taxa_list[taxa_list$ASV != "< 0.01 abund.", ])

perform_bbdml(taxa_list, ps.Pantoea.reduced, "Pantoea")



# Sphingomonas ------------------------------------------------------------

ps.Sphingomonas <- subset_taxa(ps.rarefied, Genus == "Sphingomonas")
data_bars.Sphingomonas <- MakeDataBarsOfTaxon(ps.Sphingomonas, grouping = "CPS")
number_of_phyla <- length(unique(data_bars.Sphingomonas$OTU))
colorPalette1 <- colorRampPalette(brewer.pal(9, "Set1"))(number_of_phyla)
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

my_palette <- setNames(getPalette(number_of_phyla), unique(data_bars.Sphingomonas$OTU))

p.Sphingomonas <- bar_plot_of_taxa(data_bars.Sphingomonas, positions = positions1, colorPalette = colorPalette1) +
  ggtitle("Sphingomonas")
#+
# ggplot2::theme(legend.position = "none")
print(p.Sphingomonas)
p.Sphingomonas.counts <- bar_plot_of_taxa(data_bars.Sphingomonas, positions = positions1, colorPalette = colorPalette1, abundance_type = "Counts") +
  ggtitle("Sphingomonas ") +
  ylab("Number of reads")
print(p.Sphingomonas.counts)

pp.Sphingomonas <- ggpubr::ggarrange(p.Sphingomonas, NULL, p.Sphingomonas.counts,
  ncol = 3, widths = c(30, 2, 30)
)

# save Sphingomonas bar plot
png(filename = "images/Sphingomonas ASV barplot.png", height = 3000, width = 4500, res = 450)
print(p.Sphingomonas.counts)
dev.off()



p.a.s <- plot_alpha_diversity(ps = ps.Sphingomonas) +
  ggtitle("Sphingomonas")

png(file = "images/Sphingomonas alpha diversity.png", height = 3000, width = 3500, res = 450)
print(p.a.s)
dev.off()
hist(sample_sums(ps.Sphingomonas), main = "Sphingomonas", breaks = 30)
hist(taxa_sums(ps.Sphingomonas), main = "Sphingomonas", breaks = 100)


ps.Sphingomonas <- change_metadata(ps.Sphingomonas)
ps.Sphingomonas.reduced <- subset_samples(ps.Sphingomonas, Flower.or.leaf == "Flower") %>%
  subset_samples(Cultivar == "Coscia")
ps.Sphingomonas.reduced
ntaxa <- length(unique(data_bars.Sphingomonas$OTU))
taxa_list <- data.frame(ASV = unique(data_bars.Sphingomonas$OTU), stringsAsFactors = FALSE)
# remove the bottom row which is the abundance < 0.01
taxa_list <- data.frame(ASV = taxa_list[taxa_list$ASV != "< 0.01 abund.", ])

perform_bbdml(taxa_list, ps.Sphingomonas.reduced, "Sphingomonas")



# Erwiniaceae -------------------------------------------------------------
ps.Erwiniaceae <- ps.rarefied %>%
  subset_taxa(Genus == "Erwinia")
ps.Erwiniaceae <- filter_taxa(ps.Erwiniaceae, function(x) sum(x) > 6, TRUE)
taxa_sums(ps.Erwiniaceae)
data_bars.Erwiniaceae <- MakeDataBarsOfTaxon(ps.Erwiniaceae, grouping = "CPS", threshold = 0)
number_of_phyla <- length(unique(data_bars.Erwiniaceae$OTU))
colorPalette1 <- colorRampPalette(brewer.pal(9, "Set1"))(number_of_phyla)

p.Erwiniaceae <- bar_plot_of_taxa(data_bars.Erwiniaceae, positions = positions1, colorPalette = colorPalette1) +
  ggtitle("Erwiniaceae")
#+
# ggplot2::theme(legend.position = "none")
print(p.Erwiniaceae)
p.Erwiniaceae.counts <- bar_plot_of_taxa(data_bars.Erwiniaceae, positions = positions1, colorPalette = colorPalette1, abundance_type = "Counts") +
  ggtitle("Erwiniaceae") +
  ylab("Number of reads")
print(p.Erwiniaceae.counts)

pp.Erwiniaceae <- ggpubr::ggarrange(p.Erwiniaceae, NULL, p.Erwiniaceae.counts,
  ncol = 3, widths = c(30, 2, 30)
)
# widths = c(17, 2, 30))
# save Erwiniaceae bar plot
png(filename = "images/Erwiniaceae ASV barplot.png", height = 3000, width = 4500, res = 450)
print(p.Erwiniaceae.counts)
dev.off()
png(filename = "images/Erwiniaceae ASV barplot RA.png", height = 3000, width = 4500, res = 450)
print(p.Erwiniaceae)
dev.off()

p.a.ea <- plot_alpha_diversity(ps = ps.Erwiniaceae) +
  ggtitle("Erwiniaceae")

png(file = "images/Erwiniaceae alpha diversity.png", height = 3000, width = 3500, res = 450)
print(p.a.ea)
dev.off()
hist(sample_sums(ps.Erwiniaceae), main = "Erwiniaceae", breaks = 30)
hist(taxa_sums(ps.Erwiniaceae), main = "Erwiniaceae", breaks = 100)


ps.Erwiniaceae <- change_metadata(ps.Erwiniaceae)
ps.Erwiniaceae.reduced <- subset_samples(ps.Erwiniaceae, Flower.or.leaf == "Flower") %>%
  subset_samples(Cultivar == "Coscia")
ps.Erwiniaceae.reduced
ntaxa <- length(unique(data_bars.Erwiniaceae$OTU))
taxa_list <- data.frame(ASV = unique(data_bars.Erwiniaceae$OTU), stringsAsFactors = FALSE)
# remove the bottom row which is the abundance < 0.01
taxa_list <- data.frame(ASV = taxa_list[taxa_list$ASV != "< 0.01 abund.", ])

perform_bbdml(taxa_list, ps.Erwiniaceae.reduced, "Erwiniaceae")
erwiniaceae_da <- corncob::differentialTest(
  formula = ~corncob_groups,
  phi.formula = ~1,
  formula_null = ~1,
  phi.formula_null = ~1,
  test = "Wald", boot = FALSE,
  data = ps.Erwiniaceae.reduced,
  fdr_cutoff = 0.1
)

plot(erwiniaceae_da)



# Combine all barplot and alpha figures into one  Fig S6 and Fig 5--------------------------
p.pseudo.counts
print(p.a.p)


print(p.lacto.counts)
print(p.a.l)

print(p.Pantoea.counts)
print(p.a.pa)


print(p.Sphingomonas.counts)

print(p.a.s)

print(p.Erwiniaceae.counts)
p.a.ea <- p.a.ea + theme(legend.position = "none")
print(p.a.ea)

## For figure S6

png(
  file = "images/Fig S6 Pantoea lactobacillus and sphingomonas.png",
  height = 4500, width = 5500, res = 450
)
ggpubr::ggarrange(NULL, p.Sphingomonas.counts, NULL, p.Pantoea.counts, NULL, p.lacto.counts,
  NULL, p.a.s, NULL, p.a.pa, NULL, p.a.l,
  ncol = 6, nrow = 2,
  widths = c(2, 10, 1, 10, 1, 10),
  labels = c("", "A", "", "B", "", "C", "", "D", "", "E", "", "F")
)
dev.off()

# fig S7
ggarrange(NULL, p.a.ea, NULL, p.a.p,
  ncol = 4
)

png(
  file = "images/FigS7 pseudomonas and erwiniaceae alpha diversities",
  height = 2500, width = 3500, res = 450
)
# gridExtra::grid.arrange(empty, p.a.ea, empty, p.a.p,
#                        ncol = 4,
#                        widths = c(2,10,1,10))
ggpubr::ggarrange(NULL, p.a.ea, NULL, p.a.p,
  ncol = 4,
  widths = c(2, 10, 1, 10), labels = c("", "A", "", "B")
)

dev.off()




# Figure 5 barplots of Erwiniaceae and Pseudomonas ------------------------
p.Erwiniaceae
p.pseudo
erwiniaceae_heatmap_grob <- ggplotify::as.ggplot(erwiniaceae_heatmap)

pseudo_erwinia_plots <- ggpubr::ggarrange(NULL, p.Erwiniaceae, NULL, p.pseudo,
  ncol = 4, nrow = 1,
  labels = c("", "B", "", "C"),
  widths = c(2, 10, 1, 10)
)
print(pseudo_erwinia_plots)
erwinia_plots <- ggpubr::ggarrange(NULL, NULL, erwiniaceae_heatmap_grob, NULL,
  ncol = 4, nrow = 1,
  labels = c("", "A", "", ""),
  widths = c(1, 0.5, 10, 1)
)

print(erwinia_plots)

png(filename = "images/Fig 5.png", width = 3500, height = 3500, res = 450)
ggpubr::ggarrange(erwinia_plots, NULL, pseudo_erwinia_plots,
  nrow = 3,
  heights = c(15, 2, 15)
)
dev.off()

png(filename = "images/Fig 5 B C.png", width = 3500, height = 3500, res = 450)
print(pseudo_erwinia_plots)
dev.off()
png(filename = "images/Fig 5 B C wide.png", width = 4500, height = 3500, res = 450)
print(pseudo_erwinia_plots)
dev.off()

png(filename = "images/Fig 5 A.png", width = 4500, height = 3500, res = 450)
print(erwinia_plots)
dev.off()


## notes: add the number of reads from each bar to a legend or at the top of the bars
## use the pseudo.bars with the absolute counts.

# correlation between shannon and the number of reads ---------------------

alpha_d.pseudo <- estimate_richness(ps.pseudo, measures = "Shannon")
alpha_d.pseudo$Sample <- rownames(alpha_d.pseudo)
alpha_d.pseudo <- cbind(data.frame(reads = sample_sums(ps.pseudo)), alpha_d.pseudo)
# check to see if the rownames of both dataframes match
identical(rownames(alpha_d.pseudo), alpha_d.pseudo$Sample)
corr <- cor.test(alpha_d.pseudo$reads, alpha_d.pseudo$Shannon, methos = "pearson")

p1.p <- ggplot(alpha_d.pseudo, aes(x = reads, y = Shannon)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    title = "Pseudomonas",
    subtitle = paste("Pearson R =", round(corr$estimate, digits = 3), "\np-value =", round(corr$p.value, digits = 3))
  )


alpha_d.lacto <- estimate_richness(ps.Lactobacillus, measures = "Shannon")
alpha_d.lacto$Sample <- rownames(alpha_d.lacto)
alpha_d.lacto <- cbind(data.frame(reads = sample_sums(ps.Lactobacillus)), alpha_d.lacto)
# check to see if the rownames of both dataframes match
identical(rownames(alpha_d.lacto), alpha_d.lacto$Sample)
corr <- cor.test(alpha_d.lacto$reads, alpha_d.lacto$Shannon, methos = "pearson")

p1.l <- ggplot(alpha_d.lacto, aes(x = reads, y = Shannon)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    title = "Lactobacillus",
    subtitle = paste("Pearson R =", round(corr$estimate, digits = 3), "\np-value =", format(corr$p.value, digits = 3))
  )
print(p1.l)


alpha_d.Pantoea <- estimate_richness(ps.Pantoea, measures = "Shannon")
alpha_d.Pantoea$Sample <- rownames(alpha_d.Pantoea)
alpha_d.Pantoea <- cbind(data.frame(reads = sample_sums(ps.Pantoea)), alpha_d.Pantoea)
# check to see if the rownames of both dataframes match
identical(rownames(alpha_d.Pantoea), alpha_d.Pantoea$Sample)
corr <- cor.test(alpha_d.Pantoea$reads, alpha_d.Pantoea$Shannon, methos = "pearson")

p1.pa <- ggplot(alpha_d.Pantoea, aes(x = reads, y = Shannon)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    title = "Pantoea",
    subtitle = paste("Pearson R =", round(corr$estimate, digits = 3), "\np-value =", format(corr$p.value, digits = 3))
  )
print(p1.pa)


alpha_d.Sphingomonas <- estimate_richness(ps.Sphingomonas, measures = "Shannon")
alpha_d.Sphingomonas$Sample <- rownames(alpha_d.Sphingomonas)
alpha_d.Sphingomonas <- cbind(data.frame(reads = sample_sums(ps.Sphingomonas)), alpha_d.Sphingomonas)
# check to see if the rownames of both dataframes match
identical(rownames(alpha_d.Sphingomonas), alpha_d.Sphingomonas$Sample)
corr <- cor.test(alpha_d.Sphingomonas$reads, alpha_d.Sphingomonas$Shannon, methos = "pearson")

p1.s <- ggplot(alpha_d.Sphingomonas, aes(x = reads, y = Shannon)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    title = "Sphingomonas",
    subtitle = paste("Pearson R =", round(corr$estimate, digits = 3), "\np-value =", format(corr$p.value, digits = 3))
  )
print(p1.s)


p.pp <- ggpubr::ggarrange(p1.p, p1.s, p1.pa, p1.l,
  ncol = 2, nrow = 2
)
png(file = "images/Correlation between the number of reads and shannon diversity.png", height = 3000, width = 3500, res = 450)
print(p.pp)
dev.off()


# Export ASV tables for python venn diagrams ----------------------------------------
ps.rarefied <- readRDS(file = "Data/phyloseq_all_data_rarefied_noB")
ps.Erwiniaceae <- ps.rarefied %>%
  subset_taxa(Genus == "Erwinia")
ASV.e <- data.frame(otu_table(ps.Erwiniaceae), check.names = FALSE)

ps.Sphingomonas <- subset_taxa(ps.rarefied, Genus == "Sphingomonas")
ASV.s <- data.frame(otu_table(ps.Sphingomonas), check.names = FALSE)

ps.Pantoea <- subset_taxa(ps.rarefied, Genus == "Pantoea")
ASV.p <- data.frame(otu_table(ps.Pantoea), check.names = FALSE)

ps.lacto <- subset_taxa(ps.rarefied, Genus == "Lactobacillus")
ASV.l <- data.frame(otu_table(ps.lacto), check.names = FALSE)
ps.pseudo <- subset_taxa(ps.rarefied, Genus == "Pseudomonas")
ASV.pseudo <- data.frame(otu_table(ps.pseudo), check.names = FALSE)

write.table(ASV.pseudo, file = "Data/Pseudomonas_abundance-table.csv", sep = ",", col.names = TRUE, row.names = TRUE)
write.table(ASV.s, file = "Data/Sphingomonas_abundance-table.csv", sep = ",", col.names = TRUE, row.names = TRUE)
write.table(ASV.p, file = "Data/Pantoea_abundance-table.csv", sep = ",", col.names = TRUE, row.names = TRUE)
write.table(ASV.l, file = "Data/Lactobacillus_abundance-table.csv", sep = ",", col.names = TRUE, row.names = TRUE)
write.table(ASV.e, file = "Data/Erwiniacea_abundance-table.csv", sep = ",", col.names = TRUE, row.names = TRUE)
