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
library(patchwork)
library(ggvenn)


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
      annotation_text <- paste("p-value =", format(kruskal_result$p.value, digits = 3))
    } else {
      annotation_text <- paste("p-value =", format(kruskal_result$p.value, digits = 3))
    }

    p <- ggplot(df, aes(x = Group, y = value, fill = Group)) +
      geom_boxplot() +
      labs(y = ylabel, x = "") +
      theme_minimal() +
      scale_x_discrete(limits = positions1) +
      scale_fill_manual(values = getPalette(length(unique(df$Group))), breaks = positions1) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), legend.position = "none") +
      annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, label = annotation_text, size = 5)

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

## the venn diagrams in the paper were made in python. these are the diagrams used for the 
## ISM 2024 poster. 
create_asv_venn <- function(physeq, date1, date2, date3) {
  # Function to get ASVs for a specific date
  get_asvs_for_date <- function(physeq, date) {
    # Get the sample data
    sample_data <- data.frame(sample_data(physeq))
    
    # Find samples matching the date
    samples_for_date <- rownames(sample_data)[sample_data$Short.collection.date == date]
    
    # Subset the OTU table
    otu_table_filtered <- otu_table(physeq)[, samples_for_date, drop = FALSE]
    
    # Get ASVs with non-zero abundance
    asvs <- rownames(otu_table_filtered)[rowSums(otu_table_filtered) > 0]
    
    return(asvs)
  }
  
  # Get ASVs for each date
  asvs_date1 <- get_asvs_for_date(physeq, date1)
  asvs_date2 <- get_asvs_for_date(physeq, date2)
  asvs_date3 <- get_asvs_for_date(physeq, date3)
  
  # Create a list for the Venn diagram
  asv_list <- list(
    asvs_date1,
    asvs_date2,
    asvs_date3
  )
  
  # Name the list elements with the actual dates
  names(asv_list) <- c(date1, date2, date3)
  
  # Create the Venn diagram
  venn_plot <- ggvenn(
    asv_list,
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
    stroke_size = 0.5,
    set_name_size = 4, 
    show_percentage = FALSE
  )
  
  # Print some statistics
  cat("Number of ASVs in", date1, ":", length(asvs_date1), "\n")
  cat("Number of ASVs in", date2, ":", length(asvs_date2), "\n")
  cat("Number of ASVs in", date3, ":", length(asvs_date3), "\n")
  
  # Return the plot
  return(venn_plot)
}

# check 
venn_plot <- create_asv_venn(ps.Sphingomonas, "21.04", "11.04", "24.03")

# Display the plot
print(venn_plot)


# define all phyloseq objects and place them in a list --------------------

ps.pseudo <- subset_taxa(ps.rarefied, Genus == "Pseudomonas")
ps.lacto <- subset_taxa(ps.rarefied, Genus == "Lactobacillus")
ps.Pantoea <- subset_taxa(ps.rarefied, Genus == "Pantoea")
ps.Sphingomonas <- subset_taxa(ps.rarefied, Genus == "Sphingomonas")
ps.Erwiniaceae <- ps.rarefied %>%
  subset_taxa(Genus == "Erwinia")
ps_list_taxons <- list(ps.Erwiniaceae, ps.pseudo, ps.lacto, ps.Pantoea, ps.Sphingomonas)
tag_list <- c("A", "A", "B", "B", "A", "B", "D", "E", "G", "H")
title_list <- c("Erwinia amylovora", "Pseudomonas", "Lactobacillus", "Pantoea", "Sphingomonas")


# loop over them and make the plots ---------------------------------------
# we are making two plots per phyloseq object, one for the bar plot and the other for alpha diversity (shannon)

plot_list_all <- list()
j = 1
for(i in 1:length(ps_list_taxons)){
  ps_current <- ps_list_taxons[[i]]
  current_title <- title_list[[i]]
  data_bars.current <- MakeDataBarsOfTaxon(ps_current, grouping = "CPS")
  number_of_phyla <- length(unique(data_bars.current$OTU))
  colorPalette1 <- colorRampPalette(brewer.pal(9, "Set1"))(number_of_phyla)
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  
  my_palette <- setNames(getPalette(number_of_phyla), unique(data_bars.current$OTU))
  p.current.counts <- bar_plot_of_taxa(data_bars.current, positions = positions1, colorPalette = colorPalette1, abundance_type = "Counts") +
    ggtitle(current_title) +
    ylab("Number of reads") +
    labs(tag = tag_list[[j]])
  plot_list_all[[j]] <- p.current.counts
  j = j + 1
  p.alpha <- plot_alpha_diversity(ps = ps_current) +
    ggtitle(current_title) +
    labs(tag = tag_list[[j]])
  plot_list_all[[j]] <- p.alpha
  j = j + 1
  
}
# check
wrap_plots(plot_list_all)

# Combine all barplot and alpha figures into one  Fig S7 -------------------------

plot_list_s7 <- list(plot_list_all[[5]], plot_list_all[[6]], plot_list_all[[7]], plot_list_all[[8]], plot_list_all[[9]], plot_list_all[[10]])
plot_s7 <- wrap_plots(plot_list_s7, ncol = 2) &
  theme(plot.margin = margin(5, 5, 5, 25))  # Left padding only for all plots
plot_list_s8 <- list(plot_list_all[[2]], plot_list_all[[4]])
plot_s8 <- wrap_plots(plot_list_s8, ncol = 2) &
  theme(plot.margin = margin(5, 5, 5, 25))  # Left padding only for all plots
print(plot_s8)

png(file = "images/Fig S7 Pantoea lactobacillus and sphingomonas.png",
  height = 5500, width = 3500, res = 450)
print(plot_s7)
dev.off()

png(
  file = "images/FigS8 pseudomonas and erwiniaceae alpha diversities.png",
  height = 2500, width = 3500, res = 400)
print(plot_s8) 
dev.off()




# Combine the bar and alphs plots for Erwinia, pseudomonas and sph --------
# for poster
# Assuming 'my_phyloseq' is your phyloseq object
venn_plot.p <- create_asv_venn(ps.pseudo,  "21.04", "11.04", "24.03")
venn_plot.e <- create_asv_venn(ps.Erwiniaceae,  "21.04", "11.04", "24.03")
venn_plot.s <- create_asv_venn(ps.Sphingomonas,  "21.04", "11.04", "24.03")
venn_plot <- create_asv_venn(ps.Sphingomonas, "21.04", "11.04", "24.03")
# Display the plot
print(venn_plot)


# plot_list_all[[1]] plot_list_all[[2]]
# p.pseudo p.a.p
# p.Sphingomonas p.a.s
p.a.p <- plot_alpha_diversity(ps = ps.pseudo) 
print(p.a.p)
p.a.s <- plot_alpha_diversity(ps = ps.Sphingomonas)
p.a.ea <- plot_alpha_diversity(ps = ps.Erwiniaceae)

# Remove x-axis tick labels from the top four plots
p.Erwiniaceae1 <- plot_list_all[[1]] + theme(axis.text.x = element_blank())
p.a.ea1 <- p.a.ea + theme(axis.text.x = element_blank())
p.pseudo1 <- plot_list_all[[3]] + theme(axis.text.x = element_blank())
p.a.p1 <- p.a.p + theme(axis.text.x = element_blank())

# Arrange the plots using patchwork
combined_plot <- (p.Erwiniaceae1 + p.a.ea1) / 
  (p.pseudo1 + p.a.p1) / 
  (p.Sphingomonas + p.a.s1) +
  plot_layout(ncol = 1)


combined_plot <- wrap_plots(p.Erwiniaceae1, p.a.ea1, venn_plot.e,
                            p.pseudo1, p.a.p1, venn_plot.p,
                            p.Sphingomonas, p.a.s, venn_plot.s, 
                            ncol = 3)

# Display the combined plot
print(combined_plot)

png("images/strain analysis for poster with venn.png", height = 6000, width = 5200, res = 550)
print(combined_plot)
dev.off()




# Figure 4 barplots of Erwiniaceae and Pseudomonas ------------------------
p.Erwiniaceae <- plot_list_all[[1]]
p.pseudo <- plot_list_all[[3]]
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
