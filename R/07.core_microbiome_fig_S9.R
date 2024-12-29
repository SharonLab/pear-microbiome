
# this is the script for the core microbiome analysis. 
# It generates the core microbiome heatmaps (Fig S9), 
# as well as other figures that did not go into the manuscript: the heatmaps of all the taxonomic levels, and barplots for the core microbiome. 

library(microbiome)
library(phyloseq)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(patchwork)

source("R/02.general_definitions.R")
source("R/02.helper_functions.R")

ps.all <- readRDS("Data/phyloseq_no_buchnera") %>%
subset_taxa(Genus != "Erwinia") %>%
  subset_samples(Flower.or.leaf == "Flower") %>%
  subset_samples(Cultivar == "Coscia")


## we will run these first to save on running time
ps.all.species <- ps.all %>%
  tax_glom(taxrank = "Species")
ps.all.genera <- ps.all %>%
  tax_glom(taxrank = "Genus")
  
ps.all.family <- ps.all %>%
  tax_glom(taxrank = "Family") 
ps.all.order <- ps.all %>% 
  tax_glom(taxrank = "Order") 
ps.all.class <- ps.all %>%
  tax_glom(taxrank = "Class")
ps.all.phylum <- ps.all %>%
  tax_glom(taxrank = "Phylum") 

ps.list <- list(ps.all.phylum, ps.all.class, ps.all.order, ps.all.family, ps.all.genera)
names(ps.list) <- c("Phylum","Class", "Order","Family", "Genus")

## add asvs to the analysis of the number of samples in which each taxonomic rank etc. 
## make a figure of the six taxonomic ranks
## show only the top 200 taxons for each one
## 8 panels: ASV, genus, family
## order class phylum
## order_zoom 
## the eigths one is the numer of reads for the top genera and add labels of in how many samples does each point appear. 
## the ninth - the same as 8 for order
# Functions ---------------------------------------------------------------

get_top_order_sums <- function(physeq_obj) {
  # Calculate the sum of taxa in the phyloseq object
  top_order_sums <- data.frame(taxa_sums(physeq_obj))
  top_order_sums <- top_order_sums %>% rename(Num.reads = `taxa_sums.ps.order.top.`)
  # Add hash as a column from the row names
  top_order_sums$hash <- rownames(top_order_sums)
  
  # Extract the taxonomy table and add hash as a column from the row names
  top_order_taxa <- data.frame(tax_table(physeq_obj))
  top_order_taxa$hash <- rownames(top_order_taxa)
  
  # Merge the sums and taxonomy tables by hash
  top_order_sums <- merge(top_order_sums, top_order_taxa, by = "hash")
  
  # Return the resulting dataframe
  return(top_order_sums)
}

get_sample_counts <- function(physeq_obj) {
  # Ensure the input is a phyloseq object
  if (!inherits(physeq_obj, "phyloseq")) {
    stop("Input must be a phyloseq object.")
  }
  # find out what taxonomy we are using
  taxrank <- get_taxrank(physeq_obj)
  # Extract the taxonomy table and convert it to a dataframe
  genera_taxonomy <- as.data.frame(tax_table(physeq_obj))
  
  # Add a 'hash' column using row names
  genera_taxonomy$hash <- rownames(genera_taxonomy)
  
  
  # Extract the OTU table, convert to dataframe, and create a prevalence (presence/absence) table
  genera_prevalence <- as.data.frame(otu_table(physeq_obj)) > 0
  
  # Calculate the number of samples each taxon appears in
  taxon_sample_counts <- as.data.frame(rowSums(genera_prevalence))
  
  # Add a 'hash' column using row names
  taxon_sample_counts$hash <- rownames(taxon_sample_counts)
  
  # Merge the prevalence counts with the taxonomy information using the 'hash' column
  taxon_sample_counts <- merge(taxon_sample_counts, genera_taxonomy, by = "hash")
  
  # Rename the prevalence count column to 'Num.reads'
  taxon_sample_counts <- taxon_sample_counts %>%
    rename(Num.reads = `rowSums(genera_prevalence)`)
  
  # Order the results by 'Num.reads'
  taxon_sample_counts <- taxon_sample_counts[order(taxon_sample_counts$Num.reads, decreasing = TRUE), ]
  
  return(taxon_sample_counts)
}

make_taxon_count_plot <- function(sample_counts_df, frac = "all"){
  # frac is either "all", which puts all of the taxons in one plot without names, 
  # or a number, which indicates how many taxons to put on the plot (too many are unintelligible)
  # Find out the taxonomic level of the sample_counts_df
  # Find the index of the last column that contains at least one non-NA value
  last_non_na_col <- max(which(colSums(!is.na(sample_counts_df)) > 0))
  # Retrieve the name of that column
  current_taxon <- names(sample_counts_df)[last_non_na_col]
  if(frac != "all"){
    sample_counts_df <- sample_counts_df[order(sample_counts_df$Num.reads, decreasing = TRUE),] %>%
      head(n=frac)
  }
  sample_counts_df$fraction <- sample_counts_df$Num.reads/141
  
  p <- ggplot(sample_counts_df, aes(x = reorder(hash, -Num.reads), y = fraction)) + 
    geom_point() +
    theme_minimal() +
    xlab("\n") +
    ylab("") +
    ggtitle(current_taxon) +
    scale_x_discrete(limits = sample_counts_df$hash,
                      labels = sample_counts_df[[current_taxon]]) 
  
    if(frac == "all"){
      p <- p + 
        theme(axis.text.x = element_blank(),  # Remove x-axis tick labels
            axis.ticks.x = element_blank(), 
            panel.grid.major = element_blank(),  # Remove major grid lines
            panel.grid.minor = element_blank()   # Remove minor grid lines
            ) 
    } else {
      p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis tick labels 45 degrees
                 axis.ticks.x = element_blank())   # Optionally remove x-axis ticks
    }
   print(p) 
   return(p)
}

# test
taxon_sample_counts <- get_sample_counts(ps.all.genera)
make_taxon_count_plot(taxon_sample_counts)

#test with taxon names on the x axis
taxon_sample_counts <- get_sample_counts(ps.all.order)
make_taxon_count_plot(taxon_sample_counts, 15)


#  loop over the phyloseq list --------------------------------------------

plot_list <- list()
j = 1
for(i in seq_along(ps.list)){
  temp.ps <- ps.list[[i]]
  taxon_sample_counts <- get_sample_counts(temp.ps)
  p <- make_taxon_count_plot(taxon_sample_counts) 
  plot_list[[j]] <- p
  j <- j+1
  p1 <- make_taxon_count_plot(taxon_sample_counts, 15) 
  plot_list[[j]] <- p1
  j <- j+1
}

## the ASVs need a separate plot, since the names repeat themselves so we will manipulate the data a bit
taxon_sample_counts <- get_sample_counts(ps.all)
taxon_sample_counts$Species <- taxon_sample_counts$hash
p <- make_taxon_count_plot(taxon_sample_counts)

plot_list[[11]] <- p 


taxon_sample_counts <- get_sample_counts(ps.all.species)
p <- make_taxon_count_plot(taxon_sample_counts, 15)
plot_list[[12]] <- p

# ## add plots of the top 15 orders and top 15 phyla with names
# ps_list_short <- list(ps.all.phylum, ps.all.order)
# plot_list_2 <- list()
# for(i in seq_along(ps_list_short)){
#   temp.ps <- ps_list_short[[i]]
#   taxon_sample_counts <- get_sample_counts(temp.ps)
#   p <- make_taxon_count_plot(taxon_sample_counts, 15) 
#   plot_list_2[[i]] <- p
# }
# 
# plot_list <- c(plot_list_2, plot_list)



# combined_plot <-  wrap_plots(rev(plot_list), ncol = 2)

combined_plot <-  wrap_plots(plot_list[[11]], plot_list[[9]], plot_list[[7]],
           plot_list[[12]], plot_list[[10]], plot_list[[8]],
           plot_list[[5]], plot_list[[3]], plot_list[[1]],
           plot_list[[6]], plot_list[[4]], plot_list[[2]],
           ncol = 3) +
  plot_layout() +
  plot_annotation()
print(combined_plot)
library(patchwork)
library(ggplot2)

# Create the combined plot with specific tags
combined_plot <- wrap_plots(
  plot_list[[11]] + labs(tag = "A"),
  plot_list[[9]] + labs(tag = "B"),
  plot_list[[7]] + labs(tag = "C"),
  plot_list[[12]] + labs(tag = "D"),
  plot_list[[10]] + labs(tag = "E"),
  plot_list[[8]] + labs(tag = "F"),
  plot_list[[5]] + labs(tag = "G"),
  plot_list[[3]] + labs(tag = "H"),
  plot_list[[1]] + labs(tag = "I"),
  plot_list[[6]] + labs(tag = "J"),
  plot_list[[4]] + labs(tag = "K"),
  plot_list[[2]] + labs(tag = "L"),
  ncol = 3
) + 
  plot_annotation(
    title = "Core microbiome of Coscia flowers"
  ) & 
  theme(plot.margin = margin(5, 5, 5, 25))  # Left padding only for all plots

print(combined_plot)


## keep the ratio to fit an a4 paper
# png("images/core microbiome number of samples for each taxon.png", height = 5000, width = 3750, res = 450)
# print(combined_plot)
# dev.off()


# png("images/core microbiome number of samples for each taxon coscia only.png", height = 5000, width = 3750, res = 450)
# print(combined_plot)
# dev.off()


png("images/s9 core microbiome number of samples for each taxon coscia only.png", height = 6000, width = 5100, res = 450)
print(combined_plot)
dev.off()



# core heatmap ------------------------------------------------------------
# code taken from: https://microbiome.github.io/tutorials/Core.html
# Core with compositionals:
library(RColorBrewer)
library(reshape)
library(microbiome)
library(patchwork)

# Also define gray color palette
gray <- gray(seq(0,1,length=5))

phyloseq_objects <- list(ps.all.species, ps.all.genera, ps.all.family, ps.all.order, ps.all.class, ps.all.phylum)

# Define list of names for titles, extracted from the names of the objects
taxonomic_levels <- c("Species", "Genus", "Family", "Order", "Class", "Phylum")

# Loop through each phyloseq object and replace hash names with taxon names
for (i in seq_along(phyloseq_objects)) {
  # Current phyloseq object and corresponding taxonomic level
  pseq <- phyloseq_objects[[i]]
  tax_level <- taxonomic_levels[i]
  
  # Get OTU names (hashes) from the `otu_table`
  otu_names <- taxa_names(pseq)
  
  # Extract taxonomic names at the specified level and align with OTU names
  tax_table_df <- as.data.frame(tax_table(pseq))
  tax_names <- tax_table_df[otu_names, tax_level]  # Align tax names with OTUs
  
  # Replace NA values with a placeholder
  tax_names[is.na(tax_names)] <- paste0("Unclassified_", tax_level)
  
  # Create a new vector for unique tax names
  new_tax_names <- character(length(otu_names))
  name_counts <- list()  # To track occurrences of each name
  
  # Replace each hash with the taxonomic name and handle duplicates
  for (j in seq_along(otu_names)) {
    tax_name <- tax_names[j]
    
    # Check if this tax_name has been used before
    if (tax_name %in% names(name_counts)) {
      # Increment the count and create a unique name
      name_counts[[tax_name]] <- name_counts[[tax_name]] + 1
      unique_name <- paste0(tax_name, "_", name_counts[[tax_name]])
    } else {
      # First occurrence of this name
      name_counts[[tax_name]] <- 1
      unique_name <- tax_name
    }
    
    # Assign the unique name to the new name vector
    new_tax_names[j] <- unique_name
  }
  
  # Ensure unique names for taxa
  if (any(duplicated(new_tax_names))) {
    stop("Error: Duplicate taxa names found after processing.")
  }
  
  # Set the new names as OTU names in the phyloseq object
  taxa_names(pseq) <- new_tax_names
  
  # Update the phyloseq object in the list
  phyloseq_objects[[i]] <- pseq
}

# After the loop, `phyloseq_objects` now contains updated phyloseq objects
# with hash IDs replaced by taxonomic names and serial numbers for duplicates.



prevalences <- seq(.05, 1, .05)

detections <- round(10^seq(log10(0.01), log10(.2), length = 9), 3)

# Initialize an empty list to store the plots
plot_list <- list()

# Loop through each phyloseq object and create a plot
for (i in seq_along(phyloseq_objects)) {
  # Current phyloseq object and corresponding taxonomic level for the title
  pseq <- phyloseq_objects[[i]] %>%
    transform_sample_counts(function(x) x/sum(x))
  title <- taxonomic_levels[i]
  
  # Calculate dynamic detections sequence based on max abundance in pseq
  # max_abundance <- max(abundances(pseq), na.rm = TRUE)
  # detections <- seq(from = 50, to = round(max_abundance / 10, -1), by = 100)
  # this is for the low resolution image
  prevalences <- seq(.05, 1, .05)
  # 
  detections <- round(10^seq(log10(0.01), log10(.2), length = 9), 3)
  # 
  
  # Generate the plot with plot_core
  p <- plot_core(pseq, plot.type = "heatmap",
                 prevalences = prevalences,
                 detections = detections,
                 colours = rev(brewer.pal(5, "Spectral")),
                 min.prevalence = .2, horizontal = FALSE) +
    labs(title = paste("Core Microbiome -", title)) +  # Add title dynamically
    theme(axis.text.x = element_text(size=8,angle = -45, vjust = 0.5, hjust = -0.1) ,
          axis.text.y = element_text(size=8),
          axis.title = element_text(size=10),
          legend.text = element_text(size=8),
          legend.title = element_text(size=10))
  
  # Store the plot in the list with the taxonomic level as the name
  plot_list[[title]] <- p
}

wrap_plots(plot_list)
png("images/core_microbiome_heatmaps.png", height = 4500, width = 4500, res = 450)
wrap_plots(plot_list, 
           ncol = 2)
dev.off()



# Initialize an empty list to store the plots
plot_list2 <- list()
tags <- c("A", "A", "B", "B", "C", "C")
# Loop through each phyloseq object and create a plot
for (i in c(2,4,6)) {
  # Current phyloseq object and corresponding taxonomic level for the title
  pseq <- phyloseq_objects[[i]] %>%
    transform_sample_counts(function(x) x/sum(x))
  title <- taxonomic_levels[i]
  
  # Calculate dynamic detections sequence based on max abundance in pseq

  prevalences <- seq(.07, 1, .07)
  
  detections <- c(0.001, 0.002,0.005, 0.01, 0.02, 0.05, 0.1 )
  # Generate the plot with plot_core
  p <- plot_core(pseq, plot.type = "heatmap",
                 prevalences = prevalences,
                 detections = detections,
                 colours = rev(brewer.pal(5, "Spectral")),
                 min.prevalence = .2, horizontal = FALSE) +
    labs(title = title, tag = tags[[i]]) +  # Add title dynamically
    theme(axis.text.x = element_text(size=8,angle = -45, vjust =1, hjust = -0.1) ,
          axis.text.y = element_text(size=8),
          axis.title = element_text(size=10),
          legend.text = element_text(size=8),
          legend.title = element_text(size=10))
  
  # Store the plot in the list with the taxonomic level as the name
  plot_list2[[title]] <- p
}

wrap_plots(plot_list2, ncol = 3)


png("images/core_microbiome_heatmaps_low abundances genus order phyla.png", height =3500 , width = 5500, res = 450)
wrap_plots(plot_list2, 
           ncol = 3)
dev.off()

png("images/core_microbiome_heatmaps_low abundances.png", height = 5500, width = 8500, res = 550)
wrap_plots(plot_list2, 
           ncol = 3)
dev.off()

# core microbiome using the microbiome::core function ---------------------


## the parameters used are after this paper on citrus tree rhizospheres https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6244077/
# the data we use is non rarefied after the recommendations of https://www.pnas.org/doi/full/10.1073/pnas.2104429118#sec-1
# to not use rarefied data because the data is normalized to relative abundances, and rarefaction makes it likely to miss persistent rare species. 
core_microbiome <- core(ps.all.genera, detection = 0.001, prevalence = 0.75, include.lowest = TRUE)
core_microbiome.species <- core(ps.all, detection = 0.001, prevalence = 0.75, include.lowest = TRUE)
print(core_microbiome)
# Transform data to relative abundance
ps_relabund <- transform_sample_counts(ps.all, function(x) x / sum(x))


# Aggregate at the desired taxonomic level, e.g., Genus
ps_relabund_genus <- tax_glom(ps_relabund, "Genus")

sample_sums(core_microbiome)
# Aggregate data at the CPS level to get mean abundances of core taxa
ps_core_agg <- speedyseq::merge_samples2(core_microbiome, "CPS", fun_otu = mean)
sample_sums(ps_core_agg)

ps_core_genera <- tax_glom(core_microbiome, taxrank = "Genus")
temp.ps <- ps_core_genera %>%
  speedyseq::merge_samples2(group = "CPS", fun_otu = mean) 
sample_sums(temp.ps)
sample_data(temp.ps)$Sample <- rownames(sample_data(temp.ps))  # fix a phyloseq bug. it works bad also with other packages
core_microbiome_data_bars <-  temp.ps %>%
  psmelt()
# core_microbiome_data_bars <- make_data_bars1(ps_core_genera, grouping_var = "CPS", RA = FALSE)
p_bar_plot <- taxa_bar_plot(core_microbiome_data_bars) + 
  scale_x_discrete(limits = positions_cps) +
  ggtitle("Core microbiome\n Abundances\n") + 
  ylab("Abundances") +
  guides(fill=guide_legend(ncol=4) ) + theme(legend.position = "bottom") +
   ylim(0,30000)
p_bar_plot

# Relative abundances core microbiome bar plot ----------------------------

## temp.ps <- merge_samples(core_microbiome,  "CPS", fun = mean) # for some obscure reason this simple solurtion does not work. 
temp.ps <- speedyseq::merge_samples2(core_microbiome,  "CPS", fun_otu = mean)
temp.ps_RA <- transform_sample_counts(temp.ps, function(x) x/sum(x))
sample_data(temp.ps_RA)$Sample <- rownames(sample_data(temp.ps_RA))  # fix a speedyseq bug. it works bad also with other packages
core_microbiome_data_bars_RA <-  temp.ps_RA %>%
  psmelt()
abundances <- sample_sums(temp.ps) 
abundances.df <- data.frame(Sample = names(abundances), Abundance = abundances)
core_microbiome_data_bars_RA1 <- merge(abundances.df, core_microbiome_data_bars_RA, by = "Sample" )
core_microbiome_data_bars_RA1$Genus[core_microbiome_data_bars_RA1$Genus == "d__Bacteria Kingdom"] <- "Unknown Bacteria"
number_of_phyla <- length(unique(core_microbiome_data_bars_RA1$Genus))
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))(number_of_phyla)
#p_bar_plot_RA <-  
p_bar_plot_RA <- ggplot(data = core_microbiome_data_bars_RA1, aes(x = Sample, y = Abundance.y, fill = fct_reorder(Genus, Abundance.y)) ) +
  geom_col()+
  scale_fill_manual(values = rev(getPalette)) +
  theme_minimal() +
  #labs(y = "Relative Abundance\n", x = "\n", fill = taxrank) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(1, "lines"), 
    panel.grid.major = element_blank(),## remove the grid lines
    panel.grid.minor = element_blank()
  ) + 
  expand_limits(y = 1.4) +
  scale_x_discrete(limits = positions_cps) +
  ggtitle("Core microbiome\n Abundances\n") + 
  ylab("Relative abundances") + xlab(" ") +
  guides(fill=guide_legend(ncol=5, title = "Genus") ) + 
  theme(legend.position = "bottom",  axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = -0.1), legend.title = element_text(size = 10)) +
  geom_text(aes(label = ifelse(Genus == "Erwinia", round(Abundance.x), " ")), y = 1.20) 

    
p_bar_plot_RA
png("images/The genera in the core microbiome.png", height = 2500, width = 3000, res = 350)
print(p_bar_plot_RA)
dev.off()


taxa_sums(core_microbiome)
taxa_sums(ps_core_agg)
sample_sums(ps_core_agg)



# without erwinia ---------------------------------------------------------

ps.all <- readRDS("Data/phyloseq_no_buchnera")
ps.no.erwinia.core <- subset_taxa(ps_core_genera, Genus != "Erwinia")
temp.ps.abs <- ps.no.erwinia.core %>%
  speedyseq::merge_samples2(group = "CPS", fun_otu = mean)
# ps.no.erwinia.genera <- tax_glom(ps.no.erwinia, taxrank = "Genus")
temp.ps <- ps.no.erwinia.core %>%
  speedyseq::merge_samples2(group = "CPS", fun_otu = mean) %>%
  transform_sample_counts(function(x) x/sum(x))
sample_data(temp.ps)$Sample <- rownames(sample_data(temp.ps))  # fix a phyloseq bug. it works bad also with other packages
no_erwinia_data_bars <-  temp.ps %>%
  psmelt()
p_erwinia_tax_plot <- taxa_bar_plot(no_erwinia_data_bars) + 
  scale_x_discrete(limits = positions_cps) +
  ggtitle("Absolute abundances without Erwinia")
# adjust the colors to fit the figure with the erwinia (red is reserved for EA)
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))(25)
 getPalette <- getPalette[2:25]
 # plot
 p_erwinia_tax_plot <- p_erwinia_tax_plot +  ggtitle("Core microbiome abundancess,\n no Erwinia amylovora\n") + 
  ylab("Absolute abundances") +
  scale_fill_manual(values = rev(getPalette)) 
 # combine the two plots 
 

sample_sums_labels <- sample_sums(merge_samples(ps.no.erwinia.core, group = "CPS", fun = mean))
p_bar_plot +
  geom_text(aes(label = sample_sums_labels[as.character(Sample)]), vjust = -0.5)

no_erwinia_data_bars$Genus[no_erwinia_data_bars$Genus == "d__Bacteria Kingdom"] <- "Unknown Bacterium"
number_of_phyla <- length(unique(no_erwinia_data_bars$Genus))
getPalette <- colorRampPalette(brewer.pal(9, c("Set1")))(number_of_phyla)
fills <- c(sapply(c( "Pastel2", "Paired", "Set1"),
                 function(x) brewer.pal(8, x)))

fills<-c("#E41A1C","#705C83","#377EB8","lightskyblue1","palegreen2","darkgreen","darkolivegreen3","mediumblue",
         "#AA3B50",  "#984EA3",
         "#FFFF33", "#BA5E6C","#DC6E36","#FF7F00","#FFA910","#FFD421","#E1C62F","#C38E2B", 
        "#A65628","mediumorchid1","#DB728C","orchid1","bisque","#B891A5","#999999") 

fills <- c(
  "darkred", "lightcoral",
  "darkblue", "lightskyblue1",
  "mediumblue", "turquoise1",
  "royalblue", "skyblue",
  "green4", "palegreen2",
  "orange", "yellow1",
  "black", "gray70",
  "brown", "bisque2",
  "mediumorchid4", "mediumorchid1",
  "#6A3D9A", "orchid1",
  "darkorange2", "peachpuff",
  "darkorange4", "gold1",
  "olivedrab1", "darkgreen",
  "deeppink4", "pink",
  "darkviolet", "plum1",
  "navy", "lightblue",
  "darkslategray", "lightgray",
  "firebrick", "salmon",
  "seagreen", "mediumseagreen",
  "chocolate", "sandybrown",
  "maroon", "mistyrose",
  "purple4", "thistle",
  "steelblue", "lightsteelblue",
  "darkgoldenrod", "goldenrod1",
  "chartreuse4", "chartreuse",
  "midnightblue", "cadetblue",
  "sienna", "burlywood",
  "blue4", "lightcyan1",
  "darkolivegreen", "yellowgreen"
)

fills <- c(
  "#E41A1C", "#FF9999",  # Red pair
  "#4DAF4A", "#B2DF8A",  # Green pair
  "#984EA3", "#CAB2D6",  # Purple pair
  "#FF7F00", "#FDBF6F",  # Orange pair
  "#FFFF33", "#FFEDA0",  # Yellow pair
  "#A65628", "burlywood3",  # Brown pair
  "#F781BF", "#FCCDE5",  # Pink pair
  "#999999", "#CCCCCC",  # Gray pair
  "#377EB8", "#A6CEE3",  # Blue pair
  "#FC8D62", "#FFD3B6",  # Salmon pair
  "#8DA0CB", "#D4B9DA",  # Light purple pair
  "#A6D854", "#E6F5D0",  # Light green pair
  "#FFD92F", "#FFFFCC",  # Light yellow pair
  "#B3B3B3", "#D3D3D3",  # Light gray pair
  "#BC80BD", "#E6CCE5",   # Light magenta pair
  "#80B1D3", "#CCE5FF",  # Light blue pair
  "#66A61E", "#B2DF8A",  # Dark green pair
)           
     

core_microbiome_data_bars_RA1_noe <- merge(abundances.df, no_erwinia_data_bars, by = "Sample" )
# change the "Sphingomonadaceae Family" to "Sphingomonas"


#p_bar_plot_RA <-  
p_bar_plot_RA <- ggplot(data = core_microbiome_data_bars_RA1_noe, aes(x = Sample, y = Abundance.y, fill = Genus) ) +
  geom_col()+
  #scale_fill_manual(values = rev(getPalette)) +
  scale_fill_manual(values = rev(fills))+
  theme_minimal() +
  labs(y = "Relative Abundance\n", x = "\n" , fill = "Genus") +
  guides(fill = guide_legend(ncol = 1)) +
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(1, "lines"), 
    panel.grid.major = element_blank(),## remove the grid lines
    panel.grid.minor = element_blank()
  ) + 
  expand_limits(y = 1.4) +
  scale_x_discrete(limits = positions_cps) +
  ggtitle("Core microbiome\n Abundances\n") + 
  ylab("Relative abundances") +
  guides(fill=guide_legend(ncol=5) ) + 
  theme(legend.position = "bottom",  axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = -0.1)) +
  geom_text(aes(label = ifelse(Genus == "Rhodococcus", round(Abundance.x), " ")), y = 1.20) 
p_bar_plot_RA
png("images/Core microbiome without erwinia.png", height = 3000, width = 3500, res = 450)
print(p_bar_plot_RA)
dev.off()


# Erwinia reads counts in dotplot -----------------------------------------


core_microbiome_counts <- data.frame(CPS = unique(core_microbiome_data_bars_RA1_noe$Sample), 
                                     Abundance = unique(core_microbiome_data_bars_RA1_noe$Abundance.x))

No_EA_counts <- as.data.frame(sample_sums(temp.ps.abs) )
core_microbiome_counts$No_EA_counts <- No_EA_counts$`sample_sums(temp.ps.abs)`[match(core_microbiome_counts$CPS, rownames(No_EA_counts))]
core_microbiome_counts$EA <- round(core_microbiome_counts$Abundance - core_microbiome_counts$No_EA_counts)
core_microbiome_counts_long <- reshape2::melt(core_microbiome_counts, id.vars = "CPS", measure.vars = c( "EA", "Abundance"))
ggplot(core_microbiome_counts_long, aes(x = CPS, y = round(value), color = variable)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_x_discrete(limits = positions_cps) +
  theme(legend.position = "bottom",  axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = -0.1)) +
  ylab("Number of reads") +
  guides(fill=guide_legend(title=" "))
 
  

# erwinia read counts boxplots --------------------------------------------
ps.erwinia <- subset_taxa(core_microbiome, Genus == "Erwinia")
metadata <- data.frame(sample_data(core_microbiome))  
metadata$number_of_reads<- sample_sums(core_microbiome)
metadata$EA_reads <- sample_sums(ps.erwinia)
counts_long <- reshape2::melt(metadata, id.vars = "CPS", measure.vars = c("EA_reads", "number_of_reads"))

ggplot(counts_long, aes(x = CPS, y = value, color = variable)) +
  geom_boxplot() +
  #geom_jitter() +
  theme_bw() +
  scale_x_discrete(limits = positions_cps) +
  theme(legend.position = "bottom",  axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = -0.1)) +
  ylab("Number of reads") +
  guides(fill=guide_legend(title=" "))




# Check abundances of common genera - core --------------------------------
core_microbiome <- subset_samples(core_microbiome, Flower.or.leaf == "Flower")
ps_core_agg <- subset_samples(ps_core_agg, Flower.or.leaf == "Flower")
ps.pseudo.core <- subset_taxa(core_microbiome, Genus == "Pseudomonas")%>%
speedyseq::merge_samples2("CPS", fun_otu = mean)
ps.pseudo.core
sample_sums(ps.pseudo.core)
sample_sums(ps_core_agg)
plot_bar(ps.pseudo.core)
tax_table(ps.pseudo.core)

ps.pantoea.core <- subset_taxa(core_microbiome, Genus == "Pantoea")%>%
  speedyseq::merge_samples2("CPS", fun_otu = mean)
ps.pantoea.core
sample_sums(ps.pantoea.core)
plot_bar(ps.pantoea.core)
tax_table(ps.pantoea.core)


ps.sphingo.core <- subset_taxa(core_microbiome, Genus == "Sphingomonas")%>%
  speedyseq::merge_samples2("CPS", fun_otu = mean)
ps.sphingo.core
sample_sums(ps.sphingo.core)
plot_bar(ps.sphingo.core)
tax_table(ps.sphingo.core)

ps.lacto.core <- subset_taxa(core_microbiome, Genus == "Lactobacillus")%>%
  speedyseq::merge_samples2("CPS", fun_otu = mean)
ps.lacto.core
sample_sums(ps.lacto.core)

ps.pseudomonadota.core <- subset_taxa(core_microbiome, Phylum == "Pseudomonadota")%>%
  speedyseq::merge_samples2("CPS", fun_otu = mean)
ps.pseudomonadota.core
sample_sums(ps.pseudomonadota.core)
# unify all of the core icrobiome data into one data frame
sample_sums_df <- data.frame(CPS = names(sample_sums(ps.pseudomonadota.core)), 
                             pseudomonadota  = sample_sums(ps.pseudomonadota.core), 
                             Sphingomonas = sample_sums(ps.sphingo.core), 
                             Lactobacillus = sample_sums(ps.lacto.core), 
                             Pantoea = sample_sums(ps.pantoea.core), 
                             Pseudomonas = sample_sums(ps.pseudo.core), 
                             total = sample_sums(ps_core_agg))
# calculate the fractions of the key taxa out of the core microbiome
sample_sums_df_frac <- sample_sums_df %>%
  mutate(across(pseudomonadota:Pseudomonas, ~ round((. / total) *100)))


#  Do the same analysis on rarefied data-----------------------------------

ps.all.rare <- readRDS("Data/phyloseq_all_data_rarefied_noB")
core_microbiome <- core(ps.all.rare, detection = 0.01, prevalence = 0.50, include.lowest = TRUE)
print(core_microbiome)
# Transform data to relative abundance
ps_relabund <- transform_sample_counts(ps.all.rare, function(x) x / sum(x))

# Aggregate at the desired taxonomic level, e.g., Genus
ps_relabund_genus <- tax_glom(ps_relabund, "Genus")

sample_sums(core_microbiome)
# Aggregate data at the CPS level to get mean abundances of core taxa
ps_core_agg <- speedyseq::merge_samples2(core_microbiome, "CPS", fun_otu = mean)
sample_sums(ps_core_agg)

ps_core_genera <- tax_glom(core_microbiome, taxrank = "Genus")
temp.ps <- ps_core_genera %>%
  speedyseq::merge_samples2(group = "CPS", fun_otu = mean) 
sample_sums(temp.ps)
sample_data(temp.ps)$Sample <- rownames(sample_data(temp.ps))  # fix a phyloseq bug. it works bad also with other packages
core_microbiome_data_bars <-  temp.ps %>%
  psmelt()
# core_microbiome_data_bars <- make_data_bars1(ps_core_genera, grouping_var = "CPS", RA = FALSE)
p_bar_plot <- taxa_bar_plot(core_microbiome_data_bars) + 
  scale_x_discrete(limits = positions_cps) +
  ggtitle("Core microbiome\n Abundances\n") + 
  ylab("Abundances") +
  guides(fill=guide_legend(ncol=4) ) + theme(legend.position = "bottom") +
  ylim(0,2250)
p_bar_plot

# rarefied data no erwinia ------------------------------------------------


ps.no.erwinia.core <- subset_taxa(ps_core_genera, Genus != "Erwinia")
temp.ps.abs <- ps.no.erwinia.core %>%
  speedyseq::merge_samples2(group = "CPS", fun_otu = mean)
# ps.no.erwinia.genera <- tax_glom(ps.no.erwinia, taxrank = "Genus")
temp.ps <- ps.no.erwinia.core %>%
  speedyseq::merge_samples2(group = "CPS", fun_otu = mean) %>%
  transform_sample_counts(function(x) x/sum(x))
sample_data(temp.ps)$Sample <- rownames(sample_data(temp.ps))  # fix a phyloseq bug. it works bad also with other packages
no_erwinia_data_bars <-  temp.ps %>%
  psmelt()
p_erwinia_tax_plot <- taxa_bar_plot(no_erwinia_data_bars) + 
  scale_x_discrete(limits = positions_cps) +
  ggtitle("Absolute abundances without Erwinia")
# adjust the colors to fit the figure with the erwinia (red is reserved for EA)
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))(25)
getPalette <- getPalette[2:25]
# plot
p_erwinia_tax_plot <- p_erwinia_tax_plot +  ggtitle("Core microbiome abundancess,\n no Erwinia amylovora\n") + 
  ylab("Absolute abundances") +
  scale_fill_manual(values = rev(getPalette)) 
# combine the two plots 


sample_sums_labels <- sample_sums(merge_samples(ps.no.erwinia.core, group = "CPS", fun = mean))
p_bar_plot +
  geom_text(aes(label = sample_sums_labels[as.character(Sample)]), vjust = -0.5)


number_of_phyla <- length(unique(no_erwinia_data_bars$Genus))
getPalette <- colorRampPalette(brewer.pal(8, c("Set1"))[-1])(number_of_phyla)


fills <- c(sapply(c( "Pastel2", "Paired", "Set1"),
                  function(x) brewer.pal(8, x)))

fills<-c("#E41A1C","#705C83","#377EB8","lightskyblue1","palegreen2","darkgreen","darkolivegreen3","mediumblue",
         "#AA3B50",  "#984EA3",
         "#FFFF33", "#BA5E6C","#DC6E36","#FF7F00","#FFA910","#FFD421","#E1C62F","#C38E2B", 
         "#A65628","mediumorchid1","#DB728C","orchid1","bisque","#B891A5","#999999") 


fills <- c("darkred", "lightcoral",
           "darkblue", "lightskyblue1",
           "mediumblue", "turquoise1",
           "royalblue", "skyblue", 
           "green4", "palegreen2", 
           "orange", "yellow1",
           "black", "gray70",
           "brown", "bisque2",
           "mediumorchid4", "mediumorchid1", 
           "#6A3D9A","orchid1",
           "darkorange2", "peachpuff", 
           "darkorange4", "gold1", 
           "olivedrab1")




core_microbiome_data_bars_RA1_noe <- merge(abundances.df, no_erwinia_data_bars, by = "Sample" )
# change the "Sphingomonadaceae Family" to "Sphingomonas"


#p_bar_plot_RA <-  
p_bar_plot_RA <- ggplot(data = core_microbiome_data_bars_RA1_noe, aes(x = Sample, y = Abundance.y, fill = Genus) ) +
  geom_col()+
  scale_fill_manual(values = rev(getPalette)) +
  #scale_fill_manual(values = rev(fills))+
  theme_minimal() +
  labs(y = "Relative Abundance\n", x = "\n" , fill = "Genus") +
  guides(fill = guide_legend(ncol = 1)) +
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(1, "lines"), 
    panel.grid.major = element_blank(),## remove the grid lines
    panel.grid.minor = element_blank()
  ) + 
  expand_limits(y = 1.4) +
  scale_x_discrete(limits = positions_cps) +
  ggtitle("Core microbiome\n Abundances\n") + 
  ylab("Relative abundances") +
  guides(fill=guide_legend(ncol=5) ) + 
  theme(legend.position = "bottom",  axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = -0.1)) +
  geom_text(aes(label = ifelse(Genus == "Rhodococcus", round(Abundance.x), " ")), y = 1.20) 
p_bar_plot_RA
png("images/Core microbiome without erwinia rarefied.png", height = 3000, width = 3500, res = 450)
print(p_bar_plot_RA)
dev.off()



# Erwinia_counts ----------------------------------------------------------
ps.erwiniaceae <- subset_taxa(ps.all, Family == "Erwiniaceae")
erwinia_counts <- read.delim("data/Erwinia/Erwinia_taxonomy_table_manual.tsv")
true_erwinia <- erwinia_counts[erwinia_counts$Genus == "Erwinia",]
true_erwinia
erwinia_otus_counts <- as.data.frame(taxa_sums(ps.erwiniaceae)) %>%
  tibble::rownames_to_column(var = "hash") %>%
  rename(Amount = 2) %>%
  filter(hash %in% true_erwinia$hash)
true_erwinia <- inner_join(erwinia_otus_counts, true_erwinia, by = "hash")

true_erwinia %>%
  summarise(
    sum_amylovora = sum(Amount[Species == "amylovora"]),
    sum_others = sum(Amount[Species != "amylovora"])
  )
#sum_amylovora sum_others
#1       1907353       3971

