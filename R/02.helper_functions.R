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


# Ordination --------------------------------------------------------------


plot_ordination_pears <- function(phy, OrdMeth, metadataColumn, titleAddition, shapeBy = NULL) {
  library(vegan)
  library(ggplot2)
  library(RColorBrewer)
  # this is a wrapper specifically for the pear project
  # It colors the phyloseq by collection date, arranges the grid etc.
  # OrdMeth is the distance metric
  # methods are: unifrac, wunifrac, dpcoa, jaccard, jsd, bray and anything on the phyloseq::distanceMethodList
  # titleAddition is a string
  if (missing(titleAddition)) {
    titleAddition <- " "
  }
  if (missing(OrdMeth)) {
    OrdMeth <- "bray"
    print("Default method chosen: bray-curtis")
  }
  if (metadataColumn == "Collection.date") {
    existing_positions <- unique(as.character(sample_data(phy)[[metadataColumn]]))
    positions1 <- positions_cd
    positions <- positions1[positions1 %in% existing_positions]
    colors1 <- color_values_cd[positions]
  } else if (metadataColumn == "Short.collection.date") {
    existing_positions <- unique(as.character(sample_data(phy)[[metadataColumn]]))
    positions1 <- positions_cd_s
    positions <- positions1[positions1 %in% existing_positions]
    colors1 <- color_values_cd_s[positions]
  } else if (metadataColumn == "CPS") {
    existing_positions <- unique(as.character(sample_data(phy)[[metadataColumn]]))
    positions1 <- positions_cps
    positions <- positions1[positions1 %in% existing_positions]
    colors1 <- color_values_cps[positions]
  } else if (metadataColumn == "Phenological.stage") {
    existing_positions <- unique(as.character(sample_data(phy)[[metadataColumn]]))
    positions1 <- positions_ps
    positions <- positions1[positions1 %in% existing_positions]
    colors1 <- color_values_ps[positions]
  } else {
    positions <- levels(sample_data(phy)[[metadataColumn]])
    numColors <- length(positions)
    colors1 <- colorRampPalette(colors = colors_of_set1)(numColors)
  }
  ## calculate the distance matrix
  # note to self: phyloseq uses the vegdist function of vegan to calculate the distances.
  # except for the unifracs that it calculates using its own function
  dist <- phyloseq::distance(phy, method = OrdMeth)
  # calculate the eignevectors and eigenvalues. k=2 saves 2 dimentions
  PCoA <- cmdscale(dist, eig = TRUE, k = 2)
  explainedvar1 <- round(PCoA$eig[1] / sum(PCoA$eig), 2) * 100
  explainedvar2 <- round(PCoA$eig[2] / sum(PCoA$eig), 2) * 100
  ordi_points <- as.data.frame(PCoA$points)
  ## add the metadatacolumn to the ordipoints dataframe
  metadata <- data.frame(sample_data(phy))
  ordi_points$sample <- rownames(ordi_points)
  if (!is.null(metadata$Sample.id)) {
    ordi_points$metadataCol <- metadata[[metadataColumn]][match(metadata$Sample.id, ordi_points$sample)]
  } else {
    ordi_points$metadataCol <- metadata[[metadataColumn]][match(metadata$Sample, ordi_points$sample)]
  }

  ordi_points$metadataCol <- factor(ordi_points$metadataCol, levels = positions)
  # add the shape to the ordipoints dataframe
  if (!is.null(shapeBy)) {
    if (!is.null(metadata$Sample.id)) {
      ordi_points$shapeBy <- metadata[[shapeBy]][match(metadata$Sample.id, ordi_points$sample)]
    } else {
      ordi_points$shapeBy <- metadata[[shapeBy]][match(metadata$Sample, ordi_points$sample)]
    }
  }

  p <- pear_ordiplot(ordi_points, explainedvar1, explainedvar2, metadataColumn, positions, colors1, shapeBy)
  print(p)
  return(p)
}

pear_ordiplot <- function(ordi_points, explainedvar1, explainedvar2, metadataColumn, positions, colors1, shapeBy) {
  ## ordi_points is a dataframe created by the
  ## ordi_pear function
  if (is.null(shapeBy)) {
    p <- ggplot(ordi_points, aes(x = V1, y = V2, color = metadataCol)) +
      geom_point(size = 2) +
      xlab(paste0(explainedvar1, "%")) +
      ylab(paste0(explainedvar2, "%")) +
      labs(color = metadataColumn, title = gsub(metadataColumn, ".", " ")) +
      scale_color_manual(values = colors1) +
      theme_minimal() +
      # stat_ellipse() +
      #    xlim(-1, 1) +
      #    ylim(-0.5, 0.5) +
      scale_fill_manual(breaks = positions)
  } else {
    p <- ggplot(ordi_points, aes(x = V1, y = V2, color = metadataCol, shape = shapeBy)) +
      geom_point(size = 2) +
      xlab(paste0(explainedvar1, "%")) +
      ylab(paste0(explainedvar2, "%")) +
      labs(color = metadataColumn, title = gsub(metadataColumn, ".", " ")) +
      scale_color_manual(values = colors1) +
      theme_minimal() +
      # stat_ellipse() +
      #    xlim(-1, 1) +
      #    ylim(-0.5, 0.5) +
      scale_fill_manual(breaks = positions)
  }
  return(p)
}


# bar plot functions ------------------------------------------------------


# Helper function to get taxonomic rank
get_taxrank <- function(ps) {
  tax_table <- as.data.frame(tax_table(ps))
  # The rank can be inferred from the last non-NA column in the taxonomy table
  last_non_na_column <- as.numeric(min(which(is.na(tax_table[1, ]))))
  # Get the rank name
  tax_rank <- colnames(tax_table)[last_non_na_column - 1]
  return(tax_rank)
}

make_data_bars1 <- function(your.ps, threshold, grouping_var) {
  ## this function prepares the dataframe for the taxonomic bar plot
  if (missing(grouping_var)) {
    dataBars1 <- your.ps %>%
      transform_sample_counts(function(x) {
        x / sum(x)
      }) %>% # Transform to rel. abundance
      psmelt()
  } else {
    dataBars1 <- your.ps %>%
      merge_samples(group = grouping_var, fun = mean) %>%
      transform_sample_counts(function(x) {
        x / sum(x)
      }) %>% # Transform to rel. abundance
      psmelt()
  }
  taxrank <- get_taxrank(your.ps)
  ps.melted <- your.ps %>% psmelt() ## here the abundance is count data not RA
  sum_of_seqs <- sum(ps.melted$Abundance)
  ps.melted.grouped <- aggregate(Abundance ~ get(taxrank), ps.melted, sum)
  ps.melted.grouped$ave <- ps.melted.grouped$Abundance / sum_of_seqs

  if (missing(threshold)) {
    return(dataBars1)
  } else {
    species.list <- ps.melted.grouped[ps.melted.grouped$ave > threshold, ]$`get(taxrank)`
    `%notin%` <- Negate(`%in%`)
    dataBars1[dataBars1[[taxrank]] %notin% species.list, ][[taxrank]] <- paste0(taxrank, " < ", threshold)
    return(dataBars1)
  }
}


taxa_bar_plot <- function(data_bars, x_group = "Sample") {
  # create the taxonomic bar plot with the data bars dataframe from the
  # make_data_bars or create_data_bars1 function
  taxrank <- data_bars %>%
    select(last_col(0)) %>%
    colnames()
  number_of_phyla <- length(unique(data_bars[[taxrank]]))
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))(number_of_phyla)
  ggplot(data_bars, aes(x = get(x_group), y = Abundance, fill = fct_reorder(get(taxrank), Abundance))) +
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
}


add_threshold_to_databars <- function(dataBars, threshold = 0.01) {
  taxrank <- dataBars %>%
    select(last_col(0)) %>%
    colnames()
  dataBars[[taxrank]][dataBars$Abundance < threshold] <- paste0("< ", threshold, " abund.")
  return(dataBars)
}

basic_bar_plot <- function(dataBars, taxrank = NULL) {
  library(RColorBrewer)
  library(ggplot2)
  if (is.null(taxrank)) {
    taxrank <- dataBars %>%
      select(last_col(0)) %>%
      colnames()
    number_of_phyla <- length(unique(dataBars[[taxrank]]))
    print(paste0("the taxonomic rank is ", taxrank, " and there are ", number_of_phyla, " of those"))
  }
  number_of_phyla <- length(unique(dataBars[[taxrank]]))
  if (taxrank == "Genus") {
    print("TRUE")
    all_taxa <- unique(dataBars[[taxrank]])
    all_taxa_except_Erwinia <- all_taxa[all_taxa != "Erwinia"]

    getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

    my_palette <- setNames(getPalette(length(all_taxa_except_Erwinia)), all_taxa_except_Erwinia)
    my_palette["Buchnera"] <- "green4"
    my_palette["Erwinia"] <- "red"

    dataBars[[taxrank]] <- factor(dataBars[[taxrank]], levels = c(rev(all_taxa_except_Erwinia), "Erwinia"))
    p2 <- ggplot(dataBars, aes(x = Sample, y = Abundance, fill = get(taxrank))) +
      geom_bar(stat = "identity") +
      scale_fill_manual(name = taxrank, values = my_palette) +
      theme_minimal() +
      theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
  } else {
    getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
    p2 <- ggplot(dataBars, aes(x = Sample, y = Abundance, fill = forcats::fct_reorder(get(taxrank), Abundance))) +
      geom_bar(stat = "identity") +
      scale_fill_manual(name = taxrank, values = rev(getPalette(number_of_phyla))) +
      theme_minimal() +
      theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
  }
  return(p2)
}

# boxplots for alpha diversities ------------------------------------------

plot_alpha_boxplot <- function(alpha_diversity_melt, myColors) {
  alpha_index_name <- as.character(unique(alpha_diversity_melt$measure))
  grouping_name <- as.character(colnames(alpha_diversity_melt[1]))
  p <- ggplot(alpha_diversity_melt, aes(x = get(grouping_name), y = value, fill = get(grouping_name))) +
    geom_boxplot() +
    labs(y = paste0(alpha_index_name, " index"), x = "") +
    theme_minimal() +
    scale_fill_manual(values = myColors) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), legend.position = "off")
  return(p)
}

## alpha diversity
## this function needs to be moved to the helper functions file

make_alpha_melted_df <- function(ps, grouping_var, index = "Shannon") {
  if (missing(grouping_var)) {
    print("you need to supply a grouping_var variable that is equal to one of the metadata columns")
    return(0)
  }
  # calculate the relevant alpha diversity
  alpha_d <- estimate_richness(ps, measures = index)
  alpha_d$Sample <- sub("X", "", rownames(alpha_d))
  metadata <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
  # keep only the index as some measures calculate more than one measure
  # alpha_d <- alpha_d[[ ,index]]
  # add the grouping_var column to the alpha_d data frame, with the name Group
  alpha_d$Group <- metadata[[grouping_var]][match(alpha_d$Sample, metadata$Sample)]
  alpha_diversity_melt <- reshape2::melt(alpha_d, id.vars = "Group", variable.name = "measure")
  # remove all the lines that contain the word "Sample"
  alpha_diversity_melt <- alpha_diversity_melt[alpha_diversity_melt$measure != "Sample", ]

  alpha_diversity_melt$value <- as.numeric(alpha_diversity_melt$value)
  # make sure that we only have the index
  alpha_diversity_melt_i <- alpha_diversity_melt[alpha_diversity_melt$measure == index, ]
  return(alpha_diversity_melt_i)
}

# Fasta manipulation functions --------------------------------------------


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
