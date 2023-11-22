## this is the script for figure 4 
# and also the code for checking if Buchnera are correctly classified

source("R/02.helper_functions.R")
source("R/02.general_definitions.R")
library(readxl)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape2)

# eriwinia amylovora abundance calculations -------------------------------------------------
## the creation of the erwinia amylovora phyloseq is in
## 01.load_data_tps_object_rarefy.R

ps.e.2 <- readRDS("Data/phyloseq_erwinia_amylovora")
ps <- readRDS("Data/phyloseq_all_data")
ps.rarefied <- readRDS("Data/phyloseq_all_data_rarefied_noB")
ps.pseudo ## this phyloseq is created in "05.analysis of specific taxa"

#  heatmap functions -----------------------------------------------

generate_heatmap <- function(matrix_values, palette_name = "YlOrRd",  
                             use_fixed_range = FALSE, fixed_range = c(0.001, 0.943)) {
  # Determine breaks based on fixed range or data range
  if (use_fixed_range) {
    breaks <- seq(fixed_range[1], fixed_range[2], length.out = 101)
  } else {
    data_range <- range(matrix_values, na.rm = TRUE)
    breaks <- seq(data_range[1], data_range[2], length.out = 101)
  }
  
  # Prepare the data for display
  display_matrix <- matrix_values
  display_matrix[is.na(display_matrix)] <- min(breaks) # Replace NA with the lowest break for proper color mapping
  display_matrix <- round(display_matrix, 3)
 
  # Define color palette
  pal1 <- colorRampPalette(RColorBrewer::brewer.pal(9, palette_name))(length(breaks) - 1)
  # It receives a matrix of 6 by 8, where the column names are the phenological stages
  # and the row names are collection dates. it is very specific.
  if ("27.03" %in% colnames(matrix_values)) {
    matrix_values["Leaf", "27.03"] <- as.numeric(-0.1)
  }


  # Generate display matrix - this displays the numbers in writing inside the boxes.
  display_matrix <- matrix_values %>%
    round(3)
  display_matrix[display_matrix == -0.1] <- " "
  # multiplier_matrix <- unlist(c(display_matrix))
  # multiplier_matrix[multiplier_matrix == " "] <- NA
  # multiplier_matrix <- as.numeric(multiplier_matrix)
  # # Define color palette if unspecified
  # if (is.null(multiplier)){
  #   multiplier <- max(multiplier_matrix, na.rm = TRUE) / min(multiplier_matrix, na.rm = TRUE)
  #   if (multiplier < 10) {
  #     multiplier <- 15
  #   }
  # } 
  #   fixed_breaks <- seq(0.001, multiplier/100, length.out = length(my_pallete) + 1)
  # 
  # 
  # pal1 <- colorRampPalette(colors = RColorBrewer::brewer.pal(name = palette_name, n = 9)[1:8])(round(multiplier))

  my_pallete <- c("gray88", pal1)
 
  
  # Define annotations
  annotations_rows <- data.frame(Flower.or.Leaf = rep(c("Flower", "Leaf"), times = c(5, 1)))
  rownames(annotations_rows) <- rownames(matrix_values)
  # different annotations when we have spadona or not
  if ("27.03" %in% colnames(matrix_values)) { # check if we have spadona
    annotations <- data.frame(Cultivar = rep(c("Spadona", "Coscia"), times = c(3, 5)))
    rownames(annotations) <- colnames(matrix_values)
    # Generate heatmap
    p <- pheatmap::pheatmap(matrix_values,
      breaks = breaks, 
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = my_pallete,
      fontsize = 16,
      fontsize_number = 16,
      display_numbers = display_matrix,
      gaps_col = which(diff(as.numeric(as.factor(annotations$Cultivar))) != 0) + 0.5,
      gaps_row = which(diff(as.numeric(as.factor(annotations_rows$Flower.or.Leaf))) != 0) + 0.5
    )
  } else {
    p <- pheatmap::pheatmap(matrix_values,
      breaks = breaks, 
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = my_pallete,
      fontsize = 16,
      fontsize_number = 16,
      display_numbers = display_matrix,
      gaps_row = which(diff(as.numeric(as.factor(annotations_rows$Flower.or.Leaf))) != 0) + 0.5
    )
  }

  return(p)
}
prepare_data_for_heatmap <- function(data_bars_ra_e) {
  # this function should receive a product of the phyloseq::psmelt() function
  ## that the goes through the dcast function to create a matrix
  ## change the dataframe to fit the pheatmap specs
  # change zeros to -0.1 so that the nas can be gray
  data_bars_ra_e[data_bars_ra_e == 0] <- -0.1
  data_bars_ra_e[is.na(data_bars_ra_e)] <- -0.1
  # Reorder the collection dates, and turn them into the rownames of the matrix
  positions1 <- c("24.03", "27.03", "31.03", "11.04", "13.04", "15.04", "18.04", "21.04")
  existing_positions <- unique(as.character(data_bars_ra_e$Short.collection.date))
  positions <- positions1[positions1 %in% existing_positions]
  data_bars_ra_e$Short.collection.date <- as.character(data_bars_ra_e$Short.collection.date)
  data_bars_ra_e <- data_bars_ra_e[match(positions, data_bars_ra_e$Short.collection.date), ]
  rownames(data_bars_ra_e) <- data_bars_ra_e$Short.collection.date
  data_bars_ra_e$Short.collection.date <- NULL

  # Reorder the phenological stages

  order.ps <- c("White buds", "Initial bloom", "Full bloom", "Initial petal fall", "Complete petal fall", "Leaf")
  data_bars_ra_e <- rbind(data_bars_ra_e, colnames(data_bars_ra_e))
  data_bars_ra_e <- data_bars_ra_e[, order.ps]
  # remove the last row that is the phenological stages
  data_bars_ra_e1 <- data_bars_ra_e[-nrow(data_bars_ra_e), ]
  data_bars_ra_e1 <- dplyr::mutate_all(data_bars_ra_e1, as.numeric)

  return(data_bars_ra_e1)
}

#  a function helper function that adds a title to pheatmap.
add_title_to_pheatmap <- function(p, title) {
  # P is h pheatmap object, title is the title you want to add
  title_grob <- grid::textGrob(
    title,
    x = unit(0, "npc"),
    y = unit(1, "npc"),
    just = c("left", "top"),
    gp = gpar(fontsize = 14)
  )
  # t=1 and l = 4 are locations, they can be changed
  p$gtable <- gtable::gtable_add_grob(p$gtable, title_grob, t = 1, l = 4)
  return(p)
}


# main --------------------------------------------------------------------


# turn the erwinia phyloseq into a dataframe
erwinia_data_bars <- psmelt(ps.e.2) ## this is the phyloseq object that hold only the erwinia amylovora

## extract the relative abundances of all the data and turn into a data frame.
data_bars.ra <- ps.rarefied %>%
  speedyseq::merge_samples2(group = "cd.ps", fun_otu = mean) %>%
  transform_sample_counts(function(x) {
    x / sum(x)
  }) %>% # Transform to rel. abundance
  psmelt()
## create a new column where we only have the day and mobth without the year
data_bars.ra$Collection.date <- as.character(data_bars.ra$Collection.date)
data_bars.ra$Short.collection.date <- gsub("\\.\\d{4}$", "", data_bars.ra$Collection.date)
# extract the relative abundances of the erwinia otus and add them to the erwinia dataframe
data_bars.ra.e.1 <- data_bars.ra[data_bars.ra$OTU %in% erwinia_data_bars$OTU, ]
## sum up all the different otus
data_bars.ra.e <- reshape2::dcast(data_bars.ra.e.1, Short.collection.date ~ Phenological.stage, value.var = "Abundance", fun.aggregate = sum)

#### call the function
data_bars.ra.e1 <- prepare_data_for_heatmap(data_bars.ra.e)
# transpose the data so that the matrix is the same as other figures in the paper
transposed_data <- t(data_bars.ra.e1)
## make a display matrix for the values

erwiniaceae_heatmap <- generate_heatmap(transposed_data, "YlOrRd") 
print(erwiniaceae_heatmap)
png(
  file = "images/Erwiniaceae abundance heatmap square.png",
  width = 3866, height = 3000, res = 450
)
print(erwiniaceae_heatmap)
dev.off()
dev.set(which = 1)

rm(data_bars.ra)
rm(data_bars.ra.e)
rm(data_bars.ra.e.1)
rm(data_bars.ra.e1)

alpha_melted_df_all_samples <- estimate_richness(ps.rarefied)
alpha_melted_df_all_samples <- cbind(alpha_melted_df_all_samples, data.frame(sample_data(ps.rarefied)))
alpha_melted_df_all_samples$Collection.date <- as.character(alpha_melted_df_all_samples$Collection.date)
alpha_melted_df_all_samples$Short.collection.date <- gsub("\\.\\d{4}$", "", alpha_melted_df_all_samples$Collection.date)

## sum up all the different otus
alpha_d <- reshape2::dcast(alpha_melted_df_all_samples, Short.collection.date ~ Phenological.stage,
  value.var = "Shannon",
  fun.aggregate = mean
)

alpha_d_t <- prepare_data_for_heatmap(alpha_d)
# transpose the data so that the matrix is the same as other figures in the paper

transposed_data_a <- t(as.matrix(alpha_d_t))
## make a display matrix for the values

p1 <- generate_heatmap(transposed_data_a, "GnBu")
print(p1)

library(gridExtra)
library(grid)

# Create a null plot for the NULL space
empty <- grid::nullGrob()
png("images/heatmap of alpha square.png",  width = 3866, height = 3000, res = 450)
print(p1)
dev.off()
png(filename = "images/pheatmap of alpha and erwinia contents.png", width = 7500, height = 2500, res = 400)
gridExtra::grid.arrange(p$gtable, empty, p1$gtable,
  ncol = 3,
  widths = c(10, 1, 10)
)
dev.off()

#  figure 2 without spadona -----------------------------------------------
# turn the erwinia phyloseq into a dataframe
ps.s.coscia <- ps.e.2 %>%
  subset_samples(Cultivar == "Coscia")
ps.s.coscia
ps.rarefied <- readRDS(file = "Data/phyloseq_coscia_rarefied")
ps.rarefied

erwinia_data_bars <- psmelt(ps.s.coscia) ## this is the phyloseq object that hold only the erwinia amylovora

## extract the relative abundances of all the data and turn into a data frame.
data_bars.ra <- ps.rarefied %>%
  speedyseq::merge_samples2(group = "cd.ps", fun_otu = mean) %>%
  transform_sample_counts(function(x) {
    x / sum(x)
  }) %>% # Transform to rel. abundance
  psmelt()
## create a new column where we only have the day and mobth without the year
data_bars.ra$Collection.date <- as.character(data_bars.ra$Collection.date)
data_bars.ra$Short.collection.date <- gsub("\\.\\d{4}$", "", data_bars.ra$Collection.date)
# extract the relative abundances of the erwinia otus and add them to the erwinia dataframe
data_bars.ra.e.1 <- data_bars.ra[data_bars.ra$OTU %in% erwinia_data_bars$OTU, ]
## sum up all the different otus
data_bars.ra.e <- reshape2::dcast(data_bars.ra.e.1, Short.collection.date ~ Phenological.stage, value.var = "Abundance", fun.aggregate = sum)

#### call the function
data_bars.ra.e1 <- prepare_data_for_heatmap(data_bars.ra.e)
# transpose the data so that the matrix is the same as other figures in the paper
transposed_data <- t(data_bars.ra.e1)
## make a display matrix for the values

p <- generate_heatmap(transposed_data, "YlOrRd")
print(p)
png(
  file = "images/Erwiniaceae abundance heatmap no spadona.png",
  width = 4500, height = 3000, res = 450
)
print(p)
dev.off()
dev.set(which = 1)

rm(data_bars.ra)
rm(data_bars.ra.e)
rm(data_bars.ra.e.1)
rm(data_bars.ra.e1)

alpha_melted_df_all_samples <- estimate_richness(ps.rarefied)
alpha_melted_df_all_samples <- cbind(alpha_melted_df_all_samples, data.frame(sample_data(ps.rarefied)))
alpha_melted_df_all_samples$Collection.date <- as.character(alpha_melted_df_all_samples$Collection.date)
alpha_melted_df_all_samples$Short.collection.date <- gsub("\\.\\d{4}$", "", alpha_melted_df_all_samples$Collection.date)

## sum up all the different otus
alpha_d <- reshape2::dcast(alpha_melted_df_all_samples, Short.collection.date ~ Phenological.stage,
  value.var = "Shannon",
  fun.aggregate = mean
)

alpha_d_t <- prepare_data_for_heatmap(alpha_d)
# transpose the data so that the matrix is the same as other figures in the paper

transposed_data_a <- t(as.matrix(alpha_d_t))
## make a display matrix for the values

p1 <- generate_heatmap(transposed_data_a, "GnBu")
print(p1)

library(gridExtra)
library(grid)

# Create a null plot for the NULL space
empty <- grid::nullGrob()

png(filename = "images/pheatmap of alpha and erwinia contents without spadona.png", width = 7500, height = 2500, res = 400)
gridExtra::grid.arrange(p$gtable, empty, p1$gtable,
  ncol = 3,
  widths = c(10, 1, 10)
)
dev.off()


# figure 4 C+D without leaves------------------------------------------------------------
ps.coscia.flowers <- readRDS("Data/phyloseq_costia_flowers_rarefied")
# plot the ordination with collection date and phenological stage

p.ps <- plot_ordination_pears(ps.coscia.flowers, OrdMeth = "bray", metadataColumn = "Phenological.stage") +
   ylim(-0.4, 0.45) + stat_ellipse()  
p.cd <- plot_ordination_pears(ps.coscia.flowers, OrdMeth = "bray", metadataColumn = "Short.collection.date") +
  ylim(-0.4, 0.45) + stat_ellipse()  
png(file = "images/Fig 4 A Coscia  flowerspcoa.png", height = 2500, width = 3000, res = 450)
print(p.ps)
dev.off() 

png(file = "images/Fig 4 B Coscia  flowerspcoa.png", height = 3000, width = 3500, res = 450)
print(p.cd)
dev.off() 

png(file = "images/Coscia flowers pcoa.png", height = 2000, width = 4500, res = 450)
ggpubr::ggarrange(p.ps, NULL, p.cd,
                  ncol = 3, nrow = 1, 
                  widths = c(11, 1, 11.5))
dev.off() 



# figure 4 C+D  with leaves ------------------------------------------------------------
ps.coscia <- readRDS("Data/phyloseq_coscia_rarefied")
# plot the ordination with collection date and phenological stage

p.ps <- plot_ordination_pears(ps.coscia, OrdMeth = "bray", metadataColumn = "Phenological.stage", shapeBy = "Flower.or.leaf") +
  stat_ellipse() 
print(p.ps)
p.cd <- plot_ordination_pears(ps.coscia, OrdMeth = "bray", metadataColumn = "Short.collection.date", shapeBy = "Flower.or.leaf") +
  stat_ellipse() 
png(file = "images/Coscia flowers and leaves pcoa.png", height = 2000, width = 4500, res = 450)
ggpubr::ggarrange(p.ps, NULL, p.cd,
                  ncol = 3, nrow = 1, 
                  widths = c(11, 1, 11.5))
dev.off() 


# Pseudomonas heatmap -----------------------------------------------------

# turn the erwinia phyloseq into a dataframe
erwinia_data_bars <- psmelt(ps.pseudo) ## this is the phyloseq object that hold only the pseudomonas ASVs

## extract the relative abundances of all the data and turn into a data frame.
data_bars.ra <- ps.rarefied %>%
  speedyseq::merge_samples2(group = "cd.ps", fun_otu = mean) %>%
  transform_sample_counts(function(x) {
    x / sum(x)
  }) %>% # Transform to rel. abundance
  psmelt()
## create a new column where we only have the day and mobth without the year
data_bars.ra$Collection.date <- as.character(data_bars.ra$Collection.date)
data_bars.ra$Short.collection.date <- gsub("\\.\\d{4}$", "", data_bars.ra$Collection.date)
# extract the relative abundances of the erwinia otus and add them to the erwinia dataframe
data_bars.ra.e.1 <- data_bars.ra[data_bars.ra$OTU %in% erwinia_data_bars$OTU, ]
## sum up all the different otus
data_bars.ra.e <- reshape2::dcast(data_bars.ra.e.1, Short.collection.date ~ Phenological.stage, value.var = "Abundance", fun.aggregate = sum)

#### call the function
data_bars.ra.e1 <- prepare_data_for_heatmap(data_bars.ra.e)
# transpose the data so that the matrix is the same as other figures in the paper
transposed_data <- t(data_bars.ra.e1)
## make a display matrix for the values

pseudomonas_heatmap <- generate_heatmap(transposed_data, "YlOrRd", use_fixed_range = TRUE, fixed_range = c(0.01, 0.943)) 
print(pseudomonas_heatmap)
png(
  file = "images/Pseudomonas abundance heatmap new color scheme.png",
  width = 3866, height = 3000, res = 450
)
print(pseudomonas_heatmap)
dev.off()
dev.set(which = 1)

rm(data_bars.ra)
rm(data_bars.ra.e)
rm(data_bars.ra.e.1)
rm(data_bars.ra.e1)

# location on the tree through erwinia ------------------------------------
# test the effect of location on the tree on erwinia content
ps.e <- subset_taxa(ps.rarefied, Genus == "Erwinia")
taxa_sums(ps.e)
sample_sums(ps.e)
# take only coscia flowers because only they have a meaning in this context
ps.e.c.f <- subset_samples(ps.e, Cultivar == "Coscia") %>%
  subset_samples(Flower.or.leaf == "Flower") 
ps.e.c.f

meta.e.c.f <- data.frame(sample_data(ps.e.c.f))
meta.e.c.f$Erwinia.content <- sample_sums(ps.e.c.f)  
meta.e.c.f <- meta.e.c.f %>%
  filter(!(Phenological.stage == "White buds" & Short.collection.date == "15.04"))

positions_cd_s.c <- c( "11.04", "13.04", "15.04", "18.04", "21.04")
positions_ps <- c("White buds", "Initial bloom", "Full bloom", "Initial petal fall", "Complete petal fall", "Leaf")
plot_list_EA <- list()

create_plot <- function(data, CD, positions_ps) {
  temp.meta <- data[data$Short.collection.date == CD,]
  
  # Check which phenological stages exist for this collection date
  temp_ps <- positions_ps[positions_ps %in% unique(temp.meta$Phenological.stage)]
  
  p <- ggplot(data = temp.meta, 
              aes(x = factor(Phenological.stage, levels = temp_ps), 
                  y = Erwinia.content, 
                  color = Location.on.tree)) + 
    geom_boxplot() +
    theme_minimal() +
    ggpubr::stat_compare_means(aes(label = ..p.signif..)) +
    xlab("Phenological.stage") + 
    annotate("text",
             x = 1:length(table(temp.meta$Phenological.stage)),
             y = min(temp.meta$Erwinia.content) + max(temp.meta$Erwinia.content) / 100,
             label = table(temp.meta$Phenological.stage),
             col = "blue") +
    ggtitle(paste(CD, "\n")) 
  
  return(p)
}

plot_list_EA <- lapply(positions_cd_s.c, function(CD) {
  create_plot(data = meta.e.c.f, CD = CD, positions_ps = positions_ps)
})
png("images/initial erwinia abundances comparisons location on tree.png", )
do.call(gridExtra::grid.arrange, c(plot_list_EA, ncol=2))
dev.off()
# now do the same for the group
plot_list_EA <- lapply(positions_cd_s.c, function(CD) {
  create_plot(data = meta.e.c.f, CD = CD, positions_ps = positions_ps)
})

do.call(gridExtra::grid.arrange, c(plot_list_EA_CD, ncol=2))  
# aes(label = ..p.signif..)
