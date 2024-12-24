

PlotBarDefault <- function(data_bars, x_group = "Sample", taxrank = NULL, colorPalette = NULL, abundance_type = "RA") {
  library(RColorBrewer)
  library(ggplot2)
  library(forcats)
  # Create the Abundance column
  data_bars$Abundance <- data_bars[[abundance_type]]
  
  if(is.null(taxrank)){ 
    taxrank <- data_bars %>%
      select(last_col(0)) %>%
      colnames()
  }
  number_of_phyla <- length(unique(data_bars[[taxrank]]))
  
  if (is.null(colorPalette)) {
    my_palette = colorRampPalette(brewer.pal(9, "Set1"))(number_of_phyla)
  } else {
     my_palette <- colorPalette
   }
  
  # Create label for y-axis based on abundance type
  y_label <- if (abundance_type == "RA") "Relative Abundance\n" else "Number of reads\n"
  
  p <- ggplot(data_bars, aes(x = get(x_group), y = Abundance, fill = fct_reorder(get(taxrank), Abundance))) + 
    geom_col() + 
    scale_fill_manual(values = rev(my_palette)) +
    theme_minimal() + 
    theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) +
    labs(y = y_label, x = "\n", fill = taxrank) + 
    guides(fill = guide_legend()) +
    theme(legend.title = element_text(size = 8), 
          legend.text =  element_text(size = 8), 
          legend.key.size = unit(1, "lines")) 
  return(p)
}


PlotBarsGenus <- function(data_bars, x_group = "Sample") {
  library(RColorBrewer)
  library(ggplot2)
  taxrank <- "Genus"
  
  all_taxa <- unique(data_bars[[taxrank]])
  all_taxa_except_Erwinia <- all_taxa[all_taxa != "Erwinia" & all_taxa != "< 0.01 abund."]
  # add more colors
  combined_colors <- c(  brewer.pal(8, "Set1"), brewer.pal(12, "Set3")) 

  getPalette = colorRampPalette( combined_colors)
  
  my_palette <- setNames(getPalette(length(all_taxa_except_Erwinia)), all_taxa_except_Erwinia)
  my_palette["Pseudomonas"] <- "green4"
  my_palette["Erwinia"] <- "red"
  my_palette["< 0.01 abund."] <- "darkgray"
  
  data_bars[[taxrank]] <- factor(data_bars[[taxrank]], levels = c("< 0.01 abund.", sort(all_taxa_except_Erwinia), "Erwinia"))
  
  p2 <- ggplot(data_bars, aes(x = get(x_group), y = Abundance, fill = get(taxrank))) + 
    geom_bar(stat = "identity") +
    scale_fill_manual(name = taxrank, values = my_palette) +
    theme_minimal() + 
    theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
  
  return(p2)
}


