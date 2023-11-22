## data_bars is a dataframe that is required for creating the taxonomic bar plot


MakeDataBarsDefault <- function(ps, grouping) {
  ## ps. is a phyloseq object
  if (missing(grouping)) {
    # Absolute counts
    dataBarsCounts <- ps %>%
      psmelt() # Melt to long format

    # Relative abundance
    dataBarsRA <- ps %>%
      transform_sample_counts(function(x) {
        x / sum(x)
      }) %>% # Transform to rel. abundance
      psmelt() # Melt to long format
  } else {
    # Absolute counts
    dataBarsCounts <- ps %>%
      merge_samples(group = grouping, fun = mean) %>%
      psmelt() # Melt to long format

    # Relative abundance
    dataBarsRA <- ps %>%
      merge_samples(group = grouping, fun = mean) %>%
      transform_sample_counts(function(x) {
        x / sum(x)
      }) %>% # Transform to rel. abundance
      psmelt() # Melt to long format
  }

  # Rename columns
  names(dataBarsCounts)[names(dataBarsCounts) == "Abundance"] <- "Counts"
  names(dataBarsRA)[names(dataBarsRA) == "Abundance"] <- "RA"

  # Create a vector with names of columns that are not "Counts" or "RA"
  common_cols <- setdiff(names(dataBarsCounts), "Counts")

  # Join two dataframes together
  dataBars <- merge(dataBarsCounts, dataBarsRA, by = common_cols, all = TRUE)

  return(dataBars)
}


MakeDataBarsOfTaxon <- function(phy, threshold = 0.01, min_reads = 1, counts = FALSE, grouping) {
  ps.taxa <- phy %>%
    speedyseq::merge_samples2(group = "CPS", fun_otu = mean) %>%
    filter_taxa(function(x) sum(x) > min_reads, TRUE)
  data_bars <- MakeDataBarsDefault(ps.taxa, grouping)

  data_bars$Abundance <- data_bars$RA
  data_bars$OTU[data_bars$Abundance < threshold] <- paste0("< ", threshold, " abund.")

  return(data_bars)
}
