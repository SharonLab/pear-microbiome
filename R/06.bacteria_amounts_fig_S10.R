library(phyloseq)
library(ggplot2)
library(dplyr)
library(qiime2R)
source("R/02.general_definitions.R")
source("R/02.helper_functions.R")
table_with_mitochondria <- read_qza("../02.table-pseudo.qza")$data %>%
  as.data.frame()
metadata <- read.delim("Data/input_data/metadata-filtered-copy.tsv")
dna_concentrations <- read.csv("Data/224 samples stefan.csv")

ps_with_mitochondria <- qiime2R::qza_to_phyloseq(features = "../02.table-pseudo.qza", 
                                                 taxonomy = "../04.taxonomy-silva.qza", 
                                                 tree = "../03.unrooted-tree.qza",
                                                 metadata = "Data/input_data/metadata-filtered-copy.tsv"
                                                  )
ps_with_mitochondria

ps_mito <- subset_taxa(ps_with_mitochondria, Genus == "Mitochondria")
ps_chloro <- subset_taxa(ps_with_mitochondria, Order == "Chloroplast")
ps_mito_chloro <- subset_taxa(ps_with_mitochondria, Order == "Chloroplast" | Genus == "Mitochondria")
## test
ps_mito_chloro
tax_table(ps_mito_chloro)

## load the entire data set without the mitochondria and chloroplasts
ps.all <- readRDS(file = "Data/phyloseq_all_data")
ps.all
ps.erwinia <- subset_taxa(ps.all, Genus == "Erwinia")
ps.no.erwinia <- subset_taxa(ps.all, Genus != "Erwinia")



colnames(dna_concentrations)
rownames(dna_concentrations)
dim(metadata) # 229 * 27
# remove the empty lines and the description line from the metadata and the PCR:blanks
metadata <- metadata[metadata$Cultivar != "Control", ]
metadata <- metadata[metadata$sampleid != "#q2:types", ]
metadata <- metadata[metadata$sampleid != "",  ]
dim(metadata)

dim(dna_concentrations) # 288
# remove the rows with the NAs

dna_concentrations <- dna_concentrations[! is.na(dna_concentrations$Sample_Name),]
dna_concentrations$Sample_Name <- as.character(dna_concentrations$Sample_Name)
dim(dna_concentrations) # 223

 class(metadata$sampleid)
 class(dna_concentrations$Sample_Name) 
 dna_concentrations$nucleic_acid_concentration <- gsub("low", "0.001", dna_concentrations$nucleic_acid_concentration)
 metadata$concentration <- dna_concentrations$nucleic_acid_concentration[match(metadata$sampleid, dna_concentrations$Sample_Name)] %>%
   as.numeric()

## calculate the percentage of mitochondria and chloroplast reads 
 metadata$mito_chloro_reads <- sample_sums(ps_mito_chloro)[match(metadata$sampleid, names(sample_sums(ps_mito_chloro)))]

metadata$reads_after_filtration <- sample_sums(ps.all)[match(metadata$sampleid ,names(sample_sums(ps.all)))]
metadata$number_of_reads_with_MC <- sample_sums(ps_with_mitochondria)[match(metadata$sampleid ,names(sample_sums(ps_with_mitochondria)))]
metadata$percentage_of_MC <- ( metadata$mito_chloro_reads / metadata$number_of_reads_with_MC ) * 100
#calculate the mitochondria reads
metadata$mito_reads <- sample_sums(ps_mito)[match(metadata$sampleid, names(sample_sums(ps_mito)))]
metadata$percentage_of_mito <- ( metadata$mito_reads / metadata$number_of_reads_with_MC ) * 100

## calculate the chloroplast reads
metadata$chloro_reads <- sample_sums(ps_chloro)[match(metadata$sampleid, names(sample_sums(ps_chloro)))]
metadata$percentage_of_chloro <- ( metadata$chloro_reads / metadata$number_of_reads_with_MC ) * 100
metadata$number_of_chloro_mito <- metadata$chloro_reads + metadata$mito_reads

metadata$number_of_bacterial_read <- metadata$number_of_reads_with_MC - metadata$number_of_chloro_mito

# calculate the numbers of non-erwinia reads
metadata$non_erwinia_reads <- sample_sums(ps.no.erwinia)[match(metadata$sampleid, names(sample_sums(ps.no.erwinia)))]
metadata$non_erwinia_reads[is.na(metadata$non_erwinia_reads)] <-0
metadata$percentage_of_non_erwinia_MC <- ( metadata$non_erwinia_reads / metadata$number_of_reads_with_MC ) * 100
metadata$erwinia_reads <- sample_sums(ps.erwinia)[match(metadata$sampleid, names(sample_sums(ps.erwinia)))]

metadata$percentage_of_erwinia <- (metadata$erwinia_reads / metadata$reads_after_filtration ) * 100
metadata$percentage_of_erwinia_MC <- (metadata$erwinia_reads / metadata$number_of_reads_with_MC ) * 100
metadata$mito_to_chloro_ratio <- metadata$mito_reads /metadata$chloro_reads
metadata$percentage_of_non_erwinia_reads_with_mc <-  sample_sums(ps.no.erwinia.with.buchnera)[match(metadata$sampleid, names(sample_sums(ps.no.erwinia.with.buchnera)))]
## following are two calculations of the same thing, the first one is less intuitive
# metadata$bacteria_to_plant_cells <-  (((1/(metadata$percentage_of_MC/100) )-1)*60)/6
# this is a simpler calculation
metadata$bacteria_to_plant_cells <- (metadata$number_of_bacterial_read/6) / (metadata$number_of_chloro_mito/60)
aggregate(metadata$percentage_of_MC, list(metadata$CPS), FUN = median)
aggregate(metadata$bacteria_to_plant_cells, list(metadata$CPS), FUN = median)

size_bacterial_G <- 3000000
size_pear_genome <-  57700000
size_mito_genome <- 459000
size_of_chloroplast_genome <- 160000
total_pear_DNA <- size_pear_genome + 10*size_of_chloroplast_genome + 50 * size_mito_genome
## this is for checking the equation
#total_pear_DNA <- 501000000
metadata$fraction_bacteria_to_plant <- metadata$bacteria_to_plant_cells * size_bacterial_G/(total_pear_DNA + metadata$bacteria_to_plant_cells*size_bacterial_G)
aggregate(metadata$fraction_bacteria_to_plant, list(metadata$CPS), FUN = median)
aggregate(metadata$concentration, list(metadata$CPS), FUN = median)
metadata$DNA_weight <- metadata$concentration * metadata$fraction_bacteria_to_plant
aggregate(metadata$DNA_weight, list(metadata$CPS), FUN = median)

metadata$num_bacterial_genomes <- round((metadata$DNA_weight / 30) /3.2358e-6)
aggregate(metadata$num_bacterial_genomes, list(metadata$CPS), FUN = median)
metadata$num_bacterial_genomes
metadata$genome_ratio <- round((metadata$number_of_bacterial_read/(6*30)) / metadata$num_bacterial_genomes)
aggregate(metadata$genome_ratio, list(metadata$CPS), FUN = median)
metadata$num_erwinia_genomes <- round(metadata$num_bacterial_genomes *(metadata$percentage_of_erwinia/100))
aggregate(metadata$num_erwinia_genomes, list(metadata$CPS), FUN = median)




df <- aggregate(metadata$number_of_bacterial_read, list(metadata$CPS), FUN = mean) %>%
  as.data.frame()
df$estimated_bacterial_weight <- c(12.4, 1.95, 0.14, 3.48, 0, 0.035, 0.028, 0, 0.018)
df$sample_weight <- aggregate(metadata$concentration, list(metadata$CPS), FUN = mean)

df$num_bacterial_cells <- df$x /6

(sample_weight*6.0221*10^23)



metadata %>%
  group_by(non_erwinia_reads) %>%
  summarize_at(vars(CPS), list(name = mean))

# amount of DNA boxplot  (called p4)---------------------------------------------------
metadata$concentration
metadata1 <- metadata[metadata$concentration != 217.000, ]
metadata1$concentration

## plot the amount of DNA as boxplots
p4 <- 
  ggplot(metadata1, aes(x = CPS, y = concentration, color = CPS)) + 
  geom_boxplot() +
  stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = sprintf("%.2f", after_stat(y))), # Format the median value to 2 decimal places
    vjust = -0.5, # Adjust vertical position of the label
    color = "black" # Set the color of the text
  ) +
  scale_x_discrete(limits = positions_cps) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = color_values_cps) +
  labs( title = "The amount of extracted DNA by\nphenological stages (ng)", tag = "A")+
  theme(legend.position = 'none',
        plot.title = element_text(size = 12)) +
  ylab("Amount of DNA (ng)")+
  xlab("\n")
metadata1_no_leaves <- metadata1[metadata1$Flower.or.leaf == "Flower", ]

p4_no_leaves  <- 
  ggplot(metadata1_no_leaves, aes(x = CPS, y = concentration, color = CPS)) + 
  geom_boxplot() +
  stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = sprintf("%.2f", after_stat(y))), # Format the median value to 2 decimal places
    vjust = -0.5, # Adjust vertical position of the label
    color = "black" # Set the color of the text
  ) +
  scale_x_discrete(limits = positions_cps[1:7]) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = color_values_cps) +
  labs( title = "The amount of extracted DNA by phenological stages (ng)", tag = "B")+
  theme(legend.position = 'none') +
  ylab("Amount of DNA (ng)")+
  xlab("\n")

print(p4_no_leaves)

png("images/DNA_concentrations_phenological_stages fig 5b.png", res = 350, height = 2000, width = 2500)
print(p4_no_leaves)
dev.off()   


# mitochondria and chloroplast amounts boxplots ---------------------------


## plot the number of mitochondria and chloroplast read

p3 <- ggplot(metadata, aes(x = CPS, y = percentage_of_MC, color = CPS)) + 
  geom_boxplot() +
  stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = sprintf("%.1f", after_stat(y))), # Format the median value to 2 decimal places
    vjust = -0.5, # Adjust vertical position of the label
    color = "black" # Set the color of the text
  ) +
  scale_x_discrete(limits = positions_cps) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(size = 12) ) + 
  scale_color_manual(values = color_values_cps) +
  labs(title = "The percentage of mitochondria\nand chloroplast reads", tag = "D")+
  theme(legend.position = 'none') +
  ylab("Number of reads") +
  xlab("\n")
print(p3)
p2 <- ggplot(metadata, aes(x = Collection.date, y = percentage_of_MC, color = Collection.date)) + 
  geom_boxplot() +
  stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = sprintf("%.1f", after_stat(y))), # Format the median value to 2 decimal places
    vjust = -0.5, # Adjust vertical position of the label
    color = "black" # Set the color of the text
  ) +
  scale_x_discrete(limits = positions_cd) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(size = 12))+
  scale_color_manual(values = color_values_cd) +
  labs(title = "The percentage of mitochondria\nand chloroplast reads", tag = "C")+
  theme(legend.position = 'none') +
  ylab("% of reads") +
  xlab("\n")

print(p2)


p1 <- ggplot(metadata, aes(x = CPS, y = number_of_reads_with_MC, color = CPS)) + 
  geom_boxplot() +
  scale_x_discrete(limits = positions_cps) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(size = 12))+
  scale_color_manual(values = color_values_cps) +
  labs(title = "Total number of reads\nwith mitochondria and chloroplast", tag = "B")+
  theme(legend.position = 'none') +
  ylab("Number of reads") +
  xlab("\n")
print(p1)


ggplot(metadata, aes(x = CPS, y = reads_after_filtration, color = CPS)) + 
  geom_boxplot() +
  stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = sprintf("%.1f", after_stat(y))), # Format the median value to 2 decimal places
    vjust = -0.5, # Adjust vertical position of the label
    color = "black" # Set the color of the text
  ) +
  scale_x_discrete(limits = positions_cps) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(size = 12))+
  scale_color_manual(values = color_values_cps) + 
  ggtitle("Number of reads after filtration by phenological stages") +
  theme(legend.position = 'none')

# use pathworks::wrap_plots 
# make a list of the plots
plot_list <- list(p4, p1, p2, p3)

patchwork::wrap_plots(plot_list)

# gridExtra::grid.arrangplot_listgridExtra::grid.arrange(p4, p1, p2, p3)
# png("images/Overview of read numbers and amount of DNA concentrations by phenological stage and cultivar.png", 
#     res = 350, height = 3500, width = 5000)
# gridExtra::grid.arrange(p4, p1, p2, p3)
# dev.off()


p5 <- ggplot(metadata, aes(x = CPS, y = mito_to_chloro_ratio, color = CPS)) + 
  geom_boxplot() +
  stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = sprintf("%.1f", after_stat(y))), # Format the median value to 2 decimal places
    vjust = -0.5, # Adjust vertical position of the label
    color = "black" # Set the color of the text
  ) +
  scale_x_discrete(limits = positions_cps) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(size = 12))+
  scale_color_manual(values = color_values_cps) +
  labs(title = "Mitochondria to chloroplast ratio", tag = "E")+
  theme(legend.position = 'none') +
  ylab("\n") +
  xlab("\n")
print(p5)


p6 <- ggplot(metadata, aes(x = CPS, y = percentage_of_chloro, color = CPS)) + 
  geom_boxplot() +
  stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = sprintf("%.1f", after_stat(y))), # Format the median value to 2 decimal places
    vjust = -0.5, # Adjust vertical position of the label
    color = "black" # Set the color of the text
  ) +
  scale_x_discrete(limits = positions_cps) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(size = 12))+
  scale_color_manual(values = color_values_cps) +
  labs("The percentage of chloroplast reads by phenological stages", tag = "NaN")+
  theme(legend.position = 'none')
print(p6)

write.csv(metadata, "Data/metadata_with_chloroplast_reads.csv")



## plot the amount of DNA as boxplots
p7 <- 
  ggplot(metadata, aes(x = CPS, y = percentage_of_erwinia, color = CPS)) + 
  geom_boxplot() +
  stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = sprintf("%.2f", after_stat(y))), # Format the median value to 2 decimal places
    vjust = -0.5, # Adjust vertical position of the label
    color = "black" # Set the color of the text
  ) +
  scale_x_discrete(limits = positions_cps) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(size = 12))+
  scale_color_manual(values = color_values_cps) +
  labs(title = "Percentage of Erwinia reads\nfrom total bacterial reads", tag = "F")+
  theme(legend.position = 'none') +
  ylab("% of reads") +
  xlab("\n")
print(p7)

p8 <- 
  ggplot(metadata, aes(x = CPS, y = percentage_of_non_erwinia_reads_with_mc, color = CPS)) + 
  geom_boxplot() +
  stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = sprintf("%.2f", after_stat(y))), # Format the median value to 2 decimal places
    vjust = -0.5, # Adjust vertical position of the label
    color = "black" # Set the color of the text
  ) +
  scale_x_discrete(limits = positions_cps) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(size = 12))+
  scale_color_manual(values = color_values_cps) +
  labs(title = "number of non-erwinia reads")+
  theme(legend.position = 'none') + 
  ylab("")

print(p8)
plot_list <- list(p4, p1, p2, p3, p5, p7)

combined_plot <- patchwork::wrap_plots(plot_list, ncol = 2)

png("images/supplementary fig 9 read_number_data.png", height = 5000, width = 4000, res = 450)
print(combined_plot)
dev.off()

# glm ---------------------------------------------------------------------


plot(metadata$concentration, metadata$reads_after_filtration)
plot(metadata$concentration, metadata$input)
plot(metadata$concentration, metadata$percentage_of_MC)
plot(metadata$concentration, metadata$erwinia_reads)
glm_model<-glm(percentage_of_MC ~ concentration  , data = metadata)  
# concentration: For each unit increase in concentration, the percentage of MC is expected to decrease by 0.4871.
glm_model<-glm(percentage_of_MC ~ concentration + CPS , data = metadata)  
glm_model   
# For each unit increase in concentration, the percentage of MC is expected to decrease by 0.1041, holding CPS constant.
 cor(metadata$concentration , metadata$percentage_of_MC , method = "pearson")
library(sjPlot)    
#plot_model(glm_model, vline.color = "red")    
# png("images/glm of expected percentage of MC by dna concentrations and cps.png")
# plot_model(glm_model, show.values = TRUE, value.offset = .3)   
# dev.off()
 metadata.flowers <- metadata[!grepl("Leaf", metadata$Flower.or.leaf), ]
 metadata.flowers$CPS <- as.factor(metadata.flowers$CPS)
 metadata.flowers$CPS <- relevel(metadata.flowers$CPS, ref = "Spadona White buds")
 metadata.flowers$Collection.date <- as.factor(metadata.flowers$Collection.date)
 metadata.flowers$CPS <- relevel(metadata.flowers$Collection.date, ref = "27.03.2019")
 
 plot(metadata.flowers$concentration, metadata.flowers$erwinia_reads)
 glm_model<-glm(concentration ~  erwinia_reads + CPS , data = metadata.flowers)  
 glm_model   
 sjPlot::plot_model(glm_model, vline.color = "red") 
 glm_model<-glm(concentration ~  erwinia_reads + CPS + Collection.date , data = metadata.flowers)  
 glm_model
 # 
 # Call:  glm(formula = concentration ~ erwinia_reads + CPS + Collection.date, 
 #            data = metadata.flowers)
 # 
 # Coefficients:
 #   (Intercept)              erwinia_reads              CPS11.04.2019  
 # 1.243704                   0.001117                  -3.176989  
 # CPS13.04.2019              CPS15.04.2019              CPS18.04.2019  
 # -5.627466                   0.015777                   3.987617  
 # CPS21.04.2019              CPS24.03.2019              CPS31.03.2019  
 # 19.592397                   1.972157                   0.949268  
 # Collection.date13.04.2019  Collection.date15.04.2019  Collection.date18.04.2019  
 # NA                         NA                         NA  
 # Collection.date21.04.2019  Collection.date24.03.2019  Collection.date27.03.2019  
 # NA                         NA                         NA  
 # Collection.date31.03.2019  
 # NA  
 # 
 # Degrees of Freedom: 174 Total (i.e. Null);  166 Residual
 # (1 observation deleted due to missingness)
 # Null Deviance:	    132500 
 # Residual Deviance: 73710 	AIC: 1574
 # 
 sjPlot::plot_model(glm_model, vline.color = "red") 
 # For each unit increase in concentration, the percentage of MC is expected to decrease by 0.1041, holding CPS constant.
 cor(metadata.flowers$concentration , metadata.flowers$erwinia_reads , method = "spearman")
 ## another plot option (nicer in my opinion)
 library(effects) 
 
 plot(allEffects(glm_model))
 
 library(car)
 Anova(glm_model, type = "II")
 ## pvalue for erwinia reads is 2.669e-07 ***



# new model using the cfu estimates ---------------------------------------
# we want to use both the dna concentrations and the cfu estimates to estimate the amount of bacteria and the amount of plant DNA,
# to see if there is a difference between different conditions in the amount of bacteria, based on the amount of DNA extracted. 
 colnames(metadata)
 
 metadata <- metadata %>%
   mutate(bacteria_count_estimation = 
            case_when(Phenological.stage == "Full bloom" ~ 10000, 
                      Phenological.stage == "White buds" ~100,
                      Phenological.stage == "Initial petal fall" ~1000000))
 

# check differences in dna amounts between high samples and low on --------
# the amounts of DNA are specified in the "concentration" column
 metadata_cosc_flowers <- metadata %>%
   filter(Cultivar == "Coscia") %>%
   filter(Flower.or.leaf == "Flower")
 kruskal.test(concentration~Location.on.tree, data = metadata_cosc_flowers)
 
 # Kruskal-Wallis chi-squared = 0.85919, df = 1, p-value = 0.354

 metadata_cosc_flowers_last <- metadata %>%
   filter(Cultivar == "Coscia") %>%
   filter(Flower.or.leaf == "Flower") %>%
   filter(Phenological.stage != "White buds") %>%
   filter(Phenological.stage != "Initial bloom")
 wilcox.test(concentration~Location.on.tree, data = metadata_cosc_flowers_last)
 
# Kruskal-Wallis chi-squared = 0.23403, df = 1, p-value = 0.6286
 
kruskal.test(concentration~ Phenological.stage, data = metadata_cosc_flowers) 
 
#$Kruskal-Wallis chi-squared = 100.6, df = 4, p-value < 2.2e-16
