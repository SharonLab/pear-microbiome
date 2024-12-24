library(readxl)
library(ggplot2)
library(dplyr)
input_path <- "/Users/aiaoz/MyDocuments/science/fire blight/cfu"
Coscia_cfu <- read_xlsx(paste0(input_path, "/Fire_blight_cfu.xlsx"), sheet = "Coscia_final") %>%
  as.data.frame()
Coscia_cfu
Spadona_cfu <- read_xlsx(paste0(input_path, "/Fire_blight_cfu.xlsx"), sheet = "Spadona_final") %>%
  as.data.frame()

Coscia_cfu_long <- tidyr::pivot_longer(Coscia_cfu, cols = c(NA_plates, CCT_plates), names_to = "Count_type", values_to = "Count")
Coscia_cfu_long <- Coscia_cfu_long[Coscia_cfu_long$Count_type != "CCT_plates",]
# Specify the order of levels for Phenological_stage
# Coscia_cfu_long$Phenological_stage <- factor(Coscia_cfu_long$Phenological_stage, levels = c("White Buds", "Full Bloom", "Complete Petal Fall"), ordered = TRUE)


p <- ggplot(Coscia_cfu_long, aes(x = Phenological_stage, y = Count, fill = Count_type)) +
  geom_point(position = position_dodge(width = 0.3), size = 5, shape = 21, color = "black", stroke = 0.5) +
  facet_wrap(~Cultivar) + 
  theme_bw() + 
  scale_y_continuous(trans = "log10", 
                     breaks = scales::breaks_log(n=5),
                     labels = function(x) round(log10(x),1),
  limits = c(10, NA)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none") +
  scale_fill_manual(values = c("NA_plates" = "darkgoldenrod2", "CCT_plates" = "darkorchid2")) +
  scale_x_discrete(limits = c("White Buds", "Full Bloom", "Complete Petal Fall"), 
                   labels = c("White buds", "Full bloom", "Complete petal fall")) +
  labs(x = "\n", y = "\n\nlog CFU / flower", fill = "Media")

print(p)


png("images/Coscia_cfu_NA_no_CCT.png", res = 450, width = 3500, height = 3000 )
print(p)
dev.off()

## make a version with larger fonts for a poster
p2 <- p + 
  theme(
    axis.text = element_text(size = 14),          # Axis tick labels
    axis.title = element_text(size = 16),         # Axis titles
    strip.text = element_text(size = 16),         # Facet titles
    legend.text = element_text(size = 14),        # Legend text
    legend.title = element_text(size = 16),       # Legend title
    plot.title = element_text(size = 18, hjust = 0.5) # Plot title, centered
  )

print(p2)

png("images/Coscia_cfu for poster large fonts no cct.png", res = 450, width = 3500, height = 3000 )
print(p2)
dev.off()
