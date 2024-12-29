# create Figure S1 - weather data. 
library(readxl)
library(ggplot2)
library(reshape2)
library(dplyr)
library(grid)
library(phyloseq)
library(RColorBrewer)


weather1 <- read_excel("Data/weather1.xlsx")
ps <- readRDS("Data/pear_phyloseq_filtered")

## do not run lines 15-35 (reading the weather data and adding it to the metadata, exact line numbers may change)
# if you are using the metadata from supplementary table 1,
# or if you are using the metadata from the github as it already contains the weather data. 
metadata <- sample_data(ps)
metadata$MaxTemp <- weather1$Tmax[match(metadata$Collection.date, weather1$date)]
metadata$MinTemp <- weather1$Tmin[match(metadata$Collection.date, weather1$date)]

metadata$Dew <- weather1$WetHours[match(metadata$Collection.date, weather1$date)]
metadata$FireBlightWarning <- weather1$FireBlightWarning[match(metadata$Collection.date, weather1$date)]

# sample_data(ps) <- metadata
View(metadata)

# saveRDS(ps, "pear_phyloseq_filtered")

weather2 <- weather1[-c(1:3), ]
## or start rom here:
weatherMetadata <- subset(as.data.frame(sample_data(ps)), select = c(MaxTemp, MinTemp, Collection.date, Dew, CPS))

weather2$date <- as.Date(weather2$date, format = "%d.%m.%Y")
weather2$Tmax <- as.numeric(as.character(weather2$Tmax))
weather2$Tmin <- as.integer(as.character(weather2$Tmin))
weather2 <- weather2 %>% mutate(FB = ifelse(FireBlightWarning != "c", NA, Tmax + 2))
write.table(weather2, "Data/Weather")

FireBlightWarningDates <- as.Date(weather1$date[weather1$FireBlightWarning == "c"], format = "%d.%m.%Y", group = 1)
FireBlightWarningDates1 <- as.data.frame(as.Date(weather1$date[weather1$FireBlightWarning == "c"], format = "%d.%m.%Y", group = 1))
names(FireBlightWarningDates1) <- "dates"
FireBlightWarningDates1$temperature <- weather2$Tmax[match(FireBlightWarningDates1$dates, weather2$date)] + 2

# a data frame with the annotation data
collection.dates <- as.data.frame(unique(as.Date(weatherMetadata$Collection.date, format = "%d.%m.%Y")))
names(collection.dates) <- "date"
collection.dates$pos <- 32
collection.dates$labels <- gsub("0", "", format(collection.dates$date, format = "%d.%m"))

metadata$MaxTemp <- weather1$Tmax[match(metadata$Collection.date, weather1$date)]

cultivar <- as.data.frame(c("Spadona", "Coscia"))
colnames(cultivar) <- "x"
cultivar$pos <- as.Date(c("27.3.2019", "16.4.2019"), format = "%d.%m.%Y")

phenological <- data.frame(
  dates = sort(collection.dates$date),
  PS = c(
    "WB", "WB", "IB", "WB\nIB\nFB",
    "WB\nIB\nFB", "AS", "PPF\nCPF", "FB\nPPF\nCPF"
  )
)

# collection.dates$Phenological.stages
# this is for the geom_segment()
dates <- as.data.frame(collection.dates$date)
colnames(dates) <- "x"
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))(9)

# Filter weather2 for dates where FireBlightWarnings equals "c"
fireblight_dates <- weather2 %>%
  filter(FireBlightWarning == "c")


p <- ggplot(weather2, aes(x = as.Date(date, format = "%d.%m.%Y"), group = 1)) +
  geom_segment(data = dates, aes(x = x, y = 0, xend = x, yend = 31, color = "Collection dates"), size = 3, alpha = 0.6) +
  annotate("text", x = collection.dates$date, y = collection.dates$pos + 3, label = collection.dates$labels, size = 6) +
  annotate("text", x = cultivar$pos, y = 40, label = cultivar$x, size = 8) +
  geom_line(aes(y = as.integer(as.character(Tmax)), color = "High Temp"), size = 2) +
  geom_line(aes(y = as.integer(as.character(Tmin)), color = "Low Temp"), size = 2) +
  geom_line(aes(y = as.integer(as.character(WetHours)), color = "Wet Hours"), size = 2) +
  # Add dots for FireBlightWarnings
  geom_point(data = fireblight_dates, aes(x = as.Date(date, format = "%d.%m.%Y"), 
                                          y = as.integer(as.character(Tmax))), 
             size = 4, color = "black", shape = 21, fill = "red") +
  scale_color_manual(values = c("#984EA3", "#E41A1C", "#377EB8", "#4DAF4A")) +
  scale_fill_manual(breaks = c("Collection dates", "Max Temperature", "Wet Hours", "Min temperature")) +
  scale_y_continuous("Temperature (C)\n", sec.axis = sec_axis(~ . * 1, name = "Dew (Hours/Day)\n")) +
  theme_minimal() +
  labs(color = " ") +
  ylab("Temperature (C) \n") +
  xlab(" \nDate ") +
  scale_x_date(date_labels = "%d.%m.%y") +
  theme(axis.text.x.bottom = element_text(angle = 0),
          axis.text.x = element_text(size = 18),   # Increase x-axis tick label size
          axis.text.y = element_text(size = 18),   # Increase y-axis tick label size
          axis.title.x = element_text(size = 20),  # Increase x-axis label size
          axis.title.y = element_text(size = 20),  # Increase y-axis label size
          text = element_text(size = 16))

##
# Filter weather2 for dates where FireBlightWarnings equals "c"
fireblight_dates <- weather2 %>%
  filter(FireBlightWarning == "c")

# Plot
p <- ggplot(weather2, aes(x = as.Date(date, format = "%d.%m.%Y"), group = 1)) +
  # Collection dates as vertical lines
  geom_segment(data = dates, aes(x = x, y = 0, xend = x, yend = 31, color = "Collection dates"), size = 3, alpha = 0.6) +
  # FireBlightWarnings as black vertical lines
  geom_segment(data = fireblight_dates, aes(x = as.Date(date, format = "%d.%m.%Y"), 
                                            y = 0, xend = as.Date(date, format = "%d.%m.%Y"), yend = 31), 
               color = "black", size = 0.8, alpha = 0.7) +
  # Annotate text
  annotate("text", x = collection.dates$date, y = collection.dates$pos + 3, label = collection.dates$labels, size = 6) +
  annotate("text", x = cultivar$pos, y = 40, label = cultivar$x, size = 8) +
  # Temperature and wet hours lines
  geom_line(aes(y = as.integer(as.character(Tmax)), color = "High Temp"), size = 2) +
  geom_line(aes(y = as.integer(as.character(Tmin)), color = "Low Temp"), size = 2) +
  geom_line(aes(y = as.integer(as.character(WetHours)), color = "Wet Hours"), size = 2) +
  # Manual color scales
  scale_color_manual(values = c("#984EA3", "#E41A1C", "#377EB8", "#4DAF4A")) +
  scale_fill_manual(breaks = c("Collection dates", "Max Temperature", "Wet Hours", "Min temperature")) +
  # Y-axis settings
  scale_y_continuous("Temperature (C)\n", sec.axis = sec_axis(~ . * 1, name = "Dew (Hours/Day)\n")) +
  # Axis labels and theme adjustments
  theme_minimal() +
  labs(color = " ") +
  ylab("Temperature (C) \n") +
  xlab(" \nDate ") +
  scale_x_date(date_labels = "%d.%m.%y") +
  theme(
    axis.text.x.bottom = element_text(angle = 0),
    axis.text.x = element_text(size = 18),   # Increase x-axis tick label size
    axis.text.y = element_text(size = 18),   # Increase y-axis tick label size
    axis.title.x = element_text(size = 20),  # Increase x-axis label size
    axis.title.y = element_text(size = 20),  # Increase y-axis label size
    text = element_text(size = 16)
  )

# Print the plot
print(p)


 png("images/supplementary1 weather data.png", height = 3000, width = 5900, res = 400)
print(p)
dev.off()
