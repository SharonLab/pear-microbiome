library(readxl)
library(ggplot2)
library(reshape2)
library(dplyr)
library(grid)
library(phyloseq)
setwd("~/qiime/fire-blight/pear.flower/20200125/qiime2.aia.filter/R_pears")

weather1 <- read_excel("Data/weather1.xlsx")
ps <- readRDS("Data/pear_phyloseq_filtered")


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
saveRDS(weather2, file = "Weather")

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

pdf(file = "images/weather.pdf", width = 15, height = 7.5)
ggplot(weather2, aes(x = as.Date(date, format = "%d.%m.%Y"), group = 1)) +
  geom_segment(data = dates, aes(x = x, y = 0, xend = x, yend = 31, color = "Collection dates"), size = 3, alpha = 0.6) +
  annotate("text", x = collection.dates$date, y = collection.dates$pos + 8, label = collection.dates$labels, size = 5) +
  annotate("text", x = phenological$dates, y = collection.dates$pos + 3, label = phenological$PS, size = 5) +
  annotate("text", x = cultivar$pos, y = 44, label = cultivar$x, size = 10) +
  geom_line(aes(y = as.integer(as.character(Tmax)), color = "High Temp"), size = 2) +
  geom_line(aes(y = as.integer(as.character(Tmin)), color = "Low Temp"), size = 2) +
  geom_line(aes(y = as.integer(as.character(WetHours)), color = "Wet Hours"), size = 2) +
  scale_color_manual(values = c("#984EA3", "#E41A1C", "#377EB8", "#4DAF4A")) +
  scale_fill_manual(breaks = c("Collection dates", "Max Temperature", "Wet Hours", "Min temperature")) +
  scale_y_continuous("Temperature (C)\n", sec.axis = sec_axis(~ . * 1, name = "Dew (Hours/Day)\n")) +
  theme_minimal() +
  labs(color = " ") +
  ylab("Temperature (C) \n") +
  xlab(" \nDate ") +
  scale_x_date(date_labels = "%d.%m.%y") +
  theme(axis.text.x.bottom = element_text(angle = 0), text = element_text(size = 16, ))
dev.off()

pdf(file = "images/weather_no_ps.pdf", width = 15, height = 7.5)
ggplot(weather2, aes(x = as.Date(date, format = "%d.%m.%Y"), group = 1)) +
  geom_segment(data = dates, aes(x = x, y = 0, xend = x, yend = 31, color = "Collection dates"), size = 3, alpha = 0.6) +
  annotate("text", x = collection.dates$date, y = collection.dates$pos + 3, label = collection.dates$labels, size = 5) +
  # annotate("text",x=phenological$dates, y = collection.dates$pos + 3 , label = phenological$PS, size = 5) +
  annotate("text", x = cultivar$pos, y = 40, label = cultivar$x, size = 8) +
  geom_line(aes(y = as.integer(as.character(Tmax)), color = "High Temp"), size = 2) +
  geom_line(aes(y = as.integer(as.character(Tmin)), color = "Low Temp"), size = 2) +
  geom_line(aes(y = as.integer(as.character(WetHours)), color = "Wet Hours"), size = 2) +
  scale_color_manual(values = c("#984EA3", "#E41A1C", "#377EB8", "#4DAF4A")) +
  scale_fill_manual(breaks = c("Collection dates", "Max Temperature", "Wet Hours", "Min temperature")) +
  scale_y_continuous("Temperature (C)\n", sec.axis = sec_axis(~ . * 1, name = "Dew (Hours/Day)\n")) +
  theme_minimal() +
  labs(color = " ") +
  ylab("Temperature (C) \n") +
  xlab(" \nDate ") +
  scale_x_date(date_labels = "%d.%m.%y") +
  theme(axis.text.x.bottom = element_text(angle = 0), text = element_text(size = 16, ))
dev.off()
