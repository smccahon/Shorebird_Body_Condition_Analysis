#---------------------------------------------#
#  Lesser Yellowlegs Pectoral Shape Analysis  #
#            Created 12/11/2024               #          
#           Modified 12/11/2024               #
#---------------------------------------------#

# MODIFED TO USE FAT AS ORDINAL (LOGSTIC REGRESSION, 2 GROUPS: Low and High)

# load packages
library(ggplot2)
library(AICcmodavg)
library(tidyverse)

# ---------------------------------------------------------------------------- #

# Data Processing & Manipulation ####

# Read data
setwd("processed_data")
birds <- read.csv("Shorebird_Data_Cleaned_2024-12-9.csv")

# Make neonicotinoid detection column (Detection/Non-detection)
birds$Detection <- ifelse(birds$OverallNeonic > 0, "Detection", "Non-detection")

birds <- birds %>% 
  mutate(LogAg = log10(PercentAg))

# Make sampling event column
birds <- birds %>% 
  mutate(Event = paste(Season, Year, sep = " ")) %>% 
  mutate(Event = ifelse(Season %in% c("Spring", "Fall"),
                        paste(Season, Year, sep = " "),
                        NA))

# Convert capture time to number of minutes after sunrise
birds <- birds %>%
  mutate(
    hour_of_day = as.numeric(format(strptime(Time, "%H:%M"), "%H")),  # Extract hour (0-23)
    minute_of_day = as.numeric(format(strptime(Time, "%H:%M"), "%M")),  # Extract minute (0-59)
    time_in_minutes = hour_of_day * 60 + minute_of_day  # Total time in minutes
  ) %>% 
  mutate(
    hour_of_day_s = as.numeric(format(strptime(Sunrise, "%H:%M"), "%H")),  # Extract hour (0-23)
    minute_of_day_s = as.numeric(format(strptime(Sunrise, "%H:%M"), "%M")),  # Extract minute (0-59)
    time_in_minutes_s = hour_of_day_s * 60 + minute_of_day_s  # Total time in minutes
  ) %>% 
  mutate(
    ts.sunrise = time_in_minutes - time_in_minutes_s)

# Reorder factor variables
birds$Sex <- factor(birds$Sex,
                    levels = c("M", "F"),
                    labels = c("Male", "Female"))

birds$Detection <- as.factor(birds$Detection)

birds$AgCategory <- factor(birds$AgCategory,
                           levels = c("Low", "Moderate", "High"))

# Subset data for Lesser Yellowlegs
leye <- subset(birds, Species %in% c("LesserYellowlegs"))

leye$Event <- factor(leye$Event, 
                     levels = c("Fall 2021", "Spring 2022", "Fall 2023"),
                     labels = c("Fall 2021", "Spring 2022", "Fall 2023"))

leye$Fat <- factor(leye$Fat,
                   levels = c("0", "1", "2", "3", "4", "5"))

leye$PecScore <- factor(leye$PecScore,
                   levels = c("1", "2", "3"))


# Convert Pectoral Score to two levels (low and high)

# Let's see how correlated pectoral score and pectoral muscle size are...
ggplot(na.omit(leye), aes(x = PecScore, y = PecSizeBest)) +
  geom_boxplot() +  
  theme_classic() +
  labs(x = "Lesser Yellowlegs Pectoral Shape (1-3)", 
       y = expression("Lesser Yellowlegs Pectoral Muscle Size" ~~~ (mm[score]))) +  
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")



# Standardize continuous variables
leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

# ---------------------------------------------------------------------------- #