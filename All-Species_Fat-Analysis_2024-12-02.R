#----------------------------#
#  All Species Fat Analysis  #
#     Created 12/02/2024     #          
#    Modified 01/02/2025     #
#----------------------------#

# load packages
library(dplyr)
library(ggplot2)
library(AICcmodavg)
library(RColorBrewer)
library(trtools)
library(tidyverse)
library(lme4)
library(car)
library(MASS)
library(ordinal)
library(pROC)
library(pscl)

# ---------------------------------------------------------------------------- #

# Data Processing & Manipulation ####

# Read data
setwd("processed_data")
birds <- read.csv("Shorebird_Data_Cleaned_2024-12-9.csv")

# Make neonicotinoid detection column (Detection/Non-detection)
birds$Detection <- ifelse(birds$OverallNeonic > 0, "Detection", "Non-detection")

# Reorder factor variables
birds$Sex <- factor(birds$Sex,
                    levels = c("M", "F"),
                    labels = c("Male", "Female"))

birds$Detection <- as.factor(birds$Detection)

birds$Species <- factor(birds$Species, 
                        levels = c("LesserYellowlegs", "Killdeer", "Willet", 
                                   "PectoralSandpiper", "WIPH", "LeastSandpiper", 
                                   "LongbilledDowitcher", "AmericanAvocet",
                                   "SemipalmatedSandpiper", "MarbledGodwit",
                                   "ShortbilledDowitcher", "GreaterYellowlegs"),
                        labels = c("Lesser Yellowlegs", "Killdeer", "Willet", 
                                   "Pectoral Sandpiper", "Wilson's Phalarope", 
                                   "Least Sandpiper","Long-billed Dowitcher", "American Avocet",
                                   "Semipalmated Sandpiper", "Marbled Godwit",
                                   "Short-billed Dowitcher", "Greater Yellowlegs"))

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

# only include birds with fat
birds <- birds %>% 
  filter(!is.na(Fat))

birds$Fat <- factor(birds$Fat,
                    levels = c("0", "1", "2", "3", "4", "5", "6"))

birds <- birds %>%
  mutate(Fat.G = case_when(
    Fat == 0 | Fat == 1 ~ "Low",
    Fat == 2 | Fat == 3 | Fat == 4 | Fat == 5 | Fat == 6 ~ "High"       
  ))

birds$Fat.G <- factor(birds$Fat.G,
                      levels = c("Low", "High"))

birds$Fat.G_binary <- ifelse(birds$Fat.G == "Low", 0, 1)

# Display data to determine any collinearity issues
birds %>%
  dplyr::select(Fat.G, Event) %>%
  dplyr::filter(Event == "Spring 2022") %>%
  print()

# Subset data to omit spring (multicollinearity between fat and spring)
birds <- birds[birds$Event %in% c("Fall 2021", "Fall 2023"), ]

# Only include birds with at least 3 observations per species
birds <- birds %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

birds <- birds %>% 
  mutate(MigStatus = case_when(
    Species %in% c("Lesser Yellowlegs",
                   "Semipalmated Sandpiper", "Pectoral Sandpiper",
                   "Least Sandpiper") ~ "Migratory",
    Species %in% c("Willet", "Killdeer") ~ "Resident",
    TRUE ~ "Unknown"
  ))

# Create New Ag Category column
quartiles <- birds %>%
  summarise(
    Q1 = quantile(PercentAg, 0.25, na.rm = TRUE),
    Median = quantile(PercentAg, 0.5, na.rm = TRUE),
    Q3 = quantile(PercentAg, 0.75, na.rm = TRUE)
  )

birds <- birds %>% 
  mutate(AgIntensity = case_when(
    PercentAg <= quartiles$Q1 ~ "Low",          
    PercentAg > quartiles$Q1 & PercentAg <= quartiles$Q3 ~ "Moderate",  
    PercentAg > quartiles$Q3 ~ "High" 
  ))

# Convert variables to factors
birds$AgIntensity <- factor(birds$AgIntensity,
                           levels = c("Low", "Moderate", "High"))

birds$Event <- factor(birds$Event,
                    levels = c("Fall 2021", "Fall 2023"))

birds$MigStatus <- factor(birds$MigStatus,
                      levels = c("Resident", "Migratory"))


# Standardize continuous variables
birds.cs <- birds %>%
  mutate(across(where(is.numeric), scale))

# ---------------------------------------------------------------------------- #

# Data Visualization ####

# Creating my theme
my_theme <- theme_classic() +
  theme(text = element_text(size = 16),
        axis.title.y = element_text(margin = margin(r = 13)),
        axis.title.x = element_text(margin = margin(t = 13)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")

# Fat
ggplot(birds, aes(x = Fat.G, y = PercentAg)) + 
  geom_boxplot() +
  labs(x = "Fat Score Category", y = "Surrounding Agricultural Intensity (%)") +
  my_theme 

ggplot(birds, aes(x = Fat, y = PercentAg)) + 
  geom_boxplot() +
  labs(x = "Fat Score Category", y = "Surrounding Agricultural Intensity (%)") +
  my_theme

# Date
ggplot(birds, aes(y = DaysIntoSeason_S, x = Fat.G)) +
  geom_boxplot()

ggplot(birds, aes(x = Julian, y = Fat)) +
  geom_boxplot()

# Detection
ggplot(birds, aes(x = Detection, y = Fat.G)) +
  geom_boxplot()

# ---------------------------------------------------------------------------- #

# Note: I had issues with model convergence ("boundary (singfular) fit"). 
# Species as a random effect had a variance very close to 0. 
# I decided to treat species as a fixed effect rather than a random effect.

# All Species Modeling ####

## Univariate Analysis: Date ####

m1 <- glm(Fat.G ~ DaysIntoSeason_S, data = birds.cs, family = "binomial")
m2 <- glm(Fat.G ~ Julian, data = birds.cs, family = "binomial")

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Days into season performs better

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Agricultural Intensity ####

m1 <- glm(Fat.G ~ PercentAg, data = birds.cs, family = "binomial")
m2 <- glm(Fat.G ~ PercentAg + I(PercentAg^2), data = birds.cs, family = "binomial")
m3 <- glm(Fat.G ~ AgIntensity, data = birds.cs, family = "binomial")

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# PercentAg as a continuous variable performs the best

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Event and Capture Time ####

m1 <- glm(Fat.G ~ Event * ts.sunrise, data = birds.cs, family = "binomial")
m2 <- glm(Fat.G ~ Event + ts.sunrise + I(ts.sunrise^2), data = birds.cs, family = "binomial")
m3 <- glm(Fat.G ~ Event + ts.sunrise, data = birds.cs, family = "binomial")

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with interaction performs the best

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Date * Migratory Status ####

m1 <- glm(Fat.G ~ DaysIntoSeason_S + MigStatus, data = birds.cs, family = "binomial")
m2 <- glm(Fat.G ~ DaysIntoSeason_S * MigStatus, data = birds.cs, family = "binomial")

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with interaction performs the best

# ---------------------------------------------------------------------------- #

# INTERACTIONS/TRANSFORMATIONS INCLUDED: Event * Time, Date * Migratory Status, PercentAg

# ---------------------------------------------------------------------------- #

m.global <- glm(Fat.G ~ Sex + Event * ts.sunrise + DaysIntoSeason_S + 
                   MigStatus + PercentAg + Detection +
                   MigStatus * DaysIntoSeason_S +
                   Event * ts.sunrise + Species, family = "binomial", data = birds.cs)

m.null <- glm(Fat.G ~ 1, family = "binomial", data = birds.cs)

## Single Covariate Models ####

m1 <- glm(Fat.G ~ Sex, family = "binomial", data = birds.cs)
m2 <- glm(Fat.G ~ Event, family = "binomial", data = birds.cs)
m3 <- glm(Fat.G ~ ts.sunrise, family = "binomial", data = birds.cs)
m4 <- glm(Fat.G ~ DaysIntoSeason_S, family = "binomial", data = birds.cs)
m5 <- glm(Fat.G ~ MigStatus, family = "binomial", data = birds.cs)
m6 <- glm(Fat.G ~ PercentAg, family = "binomial", data = birds.cs)
m7 <- glm(Fat.G ~ Detection, family = "binomial", data = birds.cs)
m8 <- glm(Fat.G ~ Species, family = "binomial", data = birds.cs)

##  Additive Models ####

### Two additive combinations ####

m9 <- glm(Fat.G ~ Sex + Event, family = "binomial", data = birds.cs)
m10 <- glm(Fat.G ~ Sex + ts.sunrise, family = "binomial", data = birds.cs)
m11 <- glm(Fat.G ~ Sex + DaysIntoSeason_S, family = "binomial", data = birds.cs)
m12 <- glm(Fat.G ~ Sex + MigStatus, family = "binomial", data = birds.cs)
m13 <- glm(Fat.G ~ Sex + PercentAg, family = "binomial", data = birds.cs)
m14 <- glm(Fat.G ~ Sex + Detection, family = "binomial", data = birds.cs)
m15 <- glm(Fat.G ~ Event * ts.sunrise + Event*ts.sunrise, family = "binomial", data = birds.cs)
m16 <- glm(Fat.G ~ Event + DaysIntoSeason_S, family = "binomial", data = birds.cs)
m17 <- glm(Fat.G ~ Event + MigStatus, family = "binomial", data = birds.cs)
m18 <- glm(Fat.G ~ Event + PercentAg, family = "binomial", data = birds.cs)
m19 <- glm(Fat.G ~ Event + Detection, family = "binomial", data = birds.cs)
m20 <- glm(Fat.G ~ ts.sunrise + DaysIntoSeason_S, family = "binomial", data = birds.cs)
m21 <- glm(Fat.G ~ ts.sunrise + MigStatus, family = "binomial", data = birds.cs)
m22 <- glm(Fat.G ~ ts.sunrise + PercentAg, family = "binomial", data = birds.cs)
m23 <- glm(Fat.G ~ ts.sunrise + Detection, family = "binomial", data = birds.cs)
m24 <- glm(Fat.G ~ DaysIntoSeason_S * MigStatus + DaysIntoSeason_S * MigStatus, family = "binomial", data = birds.cs)
m25 <- glm(Fat.G ~ DaysIntoSeason_S + PercentAg, family = "binomial", data = birds.cs)
m26 <- glm(Fat.G ~ DaysIntoSeason_S + Detection, family = "binomial", data = birds.cs)
m27 <- glm(Fat.G ~ MigStatus + PercentAg, family = "binomial", data = birds.cs)
m28 <- glm(Fat.G ~ MigStatus + Detection, family = "binomial", data = birds.cs)
m29 <- glm(Fat.G ~ PercentAg + Detection, family = "binomial", data = birds.cs)
m30 <- glm(Fat.G ~ Sex + Species, family = "binomial", data = birds.cs)
m31 <- glm(Fat.G ~ Event + Species, family = "binomial", data = birds.cs)
m32 <- glm(Fat.G ~ ts.sunrise + Species, family = "binomial", data = birds.cs)
m33 <- glm(Fat.G ~ DaysIntoSeason_S + Species, family = "binomial", data = birds.cs)
m34 <- glm(Fat.G ~ MigStatus + Species, family = "binomial", data = birds.cs)
m35 <- glm(Fat.G ~ PercentAg + Species, family = "binomial", data = birds.cs)

### Three Additive Combinations ####

m36 <- glm(Fat.G ~ Sex + Event * ts.sunrise + Event * ts.sunrise, family = "binomial", data = birds.cs)
m37 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S, family = "binomial", data = birds.cs)
m38 <- glm(Fat.G ~ Sex + Event + MigStatus, family = "binomial", data = birds.cs)
m39 <- glm(Fat.G ~ Sex + Event + PercentAg, family = "binomial", data = birds.cs)
m40 <- glm(Fat.G ~ Sex + Event + Detection, family = "binomial", data = birds.cs)
m41 <- glm(Fat.G ~ Sex + ts.sunrise + DaysIntoSeason_S, family = "binomial", data = birds.cs)
m42 <- glm(Fat.G ~ Sex + ts.sunrise + MigStatus, family = "binomial", data = birds.cs)
m43 <- glm(Fat.G ~ Sex + ts.sunrise + PercentAg, family = "binomial", data = birds.cs)
m44 <- glm(Fat.G ~ Sex + ts.sunrise + Detection, family = "binomial", data = birds.cs)
m45 <- glm(Fat.G ~ Sex + DaysIntoSeason_S * MigStatus + DaysIntoSeason_S * MigStatus, family = "binomial", data = birds.cs)
m46 <- glm(Fat.G ~ Sex + DaysIntoSeason_S + PercentAg, family = "binomial", data = birds.cs)
m47 <- glm(Fat.G ~ Sex + DaysIntoSeason_S + Detection, family = "binomial", data = birds.cs)
m48 <- glm(Fat.G ~ Sex + MigStatus + PercentAg, family = "binomial", data = birds.cs)
m49 <- glm(Fat.G ~ Sex + MigStatus + Detection, family = "binomial", data = birds.cs)
m50 <- glm(Fat.G ~ Sex + PercentAg + Detection, family = "binomial", data = birds.cs)
m51 <- glm(Fat.G ~ Event * ts.sunrise + DaysIntoSeason_S, family = "binomial", data = birds.cs)
m52 <- glm(Fat.G ~ Event * ts.sunrise + MigStatus, family = "binomial", data = birds.cs)
m53 <- glm(Fat.G ~ Event * ts.sunrise + PercentAg, family = "binomial", data = birds.cs)
m54 <- glm(Fat.G ~ Event * ts.sunrise + Detection, family = "binomial", data = birds.cs)
m55 <- glm(Fat.G ~ Event + DaysIntoSeason_S * MigStatus, family = "binomial", data = birds.cs)
m56 <- glm(Fat.G ~ Event + DaysIntoSeason_S + PercentAg, family = "binomial", data = birds.cs)
m57 <- glm(Fat.G ~ Event + DaysIntoSeason_S + Detection, family = "binomial", data = birds.cs)
m58 <- glm(Fat.G ~ Event + MigStatus + PercentAg, family = "binomial", data = birds.cs)
m59 <- glm(Fat.G ~ Event + MigStatus + Detection, family = "binomial", data = birds.cs)
m60 <- glm(Fat.G ~ Event + PercentAg + Detection, family = "binomial", data = birds.cs)
m61 <- glm(Fat.G ~ ts.sunrise + DaysIntoSeason_S * MigStatus, family = "binomial", data = birds.cs)
m62 <- glm(Fat.G ~ ts.sunrise + DaysIntoSeason_S + PercentAg, family = "binomial", data = birds.cs)
m63 <- glm(Fat.G ~ ts.sunrise + DaysIntoSeason_S + Detection, family = "binomial", data = birds.cs)
m64 <- glm(Fat.G ~ ts.sunrise + MigStatus + PercentAg, family = "binomial", data = birds.cs)
m65 <- glm(Fat.G ~ ts.sunrise + MigStatus + Detection, family = "binomial", data = birds.cs)
m66 <- glm(Fat.G ~ ts.sunrise + PercentAg + Detection, family = "binomial", data = birds.cs)
m67 <- glm(Fat.G ~ DaysIntoSeason_S * MigStatus + PercentAg, family = "binomial", data = birds.cs)
m68 <- glm(Fat.G ~ DaysIntoSeason_S * MigStatus + Detection, family = "binomial", data = birds.cs)
m69 <- glm(Fat.G ~ DaysIntoSeason_S + PercentAg + Detection, family = "binomial", data = birds.cs)
m70 <- glm(Fat.G ~ MigStatus + PercentAg + Detection, family = "binomial", data = birds.cs)
m71 <- glm(Fat.G ~ Sex + Event + Species, family = "binomial", data = birds.cs)
m72 <- glm(Fat.G ~ Sex + ts.sunrise + Species, family = "binomial", data = birds.cs)
m73 <- glm(Fat.G ~ Sex + DaysIntoSeason_S + Species, family = "binomial", data = birds.cs)
m74 <- glm(Fat.G ~ Sex + MigStatus + Species, family = "binomial", data = birds.cs)
m75 <- glm(Fat.G ~ Sex + PercentAg + Species, family = "binomial", data = birds.cs)
m76 <- glm(Fat.G ~ Event * ts.sunrise + Species, family = "binomial", data = birds.cs)
m77 <- glm(Fat.G ~ Event + DaysIntoSeason_S + Species, family = "binomial", data = birds.cs)
m78 <- glm(Fat.G ~ Event + MigStatus + Species, family = "binomial", data = birds.cs)
m79 <- glm(Fat.G ~ Event + PercentAg + Species, family = "binomial", data = birds.cs)
m80 <- glm(Fat.G ~ ts.sunrise + DaysIntoSeason_S + Species, family = "binomial", data = birds.cs)
m81 <- glm(Fat.G ~ ts.sunrise + MigStatus + Species, family = "binomial", data = birds.cs)
m82 <- glm(Fat.G ~ ts.sunrise + PercentAg + Species, family = "binomial", data = birds.cs)
m83 <- glm(Fat.G ~ DaysIntoSeason_S * MigStatus + Species, family = "binomial", data = birds.cs)
m84 <- glm(Fat.G ~ DaysIntoSeason_S + PercentAg + Species, family = "binomial", data = birds.cs)
m85 <- glm(Fat.G ~ MigStatus + PercentAg + Species, family = "binomial", data = birds.cs)

### Four Additive Combinations ####

m86 <- glm(Fat.G ~ Sex + Event * ts.sunrise + DaysIntoSeason_S, family = "binomial", data = birds.cs)
m87 <- glm(Fat.G ~ Sex + Event * ts.sunrise + MigStatus, family = "binomial", data = birds.cs)
m88 <- glm(Fat.G ~ Sex + Event * ts.sunrise + PercentAg, family = "binomial", data = birds.cs)
m89 <- glm(Fat.G ~ Sex + Event * ts.sunrise + Detection, family = "binomial", data = birds.cs)
m90 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S * MigStatus, family = "binomial", data = birds.cs)
m91 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S + PercentAg, family = "binomial", data = birds.cs)
m92 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S + Detection, family = "binomial", data = birds.cs)
m93 <- glm(Fat.G ~ Sex + Event + MigStatus + PercentAg, family = "binomial", data = birds.cs)
m94 <- glm(Fat.G ~ Sex + Event + MigStatus + Detection + Event + MigStatus, family = "binomial", data = birds.cs)
m95 <- glm(Fat.G ~ Sex + Event + PercentAg + Detection, family = "binomial", data = birds.cs)
m96 <- glm(Fat.G ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus, family = "binomial", data = birds.cs)
m97 <- glm(Fat.G ~ Sex + ts.sunrise + DaysIntoSeason_S + PercentAg, family = "binomial", data = birds.cs)
m98 <- glm(Fat.G ~ Sex + ts.sunrise + DaysIntoSeason_S + Detection, family = "binomial", data = birds.cs)
m99 <- glm(Fat.G ~ Sex + ts.sunrise + MigStatus + PercentAg, family = "binomial", data = birds.cs)
m100 <- glm(Fat.G ~ Sex + ts.sunrise + MigStatus + Detection, family = "binomial", data = birds.cs)
m101 <- glm(Fat.G ~ Sex + ts.sunrise + PercentAg + Detection, family = "binomial", data = birds.cs)
m102 <- glm(Fat.G ~ Sex + DaysIntoSeason_S * MigStatus + PercentAg, family = "binomial", data = birds.cs)
m103 <- glm(Fat.G ~ Sex + DaysIntoSeason_S * MigStatus + Detection, family = "binomial", data = birds.cs)
m104 <- glm(Fat.G ~ Sex + DaysIntoSeason_S + PercentAg + Detection, family = "binomial", data = birds.cs)
m105 <- glm(Fat.G ~ Sex + MigStatus + PercentAg + Detection, family = "binomial", data = birds.cs)
m106 <- glm(Fat.G ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus, family = "binomial", data = birds.cs)
m107 <- glm(Fat.G ~ Event * ts.sunrise + DaysIntoSeason_S + PercentAg, family = "binomial", data = birds.cs)
m108 <- glm(Fat.G ~ Event * ts.sunrise + DaysIntoSeason_S + Detection, family = "binomial", data = birds.cs)
m109 <- glm(Fat.G ~ Event * ts.sunrise + MigStatus + PercentAg, family = "binomial", data = birds.cs)
m110 <- glm(Fat.G ~ Event * ts.sunrise + MigStatus + Detection, family = "binomial", data = birds.cs)
m111 <- glm(Fat.G ~ Event * ts.sunrise + PercentAg + Detection, family = "binomial", data = birds.cs)
m112 <- glm(Fat.G ~ Event + DaysIntoSeason_S * MigStatus + PercentAg, family = "binomial", data = birds.cs)
m113 <- glm(Fat.G ~ Event + DaysIntoSeason_S * MigStatus + Detection, family = "binomial", data = birds.cs)
m114 <- glm(Fat.G ~ Event + DaysIntoSeason_S + PercentAg + Detection, family = "binomial", data = birds.cs)
m115 <- glm(Fat.G ~ Event + MigStatus + PercentAg + Detection, family = "binomial", data = birds.cs)
m116 <- glm(Fat.G ~ ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg, family = "binomial", data = birds.cs)
m117 <- glm(Fat.G ~ ts.sunrise + DaysIntoSeason_S * MigStatus + Detection, family = "binomial", data = birds.cs)
m118 <- glm(Fat.G ~ ts.sunrise + DaysIntoSeason_S + PercentAg + Detection, family = "binomial", data = birds.cs)
m119 <- glm(Fat.G ~ ts.sunrise + MigStatus + PercentAg + Detection, family = "binomial", data = birds.cs)
m120 <- glm(Fat.G ~ DaysIntoSeason_S * MigStatus + PercentAg + Detection, family = "binomial", data = birds.cs)
m121 <- glm(Fat.G ~ Sex + Event * ts.sunrise + Species, family = "binomial", data = birds.cs)
m122 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S + Species, family = "binomial", data = birds.cs)
m123 <- glm(Fat.G ~ Sex + Event + MigStatus + Species + Event + MigStatus, family = "binomial", data = birds.cs)
m124 <- glm(Fat.G ~ Sex + Event + PercentAg + Species, family = "binomial", data = birds.cs)
m125 <- glm(Fat.G ~ Sex + ts.sunrise + DaysIntoSeason_S + Species, family = "binomial", data = birds.cs)
m126 <- glm(Fat.G ~ Sex + ts.sunrise + MigStatus + Species, family = "binomial", data = birds.cs)
m127 <- glm(Fat.G ~ Sex + ts.sunrise + PercentAg + Species, family = "binomial", data = birds.cs)
m128 <- glm(Fat.G ~ Sex + DaysIntoSeason_S * MigStatus + Species, family = "binomial", data = birds.cs)
m129 <- glm(Fat.G ~ Sex + DaysIntoSeason_S + PercentAg + Species, family = "binomial", data = birds.cs)
m130 <- glm(Fat.G ~ Sex + MigStatus + PercentAg + Species, family = "binomial", data = birds.cs)
m131 <- glm(Fat.G ~ Event * ts.sunrise + DaysIntoSeason_S + Species, family = "binomial", data = birds.cs)
m132 <- glm(Fat.G ~ Event * ts.sunrise + MigStatus + Species, family = "binomial", data = birds.cs)
m133 <- glm(Fat.G ~ Event * ts.sunrise + PercentAg + Species, family = "binomial", data = birds.cs)
m134 <- glm(Fat.G ~ Event + DaysIntoSeason_S * MigStatus + Species, family = "binomial", data = birds.cs)
m135 <- glm(Fat.G ~ Event + DaysIntoSeason_S + PercentAg + Species, family = "binomial", data = birds.cs)
m136 <- glm(Fat.G ~ Event + MigStatus + PercentAg + Species, family = "binomial", data = birds.cs)
m137 <- glm(Fat.G ~ ts.sunrise + DaysIntoSeason_S * MigStatus + Species, family = "binomial", data = birds.cs)
m138 <- glm(Fat.G ~ ts.sunrise + DaysIntoSeason_S + PercentAg + Species, family = "binomial", data = birds.cs)
m139 <- glm(Fat.G ~ ts.sunrise + MigStatus + PercentAg + Species, family = "binomial", data = birds.cs)
m140 <- glm(Fat.G ~ DaysIntoSeason_S * MigStatus + PercentAg + Species, family = "binomial", data = birds.cs)

### Five Additive Combinations ####

m141 <- glm(Fat.G ~ Sex + Event * ts.sunrise + DaysIntoSeason_S * MigStatus +
              DaysIntoSeason_S * MigStatus + Event * ts.sunrise, family = "binomial", data = birds.cs)
m142 <- glm(Fat.G ~ Sex + Event * ts.sunrise + DaysIntoSeason_S + PercentAg +
            Event * ts.sunrise, family = "binomial", data = birds.cs)
m143 <- glm(Fat.G ~ Sex + Event * ts.sunrise + DaysIntoSeason_S + Detection +
               Event * ts.sunrise, family = "binomial", data = birds.cs)
m144 <- glm(Fat.G ~ Sex + Event * ts.sunrise + MigStatus + PercentAg +
               Event * ts.sunrise, family = "binomial", data = birds.cs)
m145 <- glm(Fat.G ~ Sex + Event * ts.sunrise + MigStatus + Detection +
               Event * ts.sunrise, family = "binomial", data = birds.cs)
m146 <- glm(Fat.G ~ Sex + Event * ts.sunrise + PercentAg + Detection +
               Event * ts.sunrise, family = "binomial", data = birds.cs)
m147 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S * MigStatus + PercentAg +
               DaysIntoSeason_S * MigStatus, family = "binomial", data = birds.cs)
m148 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S * MigStatus + Detection +
               + DaysIntoSeason_S * MigStatus, family = "binomial", data = birds.cs)
m149 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S + PercentAg + Detection, family = "binomial", data = birds.cs)
m150 <- glm(Fat.G ~ Sex + Event + MigStatus + PercentAg + Detection, family = "binomial", data = birds.cs)
m151 <- glm(Fat.G ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg +
               DaysIntoSeason_S * MigStatus, family = "binomial", data = birds.cs)
m152 <- glm(Fat.G ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + Detection +
               DaysIntoSeason_S * MigStatus, family = "binomial", data = birds.cs)
m153 <- glm(Fat.G ~ Sex + ts.sunrise + DaysIntoSeason_S + PercentAg + Detection, family = "binomial", data = birds.cs)
m154 <- glm(Fat.G ~ Sex + ts.sunrise + MigStatus + PercentAg + Detection, family = "binomial", data = birds.cs)
m155 <- glm(Fat.G ~ Sex + DaysIntoSeason_S * MigStatus + PercentAg + Detection, family = "binomial", data = birds.cs)
m156 <- glm(Fat.G ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg, family = "binomial", data = birds.cs)
m157 <- glm(Fat.G ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + Detection, family = "binomial", data = birds.cs)
m158 <- glm(Fat.G ~ Event * ts.sunrise + DaysIntoSeason_S + PercentAg + Detection, family = "binomial", data = birds.cs)
m159 <- glm(Fat.G ~ Event * ts.sunrise + MigStatus + PercentAg + Detection, family = "binomial", data = birds.cs)
m160 <- glm(Fat.G ~ Event + DaysIntoSeason_S * MigStatus + PercentAg + Detection, family = "binomial", data = birds.cs)
m161 <- glm(Fat.G ~ ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Detection, family = "binomial", data = birds.cs)
m162 <- glm(Fat.G ~ Sex + Event * ts.sunrise + DaysIntoSeason_S + Species +
              Event * ts.sunrise, family = "binomial", data = birds.cs)
m163 <- glm(Fat.G ~ Sex + Event * ts.sunrise + MigStatus + Species +
             Event * ts.sunrise, family = "binomial", data = birds.cs)
m164 <- glm(Fat.G ~ Sex + Event * ts.sunrise + PercentAg + Species +
              Event * ts.sunrise, family = "binomial", data = birds.cs)
m165 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S * MigStatus + Species +
              + DaysIntoSeason_S * MigStatus, family = "binomial", data = birds.cs)
m166 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S + PercentAg + Species, family = "binomial", data = birds.cs)
m167 <- glm(Fat.G ~ Sex + Event + MigStatus + PercentAg + Species, family = "binomial", data = birds.cs)
m168 <- glm(Fat.G ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + Species +
              DaysIntoSeason_S * MigStatus, family = "binomial", data = birds.cs)
m169 <- glm(Fat.G ~ Sex + ts.sunrise + DaysIntoSeason_S + PercentAg + Species, family = "binomial", data = birds.cs)
m170 <- glm(Fat.G ~ Sex + ts.sunrise + MigStatus + PercentAg + Species, family = "binomial", data = birds.cs)
m171 <- glm(Fat.G ~ Sex + DaysIntoSeason_S * MigStatus + PercentAg + Species, family = "binomial", data = birds.cs)
m172 <- glm(Fat.G ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg, family = "binomial", data = birds.cs)
m173 <- glm(Fat.G ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + Species, family = "binomial", data = birds.cs)
m174 <- glm(Fat.G ~ Event * ts.sunrise + DaysIntoSeason_S + PercentAg + Species, family = "binomial", data = birds.cs)
m175 <- glm(Fat.G ~ Event * ts.sunrise + MigStatus + PercentAg + Species, family = "binomial", data = birds.cs)
m176 <- glm(Fat.G ~ Event + DaysIntoSeason_S * MigStatus + PercentAg + Species, family = "binomial", data = birds.cs)
m177 <- glm(Fat.G ~ ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Species, family = "binomial", data = birds.cs)

### Six Additive Combinations ####

m178 <- glm(Fat.G ~ Sex + Event * ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg +
               Event * ts.sunrise + DaysIntoSeason_S * MigStatus, family = "binomial", data = birds.cs)
m179 <- glm(Fat.G ~ Sex + Event * ts.sunrise + DaysIntoSeason_S * MigStatus + Detection, family = "binomial", data = birds.cs)
m180 <- glm(Fat.G ~ Sex + Event * ts.sunrise + DaysIntoSeason_S + PercentAg + Detection, family = "binomial", data = birds.cs)
m181 <- glm(Fat.G ~ Sex + Event * ts.sunrise + MigStatus + PercentAg + Detection, family = "binomial", data = birds.cs)
m182 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S * MigStatus + PercentAg + Detection, family = "binomial", data = birds.cs)
m183 <- glm(Fat.G ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Detection, family = "binomial", data = birds.cs)
m184 <- glm(Fat.G ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Detection, family = "binomial", data = birds.cs)
m185 <- glm(Fat.G ~ Sex + Event * ts.sunrise + DaysIntoSeason_S * MigStatus + Species, family = "binomial", data = birds.cs)
m186 <- glm(Fat.G ~ Sex + Event * ts.sunrise + DaysIntoSeason_S + PercentAg + Species, family = "binomial", data = birds.cs)
m187 <- glm(Fat.G ~ Sex + Event * ts.sunrise + MigStatus + PercentAg + Species, family = "binomial", data = birds.cs)
m188 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S * MigStatus + PercentAg + Species, family = "binomial", data = birds.cs)
m189 <- glm(Fat.G ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Species, family = "binomial", data = birds.cs)
m190 <- glm(Fat.G ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Species, family = "binomial", data = birds.cs)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:190)

models <- mget(model_names)

models$m.null <- m.null
models$m.global <- m.global

model_names <- c(model_names, "m.null", "m.global")

aictab(models, modnames = model_names)

# ---------------------------------------------------------------------------- #

# Top Model Summaries ####

summary(m159)
pR2(m159)
# positive relationship between fat and agricultural intensity
# migrants have higher fat
# interaction between time and event on fat is significant

summary(m73)
pR2(m73)
# positive relationship between fat and date

summary(m68)
pR2(m68)
# positive relationship between fat and date

summary(m129)
pR2(m129)
# positive relationship between fat and date

summary(m157)
pR2(m157)
# positive relationship between fat and date

summary(m124)
pR2(m124)
# positive relationship between fat and agricultural intensity

summary(m167)
pR2(m167)
# positive relationship between fat and agricultural intensity

summary(m120)
pR2(m120)
# no significant covariates

summary(m184)
pR2(m184)
# no significant covariates

summary(m103)
pR2(m103)
# no significant covariates

summary(m144)
pR2(m144)
# interaction between time and event on fat is significant
# positive relationship between fat and agricultural intensity
# migrants have higher fat

summary(m181)
pR2(m181)
# interaction between time and event on fat is significant
# migrants have higher fat

summary(m160)
pR2(m160)
# no significant covariates

# ---------------------------------------------------------------------------- #

## Goodness of Fit for Top Models ####

# m124
birds.s <- birds.cs %>% 
  filter(!is.na(Sex))

m124 <- glm(Fat.G ~ Sex + Event + PercentAg + Species, family = "binomial", data = birds.s)

predicted_probs <- predict(m124, type = "response")
roc_curve <- roc(birds.s$Fat.G, predicted_probs)
plot(roc_curve)
auc(roc_curve)
#0.794

# m159
birds.s <- birds.cs %>% 
  filter(!is.na(Detection))

m159 <- glm(Fat.G ~ Event * ts.sunrise + MigStatus + PercentAg + Detection, 
            family = "binomial", data = birds.s)
predicted_probs <- predict(m159, type = "response")
roc_curve <- roc(birds.s$Fat.G, predicted_probs)
plot(roc_curve)
auc(roc_curve)
#0.767

# m167
birds.s <- birds.cs %>% 
  filter(!is.na(Sex))
m167 <- glm(Fat.G ~ Sex + Event + MigStatus + PercentAg + Species, family = "binomial", data = birds.s)
predicted_probs <- predict(m167, type = "response")
roc_curve <- roc(birds.s$Fat.G, predicted_probs)
plot(roc_curve)
auc(roc_curve)
#0.794

# m144
birds.s <- birds.cs %>% 
  filter(!is.na(Sex))
m144 <- glm(Fat.G ~ Sex + Event * ts.sunrise + MigStatus + PercentAg +
              Event * ts.sunrise, family = "binomial", data = birds.s)
predicted_probs <- predict(m144, type = "response")
roc_curve <- roc(birds.s$Fat.G, predicted_probs)
plot(roc_curve)
auc(roc_curve)
#0.772


## ALL TOP MODELS HAVE A GOOD FIT

# ---------------------------------------------------------------------------- #

# Model Assumptions

m159 <- glm(Fat.G ~ Event * ts.sunrise + MigStatus + PercentAg + Detection, 
            family = "binomial", data = birds.cs)

# Box-Tidwell transformation for PercentAg
birds$log_PercentAg <- log(birds$PercentAg)

# Fit model with interaction terms (PercentAg and its log transformation)
m_bt <- glm(Fat.G ~ PercentAg * log_PercentAg + Detection + MigStatus + 
              Event * ts.sunrise, family = "binomial", data = birds)

# View the summary of the Box-Tidwell model
summary(m_bt)

# LINEARITY ASSUMPTION HOLDS FOR AGRICULTURAL INTENSITY

# ---------------------------------------------------------------------------- #

# Plot Top Models

# Fat ~ Agricultural Intensity
birds$PercentAg <- birds$PercentAg * 100
m159 <- glm(Fat.G ~ Event * ts.sunrise + MigStatus + PercentAg + Detection, 
            family = "binomial", data = birds)

d <- data.frame(
  Event = factor("Fall 2023", levels = levels(birds$Event)),
  PercentAg = seq(min(birds$PercentAg), max(birds$PercentAg), length = 1000),
  Detection = factor("Detection", levels = levels(birds$Detection)),
  MigStatus = factor("Migratory", levels = levels(birds$MigStatus)),
  ts.sunrise = mean(birds$ts.sunrise)
)

predictions <- predict(m159, newdata = d, type = "response", se.fit = TRUE) 

d$fit <- predictions$fit

d$lower_CI <- predictions$fit - 1.96 * predictions$se.fit  # Lower CI
d$upper_CI <- predictions$fit + 1.96 * predictions$se.fit  # Upper CI

ggplot(d, aes(x = PercentAg, y = fit)) +
  geom_line(size = 0.8) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), 
              alpha = 0.25, color = NA, show.legend = FALSE) +  # Confidence intervals as ribbon
  theme_classic() +
  labs(x = "Surrounding Agricultural Intensity (%)", 
       y = "Predicted Probability of Having High Fat") +
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 13)),
        axis.title.y = element_text(size = 16, margin = margin(r = 13)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none") +
  theme(legend.position = "none")


