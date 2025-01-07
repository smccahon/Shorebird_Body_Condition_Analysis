#---------------------------------------#
#  All Species Pectoral Muscle Analysis #
#         Created 11/13/2024            #          
#         Modified 01/04/2024           #
#---------------------------------------#

# load packages
library(dplyr)
library(ggplot2)
library(AICcmodavg)
library(RColorBrewer)
library(trtools)
library(tidyverse)
library(lme4)
library(car)
library(viridis)
library(MuMIn)

# ---------------------------------------------------------------------------- #

# Data Processing & Manipulation ####

# Read data
setwd("processed_data")
birds <- read.csv("Shorebird_Data_Cleaned_2024-12-9.csv")

# Make neonicotinoid detection column (Detection/Non-detection)
birds$Detection <- ifelse(birds$OverallNeonic > 0, "Detection", "Non-detection")

# Make sampling event column
birds <- birds %>% 
  mutate(Event = paste(Season, Year, sep = " ")) %>% 
  mutate(Event = ifelse(Season %in% c("Spring", "Fall"),
                        paste(Season, Year, sep = " "),
                        NA))

# Logarithmic transformation of neonics
birds <- birds %>% 
  mutate(LogNeonic = log10(OverallNeonic + 0.0001))

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

# Filter out seasons that did not have any pectoral scores measured
birds <- birds %>% 
  group_by(Event) %>% 
  filter(any(!is.na(PecSizeBest))) %>% 
  ungroup()

# Filter data that only contain pectoral muscle scores
birds <- birds %>% 
  filter(!is.na(PecSizeBest))

# Only include species with at least three individuals
birds <- birds %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

# Reorder factor variables
birds$Sex <- factor(birds$Sex,
                    levels = c("M", "F"),
                    labels = c("Male", "Female"))

birds$Detection <- as.factor(birds$Detection)

birds$Species <- factor(birds$Species, 
                        levels = c("LesserYellowlegs", "Killdeer", "Willet", 
                                   "PectoralSandpiper", "WIPH", "LeastSandpiper", 
                                   "AmericanAvocet", "LongbilledDowitcher", 
                                   "SemipalmatedSandpiper"),
                        labels = c("Lesser Yellowlegs", "Killdeer", "Willet", 
                                   "Pectoral Sandpiper", "Wilson's Phalarope", 
                                   "Least Sandpiper", "American Avocet", 
                                   "Long-billed Dowitcher", "Semipalmated Sandpiper"))

birds$AgCategory <- factor(birds$AgCategory,
                           levels = c("Low", "Moderate", "High"))

# Create migratory status column
birds <- birds %>% 
  mutate(MigStatus = case_when(
    Species %in% c("Lesser Yellowlegs", "Long-billed Dowitcher",
                   "Semipalmated Sandpiper", "Pectoral Sandpiper",
                   "Least Sandpiper") ~ "Migratory",
    Species %in% c("Willet", "Killdeer", "American Avocet",
                   "Wilson's Phalarope") ~ "Resident",
    TRUE ~ "Unknown"
  ))

# Standardize continuous variables
birds.cs <- birds %>%
  mutate(across(where(is.numeric), scale))

# ---------------------------------------------------------------------------- #

# Data Visualization ####

# Pectoral Muscle ~ Concentration
ggplot(data = birds, aes(x = OverallNeonic, y = PecSizeBest, col = Species)) + 
  theme_classic() + 
  geom_point() + 
  labs(x = "Neonicotinoid Concentration (ug/L)", 
       y = "All Species Breast Muscle Size")

# PecSizeBest ~ Log(Concentration)
ggplot(data = birds, aes(x = LogNeonic, y = PecSizeBest, col = Species)) +
  theme_classic() +
  geom_point() + 
  labs(x = "Log(Concentration)",
       y = "All Species Breast Muscle Size")

# PecSizeBest ~ Detection
ggplot(data = birds %>% filter(!is.na(Detection)), 
       aes(x = Detection, y = PecSizeBest, col = Species)) + 
  theme_classic() + 
  geom_boxplot() + 
  geom_jitter() +
  labs(x = "Neonicotinoid Detection",
       y = "All Species Breast Muscle Size")

# PecSizeBest ~ Days into trapping season
ggplot(data = birds, aes(x = DaysIntoSeason_S, 
                         y = PecSizeBest, 
                         color = Species)) +
  theme_classic() +
  geom_point() +
  labs(x = "Date into Migration Season",
       y = "All Species Breast Muscle Size",
       color = "Species")

# PecSizeBest ~ Julian
ggplot(data = birds, aes(x = Julian, y = PecSizeBest, color = Species)) +
  theme_classic() +
  geom_point() +
  labs(x = "Julian Day",
       y = "All Species Breast Muscle Size",
       color = "Species")

# PecSizeBest ~ Sex
ggplot(data = birds %>% filter(!is.na(Sex)),
       aes(x = Sex, y = PecSizeBest, color = Species)) + 
  theme_classic() + 
  geom_boxplot() + 
  geom_jitter() +
  labs(x = "Sex",
       y = "All Species Breast Muscle Size")

# PecSizeBest ~ Sampling Event
ggplot(data = birds, aes(x = Event, y = PecSizeBest, color = Species)) +
  theme_classic() +
  geom_boxplot() +
  geom_jitter() + 
  labs(x = "Sampling Event",
       y = "All Species Breast Muscle Size",
       color = "Species")


# PecSizeBest ~ Capture Time
ggplot(data = birds, aes(x = ts.sunrise, y = PecSizeBest, color = Event)) +
  theme_classic() +
  geom_point() +
  labs(x = "Capture Time (min.)",
       y = "All Species Breast Muscle Size")

# PecSizeBest ~ AgIntensity
ggplot(data = birds, aes(x = AgCategory, y = PecSizeBest)) +
  theme_classic() +
  geom_boxplot()

ggplot(data = birds, aes(x = PercentAg, y = PecSizeBest)) +
  theme_classic() +
  geom_point()

# PecSizeBest ~ PecScore
birds$PecScore <- as.factor(birds$PecScore)
ggplot(data = birds, aes(x = PecScore, y = PecSizeBest)) +
  theme_classic() +
  geom_boxplot()

kruskal.test(PecSizeBest ~ PecScore, data = birds)

# PecSizeBest ~ MigStatus
ggplot(data = birds, aes(x = MigStatus, y = PecSizeBest)) +
  theme_classic() +
  geom_boxplot()

# ---------------------------------------------------------------------------- #

# Identify Best Random Structure: Detection ####

# Random Intercept Models
m1 <- lmer(PecSizeBest ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + 
             Detection + PercentAg + MigStatus + (1 | Species), data = birds.cs)

# Random Slope Models
# All other random slope models failed to converge
m2 <- lmer(PecSizeBest ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + 
             Detection + PercentAg + MigStatus + (Sex | Species), data = birds.cs)

m3 <- lmer(PecSizeBest ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + 
             Detection + PercentAg + MigStatus + (DaysIntoSeason_S | Species), data = birds.cs)

m4 <- lmer(PecSizeBest ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + 
             Detection + PercentAg + MigStatus + (Detection | Species), data = birds.cs)

m5 <- lmer(PecSizeBest ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + 
             Detection + PercentAg + MigStatus + (MigStatus | Species), data = birds.cs)

m6 <- lmer(PecSizeBest ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + 
             Detection + PercentAg + MigStatus + (PercentAg | Species), data = birds.cs)

model_names <- paste0("m", 1:6)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Random slope with days into season is the best model fit
# However, there are several issues with model convergence.
# The overall effect of date is weark and the random slope variance is not large.
# Random slope may not add much explanatory power to the model.
# Default to simple random intercept model.

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Neonicotinoids ####

m1 <- lmer(PecSizeBest ~ Detection + (1 | Species), 
           data = birds.cs, REML = FALSE)
m2 <- lmer(PecSizeBest ~ OverallNeonic + (1 | Species), 
           data = birds.cs, REML = FALSE)
m3 <- lmer(PecSizeBest ~ LogNeonic + (1 | Species), 
           data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# All perform similarly. Use log transformation and detection for consistency with other models.

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Agricultural Intensity ####

m1 <- lm(PecSizeBest ~ AgCategory, data = birds.cs)
m2 <- lm(PecSizeBest ~ PercentAg, data = birds.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Categorical predictor performs better

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Event and Capture Time ####

m1 <- lmer(PecSizeBest ~ Event * ts.sunrise + (1  | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(PecSizeBest ~ Event + I(ts.sunrise^2) + (1  | Species), data = birds.cs, REML = FALSE)
m3 <- lmer(PecSizeBest ~ Event + ts.sunrise + (1  | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with interaction by far performs the best.

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Date * Migratory Status ####

m1 <- lmer(PecSizeBest ~ DaysIntoSeason_S + MigStatus + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(PecSizeBest ~ DaysIntoSeason_S * MigStatus + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with interaction performs significantly better.

# ---------------------------------------------------------------------------- #

m.global <- lmer(PecSizeBest ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + 
                   MigStatus + AgCategory + Detection +
                   MigStatus * DaysIntoSeason_S +
                   Event * ts.sunrise + (1 | Species), REML = FALSE, data = birds.cs)

m.null <- lmer(PecSizeBest ~ 1 + (1|Species), data = birds.cs, REML = FALSE)

## Single Covariate Models ####

m1 <- lmer(PecSizeBest ~ Sex + (1|Species), data = birds.cs, REML = FALSE)
m2 <- lmer(PecSizeBest ~ Event + (1|Species), data = birds.cs, REML = FALSE)
m3 <- lmer(PecSizeBest ~ ts.sunrise + (1|Species), data = birds.cs, REML = FALSE)
m4 <- lmer(PecSizeBest ~ DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m5 <- lmer(PecSizeBest ~ MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m6 <- lmer(PecSizeBest ~ AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m7 <- lmer(PecSizeBest ~ Detection + (1|Species), data = birds.cs, REML = FALSE)

##  Additive Models ####

### Two additive combinations ####

m8 <- lmer(PecSizeBest ~ Sex + Event + (1|Species), data = birds.cs, REML = FALSE)
m9 <- lmer(PecSizeBest ~ Sex + ts.sunrise + (1|Species), data = birds.cs, REML = FALSE)
m10 <- lmer(PecSizeBest ~ Sex + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m11 <- lmer(PecSizeBest ~ Sex + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m12 <- lmer(PecSizeBest ~ Sex + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m13 <- lmer(PecSizeBest ~ Sex + Detection + (1|Species), data = birds.cs, REML = FALSE)
m14 <- lmer(PecSizeBest ~ Event + ts.sunrise + (1|Species) + Event*ts.sunrise, data = birds.cs, REML = FALSE)
m15 <- lmer(PecSizeBest ~ Event + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m16 <- lmer(PecSizeBest ~ Event + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m17 <- lmer(PecSizeBest ~ Event + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m18 <- lmer(PecSizeBest ~ Event + Detection + (1|Species), data = birds.cs, REML = FALSE)
m19 <- lmer(PecSizeBest ~ ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m20 <- lmer(PecSizeBest ~ ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m21 <- lmer(PecSizeBest ~ ts.sunrise + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m22 <- lmer(PecSizeBest ~ ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m23 <- lmer(PecSizeBest ~ DaysIntoSeason_S + MigStatus + (1|Species) + DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m24 <- lmer(PecSizeBest ~ DaysIntoSeason_S + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m25 <- lmer(PecSizeBest ~ DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m26 <- lmer(PecSizeBest ~ MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m27 <- lmer(PecSizeBest ~ MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m28 <- lmer(PecSizeBest ~ AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)

### Three Additive Combinations ####

m29 <- lmer(PecSizeBest ~ Sex + Event + ts.sunrise + (1|Species) + Event * ts.sunrise, data = birds.cs, REML = FALSE)
m30 <- lmer(PecSizeBest ~ Sex + Event + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m31 <- lmer(PecSizeBest ~ Sex + Event + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m32 <- lmer(PecSizeBest ~ Sex + Event + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m33 <- lmer(PecSizeBest ~ Sex + Event + Detection + (1|Species), data = birds.cs, REML = FALSE)
m34 <- lmer(PecSizeBest ~ Sex + ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m35 <- lmer(PecSizeBest ~ Sex + ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m36 <- lmer(PecSizeBest ~ Sex + ts.sunrise + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m37 <- lmer(PecSizeBest ~ Sex + ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m38 <- lmer(PecSizeBest ~ Sex + DaysIntoSeason_S + MigStatus + (1|Species) + DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m39 <- lmer(PecSizeBest ~ Sex + DaysIntoSeason_S + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m40 <- lmer(PecSizeBest ~ Sex + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m41 <- lmer(PecSizeBest ~ Sex + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m42 <- lmer(PecSizeBest ~ Sex + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m43 <- lmer(PecSizeBest ~ Sex + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m44 <- lmer(PecSizeBest ~ Event * ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m45 <- lmer(PecSizeBest ~ Event * ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m46 <- lmer(PecSizeBest ~ Event * ts.sunrise + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m47 <- lmer(PecSizeBest ~ Event * ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m48 <- lmer(PecSizeBest ~ Event + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m49 <- lmer(PecSizeBest ~ Event + DaysIntoSeason_S + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m50 <- lmer(PecSizeBest ~ Event + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m51 <- lmer(PecSizeBest ~ Event + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m52 <- lmer(PecSizeBest ~ Event + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m53 <- lmer(PecSizeBest ~ Event + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m54 <- lmer(PecSizeBest ~ ts.sunrise + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m55 <- lmer(PecSizeBest ~ ts.sunrise + DaysIntoSeason_S + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m56 <- lmer(PecSizeBest ~ ts.sunrise + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m57 <- lmer(PecSizeBest ~ ts.sunrise + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m58 <- lmer(PecSizeBest ~ ts.sunrise + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m59 <- lmer(PecSizeBest ~ ts.sunrise + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m60 <- lmer(PecSizeBest ~ DaysIntoSeason_S * MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m61 <- lmer(PecSizeBest ~ DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m62 <- lmer(PecSizeBest ~ DaysIntoSeason_S + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m63 <- lmer(PecSizeBest ~ MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)

### Four Additive Combinations ####

m64 <- lmer(PecSizeBest ~ Sex + Event * ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m65 <- lmer(PecSizeBest ~ Sex + Event * ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m66 <- lmer(PecSizeBest ~ Sex + Event * ts.sunrise + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m67 <- lmer(PecSizeBest ~ Sex + Event * ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m68 <- lmer(PecSizeBest ~ Sex + Event + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m69 <- lmer(PecSizeBest ~ Sex + Event + DaysIntoSeason_S + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m70 <- lmer(PecSizeBest ~ Sex + Event + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m71 <- lmer(PecSizeBest ~ Sex + Event + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m72 <- lmer(PecSizeBest ~ Sex + Event + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m73 <- lmer(PecSizeBest ~ Sex + Event + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m74 <- lmer(PecSizeBest ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m75 <- lmer(PecSizeBest ~ Sex + ts.sunrise + DaysIntoSeason_S + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m76 <- lmer(PecSizeBest ~ Sex + ts.sunrise + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m77 <- lmer(PecSizeBest ~ Sex + ts.sunrise + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m78 <- lmer(PecSizeBest ~ Sex + ts.sunrise + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m79 <- lmer(PecSizeBest ~ Sex + ts.sunrise + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m80 <- lmer(PecSizeBest ~ Sex + DaysIntoSeason_S * MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m81 <- lmer(PecSizeBest ~ Sex + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m82 <- lmer(PecSizeBest ~ Sex + DaysIntoSeason_S + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m83 <- lmer(PecSizeBest ~ Sex + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m84 <- lmer(PecSizeBest ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m85 <- lmer(PecSizeBest ~ Event * ts.sunrise + DaysIntoSeason_S + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m86 <- lmer(PecSizeBest ~ Event * ts.sunrise + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m87 <- lmer(PecSizeBest ~ Event * ts.sunrise + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m88 <- lmer(PecSizeBest ~ Event * ts.sunrise + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m89 <- lmer(PecSizeBest ~ Event * ts.sunrise + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m90 <- lmer(PecSizeBest ~ Event + DaysIntoSeason_S * MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m91 <- lmer(PecSizeBest ~ Event + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m92 <- lmer(PecSizeBest ~ Event + DaysIntoSeason_S + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m93 <- lmer(PecSizeBest ~ Event + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m94 <- lmer(PecSizeBest ~ ts.sunrise + DaysIntoSeason_S * MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m95 <- lmer(PecSizeBest ~ ts.sunrise + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m96 <- lmer(PecSizeBest ~ ts.sunrise + DaysIntoSeason_S + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m97 <- lmer(PecSizeBest ~ ts.sunrise + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m98 <- lmer(PecSizeBest ~ DaysIntoSeason_S * MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)

### Five Additive Combinations ####

m99 <- lmer(PecSizeBest ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + MigStatus + (1|Species) +
              DaysIntoSeason_S * MigStatus + Event * ts.sunrise, data = birds.cs, REML = FALSE)
m100 <- lmer(PecSizeBest ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + AgCategory + (1|Species) +
               Event * ts.sunrise, data = birds.cs, REML = FALSE)
m101 <- lmer(PecSizeBest ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + Detection + (1|Species) +
               Event * ts.sunrise, data = birds.cs, REML = FALSE)
m102 <- lmer(PecSizeBest ~ Sex + Event + ts.sunrise + MigStatus + AgCategory + (1|Species) +
               Event * ts.sunrise, data = birds.cs, REML = FALSE)
m103 <- lmer(PecSizeBest ~ Sex + Event + ts.sunrise + MigStatus + Detection + (1|Species) +
               Event * ts.sunrise, data = birds.cs, REML = FALSE)
m104 <- lmer(PecSizeBest ~ Sex + Event + ts.sunrise + AgCategory + Detection + (1|Species) +
               Event * ts.sunrise, data = birds.cs, REML = FALSE)
m105 <- lmer(PecSizeBest ~ Sex + Event + DaysIntoSeason_S + MigStatus + AgCategory + (1|Species) +
               DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m106 <- lmer(PecSizeBest ~ Sex + Event + DaysIntoSeason_S + MigStatus + Detection + (1|Species) +
               + DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m107 <- lmer(PecSizeBest ~ Sex + Event + DaysIntoSeason_S + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m108 <- lmer(PecSizeBest ~ Sex + Event + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m109 <- lmer(PecSizeBest ~ Sex + ts.sunrise + DaysIntoSeason_S + MigStatus + AgCategory + (1|Species) +
               DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m110 <- lmer(PecSizeBest ~ Sex + ts.sunrise + DaysIntoSeason_S + MigStatus + Detection + (1|Species) +
               DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m111 <- lmer(PecSizeBest ~ Sex + ts.sunrise + DaysIntoSeason_S + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m112 <- lmer(PecSizeBest ~ Sex + ts.sunrise + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m113 <- lmer(PecSizeBest ~ Sex + DaysIntoSeason_S * MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m114 <- lmer(PecSizeBest ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m115 <- lmer(PecSizeBest ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m116 <- lmer(PecSizeBest ~ Event * ts.sunrise + DaysIntoSeason_S + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m117 <- lmer(PecSizeBest ~ Event * ts.sunrise + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m118 <- lmer(PecSizeBest ~ Event + DaysIntoSeason_S * MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m119 <- lmer(PecSizeBest ~ ts.sunrise + DaysIntoSeason_S * MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)

### Six Additive Combinations ####

m120 <- lmer(PecSizeBest ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + MigStatus + AgCategory + (1|Species) +
               Event * ts.sunrise + DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m121 <- lmer(PecSizeBest ~ Sex + Event * ts.sunrise + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m122 <- lmer(PecSizeBest ~ Sex + Event * ts.sunrise + DaysIntoSeason_S + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m123 <- lmer(PecSizeBest ~ Sex + Event * ts.sunrise + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m124 <- lmer(PecSizeBest ~ Sex + Event + DaysIntoSeason_S * MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m125 <- lmer(PecSizeBest ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m126 <- lmer(PecSizeBest ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:126)

models <- mget(model_names)

models$m.null <- m.null
models$m.global <- m.global

model_names <- c(model_names, "m.null", "m.global")

aictab(models, modnames = model_names)

# ---------------------------------------------------------------------------- #

# Top Model Summaries ####

summary(m88)
confint(m88)

summary(m115)
confint(m115)
# ---------------------------------------------------------------------------- #

# Graph top models ####

## LAB MEETING DETECTION GRAPH ####
m <- lmer(PecSizeBest ~ Detection + Event * ts.sunrise + MigStatus + (1|Species), 
          data = birds, REML = FALSE)

d <- expand.grid(Event = c("Fall 2023"),
                 ts.sunrise = seq(min(birds$ts.sunrise), max(birds$ts.sunrise), 
                                  length = 1000),
                 Detection = c("Detection", "Non-detection"),
                 MigStatus = c("Migratory"),
                 Species = unique(birds$Species))

predictions <- predict(m, newdata = d, se.fit = TRUE, re.form = NULL)

d$fit <- predictions$fit

d$lwr <- d$fit - 1.96 * predictions$se.fit
d$upr <- d$fit + 1.96 * predictions$se.fit

d$Species <- factor(d$Species, levels = c("Lesser Yellowlegs", 
                                          "Pectoral Sandpiper", 
                                          "Killdeer", 
                                          "American Avocet", 
                                          "Long-billed Dowitcher", 
                                          "Willet", 
                                          "Wilson's Phalarope", 
                                          "Semipalmated Sandpiper",
                                          "Least Sandpiper"))

ggplot(d, aes(x = Detection, y = fit)) +
  geom_point(size = 0.1, col = "black") +  # Points showing predicted values
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1,
                col = "black",
                size = 0.2) +  # Add confidence intervals
  theme_light() +
  labs(x = NULL, 
       y = "Predicted Pectoral Muscle Size") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none",
        strip.text = element_text(size = 18)) +
  theme(legend.position = "none") + 
  facet_wrap(~ Species)


# Pectoral Muscle ~ Event * ts.sunrise (holding all else constant)
m110 <- lmer(PecSizeBest ~ Event * ts.sunrise + MigStatus + Detection + (1|Species), data = birds, REML = FALSE)

d <- expand.grid(Event = c("Spring 2022", "Spring 2023", "Fall 2023"),
                 ts.sunrise = seq(min(birds$ts.sunrise), max(birds$ts.sunrise), 
                                  length = 1000),
                 Detection = c("Detection"),
                 MigStatus = c("Migratory"),
                 Species = unique(birds$Species))

d$fit <- predict(m110, newdata = d, re.form = NULL)

ggplot(d, aes(x = ts.sunrise, y = fit, color = Event)) +
  geom_line(size = 0.9) +
  theme_light() +
  labs(x = "Time of Capture Since Sunrise (min)", 
       y = expression("Breast Muscle Size" ~~~ (mm[score])), 
       color = "Sampling Occasion") +
  facet_wrap(~ Species, scales = "free_y") +
  theme(legend.position = "right") +
  theme(legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12)) +
  scale_color_viridis_d(begin = 0.2, 
                        end = 0.8, 
                        option = "inferno") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 0, y = Inf, size = 3.5,
           label = "Sunrise", angle = 90, vjust = -0.5, hjust = 1) +
  geom_point(data = birds, aes(x = ts.sunrise, y = PecSizeBest))

# ACROSS ALL SPECIES
m <- lm(PecSizeBest ~ ts.sunrise * Event, data = birds)
d <- expand.grid(ts.sunrise = seq(min(birds$ts.sunrise),
                                  max(birds$ts.sunrise),
                                  length = 1000),
                 Event = c("Spring 2022", "Spring 2023", "Fall 2023"))
d <- cbind(d, predict(m, newdata = d, interval = "confidence"))

ggplot(d, aes(x = ts.sunrise, y = fit, col = Event)) + 
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Event), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_light() +
  labs(x = "Time of Capture since Sunrise (min)",
       y = expression("Breast Muscle Size" ~~~ (mm[score])),
       color = "Sampling Occasion") + 
  theme(legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        legend.position = c(0.98, 0.05),
        legend.justification = c(0.98, 0.05)) +
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_color_viridis_d(begin = 0.2, 
                        end = 0.8, 
                        option = "inferno") +
  scale_fill_viridis_d(begin = 0.2, 
                       end = 0.8, 
                       option = "inferno") +
  geom_point(data = birds, aes(x = ts.sunrise, y = PecSizeBest)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 0, y = max(d$fit) + 2, 
           label = "Sunrise", angle = 90, 
           vjust = -0.5, hjust = 0.5,
           size = 5)


# Pectoral Muscle ~ Detection
# Relationship not significant
ggplot(birds %>% filter(!is.na(Detection)), 
       aes(x = Detection, y = PecSizeBest, fill = Detection)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(color = Detection),
               position = position_jitterdodge(jitter.width = 0.2),  
              size = 2, 
              alpha = 0.8) + 
  theme_light() +
  facet_wrap(~ Species, scales = "free_y") +
  labs(x = "Neonicotinoid Detection Status", 
       y = expression("Breast Muscle Size" ~~~ (mm[score]))) +  
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none") +
  scale_fill_viridis_d(begin = 0.2, 
                        end = 0.8, 
                        option = "inferno",
                       alpha = 0.8) +
  scale_color_viridis_d(begin = 0.2, 
                        end = 0.8, 
                        option = "inferno", 
                        alpha = 0.8) 


# Pectoral Muscle ~ Detection
# Relationship not significant
ggplot(birds, 
       aes(x = MigStatus, y = PecSizeBest, fill = MigStatus)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(color = MigStatus),
              position = position_jitterdodge(jitter.width = 0.2),  
              size = 2, 
              alpha = 0.8) + 
  theme_light() +
  labs(x = "Migratory Status", 
       y = expression("Breast Muscle Size" ~~~ (mm[score]))) +  
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none") +
  scale_fill_viridis_d(begin = 0.2, 
                       end = 0.8, 
                       option = "inferno",
                       alpha = 0.8) +
  scale_color_viridis_d(begin = 0.2, 
                        end = 0.8, 
                        option = "inferno", 
                        alpha = 0.8) 


# Pectoral Size ~ LogNeonic (holding all else constant)

m132 <- lmer(PecSizeBest ~ Event * ts.sunrise + MigStatus + LogNeonic + (1|Species), 
             data = birds, REML = FALSE)

min_logneonic <- min(birds$LogNeonic, na.rm = TRUE)
max_logneonic <- max(birds$LogNeonic, na.rm = TRUE)

d <- expand.grid(Event = c("Fall 2023"),
                 LogNeonic = seq(min_logneonic, 
                                 max_logneonic, 
                                  length = 1000),
                 ts.sunrise = mean(birds$ts.sunrise),
                 MigStatus = c("Migratory"),
                 Species = unique(birds$Species))

d$fit <- predict(m132, newdata = d, re.form = NULL)

ggplot(d, aes(x = LogNeonic, y = fit)) +
  geom_line(size = 0.9) +
  theme_light() +
  labs(x = "Log(Neonicotinoid Concentration [ug/L])", 
       y = expression("Breast Muscle Size" ~~~ (mm[score]))) +
  facet_wrap(~ Species, scales = "free_y") +
  theme(legend.position = "right") +
  theme(legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12)) +
  scale_color_viridis_c(begin = 0.2, 
                        end = 0.8, 
                        option = "inferno") +
  geom_point(data = birds, aes(x = LogNeonic, y = PecSizeBest))

# ---------------------------------------------------------------------------- #

#OLD CODE FROM EARLY NOVEMBER

# Main finding: Birds in 2023 had significantly larger pectoral muscles...

# Pectoral Muscle ~ Event + Detection
m <- lmer(PecSizeBest ~ Event + Detection + (1 | Species), data = birds, REML = FALSE)

## All species ####
ggplot(birds %>% filter(!is.na(Detection)), aes(x = Detection, 
                                                y = PecSizeBest, 
                                                color = Event)) +
  geom_boxplot() +  
  theme_light() +
  labs(x = "Neonicotinoid Detection", 
       y = "Breast Muscle Size" ~~~ (mm[score]), 
       color = "Sampling Event") +
  facet_wrap(~ Species, scales = "free_y") +  
  theme(legend.position = "right",  
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12)) +
  scale_color_viridis_d(begin = 0.2, 
                        end = 0.8, 
                        option = "inferno")


# Pectoral Muscle ~ Event + LogNeonic
m <- lmer(PecSizeBest ~ Event + LogNeonic + (1 | Species), data = birds, REML = FALSE)

birds <- birds %>%
  filter(is.finite(LogNeonic)) 

d <- expand.grid(Event = c("Spring 2022", "Spring 2023", "Fall 2023"),
                 LogNeonic = seq(min(birds$LogNeonic), 
                                 max(birds$LogNeonic), 
                                 length = 100),
                 Species = unique(birds$Species))

d$PecSizeBest <- predict(m, newdata = d, re.form = NULL)


ggplot(d, aes(x = LogNeonic, y = PecSizeBest, color = Event)) +
  geom_line(size = 0.9) +
  theme_light() +
  labs(x = "Log(Neonicotinoid Concentration [ug/L])", 
       y = "Breast Muscle Size" ~~~ (mm[score]), 
       color = "Event") +
  facet_wrap(~ Species, scales = "free_y") +
  theme(legend.position = "right") +
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10)) +
  scale_color_viridis_d(begin = 0.2, 
                        end = 0.8, 
                        option = "inferno") +
  geom_point(data = birds, aes(x = LogNeonic, y = PecSizeBest))

# ---------------------------------------------------------------------------- #

# Model Assumptions

par(mfrow = c(2,2))
plot(m72)

# some violation of heterogeneity but not bad
fitted_vals <- fitted(m72)  
residuals_vals <- resid(m72)  
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

birds.c <- birds.cs[complete.cases(birds.cs$Sex, birds.cs$Event, birds.cs$ts.sunrise, 
                                   birds.cs$DaysIntoSeason_S, birds.cs$MigStatus, 
                                   birds.cs$PercentAg, birds.cs$Detection), ]

fitted_vals <- fitted(m72)  
residuals_vals <- resid(m72)  

# Plot residuals vs Event (okay, less variation in spring 2023)
boxplot(residuals_vals ~ birds.c$Event, xlab = "Event", ylab = "Residuals")

# Plot residuals vs MigStatus (pretty good)
boxplot(residuals_vals ~ birds.c$MigStatus, xlab = "MigStatus", ylab = "Residuals")

# Plot residuals vs Sex (pretty good)
boxplot(residuals_vals ~ birds.c$Sex, xlab = "Sex", ylab = "Residuals")

# Plot residuals vs Detection (pretty good)
boxplot(residuals_vals ~ birds.c$Detection, xlab = "Detection", ylab = "Residuals")

# Plot residuals vs Capture Time (pretty good)
plot(birds.c$ts.sunrise, residuals_vals, xlab = "Capture Time", ylab = "Residuals")

# Plot residuals vs Date (pretty good)
plot(birds.c$DaysIntoSeason_S, residuals_vals, xlab = "Date", ylab = "Residuals")

# Plot residuals vs Agricultural Intensity (slight violation)
plot(birds.c$PercentAg, residuals_vals, xlab = "Agricultural Intensity", ylab = "Residuals")

# Plot histogram of residuals (satisfied)
hist(residuals_vals, xlab = "Residuals", main = "")

# no violations of heterogeneity
m88 <- lmer(PecSizeBest ~ Event * ts.sunrise + MigStatus + Detection + (1|Species), 
            data = birds.c, REML = FALSE)

fitted_vals <- fitted(m88)  
residuals_vals <- resid(m88)  
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# Plot histogram of residuals (satisfied)
hist(residuals_vals, xlab = "Residuals", main = "")

# ---------------------------------------------------------------------------- #

# Goodness of Fit Test ####

# Marginal R²: The proportion of variance explained by the fixed effects alone.
# Conditional R²: The proportion of variance explained by both fixed and random effects.
r.squaredGLMM(m88)

r.squaredGLMM(m72)
