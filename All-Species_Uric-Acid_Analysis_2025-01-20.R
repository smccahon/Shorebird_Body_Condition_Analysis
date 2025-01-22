#-----------------------------------#
#   All Species Uric Acid Analysis  #
#           Created 01/20/2025      #          
#         Modified 01/22/2025       #
#-----------------------------------#

# load packages
library(dplyr)
library(ggplot2)
library(AICcmodavg)
library(tidyverse)
library(lme4)
library(car)

# ---------------------------------------------------------------------------- #

# Read data
setwd("processed_data")
birds <- read.csv("Shorebird_Data_Cleaned_2025-01-20.csv")

# Make neonicotinoid detection column (Detection/Non-detection)
birds$Detection <- ifelse(birds$OverallNeonic > 0, "Detection", "Non-detection")

# Reorder factor variables
birds$Sex <- factor(birds$Sex,
                    levels = c("M", "F"),
                    labels = c("Male", "Female"))

birds$Detection <- as.factor(birds$Detection)

birds$AgCategory <- factor(birds$AgCategory,
                           levels = c("Low", "Moderate", "High"))

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

# Only include birds with uric acid data
birds <- birds %>% 
  filter(!is.na(Uric))

# Only include species with at least three individuals
birds <- birds %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

table(birds$Species)

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

## Uric Acid ~ Detection ####
# No clear relationship
ggplot(birds, 
       aes(x = Detection, y = Uric, fill = Detection)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(color = Detection),
              position = position_jitterdodge(jitter.width = 0.2),  
              size = 2, 
              alpha = 0.8) + 
  theme_light() +
  labs(x = "Neonicotinoid Detection", 
       y = "Uric Acid Levels (umol)") +  
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

# Uric ~ PercentAg ####
# No clear relationship
ggplot(birds, 
       aes(x = PercentAg, y = Uric)) +
  geom_point() +
  theme_light() +
  labs(x = "Agricultural Intensity", 
       y = "All Species Uric Acid Levels (umol)") +  
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

# Uric ~ Time ####
# linear relationship looks evident
# Migration status likely has a big impact on uric acid levels
ggplot(birds, 
       aes(x = ts.sunrise, y = Uric, col = MigStatus)) +
  geom_point() +
  theme_light() +
  labs(x = "Time Since Sunrise", 
       y = "All Species Uric Acid Levels (umol)") +  
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

# Uric ~ Migration Status ####
# Migratory birds have more protein loss (makes sense)
ggplot(birds, 
       aes(x = MigStatus, y = Uric)) +
  geom_boxplot() +
  theme_light() +
  labs(x = "Migratory Status", 
       y = "All Species Uric Acid Levels (umol)") +  
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

m <- lm(Uric ~ MigStatus, data = birds)
summary(m)
confint(m)

# Uric ~ Date * MigStatus ####
# Migratory birds have more protein loss over time than residents (makes sense)
ggplot(birds, 
       aes(x = DaysIntoSeason_S, y = Uric, col = MigStatus)) +
  geom_point() +
  theme_light() +
  labs(x = "Date into Trapping Season", 
       y = "All Species Uric Acid Levels (umol)") +  
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))


# ---------------------------------------------------------------------------- #

## Univariate Analysis: Date ####

m1 <- lmer(Uric ~ DaysIntoSeason_S + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(Uric ~ Julian + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Both perform similarily; use date into season for consistency

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Event and Capture Time ####

m1 <- lmer(Uric ~ Event * ts.sunrise + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(Uric ~ Event + ts.sunrise + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Simple model without interaction by far performs the best.

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Date * Migratory Status ####

m1 <- lmer(Uric ~ DaysIntoSeason_S + MigStatus + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(Uric ~ DaysIntoSeason_S * MigStatus + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with interaction performs significantly better (wt = 90%)

# ---------------------------------------------------------------------------- #

# INTERACTIONS INCLUDED: Date * Migratory Status; PercentAg

# ---------------------------------------------------------------------------- #

# ISSUES WITH LOW VARIANCE IN SPECIES PARTICULARLY IN MODELS WITH MIGSTATUS
# MAYBE JUST USE SPECIES AS A FIXED EFFECT

m.global <- lmer(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S * MigStatus + 
                   PercentAg + Detection + (1 | Species), 
                 REML = FALSE, data = birds.cs)

m.null <- lmer(Uric ~ 1 + (1|Species), data = birds.cs, REML = FALSE)

## Single Covariate Models ####

m1 <- lmer(Uric ~ Sex + (1|Species), data = birds.cs, REML = FALSE)
m2 <- lmer(Uric ~ Event + (1|Species), data = birds.cs, REML = FALSE)
m3 <- lmer(Uric ~ ts.sunrise + (1|Species), data = birds.cs, REML = FALSE)
m4 <- lmer(Uric ~ DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m5 <- lmer(Uric ~ MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m6 <- lmer(Uric ~ PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m7 <- lmer(Uric ~ Detection + (1|Species), data = birds.cs, REML = FALSE)

##  Additive Models ####

### Two additive combinations ####

m8 <- lmer(Uric ~ Sex + Event + (1|Species), data = birds.cs, REML = FALSE)
m9 <- lmer(Uric ~ Sex + ts.sunrise + (1|Species), data = birds.cs, REML = FALSE)
m10 <- lmer(Uric ~ Sex + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m11 <- lmer(Uric ~ Sex + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m12 <- lmer(Uric ~ Sex + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m13 <- lmer(Uric ~ Sex + Detection + (1|Species), data = birds.cs, REML = FALSE)
m14 <- lmer(Uric ~ Event + ts.sunrise + (1|Species), data = birds.cs, REML = FALSE)
m15 <- lmer(Uric ~ Event + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m16 <- lmer(Uric ~ Event + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m17 <- lmer(Uric ~ Event + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m18 <- lmer(Uric ~ Event + Detection + (1|Species), data = birds.cs, REML = FALSE)
m19 <- lmer(Uric ~ ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m20 <- lmer(Uric ~ ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m21 <- lmer(Uric ~ ts.sunrise + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m22 <- lmer(Uric ~ ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m23 <- lmer(Uric ~ DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m24 <- lmer(Uric ~ DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m25 <- lmer(Uric ~ DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m26 <- lmer(Uric ~ MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m27 <- lmer(Uric ~ MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m28 <- lmer(Uric ~ PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)

### Three Additive Combinations ####

m29 <- lmer(Uric ~ Sex + Event + ts.sunrise + (1|Species), data = birds.cs, REML = FALSE)
m30 <- lmer(Uric ~ Sex + Event + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m31 <- lmer(Uric ~ Sex + Event + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m32 <- lmer(Uric ~ Sex + Event + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m33 <- lmer(Uric ~ Sex + Event + Detection + (1|Species), data = birds.cs, REML = FALSE)
m34 <- lmer(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m35 <- lmer(Uric ~ Sex + ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m36 <- lmer(Uric ~ Sex + ts.sunrise + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m37 <- lmer(Uric ~ Sex + ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m38 <- lmer(Uric ~ Sex + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m39 <- lmer(Uric ~ Sex + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m40 <- lmer(Uric ~ Sex + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m41 <- lmer(Uric ~ Sex + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m42 <- lmer(Uric ~ Sex + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m43 <- lmer(Uric ~ Sex + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m44 <- lmer(Uric ~ Event + ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m45 <- lmer(Uric ~ Event + ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m46 <- lmer(Uric ~ Event + ts.sunrise + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m47 <- lmer(Uric ~ Event + ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m48 <- lmer(Uric ~ Event + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m49 <- lmer(Uric ~ Event + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m50 <- lmer(Uric ~ Event + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m51 <- lmer(Uric ~ Event + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m52 <- lmer(Uric ~ Event + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m53 <- lmer(Uric ~ Event + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m54 <- lmer(Uric ~ ts.sunrise + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m55 <- lmer(Uric ~ ts.sunrise + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m56 <- lmer(Uric ~ ts.sunrise + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m57 <- lmer(Uric ~ ts.sunrise + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m58 <- lmer(Uric ~ ts.sunrise + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m59 <- lmer(Uric ~ ts.sunrise + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m60 <- lmer(Uric ~ DaysIntoSeason_S * MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m61 <- lmer(Uric ~ DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m62 <- lmer(Uric ~ DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m63 <- lmer(Uric ~ MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)

### Four Additive Combinations ####

m64 <- lmer(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m65 <- lmer(Uric ~ Sex + Event + ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m66 <- lmer(Uric ~ Sex + Event + ts.sunrise + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m67 <- lmer(Uric ~ Sex + Event + ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m68 <- lmer(Uric ~ Sex + Event + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m69 <- lmer(Uric ~ Sex + Event + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m70 <- lmer(Uric ~ Sex + Event + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m71 <- lmer(Uric ~ Sex + Event + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m72 <- lmer(Uric ~ Sex + Event + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m73 <- lmer(Uric ~ Sex + Event + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m74 <- lmer(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m75 <- lmer(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m76 <- lmer(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m77 <- lmer(Uric ~ Sex + ts.sunrise + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m78 <- lmer(Uric ~ Sex + ts.sunrise + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m79 <- lmer(Uric ~ Sex + ts.sunrise + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m80 <- lmer(Uric ~ Sex + DaysIntoSeason_S * MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m81 <- lmer(Uric ~ Sex + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m82 <- lmer(Uric ~ Sex + DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m83 <- lmer(Uric ~ Sex + MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m84 <- lmer(Uric ~ Event + ts.sunrise + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m85 <- lmer(Uric ~ Event + ts.sunrise + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m86 <- lmer(Uric ~ Event + ts.sunrise + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m87 <- lmer(Uric ~ Event + ts.sunrise + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m88 <- lmer(Uric ~ Event + ts.sunrise + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m89 <- lmer(Uric ~ Event + ts.sunrise + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m90 <- lmer(Uric ~ Event + DaysIntoSeason_S * MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m91 <- lmer(Uric ~ Event + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m92 <- lmer(Uric ~ Event + DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m93 <- lmer(Uric ~ Event + MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m94 <- lmer(Uric ~ ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m95 <- lmer(Uric ~ ts.sunrise + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m96 <- lmer(Uric ~ ts.sunrise + DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m97 <- lmer(Uric ~ ts.sunrise + MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m98 <- lmer(Uric ~ DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)

### Five Additive Combinations ####

m99 <- lmer(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m100 <- lmer(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m101 <- lmer(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m102 <- lmer(Uric ~ Sex + Event + ts.sunrise + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m103 <- lmer(Uric ~ Sex + Event + ts.sunrise + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m104 <- lmer(Uric ~ Sex + Event + ts.sunrise + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m105 <- lmer(Uric ~ Sex + Event + DaysIntoSeason_S * MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m106 <- lmer(Uric ~ Sex + Event + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m107 <- lmer(Uric ~ Sex + Event + DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m108 <- lmer(Uric ~ Sex + Event + MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m109 <- lmer(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m110 <- lmer(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m111 <- lmer(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m112 <- lmer(Uric ~ Sex + ts.sunrise + MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m113 <- lmer(Uric ~ Sex + DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m114 <- lmer(Uric ~ Event + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m115 <- lmer(Uric ~ Event + ts.sunrise + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m116 <- lmer(Uric ~ Event + ts.sunrise + DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m117 <- lmer(Uric ~ Event + ts.sunrise + MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m118 <- lmer(Uric ~ Event + DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m119 <- lmer(Uric ~ ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)

### Six Additive Combinations ####

m120 <- lmer(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m121 <- lmer(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m122 <- lmer(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m123 <- lmer(Uric ~ Sex + Event + ts.sunrise + MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m124 <- lmer(Uric ~ Sex + Event + DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m125 <- lmer(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m126 <- lmer(Uric ~ Event + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:126)

models <- mget(model_names)

models$m.null <- m.null
models$m.global <- m.global

model_names <- c(model_names, "m.null", "m.global")

aictab(models, modnames = model_names)

# ---------------------------------------------------------------------------- #

# Top model summaries ####
summary(m54)
confint(m54) # time and date * migstatus are significant

summary(m94)
confint(m94) # date * migstatus is significant

# ---------------------------------------------------------------------------- #

# TRY WITH SPECIES AS A FIXED EFFECT ONLY ####


## Interaction Analysis: Event and Capture Time ####

m1 <- lm(Uric ~ Event * ts.sunrise, data = birds.cs)
m2 <- lm(Uric ~ Event + ts.sunrise, data = birds.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Simple model without interaction by far performs the best.

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Date * Migratory Status ####

m1 <- lm(Uric ~ DaysIntoSeason_S + MigStatus, data = birds.cs)
m2 <- lm(Uric ~ DaysIntoSeason_S * MigStatus, data = birds.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with interaction performs significantly better (wt = 91%)

# ---------------------------------------------------------------------------- #

# INTERACTIONS INCLUDED: Date * Migratory Status; PercentAg

# ---------------------------------------------------------------------------- #


m.global <- lm(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + 
                  MigStatus + PercentAg + Detection +
                  MigStatus * DaysIntoSeason_S + Species,  data = birds.cs)

m.null <- lm(Uric ~ 1,  data = birds.cs)

## Single Covariate Models ####

m1 <- lm(Uric ~ Sex,  data = birds.cs)
m2 <- lm(Uric ~ Event,  data = birds.cs)
m3 <- lm(Uric ~ ts.sunrise,  data = birds.cs)
m4 <- lm(Uric ~ DaysIntoSeason_S,  data = birds.cs)
m5 <- lm(Uric ~ MigStatus,  data = birds.cs)
m6 <- lm(Uric ~ PercentAg,  data = birds.cs)
m7 <- lm(Uric ~ Detection,  data = birds.cs)
m8 <- lm(Uric ~ Species,  data = birds.cs)

##  Additive Models ####

### Two additive combinations ####

m9 <- lm(Uric ~ Sex + Event,  data = birds.cs)
m10 <- lm(Uric ~ Sex + ts.sunrise,  data = birds.cs)
m11 <- lm(Uric ~ Sex + DaysIntoSeason_S,  data = birds.cs)
m12 <- lm(Uric ~ Sex + MigStatus,  data = birds.cs)
m13 <- lm(Uric ~ Sex + PercentAg,  data = birds.cs)
m14 <- lm(Uric ~ Sex + Detection,  data = birds.cs)
m15 <- lm(Uric ~ Event + ts.sunrise + Event*ts.sunrise,  data = birds.cs)
m16 <- lm(Uric ~ Event + DaysIntoSeason_S,  data = birds.cs)
m17 <- lm(Uric ~ Event + MigStatus,  data = birds.cs)
m18 <- lm(Uric ~ Event + PercentAg,  data = birds.cs)
m19 <- lm(Uric ~ Event + Detection,  data = birds.cs)
m20 <- lm(Uric ~ ts.sunrise + DaysIntoSeason_S,  data = birds.cs)
m21 <- lm(Uric ~ ts.sunrise + MigStatus,  data = birds.cs)
m22 <- lm(Uric ~ ts.sunrise + PercentAg,  data = birds.cs)
m23 <- lm(Uric ~ ts.sunrise + Detection,  data = birds.cs)
m24 <- lm(Uric ~ DaysIntoSeason_S * MigStatus + DaysIntoSeason_S * MigStatus,  data = birds.cs)
m25 <- lm(Uric ~ DaysIntoSeason_S + PercentAg,  data = birds.cs)
m26 <- lm(Uric ~ DaysIntoSeason_S + Detection,  data = birds.cs)
m27 <- lm(Uric ~ MigStatus + PercentAg,  data = birds.cs)
m28 <- lm(Uric ~ MigStatus + Detection,  data = birds.cs)
m29 <- lm(Uric ~ PercentAg + Detection,  data = birds.cs)
m30 <- lm(Uric ~ Sex + Species,  data = birds.cs)
m31 <- lm(Uric ~ Event + Species,  data = birds.cs)
m32 <- lm(Uric ~ ts.sunrise + Species,  data = birds.cs)
m33 <- lm(Uric ~ DaysIntoSeason_S + Species,  data = birds.cs)
m34 <- lm(Uric ~ MigStatus + Species,  data = birds.cs)
m35 <- lm(Uric ~ PercentAg + Species,  data = birds.cs)

### Three Additive Combinations ####

m36 <- lm(Uric ~ Sex + Event + ts.sunrise + Event + ts.sunrise,  data = birds.cs)
m37 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S,  data = birds.cs)
m38 <- lm(Uric ~ Sex + Event + MigStatus,  data = birds.cs)
m39 <- lm(Uric ~ Sex + Event + PercentAg,  data = birds.cs)
m40 <- lm(Uric ~ Sex + Event + Detection,  data = birds.cs)
m41 <- lm(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S,  data = birds.cs)
m42 <- lm(Uric ~ Sex + ts.sunrise + MigStatus,  data = birds.cs)
m43 <- lm(Uric ~ Sex + ts.sunrise + PercentAg,  data = birds.cs)
m44 <- lm(Uric ~ Sex + ts.sunrise + Detection,  data = birds.cs)
m45 <- lm(Uric ~ Sex + DaysIntoSeason_S * MigStatus + DaysIntoSeason_S * MigStatus,  data = birds.cs)
m46 <- lm(Uric ~ Sex + DaysIntoSeason_S + PercentAg,  data = birds.cs)
m47 <- lm(Uric ~ Sex + DaysIntoSeason_S + Detection,  data = birds.cs)
m48 <- lm(Uric ~ Sex + MigStatus + PercentAg,  data = birds.cs)
m49 <- lm(Uric ~ Sex + MigStatus + Detection,  data = birds.cs)
m50 <- lm(Uric ~ Sex + PercentAg + Detection,  data = birds.cs)
m51 <- lm(Uric ~ Event + ts.sunrise + DaysIntoSeason_S,  data = birds.cs)
m52 <- lm(Uric ~ Event + ts.sunrise + MigStatus,  data = birds.cs)
m53 <- lm(Uric ~ Event + ts.sunrise + PercentAg,  data = birds.cs)
m54 <- lm(Uric ~ Event + ts.sunrise + Detection,  data = birds.cs)
m55 <- lm(Uric ~ Event + DaysIntoSeason_S * MigStatus,  data = birds.cs)
m56 <- lm(Uric ~ Event + DaysIntoSeason_S + PercentAg,  data = birds.cs)
m57 <- lm(Uric ~ Event + DaysIntoSeason_S + Detection,  data = birds.cs)
m58 <- lm(Uric ~ Event + MigStatus + PercentAg,  data = birds.cs)
m59 <- lm(Uric ~ Event + MigStatus + Detection,  data = birds.cs)
m60 <- lm(Uric ~ Event + PercentAg + Detection,  data = birds.cs)
m61 <- lm(Uric ~ ts.sunrise + DaysIntoSeason_S * MigStatus,  data = birds.cs)
m62 <- lm(Uric ~ ts.sunrise + DaysIntoSeason_S + PercentAg,  data = birds.cs)
m63 <- lm(Uric ~ ts.sunrise + DaysIntoSeason_S + Detection,  data = birds.cs)
m64 <- lm(Uric ~ ts.sunrise + MigStatus + PercentAg,  data = birds.cs)
m65 <- lm(Uric ~ ts.sunrise + MigStatus + Detection,  data = birds.cs)
m66 <- lm(Uric ~ ts.sunrise + PercentAg + Detection,  data = birds.cs)
m67 <- lm(Uric ~ DaysIntoSeason_S * MigStatus + PercentAg,  data = birds.cs)
m68 <- lm(Uric ~ DaysIntoSeason_S * MigStatus + Detection,  data = birds.cs)
m69 <- lm(Uric ~ DaysIntoSeason_S + PercentAg + Detection,  data = birds.cs)
m70 <- lm(Uric ~ MigStatus + PercentAg + Detection,  data = birds.cs)
m71 <- lm(Uric ~ Sex + Event + Species,  data = birds.cs)
m72 <- lm(Uric ~ Sex + ts.sunrise + Species,  data = birds.cs)
m73 <- lm(Uric ~ Sex + DaysIntoSeason_S + Species,  data = birds.cs)
m74 <- lm(Uric ~ Sex + MigStatus + Species,  data = birds.cs)
m75 <- lm(Uric ~ Sex + PercentAg + Species,  data = birds.cs)
m76 <- lm(Uric ~ Event + ts.sunrise + Species,  data = birds.cs)
m77 <- lm(Uric ~ Event + DaysIntoSeason_S + Species,  data = birds.cs)
m78 <- lm(Uric ~ Event + MigStatus + Species,  data = birds.cs)
m79 <- lm(Uric ~ Event + PercentAg + Species,  data = birds.cs)
m80 <- lm(Uric ~ ts.sunrise + DaysIntoSeason_S + Species,  data = birds.cs)
m81 <- lm(Uric ~ ts.sunrise + MigStatus + Species,  data = birds.cs)
m82 <- lm(Uric ~ ts.sunrise + PercentAg + Species,  data = birds.cs)
m83 <- lm(Uric ~ DaysIntoSeason_S * MigStatus + Species,  data = birds.cs)
m84 <- lm(Uric ~ DaysIntoSeason_S + PercentAg + Species,  data = birds.cs)
m85 <- lm(Uric ~ MigStatus + PercentAg + Species,  data = birds.cs)

### Four Additive Combinations ####

m86 <- lm(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S,  data = birds.cs)
m87 <- lm(Uric ~ Sex + Event + ts.sunrise + MigStatus,  data = birds.cs)
m88 <- lm(Uric ~ Sex + Event + ts.sunrise + PercentAg,  data = birds.cs)
m89 <- lm(Uric ~ Sex + Event + ts.sunrise + Detection,  data = birds.cs)
m90 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S * MigStatus,  data = birds.cs)
m91 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S + PercentAg,  data = birds.cs)
m92 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S + Detection,  data = birds.cs)
m93 <- lm(Uric ~ Sex + Event + MigStatus + PercentAg,  data = birds.cs)
m94 <- lm(Uric ~ Sex + Event + MigStatus + Detection + Event + MigStatus,  data = birds.cs)
m95 <- lm(Uric ~ Sex + Event + PercentAg + Detection,  data = birds.cs)
m96 <- lm(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus,  data = birds.cs)
m97 <- lm(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S + PercentAg,  data = birds.cs)
m98 <- lm(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S + Detection,  data = birds.cs)
m99 <- lm(Uric ~ Sex + ts.sunrise + MigStatus + PercentAg,  data = birds.cs)
m100 <- lm(Uric ~ Sex + ts.sunrise + MigStatus + Detection,  data = birds.cs)
m101 <- lm(Uric ~ Sex + ts.sunrise + PercentAg + Detection,  data = birds.cs)
m102 <- lm(Uric ~ Sex + DaysIntoSeason_S * MigStatus + PercentAg,  data = birds.cs)
m103 <- lm(Uric ~ Sex + DaysIntoSeason_S * MigStatus + Detection,  data = birds.cs)
m104 <- lm(Uric ~ Sex + DaysIntoSeason_S + PercentAg + Detection,  data = birds.cs)
m105 <- lm(Uric ~ Sex + MigStatus + PercentAg + Detection,  data = birds.cs)
m106 <- lm(Uric ~ Event + ts.sunrise + DaysIntoSeason_S * MigStatus,  data = birds.cs)
m107 <- lm(Uric ~ Event + ts.sunrise + DaysIntoSeason_S + PercentAg,  data = birds.cs)
m108 <- lm(Uric ~ Event + ts.sunrise + DaysIntoSeason_S + Detection,  data = birds.cs)
m109 <- lm(Uric ~ Event + ts.sunrise + MigStatus + PercentAg,  data = birds.cs)
m110 <- lm(Uric ~ Event + ts.sunrise + MigStatus + Detection,  data = birds.cs)
m111 <- lm(Uric ~ Event + ts.sunrise + PercentAg + Detection,  data = birds.cs)
m112 <- lm(Uric ~ Event + DaysIntoSeason_S * MigStatus + PercentAg,  data = birds.cs)
m113 <- lm(Uric ~ Event + DaysIntoSeason_S * MigStatus + Detection,  data = birds.cs)
m114 <- lm(Uric ~ Event + DaysIntoSeason_S + PercentAg + Detection,  data = birds.cs)
m115 <- lm(Uric ~ Event + MigStatus + PercentAg + Detection,  data = birds.cs)
m116 <- lm(Uric ~ ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg,  data = birds.cs)
m117 <- lm(Uric ~ ts.sunrise + DaysIntoSeason_S * MigStatus + Detection,  data = birds.cs)
m118 <- lm(Uric ~ ts.sunrise + DaysIntoSeason_S + PercentAg + Detection,  data = birds.cs)
m119 <- lm(Uric ~ ts.sunrise + MigStatus + PercentAg + Detection,  data = birds.cs)
m120 <- lm(Uric ~ DaysIntoSeason_S * MigStatus + PercentAg + Detection,  data = birds.cs)
m121 <- lm(Uric ~ Sex + Event + ts.sunrise + Species,  data = birds.cs)
m122 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S + Species,  data = birds.cs)
m123 <- lm(Uric ~ Sex + Event + MigStatus + Species + Event + MigStatus,  data = birds.cs)
m124 <- lm(Uric ~ Sex + Event + PercentAg + Species,  data = birds.cs)
m125 <- lm(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S + Species,  data = birds.cs)
m126 <- lm(Uric ~ Sex + ts.sunrise + MigStatus + Species,  data = birds.cs)
m127 <- lm(Uric ~ Sex + ts.sunrise + PercentAg + Species,  data = birds.cs)
m128 <- lm(Uric ~ Sex + DaysIntoSeason_S * MigStatus + Species,  data = birds.cs)
m129 <- lm(Uric ~ Sex + DaysIntoSeason_S + PercentAg + Species,  data = birds.cs)
m130 <- lm(Uric ~ Sex + MigStatus + PercentAg + Species,  data = birds.cs)
m131 <- lm(Uric ~ Event + ts.sunrise + DaysIntoSeason_S + Species,  data = birds.cs)
m132 <- lm(Uric ~ Event + ts.sunrise + MigStatus + Species,  data = birds.cs)
m133 <- lm(Uric ~ Event + ts.sunrise + PercentAg + Species,  data = birds.cs)
m134 <- lm(Uric ~ Event + DaysIntoSeason_S * MigStatus + Species,  data = birds.cs)
m135 <- lm(Uric ~ Event + DaysIntoSeason_S + PercentAg + Species,  data = birds.cs)
m136 <- lm(Uric ~ Event + MigStatus + PercentAg + Species,  data = birds.cs)
m137 <- lm(Uric ~ ts.sunrise + DaysIntoSeason_S * MigStatus + Species,  data = birds.cs)
m138 <- lm(Uric ~ ts.sunrise + DaysIntoSeason_S + PercentAg + Species,  data = birds.cs)
m139 <- lm(Uric ~ ts.sunrise + MigStatus + PercentAg + Species,  data = birds.cs)
m140 <- lm(Uric ~ DaysIntoSeason_S * MigStatus + PercentAg + Species,  data = birds.cs)

### Five Additive Combinations ####

m141 <- lm(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S * MigStatus +
              DaysIntoSeason_S * MigStatus + Event + ts.sunrise,  data = birds.cs)
m142 <- lm(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + PercentAg +
              Event + ts.sunrise,  data = birds.cs)
m143 <- lm(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + Detection +
              Event + ts.sunrise,  data = birds.cs)
m144 <- lm(Uric ~ Sex + Event + ts.sunrise + MigStatus + PercentAg +
              Event + ts.sunrise,  data = birds.cs)
m145 <- lm(Uric ~ Sex + Event + ts.sunrise + MigStatus + Detection +
              Event + ts.sunrise,  data = birds.cs)
m146 <- lm(Uric ~ Sex + Event + ts.sunrise + PercentAg + Detection +
              Event + ts.sunrise,  data = birds.cs)
m147 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S * MigStatus + PercentAg +
              DaysIntoSeason_S * MigStatus,  data = birds.cs)
m148 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S * MigStatus + Detection +
              + DaysIntoSeason_S * MigStatus,  data = birds.cs)
m149 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S + PercentAg + Detection,  data = birds.cs)
m150 <- lm(Uric ~ Sex + Event + MigStatus + PercentAg + Detection,  data = birds.cs)
m151 <- lm(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg +
              DaysIntoSeason_S * MigStatus,  data = birds.cs)
m152 <- lm(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + Detection +
              DaysIntoSeason_S * MigStatus,  data = birds.cs)
m153 <- lm(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S + PercentAg + Detection,  data = birds.cs)
m154 <- lm(Uric ~ Sex + ts.sunrise + MigStatus + PercentAg + Detection,  data = birds.cs)
m155 <- lm(Uric ~ Sex + DaysIntoSeason_S * MigStatus + PercentAg + Detection,  data = birds.cs)
m156 <- lm(Uric ~ Event + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg,  data = birds.cs)
m157 <- lm(Uric ~ Event + ts.sunrise + DaysIntoSeason_S * MigStatus + Detection,  data = birds.cs)
m158 <- lm(Uric ~ Event + ts.sunrise + DaysIntoSeason_S + PercentAg + Detection,  data = birds.cs)
m159 <- lm(Uric ~ Event + ts.sunrise + MigStatus + PercentAg + Detection,  data = birds.cs)
m160 <- lm(Uric ~ Event + DaysIntoSeason_S * MigStatus + PercentAg + Detection,  data = birds.cs)
m161 <- lm(Uric ~ ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Detection,  data = birds.cs)
m162 <- lm(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + Species +
              Event + ts.sunrise,  data = birds.cs)
m163 <- lm(Uric ~ Sex + Event + ts.sunrise + MigStatus + Species +
              Event + ts.sunrise,  data = birds.cs)
m164 <- lm(Uric ~ Sex + Event + ts.sunrise + PercentAg + Species +
              Event + ts.sunrise,  data = birds.cs)
m165 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S * MigStatus + Species +
              + DaysIntoSeason_S * MigStatus,  data = birds.cs)
m166 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S + PercentAg + Species,  data = birds.cs)
m167 <- lm(Uric ~ Sex + Event + MigStatus + PercentAg + Species,  data = birds.cs)
m168 <- lm(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + Species +
              DaysIntoSeason_S * MigStatus,  data = birds.cs)
m169 <- lm(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S + PercentAg + Species,  data = birds.cs)
m170 <- lm(Uric ~ Sex + ts.sunrise + MigStatus + PercentAg + Species,  data = birds.cs)
m171 <- lm(Uric ~ Sex + DaysIntoSeason_S * MigStatus + PercentAg + Species,  data = birds.cs)
m172 <- lm(Uric ~ Event + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg,  data = birds.cs)
m173 <- lm(Uric ~ Event + ts.sunrise + DaysIntoSeason_S * MigStatus + Species,  data = birds.cs)
m174 <- lm(Uric ~ Event + ts.sunrise + DaysIntoSeason_S + PercentAg + Species,  data = birds.cs)
m175 <- lm(Uric ~ Event + ts.sunrise + MigStatus + PercentAg + Species,  data = birds.cs)
m176 <- lm(Uric ~ Event + DaysIntoSeason_S * MigStatus + PercentAg + Species,  data = birds.cs)
m177 <- lm(Uric ~ ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Species,  data = birds.cs)

### Six Additive Combinations ####

m178 <- lm(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg +
              Event + ts.sunrise + DaysIntoSeason_S * MigStatus,  data = birds.cs)
m179 <- lm(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S * MigStatus + Detection,  data = birds.cs)
m180 <- lm(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + PercentAg + Detection,  data = birds.cs)
m181 <- lm(Uric ~ Sex + Event + ts.sunrise + MigStatus + PercentAg + Detection,  data = birds.cs)
m182 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S * MigStatus + PercentAg + Detection,  data = birds.cs)
m183 <- lm(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Detection,  data = birds.cs)
m184 <- lm(Uric ~ Event + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Detection,  data = birds.cs)
m185 <- lm(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S * MigStatus + Species,  data = birds.cs)
m186 <- lm(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + PercentAg + Species,  data = birds.cs)
m187 <- lm(Uric ~ Sex + Event + ts.sunrise + MigStatus + PercentAg + Species,  data = birds.cs)
m188 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S * MigStatus + PercentAg + Species,  data = birds.cs)
m189 <- lm(Uric ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Species,  data = birds.cs)
m190 <- lm(Uric ~ Event + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Species,  data = birds.cs)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:190)

models <- mget(model_names)

models$m.null <- m.null
models$m.global <- m.global

model_names <- c(model_names, "m.null", "m.global")

aictab(models, modnames = model_names)

# ---------------------------------------------------------------------------- #

# Top model summaries ####
summary(m61)
confint(m61) # time and date * migstatus are significant
plot(m61)
# heterogeneity issues

summary(m116)
confint(m116) # date * migstatus is significant
plot(m116)
# heterogeneity issues

summary(m117)
confint(m117) # time and date * migstatus are significant
plot(m117) # heterogeneity issues

# ---------------------------------------------------------------------------- #

# No effect of detection or ag intensity on uric acid levels

# Model Assumptions ####
hist(birds$Uric) # normally distributed data, no errors

# Set up the 2x2 plot layout
op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))

# Plot residuals vs fitted values (CLEAR VIOLATION)
fitted_vals <- fitted(m61)  # Get fitted values from the model
residuals_vals <- resid(m61)  # Get residuals from the model
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")  # Add a horizontal line at zero for reference

# Plot histogram of residuals (satisfied)
hist(residuals_vals, xlab = "Residuals", main = "")

# Plot residuals vs MigStatus (no violation)
boxplot(residuals_vals ~ birds.cs$MigStatus, xlab = "MigStatus", ylab = "Residuals")

# Plot residuals vs Capture Time (no violation)
plot(birds.cs$ts.sunrise, residuals_vals, xlab = "Capture Time", ylab = "Residuals")

# Plot residuals vs Date (no violation)
plot(birds.cs$DaysIntoSeason_S, residuals_vals, xlab = "Date", ylab = "Residuals")

# Restore original plot settings
par(op)

# VIOLATION IS PROBABLY BECAUSE A VARIABLE HAS A NONLINEAR RELATIONSHIP?

# ---------------------------------------------------------------------------- #
# Uric ~ MigStatus ####
m <- lm(Uric ~ ts.sunrise + DaysIntoSeason_S * MigStatus, data = birds)

d <- expand.grid(ts.sunrise = mean(birds$ts.sunrise),
                 DaysIntoSeason_S = mean(birds$DaysIntoSeason_S),
                 Detection = c("Non-detection"),
                 MigStatus = c("Migratory", "Resident")) 

predictions <- predict(m, newdata = d, se.fit = TRUE)

d$predicted_Mass <- predictions$fit

d$lower_CI <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upper_CI <- d$predicted_Mass + 1.96 * predictions$se.fit

ggplot(d, aes(x = MigStatus, y = predicted_Mass)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.1,
                col = "black",
                size = 1) +
  theme_light() +
  labs(x = "Migratory Status", 
       y = "All Species Uric Acid Levels (umol)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none",
        strip.text = element_text(size = 18)) +
  theme(legend.position = "none")



# Uric ~ Time ####
m <- lm(Uric ~ ts.sunrise + DaysIntoSeason_S * MigStatus, data = birds)

d <- expand.grid(ts.sunrise = seq(min(birds$ts.sunrise),
                                  max(birds$ts.sunrise),
                                  length = 1000),
                 DaysIntoSeason_S = mean(birds$DaysIntoSeason_S),
                 Detection = c("Non-detection"),
                 MigStatus = c("Migratory", "Resident")) 

predictions <- predict(m, newdata = d, se.fit = TRUE, re.form = NULL)

d$predicted_Mass <- predictions$fit

d$lower_CI <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upper_CI <- d$predicted_Mass + 1.96 * predictions$se.fit

ggplot(d, aes(x = ts.sunrise, y = predicted_Mass)) + 
  geom_line(aes(color = MigStatus), linewidth = 1, show.legend = FALSE) +  
  theme_light() +
  labs(x = "Time of Capture since Sunrise (min)",
       y = "All Species Uric Acid Levels (umol)",
       color = "Migratory Status") +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI, fill = MigStatus),  
              alpha = 0.25, color = NA, show.legend = FALSE) +
  geom_point(data = birds, aes(x = ts.sunrise, y = Uric, color = MigStatus)) +
  theme(axis.title.x = element_text(size = 19,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 19,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 12, face = "italic"), 
        legend.title = element_text(size = 14, face = "bold", color = "black"),
        legend.position = c(0.875, 0.875)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 0, y = 1900, 
           label = "Sunrise", angle = 90, 
           vjust = -0.5, hjust = 0.5,
           size = 6)
