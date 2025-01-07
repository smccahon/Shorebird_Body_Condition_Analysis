#-----------------------------------#
#  All Species Fat Analysis  #
#          Created 12/11/2024       #          
#         Modified 12/12/2024       #
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
birds <- read.csv("Shorebird_Data_Cleaned_2024-12-9.csv")

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

# Only include birds with fat data
birds <- birds %>% 
  filter(!is.na(Fat))

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

# Data visualization ####

# Fat ~ Agricultural Intensity
ggplot(birds, aes(x = PercentAg, y = Fat)) + 
  geom_point() +
  labs(x = "Surrounding Agricultural Intensity", y = "All Species Fat") +
  theme_classic()

# Fat ~ Detection
ggplot(birds, aes(x = Detection, y = Fat)) + 
  geom_point() +
  labs(x = "Detection", y = "All Species Fat") +
  theme_classic()

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Event and Capture Time ####

m1 <- lmer(Fat ~ Event * ts.sunrise + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(Fat ~ Event + ts.sunrise + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with interaction performs the best.

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Time

m1 <- lmer(Fat ~ ts.sunrise + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(Fat ~ ts.sunrise + I(ts.sunrise^2) + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# quadratic relationship performs the best

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Date * Migratory Status ####

m1 <- lmer(Fat ~ DaysIntoSeason_S + MigStatus + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(Fat ~ DaysIntoSeason_S * MigStatus + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with interaction performs significantly better

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Agricultural Intensity ####

m1 <- lmer(Fat ~ PercentAg + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(Fat ~ PercentAg + I(PercentAg^2) + (1 | Species), data = birds.cs, REML = FALSE)
m3 <- lmer(Fat ~ AgCategory + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# quadratic transformation of agricultural intensity performs best

# ---------------------------------------------------------------------------- #

# INTERACTIONS INCLUDED: date * migrant status ; event * time 
# TRANSFORMATIONS: agricultural intensity and time as a quadratic relationship; 
# ISSUES WITH OVERFITTING SO TRY REMOVE QUADRATIC RELATIONSHIPS

# ---------------------------------------------------------------------------- #

# RANDOM EFFECTS STRUCTURE IS OVERFITTING MODEL; REMOVE RANDOM EFFECTS AS USE SPECIES AS A FIXED EFFECT

m.global <- lm(Fat ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + 
                   MigStatus + PercentAg + I(PercentAg^2) + Detection +
                   MigStatus * DaysIntoSeason_S +
                   Event * ts.sunrise + Species, data = birds.cs)

m.null <- lm(Fat ~ 1, data = birds.cs)

## Single Covariate Models ####

m1 <- lm(Fat ~ Sex, data = birds.cs)
m2 <- lm(Fat ~ Event, data = birds.cs)
m3 <- lm(Fat ~ ts.sunrise + I(ts.sunrise^2), data = birds.cs)
m4 <- lm(Fat ~ DaysIntoSeason_S, data = birds.cs)
m5 <- lm(Fat ~ MigStatus, data = birds.cs)
m6 <- lm(Fat ~ PercentAg + I(PercentAg^2), data = birds.cs)
m7 <- lm(Fat ~ Detection, data = birds.cs)
m8 <- lm(Fat ~ Species, data = birds.cs)

##  Additive Models ####

### Two additive combinations ####

m9 <- lm(Fat ~ Sex + Event, data = birds.cs)
m10 <- lm(Fat ~ Sex + ts.sunrise + I(ts.sunrise^2), data = birds.cs)
m11 <- lm(Fat ~ Sex + DaysIntoSeason_S, data = birds.cs)
m12 <- lm(Fat ~ Sex + MigStatus, data = birds.cs)
m13 <- lm(Fat ~ Sex + PercentAg + I(PercentAg^2), data = birds.cs)
m14 <- lm(Fat ~ Sex + Detection, data = birds.cs)
m15 <- lm(Fat ~ Sex + Species, data = birds.cs)
m16 <- lm(Fat ~ Event + ts.sunrise + I(ts.sunrise^2) + Event*ts.sunrise + I(ts.sunrise^2), data = birds.cs)
m17 <- lm(Fat ~ Event + DaysIntoSeason_S, data = birds.cs)
m18 <- lm(Fat ~ Event + MigStatus, data = birds.cs)
m19 <- lm(Fat ~ Event + PercentAg + I(PercentAg^2), data = birds.cs)
m20 <- lm(Fat ~ Event + Detection, data = birds.cs)
m21 <- lm(Fat ~ Event + Species, data = birds.cs)
m22 <- lm(Fat ~ ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S, data = birds.cs)
m23 <- lm(Fat ~ ts.sunrise + I(ts.sunrise^2) + MigStatus, data = birds.cs)
m24 <- lm(Fat ~ ts.sunrise + I(ts.sunrise^2) + PercentAg + I(PercentAg^2), data = birds.cs)
m25 <- lm(Fat ~ ts.sunrise + I(ts.sunrise^2) + Detection, data = birds.cs)
m26 <- lm(Fat ~ ts.sunrise + I(ts.sunrise^2) + Species, data = birds.cs)
m27 <- lm(Fat ~ DaysIntoSeason_S + MigStatus + DaysIntoSeason_S * MigStatus, data = birds.cs)
m28 <- lm(Fat ~ DaysIntoSeason_S + PercentAg + I(PercentAg^2), data = birds.cs)
m29 <- lm(Fat ~ DaysIntoSeason_S + Detection, data = birds.cs)
m30 <- lm(Fat ~ DaysIntoSeason_S + Species, data = birds.cs)
m31 <- lm(Fat ~ MigStatus + PercentAg + I(PercentAg^2), data = birds.cs)
m32 <- lm(Fat ~ MigStatus + Detection, data = birds.cs)
m33 <- lm(Fat ~ MigStatus + Species, data = birds.cs)
m34 <- lm(Fat ~ PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m35 <- lm(Fat ~ PercentAg + I(PercentAg^2) + Species, data = birds.cs)

### Three Additive Combinations ####

m29 <- lm(Fat ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + Event * ts.sunrise + I(ts.sunrise^2), data = birds.cs)
m30 <- lm(Fat ~ Sex + Event + DaysIntoSeason_S, data = birds.cs)
m31 <- lm(Fat ~ Sex + Event + MigStatus, data = birds.cs)
m32 <- lm(Fat ~ Sex + Event + PercentAg + I(PercentAg^2), data = birds.cs)
m33 <- lm(Fat ~ Sex + Event + Detection, data = birds.cs)
m34 <- lm(Fat ~ Sex + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S, data = birds.cs)
m35 <- lm(Fat ~ Sex + ts.sunrise + I(ts.sunrise^2) + MigStatus, data = birds.cs)
m36 <- lm(Fat ~ Sex + ts.sunrise + I(ts.sunrise^2) + PercentAg + I(PercentAg^2), data = birds.cs)
m37 <- lm(Fat ~ Sex + ts.sunrise + I(ts.sunrise^2) + Detection, data = birds.cs)
m38 <- lm(Fat ~ Sex + DaysIntoSeason_S + MigStatus + DaysIntoSeason_S * MigStatus, data = birds.cs)
m39 <- lm(Fat ~ Sex + DaysIntoSeason_S + PercentAg + I(PercentAg^2), data = birds.cs)
m40 <- lm(Fat ~ Sex + DaysIntoSeason_S + Detection, data = birds.cs)
m41 <- lm(Fat ~ Sex + MigStatus + PercentAg + I(PercentAg^2), data = birds.cs)
m42 <- lm(Fat ~ Sex + MigStatus + Detection, data = birds.cs)
m43 <- lm(Fat ~ Sex + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m44 <- lm(Fat ~ Event * ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S, data = birds.cs)
m45 <- lm(Fat ~ Event * ts.sunrise + I(ts.sunrise^2) + MigStatus, data = birds.cs)
m46 <- lm(Fat ~ Event * ts.sunrise + I(ts.sunrise^2) + PercentAg + I(PercentAg^2), data = birds.cs)
m47 <- lm(Fat ~ Event * ts.sunrise + I(ts.sunrise^2) + Detection, data = birds.cs)
m48 <- lm(Fat ~ Event + DaysIntoSeason_S * MigStatus, data = birds.cs)
m49 <- lm(Fat ~ Event + DaysIntoSeason_S + PercentAg + I(PercentAg^2), data = birds.cs)
m50 <- lm(Fat ~ Event + DaysIntoSeason_S + Detection, data = birds.cs)
m51 <- lm(Fat ~ Event + MigStatus + PercentAg + I(PercentAg^2), data = birds.cs)
m52 <- lm(Fat ~ Event + MigStatus + Detection, data = birds.cs)
m53 <- lm(Fat ~ Event + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m54 <- lm(Fat ~ ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S * MigStatus, data = birds.cs)
m55 <- lm(Fat ~ ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + PercentAg + I(PercentAg^2), data = birds.cs)
m56 <- lm(Fat ~ ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + Detection, data = birds.cs)
m57 <- lm(Fat ~ ts.sunrise + I(ts.sunrise^2) + MigStatus + PercentAg + I(PercentAg^2), data = birds.cs)
m58 <- lm(Fat ~ ts.sunrise + I(ts.sunrise^2) + MigStatus + Detection, data = birds.cs)
m59 <- lm(Fat ~ ts.sunrise + I(ts.sunrise^2) + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m60 <- lm(Fat ~ DaysIntoSeason_S * MigStatus + PercentAg + I(PercentAg^2), data = birds.cs)
m61 <- lm(Fat ~ DaysIntoSeason_S * MigStatus + Detection, data = birds.cs)
m62 <- lm(Fat ~ DaysIntoSeason_S + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m63 <- lm(Fat ~ MigStatus + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)

### Four Additive Combinations ####

m64 <- lm(Fat ~ Sex + Event * ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S, data = birds.cs)
m65 <- lm(Fat ~ Sex + Event * ts.sunrise + I(ts.sunrise^2) + MigStatus, data = birds.cs)
m66 <- lm(Fat ~ Sex + Event * ts.sunrise + I(ts.sunrise^2) + PercentAg + I(PercentAg^2), data = birds.cs)
m67 <- lm(Fat ~ Sex + Event * ts.sunrise + I(ts.sunrise^2) + Detection, data = birds.cs)
m68 <- lm(Fat ~ Sex + Event + DaysIntoSeason_S * MigStatus, data = birds.cs)
m69 <- lm(Fat ~ Sex + Event + DaysIntoSeason_S + PercentAg + I(PercentAg^2), data = birds.cs)
m70 <- lm(Fat ~ Sex + Event + DaysIntoSeason_S + Detection, data = birds.cs)
m71 <- lm(Fat ~ Sex + Event + MigStatus + PercentAg + I(PercentAg^2), data = birds.cs)
m72 <- lm(Fat ~ Sex + Event + MigStatus + Detection + Event * MigStatus, data = birds.cs)
m73 <- lm(Fat ~ Sex + Event + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m74 <- lm(Fat ~ Sex + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S * MigStatus, data = birds.cs)
m75 <- lm(Fat ~ Sex + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + PercentAg + I(PercentAg^2), data = birds.cs)
m76 <- lm(Fat ~ Sex + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + Detection, data = birds.cs)
m77 <- lm(Fat ~ Sex + ts.sunrise + I(ts.sunrise^2) + MigStatus + PercentAg + I(PercentAg^2), data = birds.cs)
m78 <- lm(Fat ~ Sex + ts.sunrise + I(ts.sunrise^2) + MigStatus + Detection, data = birds.cs)
m79 <- lm(Fat ~ Sex + ts.sunrise + I(ts.sunrise^2) + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m80 <- lm(Fat ~ Sex + DaysIntoSeason_S * MigStatus + PercentAg + I(PercentAg^2), data = birds.cs)
m81 <- lm(Fat ~ Sex + DaysIntoSeason_S * MigStatus + Detection, data = birds.cs)
m82 <- lm(Fat ~ Sex + DaysIntoSeason_S + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m83 <- lm(Fat ~ Sex + MigStatus + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m84 <- lm(Fat ~ Event * ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S * MigStatus, data = birds.cs)
m85 <- lm(Fat ~ Event * ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + PercentAg + I(PercentAg^2), data = birds.cs)
m86 <- lm(Fat ~ Event * ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + Detection, data = birds.cs)
m87 <- lm(Fat ~ Event * ts.sunrise + I(ts.sunrise^2) + MigStatus + PercentAg + I(PercentAg^2), data = birds.cs)
m88 <- lm(Fat ~ Event * ts.sunrise + I(ts.sunrise^2) + MigStatus + Detection, data = birds.cs)
m89 <- lm(Fat ~ Event * ts.sunrise + I(ts.sunrise^2) + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m90 <- lm(Fat ~ Event + DaysIntoSeason_S * MigStatus + PercentAg + I(PercentAg^2), data = birds.cs)
m91 <- lm(Fat ~ Event + DaysIntoSeason_S * MigStatus + Detection, data = birds.cs)
m92 <- lm(Fat ~ Event + DaysIntoSeason_S + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m93 <- lm(Fat ~ Event + MigStatus + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m94 <- lm(Fat ~ ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S * MigStatus + PercentAg + I(PercentAg^2), data = birds.cs)
m95 <- lm(Fat ~ ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S * MigStatus + Detection, data = birds.cs)
m96 <- lm(Fat ~ ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m97 <- lm(Fat ~ ts.sunrise + I(ts.sunrise^2) + MigStatus + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m98 <- lm(Fat ~ DaysIntoSeason_S * MigStatus + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)

### Five Additive Combinations ####

m99 <- lm(Fat ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus +
              DaysIntoSeason_S * MigStatus + Event * ts.sunrise + I(ts.sunrise^2), data = birds.cs)
m100 <- lm(Fat ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + PercentAg + I(PercentAg^2) +
               Event * ts.sunrise + I(ts.sunrise^2), data = birds.cs)
m101 <- lm(Fat ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + Detection +
               Event * ts.sunrise + I(ts.sunrise^2), data = birds.cs)
m102 <- lm(Fat ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + MigStatus + PercentAg + I(PercentAg^2) +
               Event * ts.sunrise + I(ts.sunrise^2), data = birds.cs)
m103 <- lm(Fat ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + MigStatus + Detection +
               Event * ts.sunrise + I(ts.sunrise^2), data = birds.cs)
m104 <- lm(Fat ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + PercentAg + I(PercentAg^2) + Detection +
               Event * ts.sunrise + I(ts.sunrise^2), data = birds.cs)
m105 <- lm(Fat ~ Sex + Event + DaysIntoSeason_S + MigStatus + PercentAg + I(PercentAg^2) +
               DaysIntoSeason_S * MigStatus, data = birds.cs)
m106 <- lm(Fat ~ Sex + Event + DaysIntoSeason_S + MigStatus + Detection +
               + DaysIntoSeason_S * MigStatus, data = birds.cs)
m107 <- lm(Fat ~ Sex + Event + DaysIntoSeason_S + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m108 <- lm(Fat ~ Sex + Event + MigStatus + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m109 <- lm(Fat ~ Sex + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus + PercentAg + I(PercentAg^2) +
               DaysIntoSeason_S * MigStatus, data = birds.cs)
m110 <- lm(Fat ~ Sex + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus + Detection +
               DaysIntoSeason_S * MigStatus, data = birds.cs)
m111 <- lm(Fat ~ Sex + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m112 <- lm(Fat ~ Sex + ts.sunrise + I(ts.sunrise^2) + MigStatus + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m113 <- lm(Fat ~ Sex + DaysIntoSeason_S * MigStatus + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m114 <- lm(Fat ~ Event * ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S * MigStatus + PercentAg + I(PercentAg^2), data = birds.cs)
m115 <- lm(Fat ~ Event * ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S * MigStatus + Detection, data = birds.cs)
m116 <- lm(Fat ~ Event * ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m117 <- lm(Fat ~ Event * ts.sunrise + I(ts.sunrise^2) + MigStatus + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m118 <- lm(Fat ~ Event + DaysIntoSeason_S * MigStatus + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m119 <- lm(Fat ~ ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S * MigStatus + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)

### Six Additive Combinations ####

m120 <- lm(Fat ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus + PercentAg + I(PercentAg^2) +
               Event * ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S * MigStatus, data = birds.cs)
m121 <- lm(Fat ~ Sex + Event * ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S * MigStatus + Detection, data = birds.cs)
m122 <- lm(Fat ~ Sex + Event * ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m123 <- lm(Fat ~ Sex + Event * ts.sunrise + I(ts.sunrise^2) + MigStatus + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m124 <- lm(Fat ~ Sex + Event + DaysIntoSeason_S * MigStatus + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m125 <- lm(Fat ~ Sex + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S * MigStatus + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)
m126 <- lm(Fat ~ Event * ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S * MigStatus + PercentAg + I(PercentAg^2) + Detection, data = birds.cs)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:35)

models <- mget(model_names)

models$m.null <- m.null
models$m.global <- m.global

model_names <- c(model_names, "m.null", "m.global")

aictab(models, modnames = model_names)

# ---------------------------------------------------------------------------- #

confint(m106)






confint(m.global)
















