#---------------------------------#
# All Species Body Mass Analysis  #
#         Created 11/11/2024      #          
#        Modified 12/09/2024      #
#---------------------------------#

# load packages
library(dplyr)
library(ggplot2)
library(AICcmodavg)
library(RColorBrewer)
library(trtools)
library(tidyverse)
library(lme4)
library(car)
library(tidyr)
library(viridis)

# ---------------------------------------------------------------------------- #

# Data Processing & Manipulation ####

# Read data
setwd("processed_data")
birds <- read.csv("Shorebird_Data_Cleaned_2024-12-9.csv")

# Make neonicotinoid detection column (Detection/Non-detection)
birds$Detection <- ifelse(birds$OverallNeonic > 0, "Detection", "Non-detection")

# Make sampling event column
birds <- birds %>% 
  mutate(Event = paste(Season, Year, sep = "_")) %>% 
  mutate(Event = ifelse(Season %in% c("Spring", "Fall"),
                        paste(Season, Year, sep = "_"),
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

# Only include species with mass and at least three individuals
birds <- birds %>%
  filter(!is.na(Mass)) %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

table(birds$Species)

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

# Mass ~ Concentration
ggplot(data = birds %>% drop_na(OverallNeonic, Mass), aes(x = OverallNeonic, y = Mass)) + 
  theme_classic() + 
  geom_point() + 
  labs(x = "Neonicotinoid Concentration (ug/L)", 
       y = "All Species Body Mass")

# Mass ~ Log(Concentration)
ggplot(data = birds %>% drop_na(LogNeonic, Mass), aes(x = LogNeonic, y = Mass)) +
  theme_classic() +
  geom_point() + 
  labs(x = "Log(Concentration)",
       y = "All Species Body Mass")

# Mass ~ Detection
ggplot(data = birds %>% drop_na(Detection, Mass), aes(x = Detection, y = Mass)) + 
  theme_classic() + 
  geom_boxplot() +  
  geom_jitter() +   
  labs(x = "Neonicotinoid Detection", y = "All Species Body Mass")

# Mass ~ DaysIntoSeason_S
birds$Event <- factor(birds$Event, levels = c("Fall_2021", "Spring_2022", 
                                                    "Spring_2023","Fall_2023"))
ggplot(data = birds, aes(x = DaysIntoSeason_S, y = Mass, color = Event)) +
  theme_classic() +
  geom_point() +
  labs(x = "Date into Trapping Season",
       y = "All Species Body Mass",
       color = "Sampling Event") +
  scale_color_manual(
    values = c("orange", "purple", "green", "blue"),
    labels = c("Fall 2021", "Spring 2022","Spring 2023", "Fall 2023")
  )

# Mass ~ DaysIntoSeason_MR
birds$Event <- factor(birds$Event, levels = c("Fall_2021", "Spring_2022", 
                                              "Spring_2023","Fall_2023"))
ggplot(data = birds, aes(x = DaysIntoSeason_MR, y = Mass, color = Event)) +
  theme_classic() +
  geom_point() +
  labs(x = "Date into Migration Season",
       y = "All Species Body Mass",
       color = "Sampling Event") +
  scale_color_manual(
    values = c("orange", "purple", "green", "blue"),
    labels = c("Fall 2021", "Spring 2022","Spring 2023", "Fall 2023")
  )

# Mass ~ Julian
birds$Event <- factor(birds$Event, levels = c("Fall_2021", "Spring_2022", 
                                                    "Spring_2023", "Fall_2023"))
ggplot(data = birds, aes(x = Julian, y = Mass, color = Event)) +
  theme_classic() +
  geom_point() +
  labs(x = "Julian Day",
       y = "All Species Body Mass",
       color = "Sampling Event") +
  scale_color_manual(
    values = c("orange", "purple", "green", "blue"),
    labels = c("Fall 2021", "Spring 2022", "Spring 2023", "Fall 2023")
  )

# Mass ~ Sex
birds$Sex <- factor(birds$Sex, levels = c("M", "F"))
ggplot(data = birds %>% drop_na(Sex, Mass), aes(x = Sex, y = Mass)) + 
  theme_classic() + 
  geom_boxplot() + 
  geom_jitter() +
  labs(x = "Sex",
       y = "All Species Body Mass")

# Mass ~ Sampling Event
birds$Event <- factor(birds$Event, levels = c("Fall_2021", "Spring_2022",
                                              "Spring_2023", "Fall_2023"))

ggplot(data = birds, aes(x = Event, y = Mass, color = Event)) +
  theme_classic() +
  geom_boxplot() +
  geom_jitter() + 
  labs(x = "Sampling Event",
       y = "All Species Body Mass",
       color = "Sampling Event") +
  scale_color_manual(
    values = c("orange", "purple", "green", "blue"),
    labels = c("Fall 2021", "Spring 2022", "Spring 2023", "Fall 2023")
  ) + 
  scale_x_discrete(
    labels = c("Fall_2021" = "Fall 2021", "Spring_2022" = "Spring 2022", 
               "Spring_2023" = "Spring 2023", "Fall_2023" = "Fall 2023"))
  

# Mass ~ Capture Time
ggplot(data = birds, aes(x = ts.sunrise, y = Mass)) +
  theme_classic() +
  geom_smooth(method = "lm",  color = "blue") +
  geom_point() +
  labs(x = "Capture Time Since Sunrise (min.)",
       y = "All Species Body Mass (g)")


# ---------------------------------------------------------------------------- #

# All Species Modeling ####

# Identify Best Random Structure: Detection ####

# Random Intercept Models
m1 <- lmer(Mass ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + 
             Detection + (1 | Species), data = birds.cs)

# Random Slope Models
# All other random slope models failed to converge
m2 <- lmer(Mass ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + 
             Detection + (Sex | Species), data = birds.cs)

m3 <- lmer(Mass ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + 
             Detection + (DaysIntoSeason_S | Species), data = birds.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Random slope of Sex | Species is the best 
# However, this random slope has issues of model convergence
# Default to (1 | Species) only. 

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Date ####

m1 <- lmer(Mass ~ DaysIntoSeason_S + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(Mass ~ Julian + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Days into trapping season performs best

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Agricultural Intensity ####

m1 <- lmer(Mass ~ AgCategory + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(Mass ~ PercentAg + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Agricultural intensity as a continuous variable performs significantly better.

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Neonicotinoid ####

m1 <- lmer(Mass ~ OverallNeonic + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(Mass ~ LogNeonic + (1 | Species), data = birds.cs, REML = FALSE)
m3 <- lmer(Mass ~ Detection + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# All three transformations perform very similarly
# Log transformation is necessary to address issues of homoscedasticity and normality
# Detection makes sense biologically
# Just use detection for consistency with other models

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Event and Capture Time ####

m1 <- lmer(Mass ~ Event * ts.sunrise + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(Mass ~ Event + I(ts.sunrise^2) + (1 | Species), data = birds.cs, REML = FALSE)
m3 <- lmer(Mass ~ Event + ts.sunrise + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with interaction by far performs the best.

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Date * Migratory Status ####

m1 <- lmer(Mass ~ DaysIntoSeason_S + MigStatus + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(Mass ~ DaysIntoSeason_S * MigStatus + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with interaction performs significantly

# ---------------------------------------------------------------------------- #

m.global <- lmer(Mass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + 
                     MigStatus + PercentAg + Detection +
                     MigStatus * DaysIntoSeason_S + 
                       Detection * Event +
                       Event * ts.sunrise + (1 | Species), REML = FALSE, data = birds.cs)

m.null <- lmer(Mass ~ 1 + (1|Species), data = birds.cs, REML = FALSE)

## Single Covariate Models ####

m1 <- lmer(Mass ~ Sex + (1|Species), data = birds.cs, REML = FALSE)
m2 <- lmer(Mass ~ Event + (1|Species), data = birds.cs, REML = FALSE)
m3 <- lmer(Mass ~ ts.sunrise + (1|Species), data = birds.cs, REML = FALSE)
m4 <- lmer(Mass ~ DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m5 <- lmer(Mass ~ MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m6 <- lmer(Mass ~ PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m7 <- lmer(Mass ~ Detection + (1|Species), data = birds.cs, REML = FALSE)
m8 <- lmer(Mass ~ LogNeonic + (1|Species), data = birds.cs, REML = FALSE)


##  Additive Models ####

### Two additive combinations ####

m9 <- lmer(Mass ~ Sex + Event + (1|Species), data = birds.cs, REML = FALSE)
m10 <- lmer(Mass ~ Sex + ts.sunrise + (1|Species), data = birds.cs, REML = FALSE)
m11 <- lmer(Mass ~ Sex + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m12 <- lmer(Mass ~ Sex + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m13 <- lmer(Mass ~ Sex + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m14 <- lmer(Mass ~ Sex + Detection + (1|Species), data = birds.cs, REML = FALSE)
m15 <- lmer(Mass ~ Event + ts.sunrise + (1|Species) + Event*ts.sunrise, data = birds.cs, REML = FALSE)
m16 <- lmer(Mass ~ Event + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m17 <- lmer(Mass ~ Event + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m18 <- lmer(Mass ~ Event + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m19 <- lmer(Mass ~ Event + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m20 <- lmer(Mass ~ ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m21 <- lmer(Mass ~ ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m22 <- lmer(Mass ~ ts.sunrise + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m23 <- lmer(Mass ~ ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m24 <- lmer(Mass ~ DaysIntoSeason_S + MigStatus + (1|Species) + DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m25 <- lmer(Mass ~ DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m26 <- lmer(Mass ~ DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m27 <- lmer(Mass ~ MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m28 <- lmer(Mass ~ MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m29 <- lmer(Mass ~ PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m30 <- lmer(Mass ~ Sex + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m31 <- lmer(Mass ~ Event + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m32 <- lmer(Mass ~ ts.sunrise + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m33 <- lmer(Mass ~ DaysIntoSeason_S + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m34 <- lmer(Mass ~ MigStatus + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m35 <- lmer(Mass ~ PercentAg + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)

### Three Additive Combinations ####

m36 <- lmer(Mass ~ Sex + Event + ts.sunrise + (1|Species) + Event * ts.sunrise, data = birds.cs, REML = FALSE)
m37 <- lmer(Mass ~ Sex + Event + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m38 <- lmer(Mass ~ Sex + Event + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m39 <- lmer(Mass ~ Sex + Event + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m40 <- lmer(Mass ~ Sex + Event + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m41 <- lmer(Mass ~ Sex + ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m42 <- lmer(Mass ~ Sex + ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m43 <- lmer(Mass ~ Sex + ts.sunrise + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m44 <- lmer(Mass ~ Sex + ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m45 <- lmer(Mass ~ Sex + DaysIntoSeason_S + MigStatus + (1|Species) + DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m46 <- lmer(Mass ~ Sex + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m47 <- lmer(Mass ~ Sex + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m48 <- lmer(Mass ~ Sex + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m49 <- lmer(Mass ~ Sex + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m50 <- lmer(Mass ~ Sex + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m51 <- lmer(Mass ~ Event * ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m52 <- lmer(Mass ~ Event * ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m53 <- lmer(Mass ~ Event * ts.sunrise + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m54 <- lmer(Mass ~ Event * ts.sunrise + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m55 <- lmer(Mass ~ Event + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m56 <- lmer(Mass ~ Event + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m57 <- lmer(Mass ~ Event + DaysIntoSeason_S + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m58 <- lmer(Mass ~ Event + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m59 <- lmer(Mass ~ Event + MigStatus + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m60 <- lmer(Mass ~ Event + PercentAg + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m61 <- lmer(Mass ~ ts.sunrise + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m62 <- lmer(Mass ~ ts.sunrise + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m63 <- lmer(Mass ~ ts.sunrise + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m64 <- lmer(Mass ~ ts.sunrise + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m65 <- lmer(Mass ~ ts.sunrise + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m66 <- lmer(Mass ~ ts.sunrise + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m67 <- lmer(Mass ~ DaysIntoSeason_S * MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m68 <- lmer(Mass ~ DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m69 <- lmer(Mass ~ DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m70 <- lmer(Mass ~ MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m71 <- lmer(Mass ~ Sex + Event + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m72 <- lmer(Mass ~ Sex + ts.sunrise + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m73 <- lmer(Mass ~ Sex + DaysIntoSeason_S + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m74 <- lmer(Mass ~ Sex + MigStatus + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m75 <- lmer(Mass ~ Sex + PercentAg + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m76 <- lmer(Mass ~ Event * ts.sunrise + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m77 <- lmer(Mass ~ Event + DaysIntoSeason_S + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m78 <- lmer(Mass ~ Event + MigStatus + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m79 <- lmer(Mass ~ Event + PercentAg + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m80 <- lmer(Mass ~ ts.sunrise + DaysIntoSeason_S + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m81 <- lmer(Mass ~ ts.sunrise + MigStatus + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m82 <- lmer(Mass ~ ts.sunrise + PercentAg + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m83 <- lmer(Mass ~ DaysIntoSeason_S * MigStatus + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m84 <- lmer(Mass ~ DaysIntoSeason_S + PercentAg + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m85 <- lmer(Mass ~ MigStatus + PercentAg + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)


### Four Additive Combinations ####

m86 <- lmer(Mass ~ Sex + Event * ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m87 <- lmer(Mass ~ Sex + Event * ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m88 <- lmer(Mass ~ Sex + Event * ts.sunrise + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m89 <- lmer(Mass ~ Sex + Event * ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m90 <- lmer(Mass ~ Sex + Event + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m91 <- lmer(Mass ~ Sex + Event + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m92 <- lmer(Mass ~ Sex + Event + DaysIntoSeason_S + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m93 <- lmer(Mass ~ Sex + Event + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m94 <- lmer(Mass ~ Sex + Event + MigStatus + Detection + (1|Species) + Event * MigStatus, data = birds.cs, REML = FALSE)
m95 <- lmer(Mass ~ Sex + Event + PercentAg + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m96 <- lmer(Mass ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m97 <- lmer(Mass ~ Sex + ts.sunrise + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m98 <- lmer(Mass ~ Sex + ts.sunrise + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m99 <- lmer(Mass ~ Sex + ts.sunrise + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m100 <- lmer(Mass ~ Sex + ts.sunrise + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m101 <- lmer(Mass ~ Sex + ts.sunrise + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m102 <- lmer(Mass ~ Sex + DaysIntoSeason_S * MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m103 <- lmer(Mass ~ Sex + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m104 <- lmer(Mass ~ Sex + DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m105 <- lmer(Mass ~ Sex + MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m106 <- lmer(Mass ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m107 <- lmer(Mass ~ Event * ts.sunrise + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m108 <- lmer(Mass ~ Event * ts.sunrise + DaysIntoSeason_S + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m109 <- lmer(Mass ~ Event * ts.sunrise + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m110 <- lmer(Mass ~ Event * ts.sunrise + MigStatus + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m111 <- lmer(Mass ~ Event * ts.sunrise + PercentAg + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m112 <- lmer(Mass ~ Event + DaysIntoSeason_S * MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m113 <- lmer(Mass ~ Event + DaysIntoSeason_S * MigStatus + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m114 <- lmer(Mass ~ Event + DaysIntoSeason_S + PercentAg + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m115 <- lmer(Mass ~ Event + MigStatus + PercentAg + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m116 <- lmer(Mass ~ ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m117 <- lmer(Mass ~ ts.sunrise + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m118 <- lmer(Mass ~ ts.sunrise + DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m119 <- lmer(Mass ~ ts.sunrise + MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m120 <- lmer(Mass ~ DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m121 <- lmer(Mass ~ Sex + Event * ts.sunrise + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m122 <- lmer(Mass ~ Sex + Event + DaysIntoSeason_S + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m123 <- lmer(Mass ~ Sex + Event + MigStatus + LogNeonic + (1|Species) + Event * MigStatus, data = birds.cs, REML = FALSE)
m124 <- lmer(Mass ~ Sex + Event + PercentAg + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m125 <- lmer(Mass ~ Sex + ts.sunrise + DaysIntoSeason_S + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m126 <- lmer(Mass ~ Sex + ts.sunrise + MigStatus + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m127 <- lmer(Mass ~ Sex + ts.sunrise + PercentAg + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m128 <- lmer(Mass ~ Sex + DaysIntoSeason_S * MigStatus + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m129 <- lmer(Mass ~ Sex + DaysIntoSeason_S + PercentAg + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m130 <- lmer(Mass ~ Sex + MigStatus + PercentAg + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m131 <- lmer(Mass ~ Event * ts.sunrise + DaysIntoSeason_S + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m132 <- lmer(Mass ~ Event * ts.sunrise + MigStatus + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m133 <- lmer(Mass ~ Event * ts.sunrise + PercentAg + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m134 <- lmer(Mass ~ Event + DaysIntoSeason_S * MigStatus + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m135 <- lmer(Mass ~ Event + DaysIntoSeason_S + PercentAg + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m136 <- lmer(Mass ~ Event + MigStatus + PercentAg + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m137 <- lmer(Mass ~ ts.sunrise + DaysIntoSeason_S * MigStatus + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m138 <- lmer(Mass ~ ts.sunrise + DaysIntoSeason_S + PercentAg + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m139 <- lmer(Mass ~ ts.sunrise + MigStatus + PercentAg + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m140 <- lmer(Mass ~ DaysIntoSeason_S * MigStatus + PercentAg + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)

### Five Additive Combinations ####

m141 <- lmer(Mass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + MigStatus + (1|Species) +
               DaysIntoSeason_S * MigStatus + Event * ts.sunrise, data = birds.cs, REML = FALSE)
m142 <- lmer(Mass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + PercentAg + (1|Species) +
               Event * ts.sunrise, data = birds.cs, REML = FALSE)
m143 <- lmer(Mass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + Detection + (1|Species) +
               Event * ts.sunrise + Event * Detection, data = birds.cs, REML = FALSE)
m144 <- lmer(Mass ~ Sex + Event + ts.sunrise + MigStatus + PercentAg + (1|Species) +
               Event * ts.sunrise, data = birds.cs, REML = FALSE)
m145 <- lmer(Mass ~ Sex + Event + ts.sunrise + MigStatus + Detection + (1|Species) +
               Event * ts.sunrise + Event * Detection, data = birds.cs, REML = FALSE)
m146 <- lmer(Mass ~ Sex + Event + ts.sunrise + PercentAg + Detection + (1|Species) +
               Event * ts.sunrise + Event * Detection, data = birds.cs, REML = FALSE)
m147 <- lmer(Mass ~ Sex + Event + DaysIntoSeason_S + MigStatus + PercentAg + (1|Species) +
               DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m148 <- lmer(Mass ~ Sex + Event + DaysIntoSeason_S + MigStatus + Detection + (1|Species) +
               Event * Detection + DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m149 <- lmer(Mass ~ Sex + Event + DaysIntoSeason_S + PercentAg + Detection + (1|Species) +
               Event * Detection, data = birds.cs, REML = FALSE)
m150 <- lmer(Mass ~ Sex + Event + MigStatus + PercentAg + Detection + (1|Species) +
               Event * Detection, data = birds.cs, REML = FALSE)
m151 <- lmer(Mass ~ Sex + ts.sunrise + DaysIntoSeason_S + MigStatus + PercentAg + (1|Species) +
               DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m152 <- lmer(Mass ~ Sex + ts.sunrise + DaysIntoSeason_S + MigStatus + Detection + (1|Species) +
               DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m153 <- lmer(Mass ~ Sex + ts.sunrise + DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m154 <- lmer(Mass ~ Sex + ts.sunrise + MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m155 <- lmer(Mass ~ Sex + DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m156 <- lmer(Mass ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m157 <- lmer(Mass ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m158 <- lmer(Mass ~ Event * ts.sunrise + DaysIntoSeason_S + PercentAg + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m159 <- lmer(Mass ~ Event * ts.sunrise + MigStatus + PercentAg + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m160 <- lmer(Mass ~ Event + DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m161 <- lmer(Mass ~ ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m162 <- lmer(Mass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + LogNeonic + (1|Species) +
               Event * ts.sunrise + Event * LogNeonic, data = birds.cs, REML = FALSE)
m163 <- lmer(Mass ~ Sex + Event + ts.sunrise + MigStatus + LogNeonic + (1|Species) +
               Event * ts.sunrise + Event * LogNeonic, data = birds.cs, REML = FALSE)
m164 <- lmer(Mass ~ Sex + Event + ts.sunrise + PercentAg + LogNeonic + (1|Species) +
               Event * ts.sunrise + Event * LogNeonic, data = birds.cs, REML = FALSE)
m165 <- lmer(Mass ~ Sex + Event + DaysIntoSeason_S + MigStatus + LogNeonic + (1|Species) +
               Event * LogNeonic + DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m166 <- lmer(Mass ~ Sex + Event + DaysIntoSeason_S + PercentAg + LogNeonic + (1|Species) +
               Event * LogNeonic, data = birds.cs, REML = FALSE)
m167 <- lmer(Mass ~ Sex + Event + MigStatus + PercentAg + LogNeonic + (1|Species) +
               Event * LogNeonic, data = birds.cs, REML = FALSE)
m168 <- lmer(Mass ~ Sex + ts.sunrise + DaysIntoSeason_S + MigStatus + LogNeonic + (1|Species) +
               DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m169 <- lmer(Mass ~ Sex + ts.sunrise + DaysIntoSeason_S + PercentAg + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m170 <- lmer(Mass ~ Sex + ts.sunrise + MigStatus + PercentAg + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m171 <- lmer(Mass ~ Sex + DaysIntoSeason_S * MigStatus + PercentAg + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m172 <- lmer(Mass ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m173 <- lmer(Mass ~ Event * ts.sunrise + DaysIntoSeason_S + PercentAg + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m174 <- lmer(Mass ~ Event * ts.sunrise + MigStatus + PercentAg + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m175 <- lmer(Mass ~ Event + DaysIntoSeason_S * MigStatus + PercentAg + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m176 <- lmer(Mass ~ ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)


### Six Additive Combinations ####

m177 <- lmer(Mass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + MigStatus + PercentAg + (1|Species) +
               Event * ts.sunrise + DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m178 <- lmer(Mass ~ Sex + Event * ts.sunrise + DaysIntoSeason_S * MigStatus + Detection + (1|Species) +
               Event * Detection, data = birds.cs, REML = FALSE)
m179 <- lmer(Mass ~ Sex + Event * ts.sunrise + DaysIntoSeason_S + PercentAg + Detection + (1|Species) +
               Event * Detection, data = birds.cs, REML = FALSE)
m180 <- lmer(Mass ~ Sex + Event * ts.sunrise + MigStatus + PercentAg + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m181 <- lmer(Mass ~ Sex + Event + DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species) + Event * Detection, data = birds.cs, REML = FALSE)
m182 <- lmer(Mass ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m183 <- lmer(Mass ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m184 <- lmer(Mass ~ Sex + Event * ts.sunrise + DaysIntoSeason_S * MigStatus + LogNeonic + (1|Species) +
               Event * LogNeonic, data = birds.cs, REML = FALSE)
m185 <- lmer(Mass ~ Sex + Event * ts.sunrise + DaysIntoSeason_S + PercentAg + LogNeonic + (1|Species) +
               Event * LogNeonic, data = birds.cs, REML = FALSE)
m186 <- lmer(Mass ~ Sex + Event * ts.sunrise + MigStatus + PercentAg + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m187 <- lmer(Mass ~ Sex + Event + DaysIntoSeason_S * MigStatus + PercentAg + LogNeonic + (1|Species) + Event * LogNeonic, data = birds.cs, REML = FALSE)
m188 <- lmer(Mass ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)
m189 <- lmer(Mass ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + LogNeonic + (1|Species), data = birds.cs, REML = FALSE)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:189)

models <- mget(model_names)

models$m.null <- m.null
models$m.global.d <- m.global.d
models$m.global.l <- m.global.l

model_names <- c(model_names, "m.null", "m.global.d", "m.global.l")

aictab(models, modnames = model_names)

# ---------------------------------------------------------------------------- #

# Top Model Summaries ####

m145 <- lmer(Mass ~ Sex + Event + ts.sunrise + MigStatus + Detection + (1|Species) +
               Event * ts.sunrise + Event * Detection, data = birds.cs, REML = FALSE)

summary(m145)
confint(m145)

m178 <- lmer(Mass ~ Sex + Event * ts.sunrise + DaysIntoSeason_S * MigStatus + Detection + (1|Species) +
               Event * Detection, data = birds.cs, REML = FALSE)

summary(m178)
confint(m178)

# Should I include interactions?? Does it make a difference in my output?

# not significant, slightly positive
m <- lmer(Mass ~ Event * ts.sunrise + Detection + (1 | Species), data = birds.cs, REML = FALSE)

# not significant, slightly positive
m <- lmer(Mass ~ Event + ts.sunrise + Detection + (1 | Species), data = birds.cs, REML = FALSE)

# not significant, slightly positive
m <- lmer(Mass ~ MigStatus * DaysIntoSeason_S + Detection + (1 | Species), data = birds.cs, REML = FALSE)

# not significant, slightly positive
m <- lmer(Mass ~ MigStatus + DaysIntoSeason_S + Detection + (1 | Species), data = birds.cs, REML = FALSE)

# Probably should exclude these!

# ---------------------------------------------------------------------------- #

# Graph top models ####

m <- lmer(Mass ~ Event * ts.sunrise + (1 | Species), data = birds, REML = FALSE)

## All species ####
d <- expand.grid(Event = c("Fall_2021", "Spring_2022", "Spring_2023", "Fall_2023"),
                 ts.sunrise = seq(min(birds$ts.sunrise), max(birds$ts.sunrise), 
                                  length = 1000),
                 Species = unique(birds$Species))

d$Mass <- predict(m, newdata = d, re.form = NULL)

d$Species <- factor(d$Species, 
                    levels = c("LesserYellowlegs", "Killdeer", "Willet", 
                               "PectoralSandpiper", "WIPH", "LeastSandpiper", 
                               "AmericanAvocet", "LongbilledDowitcher", 
                               "SemipalmatedSandpiper"),
                    labels = c("Lesser Yellowlegs", "Killdeer", "Willet", 
                               "Pectoral Sandpiper", "Wilson's Phalarope", 
                               "Least Sandpiper", "American Avocet", 
                               "Long-billed Dowitcher", "Semipalmated Sandpiper"))

d$Event <- factor(d$Event, 
                  levels = c("Fall_2021", "Spring_2022", "Spring_2023", "Fall_2023"),
                  labels = c("Fall 2021", "Spring 2022", "Spring 2023", "Fall 2023"))

birds$Species <- factor(birds$Species, 
                    levels = c("LesserYellowlegs", "Killdeer", "Willet", 
                               "PectoralSandpiper", "WIPH", "LeastSandpiper", 
                               "AmericanAvocet", "LongbilledDowitcher", 
                               "SemipalmatedSandpiper"),
                    labels = c("Lesser Yellowlegs", "Killdeer", "Willet", 
                               "Pectoral Sandpiper", "Wilson's Phalarope", 
                               "Least Sandpiper", "American Avocet", 
                               "Long-billed Dowitcher", "Semipalmated Sandpiper"))

birds$Event <- factor(birds$Event, 
                  levels = c("Fall_2021", "Spring_2022", "Spring_2023", "Fall_2023"),
                  labels = c("Fall 2021", "Spring 2022", "Spring 2023", "Fall 2023"))



ggplot(d, aes(x = ts.sunrise, y = Mass, color = Event)) +
  geom_line(size = 0.9) +
  theme_light() +
  labs(x = "Capture Time (min)", 
       y = "Body Mass (g)", 
       color = "Event") +
  facet_wrap(~ Species, scales = "free_y") +
  theme(legend.position = "right") +
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10)) +
  scale_color_viridis_d(begin = 0.2, 
                        end = 0.8, 
                        option = "inferno") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
    annotate("text", x = 0, y = Inf, size = 3.5,
             label = "Sunrise", angle = 90, vjust = -0.5, hjust = 1) +
  geom_point(data = birds, aes(x = ts.sunrise, y = Mass))
  
  
## Lesser Yellowlegs ####
d_leye <- d[d$Species == "Lesser Yellowlegs",]

ggplot(d_leye, aes(x = ts.sunrise, y = Mass, color = Event)) +
  geom_line(size = 0.9) +
  theme_light() +
  labs(x = "Capture Time (min)", 
       y = "Lesser Yellowlegs Body Mass (g)", 
       color = "Event") +
  theme(legend.position = "right") +
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10)) +
  scale_color_viridis_d(begin = 0.2, 
                        end = 0.8, 
                        option = "inferno") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 0, y = max(d$Mass), 
           label = "Sunrise", angle = 90, vjust = -0.5, hjust = 0.5) +
  ylim(min(d_leye$Mass), max(d_leye$Mass))


## Mass ~ Event + Detection + (1 | Species)
m <- lmer(Mass ~ Event + Detection + (1 | Species), data = birds, REML = FALSE)

birds$Species <- factor(birds$Species, 
                    levels = c("LesserYellowlegs", "Killdeer", "Willet", 
                               "PectoralSandpiper", "WIPH", "LeastSandpiper", 
                               "AmericanAvocet", "LongbilledDowitcher", 
                               "SemipalmatedSandpiper"),
                    labels = c("Lesser Yellowlegs", "Killdeer", "Willet", 
                               "Pectoral Sandpiper", "Wilson's Phalarope", 
                               "Least Sandpiper", "American Avocet", 
                               "Long-billed Dowitcher", "Semipalmated Sandpiper"))

birds$Event <- factor(birds$Event, 
                  levels = c("Fall_2021", "Spring_2022", "Spring_2023", "Fall_2023"),
                  labels = c("Fall 2021", "Spring 2022", "Spring 2023", "Fall 2023"))


ggplot(birds %>% filter(!is.na(Detection)), aes(x = factor(Detection), y = Mass, 
                                                color = Event)) +
  geom_boxplot() +  
  geom_jitter(aes(color = Event),
              outlier.size = 2,
              position = position_jitterdodge(jitter.width = 0.2),  
              size = 1, 
              alpha = 0.5) + 
  theme_light() +
  labs(x = "Neonicotinoid Detection", 
       y = "Body Mass (g)", 
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


## Mass ~ Event + DaysIntoSeason_S + (1 | Species) ####
m <- lmer(Mass ~ Event + DaysIntoSeason_S + (1 | Species), data = birds, REML= FALSE)

d <- expand.grid(Event = c("Fall_2021", "Spring_2022", "Spring_2023", "Fall_2023"),
                 DaysIntoSeason_S = seq(min(birds$DaysIntoSeason_S), 
                                        max(birds$DaysIntoSeason_S), 
                                  length = 1000),
                 Species = unique(birds$Species))

d$Mass <- predict(m, newdata = d, re.form = NULL)

d$Species <- factor(d$Species, 
                    levels = c("LesserYellowlegs", "Killdeer", "Willet", 
                               "PectoralSandpiper", "WIPH", "LeastSandpiper", 
                               "AmericanAvocet", "LongbilledDowitcher", 
                               "SemipalmatedSandpiper"),
                    labels = c("Lesser Yellowlegs", "Killdeer", "Willet", 
                               "Pectoral Sandpiper", "Wilson's Phalarope", 
                               "Least Sandpiper", "American Avocet", 
                               "Long-billed Dowitcher", "Semipalmated Sandpiper"))

d$Event <- factor(d$Event, 
                  levels = c("Fall_2021", "Spring_2022", "Spring_2023", "Fall_2023"),
                  labels = c("Fall 2021", "Spring 2022", "Spring_2023", "Fall 2023"))

birds$Species <- factor(birds$Species, 
                    levels = c("LesserYellowlegs", "Killdeer", "Willet", 
                               "PectoralSandpiper", "WIPH", "LeastSandpiper", 
                               "AmericanAvocet", "LongbilledDowitcher", 
                               "SemipalmatedSandpiper"),
                    labels = c("Lesser Yellowlegs", "Killdeer", "Willet", 
                               "Pectoral Sandpiper", "Wilson's Phalarope", 
                               "Least Sandpiper", "American Avocet", 
                               "Long-billed Dowitcher", "Semipalmated Sandpiper"))

birds$Event <- factor(birds$Event, 
                  levels = c("Fall_2021", "Spring_2022", "Spring_2023", "Fall_2023"),
                  labels = c("Fall 2021", "Spring 2022", "Spring_2023", "Fall 2023"))


ggplot(d, aes(x = DaysIntoSeason_S, y = Mass, color = Event)) +
  geom_line(size = 0.9) +
  theme_light() +
  labs(x = "Days Into Trapping Season", 
       y = "Body Mass (g)", 
       color = "Event") +
  facet_wrap(~ Species, scales = "free_y") +
  theme(legend.position = "right") +
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10)) +
  scale_color_viridis_d(begin = 0.2, 
                        end = 0.8, 
                        option = "inferno") +
  geom_point(data = birds, aes(x = DaysIntoSeason_S, y = Mass))


## Mass ~ LogNeonic + (1 | Species) ####
# Model not significant

m <- lmer(Mass ~ LogNeonic + Event + (1 | Species), data = birds, REML = FALSE)

birds <- birds %>% filter(!is.na(LogNeonic))

d <- expand.grid(Event = c("Fall_2021", "Spring_2022", "Spring_2023", "Fall_2023"),
                 LogNeonic = seq(min(birds$LogNeonic), 
                                        max(birds$LogNeonic), 
                                        length = 1000),
                 Species = unique(birds$Species))

d$Mass <- predict(m, newdata = d, re.form = NULL)

d$Species <- factor(d$Species, 
                    levels = c("LesserYellowlegs", "Killdeer", "Willet", 
                               "PectoralSandpiper", "WIPH", "LeastSandpiper", 
                               "AmericanAvocet", "LongbilledDowitcher", 
                               "SemipalmatedSandpiper"),
                    labels = c("Lesser Yellowlegs", "Killdeer", "Willet", 
                               "Pectoral Sandpiper", "Wilson's Phalarope", 
                               "Least Sandpiper", "American Avocet", 
                               "Long-billed Dowitcher", "Semipalmated Sandpiper"))

d$Event <- factor(d$Event, 
                  levels = c("Fall_2021", "Spring_2022", "Spring_2023", "Fall_2023"),
                  labels = c("Fall 2021", "Spring 2022", "Spring_2023", "Fall 2023"))

birds$Species <- factor(birds$Species, 
                        levels = c("LesserYellowlegs", "Killdeer", "Willet", 
                                   "PectoralSandpiper", "WIPH", "LeastSandpiper", 
                                   "AmericanAvocet", "LongbilledDowitcher", 
                                   "SemipalmatedSandpiper"),
                        labels = c("Lesser Yellowlegs", "Killdeer", "Willet", 
                                   "Pectoral Sandpiper", "Wilson's Phalarope", 
                                   "Least Sandpiper", "American Avocet", 
                                   "Long-billed Dowitcher", "Semipalmated Sandpiper"))

birds$Event <- factor(birds$Event, 
                      levels = c("Fall_2021", "Spring_2022", "Spring_2023", "Fall_2023"),
                      labels = c("Fall 2021", "Spring 2022", "Spring_2023", "Fall 2023"))


ggplot(d, aes(x = LogNeonic, y = Mass, color = Event)) +
  geom_line(size = 0.9) +
  theme_light() +
  labs(x = "Log(Neonicotinoid Concentration [ug/L])", 
       y = "Body Mass (g)", 
       color = "Event") +
  facet_wrap(~ Species, scales = "free_y") +
  theme(legend.position = "right") +
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10)) +
  scale_color_viridis_d(begin = 0.2, 
                        end = 0.8, 
                        option = "inferno") +
  geom_point(data = birds, aes(x = LogNeonic, y = Mass))








