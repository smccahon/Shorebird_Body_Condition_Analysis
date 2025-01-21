#----------------------------------#
#  Lesser Yellowlegs Fat Analysis  #
#          Created 12/11/2024      #          
#         Modified 01/20/2025      #
#----------------------------------#

# MODIFIED TO USE FAT AS A NUMERIC VARIABLE

# load packages
library(ggplot2)
library(AICcmodavg)
library(tidyverse)

# ---------------------------------------------------------------------------- #

# Data Processing & Manipulation ####

# Read data
setwd("processed_data")
birds <- read.csv("Shorebird_Data_Cleaned_2025-01-20.csv")

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

# Standardize continuous variables
leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

# ---------------------------------------------------------------------------- #

# Data Visualization ####

# Fat ~ Detection
ggplot(data = leye, aes(x = Detection, y = Fat)) + 
  theme_classic() + 
  geom_boxplot() + 
  labs(x = "Neonicotinoid Detection",
       y = "Lesser Yellowlegs Fat Score")

ggplot(leye, 
       aes(x = Detection, y = Fat)) +
  geom_boxplot() +  
  theme_classic() +
  labs(x = "Neonicotinoid Detection", 
       y = "Lesser Yellowlegs Fat Score") +  
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

wilcox.test(Fat ~ Detection, data = leye)

# Fat ~ PercentAg
ggplot(data = leye, aes(x = PercentAg, y = Fat)) + 
  theme_classic() + 
  geom_point() + 
  labs(x = "Surrounding Agricultural Intensity (%)",
       y = "Lesser Yellowlegs Fat Score")

ggplot(leye, 
       aes(x = PercentAg, y = Fat)) +
  geom_point(size = 2) +  
  theme_classic() +
  labs(x = "Surrounding Agricultural Intensity", 
       y = "Lesser Yellowlegs Fat Score") +  
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = scales::percent)

# ---------------------------------------------------------------------------- #

# Lesser Yellowlegs Modeling ####

## Univariate Analysis: Date ####

m1 <- lm(Mass ~ DaysIntoSeason_S, data = leye.cs)
m2 <- lm(Mass ~ Julian, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Date into trapping season performs significantly better than Julian date

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Agricultural Intensity ####

m1 <- lm(Mass ~ PercentAg, data = leye.cs)
m2 <- lm(Mass ~ LogAg, data = leye.cs)
m3 <- lm(Mass ~ PercentAg + I(PercentAg^2), data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Percent Ag as a linear, continuous variable performs the best (no transformation)

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Event and Capture Time ####

m1 <- lm(Mass ~ Event * ts.sunrise, data = leye.cs)
m2 <- lm(Mass ~ Event + I(ts.sunrise^2), data = leye.cs)
m3 <- lm(Mass ~ Event + ts.sunrise, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with interaction by far performs the best.

# ---------------------------------------------------------------------------- #

# INTERACTIONS INCLUDED: Event * Time, PercentAg

# ---------------------------------------------------------------------------- #

## Single Models ####

m.global <- lm(Fat ~ Sex + Event*ts.sunrise + DaysIntoSeason_S + Detection + PercentAg, 
               data = leye.cs)


m.null <- lm(Fat ~ 1, data = leye.cs)

m1 <- lm(Fat ~ Sex, data = leye.cs)
m2 <- lm(Fat ~ Event, data = leye.cs)
m3 <- lm(Fat ~ DaysIntoSeason_S, data = leye.cs)
m4 <- lm(Fat ~ ts.sunrise, data = leye.cs)
m5 <- lm(Fat ~ PercentAg, data = leye.cs)
m6 <- lm(Fat ~ Detection, data = leye.cs)

## Additive Models ####

### Two additive combinations ####
m7 <- lm(Fat ~ Sex + Event, data = leye.cs)
m8 <- lm(Fat ~ Sex + DaysIntoSeason_S, data = leye.cs)
m9 <- lm(Fat ~ Sex + ts.sunrise, data = leye.cs)
m10 <- lm(Fat ~ Sex + PercentAg, data = leye.cs)
m11 <- lm(Fat ~ Sex + Detection, data = leye.cs)
m12 <- lm(Fat ~ Event + DaysIntoSeason_S, data = leye.cs)
m13 <- lm(Fat ~ Event * ts.sunrise, data = leye.cs)
m14 <- lm(Fat ~ Event + PercentAg, data = leye.cs)
m15 <- lm(Fat ~ Event + Detection, data = leye.cs)
m16 <- lm(Fat ~ DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m17 <- lm(Fat ~ DaysIntoSeason_S + PercentAg, data = leye.cs)
m18 <- lm(Fat ~ DaysIntoSeason_S + Detection, data = leye.cs)
m19 <- lm(Fat ~ ts.sunrise + PercentAg, data = leye.cs)
m20 <- lm(Fat ~ ts.sunrise + Detection, data = leye.cs)
m21 <- lm(Fat ~ PercentAg + Detection, data = leye.cs)

### Three additive combinations ####
m22 <- lm(Fat ~ Sex + Event + DaysIntoSeason_S, data = leye.cs)
m23 <- lm(Fat ~ Sex + Event * ts.sunrise, data = leye.cs)
m24 <- lm(Fat ~ Sex + Event + PercentAg, data = leye.cs)
m25 <- lm(Fat ~ Sex + Event + Detection, data = leye.cs)
m26 <- lm(Fat ~ Sex + DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m27 <- lm(Fat ~ Sex + DaysIntoSeason_S + PercentAg, data = leye.cs)
m28 <- lm(Fat ~ Sex + DaysIntoSeason_S + Detection, data = leye.cs)
m29 <- lm(Fat ~ Sex + ts.sunrise + PercentAg, data = leye.cs)
m30 <- lm(Fat ~ Sex + ts.sunrise + Detection, data = leye.cs)
m31 <- lm(Fat ~ Sex + PercentAg + Detection, data = leye.cs)
m32 <- lm(Fat ~ Event + DaysIntoSeason_S + ts.sunrise + Event * ts.sunrise, data = leye.cs)
m33 <- lm(Fat ~ Event + DaysIntoSeason_S + PercentAg, data = leye.cs)
m34 <- lm(Fat ~ Event + DaysIntoSeason_S + Detection, data = leye.cs)
m35 <- lm(Fat ~ Event * ts.sunrise + PercentAg, data = leye.cs)
m36 <- lm(Fat ~ Event * ts.sunrise + Detection, data = leye.cs)
m37 <- lm(Fat ~ Event + PercentAg + Detection, data = leye.cs)
m38 <- lm(Fat ~ DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m39 <- lm(Fat ~ DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m40 <- lm(Fat ~ DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m41 <- lm(Fat ~ ts.sunrise + PercentAg + Detection, data = leye.cs)

### Four additive combinations ####
m42 <- lm(Fat ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + Event*ts.sunrise, data = leye.cs)
m43 <- lm(Fat ~ Sex + Event + DaysIntoSeason_S + PercentAg, data = leye.cs)
m44 <- lm(Fat ~ Sex + Event + DaysIntoSeason_S + Detection, data = leye.cs)
m45 <- lm(Fat ~ Sex + Event * ts.sunrise + PercentAg, data = leye.cs)
m46 <- lm(Fat ~ Sex + Event * ts.sunrise + Detection, data = leye.cs)
m47 <- lm(Fat ~ Sex + Event + PercentAg + Detection, data = leye.cs)
m48 <- lm(Fat ~ Sex + DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m49 <- lm(Fat ~ Sex + DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m50 <- lm(Fat ~ Sex + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m51 <- lm(Fat ~ Sex + ts.sunrise + PercentAg + Detection, data = leye.cs)
m52 <- lm(Fat ~ Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Event*ts.sunrise, data = leye.cs)
m53 <- lm(Fat ~ Event + DaysIntoSeason_S + ts.sunrise + Detection + Event*ts.sunrise, data = leye.cs)
m54 <- lm(Fat ~ Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m55 <- lm(Fat ~ Event * ts.sunrise + PercentAg + Detection, data = leye.cs)
m56 <- lm(Fat ~ DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)

### Five additive combinations ####

m57 <- lm(Fat ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Event * ts.sunrise, data = leye.cs)
m58 <- lm(Fat ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + Detection + Event * ts.sunrise, data = leye.cs)
m59 <- lm(Fat ~ Sex + Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m60 <- lm(Fat ~ Sex + Event * ts.sunrise + PercentAg + Detection, data = leye.cs)
m61 <- lm(Fat ~ Sex + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)
m62 <- lm(Fat ~ Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection + Event*ts.sunrise, data = leye.cs)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:62)

models <- mget(model_names)

models$m.null <- m.null
models$m.global <- m.global

model_names <- c(model_names, "m.null", "m.global")

aictab(models, modnames = model_names)

# NEW RESULTS 2025-01-20 ####
confint(m3) # date significant
confint(m18) # date significant
confint(m12) # date significant


# ---------------------------------------------------------------------------- #

# Top Model Summaries

confint(m14)
summary(m14)

# PercentAg has significant impact on fat but in a positive direction...

summary(m35)
confint(m35)

summary(m45)
confint(m45)

summary(m40)
confint(m40)

summary(m24)
confint(m24)

summary(m17)
confint(m17)

summary(m3)
confint(m3)

summary(m33)
confint(m33)

summary(m37)
confint(m37)

# Detection has no significant effect on Lesser Yellowlegs fat.
# Agricultural intensity does have a significant effect on Lesser Yellowlegs fat but it's positive.
# Interaction between event * time significant.

# ---------------------------------------------------------------------------- #

# Model Assumptions:
par(mfrow = c(2,2))
plot(m14)

# Severe violations of homoscedasticity

# ---------------------------------------------------------------------------- #

# Graph top models

