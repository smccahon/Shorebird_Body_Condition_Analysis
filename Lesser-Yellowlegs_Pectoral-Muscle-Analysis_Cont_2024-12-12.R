#---------------------------------------------#
#  Lesser Yellowlegs Pectoral Score Analysis  #
#          Created 12/12/2024                 #          
#         Modified 12/12/2024                 #
#---------------------------------------------#

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

# ---------------------------------------------------------------------------- #

# Data Processing & Manipulation ####

# Read data
setwd("processed_data")
birds <- read.csv("LM_ShorebirdsALLNeg.csv")

# Make neonicotinoid detection column (Detection/Non-detection)
birds$Detection <- ifelse(birds$OverallNeonic > 0, "Detection", "Non-detection")

# Make sampling event column
birds <- birds %>% 
  mutate(Event = paste(Season, Year, sep = " ")) %>% 
  mutate(Event = ifelse(Season %in% c("Spring", "Fall"),
                        paste(Season, Year, sep = " "),
                        NA))

birds <- birds %>% 
  mutate(LogAg = log10(PercentAg))

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


# filter birds that only contain pectoral score information
leye <- leye %>% 
  filter(!is.na(PecScore))


# Standardize continuous variables
leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

# ---------------------------------------------------------------------------- #

# Data Visualization ####

# Pectoral Shape ~ Agricultural Intensity

ggplot(leye, 
       aes(x = PercentAg, y = PecScore)) +
  geom_point(size = 2) +  
  theme_classic() +
  labs(x = "Surrounding Agricultural Intensity", 
       y = "Lesser Yellowlegs Pectoral Shape") +  
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = scales::percent) +
  scale_y_continuous(breaks = 0:3)

# Pectoral Shape ~ Detection

ggplot(leye, 
       aes(x = Detection, y = PecScore)) +
  geom_boxplot() +  
  theme_classic() +
  labs(x = "Neonicotinoid Detection", 
       y = "Lesser Yellowlegs Pectoral Shape") +  
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_y_continuous(breaks = 0:3)

# Pectoral Shape ~ Date
ggplot(leye, 
       aes(x = DaysIntoSeason_S, y = PecScore)) +
  geom_point() +  
  theme_classic() +
  labs(x = "Days Into Trapping Season", 
       y = "Lesser Yellowlegs Pectoral Shape") +  
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_y_continuous(breaks = 0:3)

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Agricultural Intensity ####

# Can't consider categorical variable because there's only 3 birds in low ag category

m1 <- lm(PecScore ~ PercentAg, data = leye.cs)
m2 <- lm(PecScore ~ LogAg, data = leye.cs)
m3 <- lm(PecScore ~ PercentAg + I(PercentAg^2), data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Linear transformation of agricultural intensity has best fit

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Event and Capture Time ####

m1 <- lm(PecScore ~ Event * ts.sunrise, data = leye.cs)
m2 <- lm(PecScore ~ Event + ts.sunrise, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model without interaction performs significantly better

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Time

m1 <- lm(PecScore ~ ts.sunrise, data = leye.cs)
m2 <- lm(PecScore ~ ts.sunrise + I(ts.sunrise^2), data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# quadratic relationship performs the best

# ---------------------------------------------------------------------------- #

# INTERACTIONS INCLUDED: None, time (quadratic), agricultural intensity (linear)

# ---------------------------------------------------------------------------- #

m.global <- lm(PecScore ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + Detection + PercentAg, 
               data = leye.cs)


m.null <- lm(PecScore ~ 1, data = leye.cs)

m1 <- lm(PecScore ~ Sex, data = leye.cs)
m2 <- lm(PecScore ~ Event, data = leye.cs)
m3 <- lm(PecScore ~ DaysIntoSeason_S, data = leye.cs)
m4 <- lm(PecScore ~ ts.sunrise + I(ts.sunrise^2), data = leye.cs)
m5 <- lm(PecScore ~ PercentAg, data = leye.cs)
m6 <- lm(PecScore ~ Detection, data = leye.cs)

## Additive Models ####

### Two additive combinations ####
m7 <- lm(PecScore ~ Sex + Event, data = leye.cs)
m8 <- lm(PecScore ~ Sex + DaysIntoSeason_S, data = leye.cs)
m9 <- lm(PecScore ~ Sex + ts.sunrise + I(ts.sunrise^2), data = leye.cs)
m10 <- lm(PecScore ~ Sex + PercentAg, data = leye.cs)
m11 <- lm(PecScore ~ Sex + Detection, data = leye.cs)
m12 <- lm(PecScore ~ Event + DaysIntoSeason_S, data = leye.cs)
m13 <- lm(PecScore ~ Event + ts.sunrise + I(ts.sunrise^2), data = leye.cs)
m14 <- lm(PecScore ~ Event + PercentAg, data = leye.cs)
m15 <- lm(PecScore ~ Event + Detection, data = leye.cs)
m16 <- lm(PecScore ~ DaysIntoSeason_S + ts.sunrise + I(ts.sunrise^2), data = leye.cs)
m17 <- lm(PecScore ~ DaysIntoSeason_S + PercentAg, data = leye.cs)
m18 <- lm(PecScore ~ DaysIntoSeason_S + Detection, data = leye.cs)
m19 <- lm(PecScore ~ ts.sunrise + I(ts.sunrise^2) + PercentAg, data = leye.cs)
m20 <- lm(PecScore ~ ts.sunrise + I(ts.sunrise^2) + Detection, data = leye.cs)
m21 <- lm(PecScore ~ PercentAg + Detection, data = leye.cs)

### Three additive combinations ####
m22 <- lm(PecScore ~ Sex + Event + DaysIntoSeason_S, data = leye.cs)
m23 <- lm(PecScore ~ Sex + Event + ts.sunrise + I(ts.sunrise^2), data = leye.cs)
m24 <- lm(PecScore ~ Sex + Event + PercentAg, data = leye.cs)
m25 <- lm(PecScore ~ Sex + Event + Detection, data = leye.cs)
m26 <- lm(PecScore ~ Sex + DaysIntoSeason_S + ts.sunrise + I(ts.sunrise^2), data = leye.cs)
m27 <- lm(PecScore ~ Sex + DaysIntoSeason_S + PercentAg, data = leye.cs)
m28 <- lm(PecScore ~ Sex + DaysIntoSeason_S + Detection, data = leye.cs)
m29 <- lm(PecScore ~ Sex + ts.sunrise + I(ts.sunrise^2) + PercentAg, data = leye.cs)
m30 <- lm(PecScore ~ Sex + ts.sunrise + I(ts.sunrise^2) + Detection, data = leye.cs)
m31 <- lm(PecScore ~ Sex + PercentAg + Detection, data = leye.cs)
m32 <- lm(PecScore ~ Event + DaysIntoSeason_S + ts.sunrise + I(ts.sunrise^2) + Event + ts.sunrise + I(ts.sunrise^2), data = leye.cs)
m33 <- lm(PecScore ~ Event + DaysIntoSeason_S + PercentAg, data = leye.cs)
m34 <- lm(PecScore ~ Event + DaysIntoSeason_S + Detection, data = leye.cs)
m35 <- lm(PecScore ~ Event + ts.sunrise + I(ts.sunrise^2) + PercentAg, data = leye.cs)
m36 <- lm(PecScore ~ Event + ts.sunrise + I(ts.sunrise^2) + Detection, data = leye.cs)
m37 <- lm(PecScore ~ Event + PercentAg + Detection, data = leye.cs)
m38 <- lm(PecScore ~ DaysIntoSeason_S + ts.sunrise + I(ts.sunrise^2) + PercentAg, data = leye.cs)
m39 <- lm(PecScore ~ DaysIntoSeason_S + ts.sunrise + I(ts.sunrise^2) + Detection, data = leye.cs)
m40 <- lm(PecScore ~ DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m41 <- lm(PecScore ~ ts.sunrise + I(ts.sunrise^2) + PercentAg + Detection, data = leye.cs)

### Four additive combinations ####
m42 <- lm(PecScore ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + I(ts.sunrise^2) + Event + ts.sunrise + I(ts.sunrise^2), data = leye.cs)
m43 <- lm(PecScore ~ Sex + Event + DaysIntoSeason_S + PercentAg, data = leye.cs)
m44 <- lm(PecScore ~ Sex + Event + DaysIntoSeason_S + Detection, data = leye.cs)
m45 <- lm(PecScore ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + PercentAg, data = leye.cs)
m46 <- lm(PecScore ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + Detection, data = leye.cs)
m47 <- lm(PecScore ~ Sex + Event + PercentAg + Detection, data = leye.cs)
m48 <- lm(PecScore ~ Sex + DaysIntoSeason_S + ts.sunrise + I(ts.sunrise^2) + PercentAg, data = leye.cs)
m49 <- lm(PecScore ~ Sex + DaysIntoSeason_S + ts.sunrise + I(ts.sunrise^2) + Detection, data = leye.cs)
m50 <- lm(PecScore ~ Sex + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m51 <- lm(PecScore ~ Sex + ts.sunrise + I(ts.sunrise^2) + PercentAg + Detection, data = leye.cs)
m52 <- lm(PecScore ~ Event + DaysIntoSeason_S + ts.sunrise + I(ts.sunrise^2) + PercentAg + Event + ts.sunrise + I(ts.sunrise^2), data = leye.cs)
m53 <- lm(PecScore ~ Event + DaysIntoSeason_S + ts.sunrise + I(ts.sunrise^2) + Detection + Event + ts.sunrise + I(ts.sunrise^2), data = leye.cs)
m54 <- lm(PecScore ~ Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m55 <- lm(PecScore ~ Event + ts.sunrise + I(ts.sunrise^2) + PercentAg + Detection, data = leye.cs)
m56 <- lm(PecScore ~ DaysIntoSeason_S + ts.sunrise + I(ts.sunrise^2) + PercentAg + Detection, data = leye.cs)

### Five additive combinations ####

m57 <- lm(PecScore ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + I(ts.sunrise^2) + PercentAg + Event + ts.sunrise + I(ts.sunrise^2), data = leye.cs)
m58 <- lm(PecScore ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + I(ts.sunrise^2) + Detection + Event + ts.sunrise + I(ts.sunrise^2), data = leye.cs)
m59 <- lm(PecScore ~ Sex + Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m60 <- lm(PecScore ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + PercentAg + Detection, data = leye.cs)
m61 <- lm(PecScore ~ Sex + DaysIntoSeason_S + ts.sunrise + I(ts.sunrise^2) + PercentAg + Detection, data = leye.cs)
m62 <- lm(PecScore ~ Event + DaysIntoSeason_S + ts.sunrise + I(ts.sunrise^2) + PercentAg + Detection + Event + ts.sunrise + I(ts.sunrise^2), data = leye.cs)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:62)

models <- mget(model_names)

models$m.null <- m.null
models$m.global <- m.global

model_names <- c(model_names, "m.null", "m.global")

aictab(models, modnames = model_names)

# ---------------------------------------------------------------------------- #

# Top Model Summaries

confint(m3)

summary(m18)
confint(m18)
confint(m8)
confint(m5)

summary(m28)
confint(m28)

# Date has a significant impact on pectoral score (positive relationship)

# ---------------------------------------------------------------------------- #

# Model Assumptions

par(mfrow = c(2,2))
plot(m3)

# some severe violations and patterns in the residuals...