#-----------------------------------------#
#   Lesser Yellowlegs Uric Acid Analysis  #
#          Created 01/20/2025             #          
#         Modified 01/20/2025             #
#-----------------------------------------#

# load packages
library(dplyr)
library(ggplot2)
library(AICcmodavg)
library(RColorBrewer)
library(trtools)
library(tidyverse)
library(lme4)
library(car)

# ---------------------------------------------------------------------------- #

# Data Processing & Manipulation ####

# Read data
setwd("processed_data")
birds <- read.csv("Shorebird_Data_Cleaned_2025-01-20.csv")

# Make neonicotinoid detection column (Detection/Non-detection)
birds$Detection <- ifelse(birds$OverallNeonic > 0, "Detection", "Non-detection")

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

# filter birds that only contain uric acid information (n = 29)
leye <- leye %>% 
  filter(!is.na(Uric))

# Standardize continuous variables
leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

# ---------------------------------------------------------------------------- #

# Data Visualization ####

# Uric Acid ~ Detection
ggplot(leye, 
       aes(x = Detection, y = Uric, fill = Detection)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(color = Detection),
              position = position_jitterdodge(jitter.width = 0.2),  
              size = 2, 
              alpha = 0.8) + 
  theme_light() +
  labs(x = "Neonicotinoid Detection", 
       y = "Uric Acid") +  
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

# Uric Acid ~ Percent Ag #### # no clear relationship
ggplot(leye, 
       aes(x = PercentAg, y = Uric)) + 
  geom_point(size = 2) + 
  labs(x = "Surrounding Agricultural Intensity", 
       y = "LEYE Uric Acid Levels (umol)") +  
  theme_classic() +
  theme(text = element_text(size = 20),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  geom_smooth(method = "lm", color = "black", size = 1)  +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = scales::percent)

# Uric Acid ~ Time
ggplot(leye, 
       aes(x = ts.sunrise, y = Uric)) + 
  geom_point(size = 2) + 
  labs(x = "Capture Time Since Sunrise", 
       y = "LEYE Uric Acid Levels (umol)") +  
  theme_classic() +
  theme(text = element_text(size = 20),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  geom_smooth(method = "lm", color = "black", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 0, y = 2000, 
           label = "Sunrise", angle = 90, 
           vjust = -0.5, hjust = 0.5,
           size = 5) 

# ---------------------------------------------------------------------------- #

# Lesser Yellowlegs Modeling ####

## Univariate Analysis: Date ####

m1 <- lm(Uric ~ DaysIntoSeason_S, data = leye.cs)
m2 <- lm(Uric ~ Julian, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Julian performs better (wt = 0.61) but use date into trapping season for consistency

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Event and Capture Time ####

m1 <- lm(Uric ~ Event * ts.sunrise, data = leye.cs)
m2 <- lm(Uric ~ Event + ts.sunrise, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model without interaction performs significantly better.

# ---------------------------------------------------------------------------- #

## Best Time Variable: Time or Time^2 ####

m1 <- lm(Uric ~ ts.sunrise, data = leye.cs)
m2 <- lm(Uric ~ ts.sunrise + I(ts.sunrise^2), data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with linear transformation performed best

# ---------------------------------------------------------------------------- #

# INTERACTIONS INCLUDED: NONE

# ---------------------------------------------------------------------------- #


## Single Models ####

m.global <- lm(Uric ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + Detection + PercentAg, 
               data = leye.cs)

m.null <- lm(Uric ~ 1, data = leye.cs)

m1 <- lm(Uric ~ Sex, data = leye.cs)
m2 <- lm(Uric ~ Event, data = leye.cs)
m3 <- lm(Uric ~ DaysIntoSeason_S, data = leye.cs)
m4 <- lm(Uric ~ ts.sunrise, data = leye.cs)
m5 <- lm(Uric ~ PercentAg, data = leye.cs)
m6 <- lm(Uric ~ Detection, data = leye.cs)

## Additive Models ####

### Two additive combinations ####
m7 <- lm(Uric ~ Sex + Event, data = leye.cs)
m8 <- lm(Uric ~ Sex + DaysIntoSeason_S, data = leye.cs)
m9 <- lm(Uric ~ Sex + ts.sunrise, data = leye.cs)
m10 <- lm(Uric ~ Sex + PercentAg, data = leye.cs)
m11 <- lm(Uric ~ Sex + Detection, data = leye.cs)
m12 <- lm(Uric ~ Event + DaysIntoSeason_S, data = leye.cs)
m13 <- lm(Uric ~ Event + ts.sunrise, data = leye.cs)
m14 <- lm(Uric ~ Event + PercentAg, data = leye.cs)
m15 <- lm(Uric ~ DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m16 <- lm(Uric ~ DaysIntoSeason_S + PercentAg, data = leye.cs)
m17 <- lm(Uric ~ DaysIntoSeason_S + Detection, data = leye.cs)
m18 <- lm(Uric ~ ts.sunrise + PercentAg, data = leye.cs)
m19 <- lm(Uric ~ ts.sunrise + Detection, data = leye.cs)
m20 <- lm(Uric ~ PercentAg + Detection, data = leye.cs)

### Three additive combinations ####
m21 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S, data = leye.cs)
m22 <- lm(Uric ~ Sex + Event + ts.sunrise, data = leye.cs)
m23 <- lm(Uric ~ Sex + Event + PercentAg, data = leye.cs)
m24 <- lm(Uric ~ Sex + Event + Detection, data = leye.cs)
m25 <- lm(Uric ~ Sex + DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m26 <- lm(Uric ~ Sex + DaysIntoSeason_S + PercentAg, data = leye.cs)
m27 <- lm(Uric ~ Sex + DaysIntoSeason_S + Detection, data = leye.cs)
m28 <- lm(Uric ~ Sex + ts.sunrise + PercentAg, data = leye.cs)
m29 <- lm(Uric ~ Sex + ts.sunrise + Detection, data = leye.cs)
m30 <- lm(Uric ~ Sex + PercentAg + Detection, data = leye.cs)
m31 <- lm(Uric ~ Event + DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m32 <- lm(Uric ~ Event + DaysIntoSeason_S + PercentAg, data = leye.cs)
m33 <- lm(Uric ~ Event + DaysIntoSeason_S + Detection, data = leye.cs)
m34 <- lm(Uric ~ Event + ts.sunrise + PercentAg, data = leye.cs)
m35 <- lm(Uric ~ Event + ts.sunrise + Detection, data = leye.cs)
m36 <- lm(Uric ~ Event + PercentAg + Detection, data = leye.cs)
m37 <- lm(Uric ~ DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m38 <- lm(Uric ~ DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m39 <- lm(Uric ~ DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m40 <- lm(Uric ~ ts.sunrise + PercentAg + Detection, data = leye.cs)

### Four additive combinations ####
m41 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m42 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S + PercentAg, data = leye.cs)
m43 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S + Detection, data = leye.cs)
m44 <- lm(Uric ~ Sex + Event + ts.sunrise + PercentAg, data = leye.cs)
m45 <- lm(Uric ~ Sex + Event + ts.sunrise + Detection, data = leye.cs)
m46 <- lm(Uric ~ Sex + Event + PercentAg + Detection, data = leye.cs)
m47 <- lm(Uric ~ Sex + DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m48 <- lm(Uric ~ Sex + DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m49 <- lm(Uric ~ Sex + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m50 <- lm(Uric ~ Sex + ts.sunrise + PercentAg + Detection, data = leye.cs)
m51 <- lm(Uric ~ Event + DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m52 <- lm(Uric ~ Event + DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m53 <- lm(Uric ~ Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m54 <- lm(Uric ~ Event + ts.sunrise + PercentAg + Detection, data = leye.cs)
m55 <- lm(Uric ~ DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)

### Five additive combinations ####

m56 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m57 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m58 <- lm(Uric ~ Sex + Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m59 <- lm(Uric ~ Sex + Event + ts.sunrise + PercentAg + Detection, data = leye.cs)
m60 <- lm(Uric ~ Sex + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)
m61 <- lm(Uric ~ Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:61)

models <- mget(model_names)

models$m.null <- m.null
models$m.global <- m.global

model_names <- c(model_names, "m.null", "m.global")

aictab(models, modnames = model_names)

# ---------------------------------------------------------------------------- #

# Top model summaries ####

summary(m4)
confint(m4)
plot(m4)

# plot top model ####
m <- lm(Uric ~ ts.sunrise, data = leye)
d <- expand.grid(ts.sunrise = seq(min(leye$ts.sunrise),
                                  max(leye$ts.sunrise),
                                  length = 1000))

d$yhat <- predict(m, newdata = d)

d$se <- predict(m, newdata = d, se.fit = TRUE)$se.fit
d$lower <- d$yhat - 2*d$se
d$upper <- d$yhat + 2*d$se

ggplot(d,
       aes(x = ts.sunrise,
           y = Uric)) +
  geom_point(color = "tan2") + geom_line(aes(y = yhat), data = d, color = "tan2") +
  geom_ribbon(aes(x = ts.sunrise, ymin = lower, ymax = upper, 
                  y = NULL), data = d, fill = "tan1", alpha = 0.25) +
  theme_light() +
  labs(x = "Time of Capture since Sunrise (min)",
       y = "Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 0, y = 2, 
           label = "Sunrise", angle = 90, 
           vjust = -0.5, hjust = 0.5,
           size = 5) 

# ---------------------------------------------------------------------------- #

library(splines)
library(mgcv)
library(MASS)

# class example
knots <- seq(2.4, 57.6, length = 27)
m <- gam(accel ~ s(times, bs = "bs", k = length(knots), sp = 0),
         data = mcycle, knots = list(x = knots))

## LEYE GAM (standardized) ####
knots <- seq(min(leye.cs$ts.sunrise), max(leye.cs$ts.sunrise), length = 6)

m <- gam(Uric ~ s(ts.sunrise, bs = "bs", k = length(knots), sp = 0), 
         data = leye.cs, knots = list(x = knots))

summary(m)
par(mfrow = c(2,2))
plot(m)

E.m <- resid(m)
Fit.m <- fitted(m)
plot(x = Fit.m, y = E.m, xlab = "Fitted values",
     ylab = "Residuals")

# Plot model
d <- expand.grid(ts.sunrise = seq(min(leye.cs$ts.sunrise),
                                  max(leye.cs$ts.sunrise),
                                  length = 1000))

d$yhat <- predict(m, newdata = d)

d$se <- predict(m, newdata = d, se.fit = TRUE)$se.fit
d$lower <- d$yhat - 2*d$se
d$upper <- d$yhat + 2*d$se

ggplot() +
  geom_point(data = leye, aes(x = ts.sunrise, y = Uric), color = "tan2") +
  geom_line(data = d, aes(x = ts.sunrise, y = yhat), color = "tan2") +
  geom_ribbon(data = d, aes(x = ts.sunrise, ymin = lower, ymax = upper), 
              fill = "tan1", alpha = 0.25) +
  theme_light() +
  labs(x = "Time of Capture since Sunrise (min)",
       y = "LEYE Uric Acid Levels (umol)") +
  theme(axis.title.x = element_text(size = 21, margin = margin(t = 13)),
        axis.title.y = element_text(size = 21, margin = margin(r = 13)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 0, y = 2000, label = "Sunrise", angle = 90, 
           vjust = -0.5, hjust = 0.5, size = 5)

## LEYE GAM (unstandardized) ####
knots <- seq(min(leye$ts.sunrise), max(leye$ts.sunrise), length = 6)

m <- gam(Uric ~ s(ts.sunrise, bs = "bs", k = length(knots), sp = 0), 
         data = leye, knots = list(x = knots))

# Plot model
d <- expand.grid(ts.sunrise = seq(min(leye$ts.sunrise),
                                  max(leye$ts.sunrise),
                                  length = 1000))

d$yhat <- predict(m, newdata = d)

d$se <- predict(m, newdata = d, se.fit = TRUE)$se.fit
d$lower <- d$yhat - 2*d$se
d$upper <- d$yhat + 2*d$se


ggplot() +
  geom_point(data = leye, aes(x = ts.sunrise, y = Uric), color = "tan2") +
  geom_line(data = d, aes(x = ts.sunrise, y = yhat), color = "tan2") +
  geom_ribbon(data = d, aes(x = ts.sunrise, ymin = lower, ymax = upper), 
              fill = "tan1", alpha = 0.25) +
  theme_light() +
  labs(x = "Time of Capture since Sunrise (min)",
       y = "LEYE Uric Acid Levels (umol)") +
  theme(axis.title.x = element_text(size = 21, margin = margin(t = 13)),
        axis.title.y = element_text(size = 21, margin = margin(r = 13)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 0, y = 2000, label = "Sunrise", angle = 90, 
           vjust = -0.5, hjust = 0.5, size = 5)


m1 <- gam(Uric ~ s(ts.sunrise), data = leye.cs)

knots <- seq(min(leye.cs$ts.sunrise), max(leye.cs$ts.sunrise), length = 6)
m2 <- gam(Uric ~ s(ts.sunrise, bs = "bs", k = length(knots), sp = 0), 
          data = leye.cs, knots = list(x = knots))

aic_m1 <- AIC(m1)
aic_m2 <- AIC(m2)

# Compare AIC values
aic_values <- c(m1 = aic_m1, m2 = aic_m2)

# linear model slightly better fit than gam but there is a pattern

