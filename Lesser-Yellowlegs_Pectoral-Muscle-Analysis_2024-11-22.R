#---------------------------------------------#
#  Lesser Yellowlegs Pectoral Muscle Analysis #
#           Created 11/11/2024                #          
#           Modified 01/06/2025               #
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


# Subset data for Lesser Yellowlegs
leye <- subset(birds, Species %in% c("LesserYellowlegs"))

# Filter data that only contain pectoral muscle scores
leye <- leye %>% 
  filter(!is.na(PecSizeBest))

# Standardize continuous variables
leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

# ---------------------------------------------------------------------------- #

# Data Visualization ####

# Pectoral Muscle ~ Concentration
ggplot(data = leye, aes(x = OverallNeonic, y = PecSizeBest)) + 
  theme_classic() + 
  geom_point() + 
  labs(x = "Neonicotinoid Concentration (ug/L)", 
       y = "Lesser Yellowlegs Body PecSizeBest")

# PecSizeBest ~ Log(Concentration)
ggplot(data = leye, aes(x = LogNeonic, y = PecSizeBest)) +
  theme_classic() +
  geom_point() + 
  labs(x = "Log(Concentration)",
       y = "Lesser Yellowlegs Body PecSizeBest")

# PecSizeBest ~ Detection
ggplot(leye, aes(x = Detection, y = PecSizeBest)) + 
  geom_boxplot(fill = "tan2") +  
  theme_classic() +
  labs(x = "Neonicotinoid Detection",
       y = expression("Lesser Yellowlegs Breast Muscle Size" ~~~ (mm[score]))) +
  theme(text = element_text(size = 20),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.title.y = element_text(margin = margin(r = 10), size = 15),
        axis.title.x = element_text(margin = margin(t = 10), size = 15))

# PecSizeBest ~ Migration Date
leye$Event <- factor(leye$Event, levels = c("Fall_2021", "Spring_2022", "Fall_2023"))
ggplot(data = leye, aes(x = DaysIntoSeason_S, y = PecSizeBest, color = Event)) +
  theme_classic() +
  geom_point() +
  labs(x = "Date into Migration Season",
       y = "Lesser Yellowlegs Body PecSizeBest",
       color = "Sampling Event") +
  scale_color_manual(
    values = c("orange", "purple", "green"),
    labels = c("Fall 2021", "Spring 2022", "Fall 2023")
  )

# PecSizeBest ~ Julian
ggplot(data = leye, aes(x = Julian, y = PecSizeBest, color = Event)) +
  theme_classic() +
  geom_point() +
  labs(x = "Julian Day",
       y = "Lesser Yellowlegs Body PecSizeBest",
       color = "Sampling Event") +
  scale_color_manual(
    values = c("orange", "purple", "green"),
    labels = c("Fall 2021", "Spring 2022", "Fall 2023")
  )

# PecSizeBest ~ Sex
ggplot(data = leye, aes(x = Sex, y = PecSizeBest)) + 
  theme_classic() + 
  geom_boxplot() + 
  geom_jitter() +
  labs(x = "Sex",
       y = "Lesser Yellowlegs Body PecSizeBest")

# PecSizeBest ~ Sampling Event
ggplot(data = leye, aes(x = Event, y = PecSizeBest, color = Event)) +
  theme_classic() +
  geom_boxplot() +
  geom_jitter() + 
  labs(x = "Sampling Event",
       y = "Lesser Yellowlegs Body PecSizeBest",
       color = "Sampling Event") +
  scale_color_manual(
    values = c("orange", "purple", "green"),
    labels = c("Fall 2021", "Spring 2022", "Fall 2023")
  ) + 
  scale_x_discrete(
    labels = c("Fall_2021" = "Fall 2021", "Spring_2022" = "Spring 2022", "Fall_2023" = "Fall 2023")
  )

# PecSizeBest ~ Capture Time
ggplot(data = leye, aes(x = ts.sunrise, y = PecSizeBest)) +
  theme_classic() +
  geom_smooth(method = "lm",  color = "blue") +
  geom_point() +
  labs(x = "Capture Time (min.)",
       y = "Lesser Yellowlegs Body PecSizeBest (g)")

# PecSizeBest ~ AgIntensity
# Negative relationship!
ggplot(data = leye, aes(x = PercentAg, y = PecSizeBest)) +
  theme_classic() +
  geom_smooth(method = "lm",  color = "blue") +
  geom_point()

# ---------------------------------------------------------------------------- #

# Lesser Yellowlegs Modeling ####

## Univariate Analysis: Date ####

m1 <- lm(PecSizeBest ~ DaysIntoSeason_S, data = leye.cs)
m2 <- lm(PecSizeBest ~ Julian, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Date into migration season performs better than Julian date but fairly similar

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Neonicotinoids ####

m1 <- lm(PecSizeBest ~ Detection, data = leye.cs)
m2 <- lm(PecSizeBest ~ OverallNeonic, data = leye.cs)
m3 <- lm(PecSizeBest ~ LogNeonic, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Detection and Log transformation of concentrations perform significantly better
# Conclusion: use just detection for consistency with other models

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Agricultural Intensity ####

# Can't do agricultural category due to low samples size in low ag class

m1 <- lm(PecSizeBest ~ PercentAg, data = leye.cs)
m2 <- lm(PecSizeBest ~ I(PercentAg^2), data = leye.cs)
m3 <- lm(PecSizeBest ~ LogAg, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Continous predictor performs best (PercentAg)

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Event and Capture Time ####

m1 <- lm(PecSizeBest ~ Event * ts.sunrise, data = leye.cs)
m2 <- lm(PecSizeBest ~ Event + I(ts.sunrise^2), data = leye.cs)
m3 <- lm(PecSizeBest ~ Event + ts.sunrise, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with interaction performs the best.

# ---------------------------------------------------------------------------- #

## Single Models ####

m.global <- lm(PecSizeBest ~ Sex + Event*ts.sunrise + DaysIntoSeason_S + Detection + PercentAg, 
                 data = leye.cs)

m.null <- lm(PecSizeBest ~ 1, data = leye.cs)

m1 <- lm(PecSizeBest ~ Sex, data = leye.cs)
m2 <- lm(PecSizeBest ~ Event, data = leye.cs)
m3 <- lm(PecSizeBest ~ DaysIntoSeason_S, data = leye.cs)
m4 <- lm(PecSizeBest ~ ts.sunrise, data = leye.cs)
m5 <- lm(PecSizeBest ~ PercentAg, data = leye.cs)
m6 <- lm(PecSizeBest ~ Detection, data = leye.cs)

## Additive Models ####

### Two additive combinations ####
m7 <- lm(PecSizeBest ~ Sex + Event, data = leye.cs)
m8 <- lm(PecSizeBest ~ Sex + DaysIntoSeason_S, data = leye.cs)
m9 <- lm(PecSizeBest ~ Sex + ts.sunrise, data = leye.cs)
m10 <- lm(PecSizeBest ~ Sex + PercentAg, data = leye.cs)
m11 <- lm(PecSizeBest ~ Sex + Detection, data = leye.cs)
m12 <- lm(PecSizeBest ~ Event + DaysIntoSeason_S, data = leye.cs)
m13 <- lm(PecSizeBest ~ Event * ts.sunrise, data = leye.cs)
m14 <- lm(PecSizeBest ~ Event + PercentAg, data = leye.cs)
m15 <- lm(PecSizeBest ~ DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m16 <- lm(PecSizeBest ~ DaysIntoSeason_S + PercentAg, data = leye.cs)
m17 <- lm(PecSizeBest ~ DaysIntoSeason_S + Detection, data = leye.cs)
m18 <- lm(PecSizeBest ~ ts.sunrise + PercentAg, data = leye.cs)
m19 <- lm(PecSizeBest ~ ts.sunrise + Detection, data = leye.cs)
m20 <- lm(PecSizeBest ~ PercentAg + Detection, data = leye.cs)

### Three additive combinations ####
m21 <- lm(PecSizeBest ~ Sex + Event + DaysIntoSeason_S, data = leye.cs)
m22 <- lm(PecSizeBest ~ Sex + Event * ts.sunrise, data = leye.cs)
m23 <- lm(PecSizeBest ~ Sex + Event + PercentAg, data = leye.cs)
m24 <- lm(PecSizeBest ~ Sex + Event + Detection, data = leye.cs)
m25 <- lm(PecSizeBest ~ Sex + DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m26 <- lm(PecSizeBest ~ Sex + DaysIntoSeason_S + PercentAg, data = leye.cs)
m27 <- lm(PecSizeBest ~ Sex + DaysIntoSeason_S + Detection, data = leye.cs)
m28 <- lm(PecSizeBest ~ Sex + ts.sunrise + PercentAg, data = leye.cs)
m29 <- lm(PecSizeBest ~ Sex + ts.sunrise + Detection, data = leye.cs)
m30 <- lm(PecSizeBest ~ Sex + PercentAg + Detection, data = leye.cs)
m31 <- lm(PecSizeBest ~ Event + DaysIntoSeason_S + ts.sunrise + Event * ts.sunrise, data = leye.cs)
m32 <- lm(PecSizeBest ~ Event + DaysIntoSeason_S + PercentAg, data = leye.cs)
m33 <- lm(PecSizeBest ~ Event + DaysIntoSeason_S + Detection, data = leye.cs)
m34 <- lm(PecSizeBest ~ Event * ts.sunrise + PercentAg, data = leye.cs)
m35 <- lm(PecSizeBest ~ Event * ts.sunrise + Detection, data = leye.cs)
m36 <- lm(PecSizeBest ~ Event + PercentAg + Detection, data = leye.cs)
m37 <- lm(PecSizeBest ~ DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m38 <- lm(PecSizeBest ~ DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m39 <- lm(PecSizeBest ~ DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m40 <- lm(PecSizeBest ~ ts.sunrise + PercentAg + Detection, data = leye.cs)

### Four additive combinations ####
m41 <- lm(PecSizeBest ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + Event*ts.sunrise, data = leye.cs)
m42 <- lm(PecSizeBest ~ Sex + Event + DaysIntoSeason_S + PercentAg, data = leye.cs)
m43 <- lm(PecSizeBest ~ Sex + Event + DaysIntoSeason_S + Detection, data = leye.cs)
m44 <- lm(PecSizeBest ~ Sex + Event * ts.sunrise + PercentAg, data = leye.cs)
m45 <- lm(PecSizeBest ~ Sex + Event * ts.sunrise + Detection, data = leye.cs)
m46 <- lm(PecSizeBest ~ Sex + Event + PercentAg + Detection, data = leye.cs)
m47 <- lm(PecSizeBest ~ Sex + DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m48 <- lm(PecSizeBest ~ Sex + DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m49 <- lm(PecSizeBest ~ Sex + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m50 <- lm(PecSizeBest ~ Sex + ts.sunrise + PercentAg + Detection, data = leye.cs)
m51 <- lm(PecSizeBest ~ Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Event*ts.sunrise, data = leye.cs)
m52 <- lm(PecSizeBest ~ Event + DaysIntoSeason_S + ts.sunrise + Detection + Event*ts.sunrise, data = leye.cs)
m53 <- lm(PecSizeBest ~ Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m54 <- lm(PecSizeBest ~ Event * ts.sunrise + PercentAg + Detection, data = leye.cs)
m55 <- lm(PecSizeBest ~ DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)

### Five additive combinations ####

m56 <- lm(PecSizeBest ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Event * ts.sunrise, data = leye.cs)
m57 <- lm(PecSizeBest ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + Detection + Event * ts.sunrise, data = leye.cs)
m58 <- lm(PecSizeBest ~ Sex + Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m59 <- lm(PecSizeBest ~ Sex + Event * ts.sunrise + PercentAg + Detection, data = leye.cs)
m60 <- lm(PecSizeBest ~ Sex + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)
m61 <- lm(PecSizeBest ~ Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection + Event*ts.sunrise, data = leye.cs)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:61)

models <- mget(model_names)

models$m.null <- m.null
models$m.global <- m.global

model_names <- c(model_names, "m.null", "m.global")

aictab(models, modnames = model_names)

# ---------------------------------------------------------------------------- #

# Top Model Summaries: All Models ####
cbind(summary(m6)$coefficients, confint(m6))

# Birds with detection had significantly larger pectoral muscles...

cbind(summary(m3)$coefficients, confint(m3))

# Conclusion: There is a negative effect of date into trapping season on pectoral muscle size
# As the date into the trapping season progresses, pectoral muscle size decreases
# Opposite of prediction

cbind(summary(m17)$coefficients, confint(m17))

# Conclusion: Neither effect is significant.

cbind(summary(m20)$coefficients, confint(m20))
cbind(summary(m19)$coefficients, confint(m19))
cbind(summary(m15)$coefficients, confint(m15))
cbind(summary(m11)$coefficients, confint(m11))

summary(m5)
confint(m5)

## Goodness of fit ####
summary(m6)
summary(m20)
summary(m19)
summary(m15)
summary(m11)
summary(m3)
summary(m17)
# ---------------------------------------------------------------------------- #

# Graph Top Models ####

## LAB MEETING DETECTION ####
m <- lm(PecSizeBest ~ Detection, data = leye)

d <- expand.grid(Detection = unique(leye$Detection))

predictions <- predict(m, newdata = d, se.fit = TRUE)

d$predicted_Mass <- predictions$fit

d$lower_CI <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upper_CI <- d$predicted_Mass + 1.96 * predictions$se.fit

ggplot(d, aes(x = Detection, y = predicted_Mass)) +
  geom_point(size = 5, col = "black") +  # Points showing predicted values
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.1,
                col = "black",
                size = 1) +  # Add confidence intervals
  theme_light() +
  labs(x = NULL, 
       y = "Predicted Pectoral Muscle Size") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  theme(legend.position = "none")


## LAB MEETING AG ####
m <- lm(PecSizeBest ~ Detection + PercentAg, data = leye)

d <- expand.grid(Detection = c("Non-detection", "Detection"),
                 PercentAg = seq(min(leye$PercentAg),
                                 max(leye$PercentAg),
                                 length = 1000)) 

d <- cbind(d, predict(m, newdata = d, interval = "confidence"))

ggplot(d, aes(x = (PercentAg * 100), y = fit, color = Detection)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Detection), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_light() +
  labs(x = "Surrounding Agricultural Intensity (%)", 
       y = "Predicted Pectoral Muscle Size") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  


# Pectoral ~ Detection
m <- lm(PecSizeBest ~ Detection, data = leye)
summary(m)

d <- expand.grid(Detection = unique(leye$Detection))

predictions <- predict(m, newdata = d, se.fit = TRUE)

d$predicted_Mass <- predictions$fit

d$lower_CI <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upper_CI <- d$predicted_Mass + 1.96 * predictions$se.fit

ggplot(d, aes(x = Detection, y = predicted_Mass, color = Detection)) +
  geom_point(size =3,
             col = "black") +  # Points showing predicted values
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.1,
                col = "black",
                size = 0.7) +  # Add confidence intervals
  theme_classic() +
  labs(x = "Neonicotinoid Detection", 
       y = "Predicted Lesser Yellowlegs Breast Muscle Size" ~~~ (mm[score])) +
  theme(axis.title.x = element_text(size = 16,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 16,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none") +
  theme(legend.position = "none")

ggplot(d, aes(x = Detection, y = fit)) + 
  geom_line(aes(color = "tan2"), size = 1, show.legend = FALSE) +  
  theme_light() +
  labs(x = "Neonicotinoid Detection",
       y = expression("Lesser Yellowlegs Breast Muscle Size" ~~~ (mm[score])))


# Pectoral ~ Detection
ggplot(leye, aes(x = Detection, y = PecSizeBest, fill = Detection)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(color = Detection),
              position = position_jitterdodge(jitter.width = 0.2),  
              size = 2, 
              alpha = 0.8) + 
  theme_light() +
  labs(x = "Neonicotinoid Detection", 
       y = expression("Lesser Yellowlegs Breast Muscle Size" ~~~ (mm[score]))) +  
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

# Pectoral ~ Detection (different colors)
ggplot(leye, aes(x = Detection, y = PecSizeBest, fill = Detection)) +
  geom_boxplot(outlier.shape = NA, fill = "tan") +
  theme_light() +
  labs(x = "Neonicotinoid Detection", 
       y = expression("Lesser Yellowlegs Breast Muscle Size" ~~~ (mm[score]))) +  
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")


# Pectoral ~ DaysIntoSeason_S
m <- lm(PecSizeBest ~ DaysIntoSeason_S, data = leye)

d <- expand.grid(DaysIntoSeason_S = seq(min(leye$DaysIntoSeason_S), 
                                 max(leye$DaysIntoSeason_S),
                                 length.out = 1000))

d <- cbind(d, predict(m, newdata = d, interval = "confidence"))

ggplot(d, aes(x = DaysIntoSeason_S, y = fit)) + 
  geom_line(aes(color = "tan2"), size = 1, show.legend = FALSE) +  
  theme_light() +
  labs(x = "Days Into Trapping Season",
       y = expression("Lesser Yellowlegs Breast Muscle Size" ~~~ (mm[score]))) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = LogNeonic),  
              alpha = 0.25, color = NA, show.legend = FALSE,
              fill = "tan1") +
  geom_point(data = leye, aes(x = DaysIntoSeason_S, y = PecSizeBest), color = "tan2")




# Pectoral ~ Detection + Agricultural Intensity

m <- lm(PecSizeBest ~ Detection + PercentAg, data = leye)
summary(m)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg), 
                                 max(leye$PercentAg),
                                 length.out = 1000),
                 Detection = c("Detection"))

d <- cbind(d, predict(m, newdata = d, interval = "confidence"))

ggplot(d, aes(x = PercentAg, y = fit)) + 
  geom_line(aes(color = "tan2"), size = 1, show.legend = FALSE) +  
  theme_light() +
  labs(x = "Surrounding Agricultural Intensity (%)",
       y = expression("Lesser Yellowlegs Breast Muscle Size" ~~~ (mm[score]))) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),  
              alpha = 0.25, color = NA, show.legend = FALSE,
              fill = "tan1") +
  geom_point(data = leye, aes(x = PercentAg, y = PecSizeBest), color = "tan2")

# Pectoral ~ Agricultural Intensity
m <- lm(PecSizeBest ~ PercentAg, data = leye)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg), 
                                 max(leye$PercentAg),
                                 length.out = 1000))

d <- cbind(d, predict(m, newdata = d, interval = "confidence"))

ggplot(d, aes(x = PercentAg, y = fit)) + 
  geom_line(aes(color = "tan2"), size = 1, show.legend = FALSE) +  
  theme_classic() +
  labs(x = "Surrounding Agricultural Intensity (%)",
       y = expression("Lesser Yellowlegs Breast Muscle Size" ~~~ (mm[score]))) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),  
              alpha = 0.25, color = NA, show.legend = FALSE,
              fill = "tan1") +
  geom_point(data = leye, aes(x = PercentAg, y = PecSizeBest), color = "tan2") +
  theme(text = element_text(size = 20),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = scales::percent)

# ---------------------------------------------------------------------------- #

# Model assumptions ####

## Pectoral ~ Detection ####
m <- lm(PecSizeBest ~ Detection, data = leye.cs)
par(mfrow = c(1, 1))

# Plot boxplots to assess linearity assumption
# Two groups should be different from one another

ggplot(leye,
       aes(x = Detection, y = PecSizeBest)) +
  geom_boxplot()

# Plot residuals by categorical variable to assess homoscedasticity
boxplot(residuals(m) ~ leye.cs$Detection, 
        main = "Residuals by Detection",
        xlab = "Detection",
        ylab = "Residuals")

# Conclusion: No severe violations of homoscedasticity

## Pectoral ~ Agricultural Intensity + Detection ####
m <- lm(PecSizeBest ~ PercentAg + Detection, data = leye.cs)
par(mfrow = c(2, 2))
plot(m)

m <- lm(PecSizeBest ~ PercentAg, data = leye.cs)
plot(m)

# mostly good...some pattern in residuals (opposite megaphone) so no big concern

ggplot(leye,
       aes(x = PercentAg, y = PecSizeBest)) +
  geom_point()

## Pectoral ~ Agricultural Intensity
m <- lm(PecSizeBest ~ PercentAg, data = leye.cs)
par(mfrow = c(2, 2))
plot(m)

# Plot residuals by categorical variable to assess homoscedasticity
boxplot(residuals(m) ~ leye.cs$Detection, 
        main = "Residuals by Detection",
        xlab = "Detection",
        ylab = "Residuals")


## Pectoral ~ Date ####
par(mfrow = c(2, 2))
m <- lm(PecSizeBest ~ DaysIntoSeason_S, data = leye)
plot(m)

# Plot 1: Residuals vs Fitted Values (Check for linearity and homoscedasticity)
plot(fitted(m), resid(m), 
     main = "Residuals vs Fitted", 
     xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")  # Add horizontal line at 0

# No violations of homoscedasticity

# Plot 2: Normal Q-Q Plot (Check for normality of residuals)
qqnorm(resid(m), main = "Normal Q-Q Plot")
qqline(resid(m), col = "red")

# No severe deviation from normality

# Plot 3: Scale-Location Plot (Check for homoscedasticity)
plot(fitted(m), sqrt(abs(resid(m))), 
     main = "Scale-Location Plot", 
     xlab = "Fitted Values", ylab = "Sqrt(|Residuals|)")
abline(h = 0, col = "red")

# No violations of homoscedasticity

# Plot 4: Cook's Distance (Check for influential points)
plot(cooks.distance(m), 
     main = "Cook's Distance", 
     ylab = "Cook's Distance", 
     xlab = "Index")
abline(h = 4/(length(resid(m)) - length(coef(m))), col = "red")  # Threshold for influential points

# 3 influential points but not too severe

## Pectoral ~ Date + Detection
m <- lm(PecSizeBest ~ DaysIntoSeason_S + Detection, data = leye.cs)
plot(m)
