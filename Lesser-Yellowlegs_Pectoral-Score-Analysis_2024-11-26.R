#---------------------------------------------#
#  Lesser Yellowlegs Pectoral Score Analysis  #
#          Created 11/26/2024                 #          
#         Modified 11/26/2024                 #
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

# Reorder factor variables
birds$Sex <- factor(birds$Sex,
                    levels = c("M", "F"),
                    labels = c("Male", "Female"))

birds$Detection <- as.factor(birds$Detection)

birds$Fat <- as.factor(birds$Fat)

# Subset data for Lesser Yellowlegs
leye <- subset(birds, Species %in% c("LesserYellowlegs"))

leye$Event <- factor(leye$Event, 
                     levels = c("Fall 2021", "Spring 2022", "Fall 2023"),
                     labels = c("Fall 2021", "Spring 2022", "Fall 2023"))

leye$PecScore <- factor(leye$PecScore,
                   levels = c("0", "1", "2", "3"))

# Create New Ag Category column
quartiles <- leye %>%
  summarise(
    Q1 = quantile(PercentAg, 0.25, na.rm = TRUE),
    Median = quantile(PercentAg, 0.5, na.rm = TRUE),
    Q3 = quantile(PercentAg, 0.75, na.rm = TRUE)
  )

leye <- leye %>% 
  mutate(AgIntensity = case_when(
    PercentAg <= quartiles$Q1 ~ "Low",          
    PercentAg > quartiles$Q1 & PercentAg <= quartiles$Q3 ~ "Moderate",  
    PercentAg > quartiles$Q3 ~ "High" 
  ))

leye$AgIntensity <- factor(leye$AgIntensity,
                           levels = c("Low", "Moderate", "High"))

# filter birds that only contain pectoral score information
leye <- leye %>% 
  filter(!is.na(PecScore))


# Standardize continuous variables
leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

# ---------------------------------------------------------------------------- #

# Data Visualization ####

ggplot(leye,
       aes(x = PecScore,
           y = PercentAg)) + 
  geom_boxplot()

plot(leye.cs$PercentAg, leye.cs$PecScore)

ggplot(leye, aes(x = PecScore, fill = AgIntensity)) +
  geom_bar(position = "dodge") +
  labs(title = "Count of PecScore by Agricultural Intensity",
       x = "PecScore Category", y = "Count") +
  theme_light()


# ---------------------------------------------------------------------------- #

# Lesser Yellowlegs Modeling ####

## Univariate Analysis: Date ####

m1 <- polr(PecScore ~ DaysIntoSeason_S, data = leye.cs)
m2 <- polr(PecScore ~ Julian, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Days into trapping season performs better (wt = 0.64)

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Neonicotinoids ####

# Cannot use concentrations due to extremely skewed data

m1 <- polr(PecScore ~ Detection, data = leye.cs)
m2 <- polr(PecScore ~ LogNeonic, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Both perform similarly (log(neonic) performs best)
# Use detection for consistency

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Agricultural Intensity ####

m1 <- polr(PecScore ~ PercentAg, data = leye.cs)
m2 <- polr(PecScore ~ PercentAg + I(PercentAg^2), data = leye.cs)
m3 <- polr(PecScore ~ AgIntensity, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# percentag with no transformation performs best

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Event and Capture Time ####

m1 <- polr(PecScore ~ Event * ts.sunrise, data = leye.cs)
m2 <- polr(PecScore ~ Event + ts.sunrise + I(ts.sunrise^2), data = leye.cs)
m3 <- polr(PecScore ~ Event + ts.sunrise, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model without interaction performs the worst.

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Detection & Sampling Event ####

# Cannot consider interaction due to perfect collinearity

m1 <- polr(PecScore ~ Event * Detection, data = leye.cs)
m2 <- polr(PecScore ~ Event + Detection, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# ---------------------------------------------------------------------------- #

# INTERACTIONS/TRANSFORMATIONS INCLUDED: None; PercentAg

# ---------------------------------------------------------------------------- #

## Single Models ####

m.global <- polr(PecScore ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + Detection + 
                   PercentAg, 
                 data = leye.cs)

m.null <- polr(PecScore ~ 1, data = leye.cs)

m1 <- polr(PecScore ~ Sex, data = leye.cs)
m2 <- polr(PecScore ~ Event, data = leye.cs)
m3 <- polr(PecScore ~ DaysIntoSeason_S, data = leye.cs)
m4 <- polr(PecScore ~ ts.sunrise, data = leye.cs)
m5 <- polr(PecScore ~ PercentAg, data = leye.cs)
m6 <- polr(PecScore ~ Detection, data = leye.cs)

## Additive Models ####

### Two additive combinations ####
m7 <- polr(PecScore ~ Sex + Event, data = leye.cs)
m8 <- polr(PecScore ~ Sex + DaysIntoSeason_S, data = leye.cs)
m9 <- polr(PecScore ~ Sex + ts.sunrise, data = leye.cs)
m10 <- polr(PecScore ~ Sex + PercentAg, data = leye.cs)
m11 <- polr(PecScore ~ Sex + Detection, data = leye.cs)
m12 <- polr(PecScore ~ Event + DaysIntoSeason_S, data = leye.cs)
m13 <- polr(PecScore ~ Event * ts.sunrise, data = leye.cs)
m14 <- polr(PecScore ~ Event + PercentAg, data = leye.cs)
m15 <- polr(PecScore ~ DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m16 <- polr(PecScore ~ DaysIntoSeason_S + PercentAg, data = leye.cs)
m17 <- polr(PecScore ~ DaysIntoSeason_S + Detection, data = leye.cs)
m18 <- polr(PecScore ~ ts.sunrise + PercentAg, data = leye.cs)
m19 <- polr(PecScore ~ ts.sunrise + Detection, data = leye.cs)
m20 <- polr(PecScore ~ PercentAg + Detection, data = leye.cs)

### Three additive combinations ####
m21 <- polr(PecScore ~ Sex + Event + DaysIntoSeason_S, data = leye.cs)
m22 <- polr(PecScore ~ Sex + Event * ts.sunrise, data = leye.cs)
m23 <- polr(PecScore ~ Sex + Event + PercentAg, data = leye.cs)
m24 <- polr(PecScore ~ Sex + Event + Detection, data = leye.cs)
m25 <- polr(PecScore ~ Sex + DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m26 <- polr(PecScore ~ Sex + DaysIntoSeason_S + PercentAg, data = leye.cs)
m27 <- polr(PecScore ~ Sex + DaysIntoSeason_S + Detection, data = leye.cs)
m28 <- polr(PecScore ~ Sex + ts.sunrise + PercentAg, data = leye.cs)
m29 <- polr(PecScore ~ Sex + ts.sunrise + Detection, data = leye.cs)
m30 <- polr(PecScore ~ Sex + PercentAg + Detection, data = leye.cs)
m31 <- polr(PecScore ~ Event + DaysIntoSeason_S + ts.sunrise + Event * ts.sunrise, data = leye.cs)
m32 <- polr(PecScore ~ Event + DaysIntoSeason_S + PercentAg, data = leye.cs)
m33 <- polr(PecScore ~ Event + DaysIntoSeason_S + Detection, data = leye.cs)
m34 <- polr(PecScore ~ Event * ts.sunrise + PercentAg, data = leye.cs)
m35 <- polr(PecScore ~ Event * ts.sunrise + Detection, data = leye.cs)
m36 <- polr(PecScore ~ Event + PercentAg + Detection, data = leye.cs)
m37 <- polr(PecScore ~ DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m38 <- polr(PecScore ~ DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m39 <- polr(PecScore ~ DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m40 <- polr(PecScore ~ ts.sunrise + PercentAg + Detection, data = leye.cs)

### Four additive combinations ####
m41 <- polr(PecScore ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + Event + ts.sunrise, data = leye.cs)
m42 <- polr(PecScore ~ Sex + Event + DaysIntoSeason_S + PercentAg, data = leye.cs)
m43 <- polr(PecScore ~ Sex + Event + DaysIntoSeason_S + Detection, data = leye.cs)
m44 <- polr(PecScore ~ Sex + Event * ts.sunrise + PercentAg, data = leye.cs)
m45 <- polr(PecScore ~ Sex + Event * ts.sunrise + Detection, data = leye.cs)
m46 <- polr(PecScore ~ Sex + Event + PercentAg + Detection, data = leye.cs)
m47 <- polr(PecScore ~ Sex + DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m48 <- polr(PecScore ~ Sex + DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m49 <- polr(PecScore ~ Sex + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m50 <- polr(PecScore ~ Sex + ts.sunrise + PercentAg + Detection, data = leye.cs)
m51 <- polr(PecScore ~ Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Event + ts.sunrise, data = leye.cs)
m52 <- polr(PecScore ~ Event + DaysIntoSeason_S + ts.sunrise + Detection + Event + ts.sunrise, data = leye.cs)
m53 <- polr(PecScore ~ Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m54 <- polr(PecScore ~ Event * ts.sunrise + PercentAg + Detection, data = leye.cs)
m55 <- polr(PecScore ~ DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)

### Five additive combinations ####

m56 <- polr(PecScore ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Event * ts.sunrise, data = leye.cs)
m57 <- polr(PecScore ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + Detection + Event * ts.sunrise, data = leye.cs)
m58 <- polr(PecScore ~ Sex + Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m59 <- polr(PecScore ~ Sex + Event * ts.sunrise + PercentAg + Detection, data = leye.cs)
m60 <- polr(PecScore ~ Sex + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)
m61 <- polr(PecScore ~ Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection + Event + ts.sunrise, data = leye.cs)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:61)

models <- mget(model_names)

models$m.null <- m.null
models$m.global <- m.global

model_names <- c(model_names, "m.null", "m.global")

aictab(models, modnames = model_names)

# ---------------------------------------------------------------------------- #

# Top Model Summaries ####

summary(m16)
confint(m16)

# Conclusion: Date is a significant predictor of pectoral score. Detection and agricultural intensity has no effect

# ---------------------------------------------------------------------------- #

# Plot Top Models ####


## PecScore ~ Julian + PercentAg
m <- polr(as.factor(PecScore) ~ Julian + PercentAg + I(PercentAg^2), data = leye)

# hold julian date constant
d <- expand.grid(PercentAg = seq(min(leye$PercentAg),
                                 max(leye$PercentAg), length = 100),
                 Julian = mean(leye$Julian))

d$yhat <- predict(m, newdata = d, type = "probs", se.fit = TRUE)

d <- cbind(d, as.data.frame(d$yhat))

colnames(d)[(ncol(d) - 5):ncol(d)] <- as.character(0:5)

# long format for easier plotting
d_long <- d %>%
  pivot_longer(
    cols = as.character(0:5), 
    names_to = "PecScore",           
    values_to = "Probability" 
  )

# Predicted Probabilities for Each Category of PecScore
ggplot(d_long, aes(x = PercentAg, y = Probability, color = PecScore)) +
  geom_line() +
  labs(x = "% Surrounding Agriculture", y = "Predicted Probability for Each PecScore Score") +
  theme_light() +
  theme(legend.title = element_blank())



## PecScore ~ Julian + AgIntensity ####

m <- polr(PecScore ~ Julian + AgIntensity, data = leye, Hess = TRUE)

d <- expand.grid(AgIntensity = c("Low", "Moderate", "High"),
                 Julian = mean(leye$Julian))

d <- cbind(d, predict(m, newdata = d, type = "probs"))

# transform data to long format
d_long <- d %>%
  pivot_longer(cols = `0`:`5`,   
               names_to = "PecScore", 
               values_to = "Probability")

d_long <- d_long %>%
  mutate(PecScore = factor(PecScore, levels = 0:5))

# plot
# cannot figure out how to extract or add standard error bars to figure...
# must use delta method or bootstrapping
ggplot(d_long, aes(x = PecScore, y = Probability)) +
  geom_point() +
  facet_wrap(~AgIntensity) +
  labs(title = "Distribution of Predicted Probabilities for Each Category of PecScore",
       x = "PecScore Score", y = "Predicted Probability") +
  theme_light() +
  theme(legend.title = element_blank())

# ---------------------------------------------------------------------------- #

## Single Models ####

m.global <- polr(PecScore ~ Sex + Event*ts.sunrise + Julian + Detection + 
                   AgIntensity, 
                 data = leye.cs)

m.null <- polr(PecScore ~ 1, data = leye.cs)

m1 <- polr(PecScore ~ Sex, data = leye.cs)
m2 <- polr(PecScore ~ Event, data = leye.cs)
m3 <- polr(PecScore ~ Julian, data = leye.cs)
m4 <- polr(PecScore ~ ts.sunrise, data = leye.cs)
m5 <- polr(PecScore ~ AgIntensity, data = leye.cs)
m6 <- polr(PecScore ~ Detection, data = leye.cs)

## Additive Models ####

### Two additive combinations ####
m7 <- polr(PecScore ~ Sex + Event, data = leye.cs)
m8 <- polr(PecScore ~ Sex + Julian, data = leye.cs)
m9 <- polr(PecScore ~ Sex + ts.sunrise, data = leye.cs)
m10 <- polr(PecScore ~ Sex + AgIntensity, data = leye.cs)
m11 <- polr(PecScore ~ Sex + Detection, data = leye.cs)
m12 <- polr(PecScore ~ Event + Julian, data = leye.cs)
m13 <- polr(PecScore ~ Event * ts.sunrise, data = leye.cs)
m14 <- polr(PecScore ~ Event + AgIntensity, data = leye.cs)
m15 <- polr(PecScore ~ Julian + ts.sunrise, data = leye.cs)
m16 <- polr(PecScore ~ Julian + AgIntensity, data = leye.cs)
m17 <- polr(PecScore ~ Julian + Detection, data = leye.cs)
m18 <- polr(PecScore ~ ts.sunrise + AgIntensity, data = leye.cs)
m19 <- polr(PecScore ~ ts.sunrise + Detection, data = leye.cs)
m20 <- polr(PecScore ~ AgIntensity + Detection, data = leye.cs)

### Three additive combinations ####
m21 <- polr(PecScore ~ Sex + Event + Julian, data = leye.cs)
m22 <- polr(PecScore ~ Sex + Event * ts.sunrise, data = leye.cs)
m23 <- polr(PecScore ~ Sex + Event + AgIntensity, data = leye.cs)
m24 <- polr(PecScore ~ Sex + Event + Detection, data = leye.cs)
m25 <- polr(PecScore ~ Sex + Julian + ts.sunrise, data = leye.cs)
m26 <- polr(PecScore ~ Sex + Julian + AgIntensity, data = leye.cs)
m27 <- polr(PecScore ~ Sex + Julian + Detection, data = leye.cs)
m28 <- polr(PecScore ~ Sex + ts.sunrise + AgIntensity, data = leye.cs)
m29 <- polr(PecScore ~ Sex + ts.sunrise + Detection, data = leye.cs)
m30 <- polr(PecScore ~ Sex + AgIntensity + Detection, data = leye.cs)
m31 <- polr(PecScore ~ Event + Julian + ts.sunrise + Event * ts.sunrise, data = leye.cs)
m32 <- polr(PecScore ~ Event + Julian + AgIntensity, data = leye.cs)
m33 <- polr(PecScore ~ Event + Julian + Detection, data = leye.cs)
m34 <- polr(PecScore ~ Event * ts.sunrise + AgIntensity, data = leye.cs)
m35 <- polr(PecScore ~ Event * ts.sunrise + Detection, data = leye.cs)
m36 <- polr(PecScore ~ Event + AgIntensity + Detection, data = leye.cs)
m37 <- polr(PecScore ~ Julian + ts.sunrise + AgIntensity, data = leye.cs)
m38 <- polr(PecScore ~ Julian + ts.sunrise + Detection, data = leye.cs)
m39 <- polr(PecScore ~ Julian + AgIntensity + Detection, data = leye.cs)
m40 <- polr(PecScore ~ ts.sunrise + AgIntensity + Detection, data = leye.cs)

### Four additive combinations ####
m41 <- polr(PecScore ~ Sex + Event + Julian + ts.sunrise + Event*ts.sunrise, data = leye.cs)
m42 <- polr(PecScore ~ Sex + Event + Julian + AgIntensity, data = leye.cs)
m43 <- polr(PecScore ~ Sex + Event + Julian + Detection, data = leye.cs)
m44 <- polr(PecScore ~ Sex + Event * ts.sunrise + AgIntensity, data = leye.cs)
m45 <- polr(PecScore ~ Sex + Event * ts.sunrise + Detection, data = leye.cs)
m46 <- polr(PecScore ~ Sex + Event + AgIntensity + Detection, data = leye.cs)
m47 <- polr(PecScore ~ Sex + Julian + ts.sunrise + AgIntensity, data = leye.cs)
m48 <- polr(PecScore ~ Sex + Julian + ts.sunrise + Detection, data = leye.cs)
m49 <- polr(PecScore ~ Sex + Julian + AgIntensity + Detection, data = leye.cs)
m50 <- polr(PecScore ~ Sex + ts.sunrise + AgIntensity + Detection, data = leye.cs)
m51 <- polr(PecScore ~ Event + Julian + ts.sunrise + AgIntensity + Event*ts.sunrise, data = leye.cs)
m52 <- polr(PecScore ~ Event + Julian + ts.sunrise + Detection + Event*ts.sunrise, data = leye.cs)
m53 <- polr(PecScore ~ Event + Julian + AgIntensity + Detection, data = leye.cs)
m54 <- polr(PecScore ~ Event * ts.sunrise + AgIntensity + Detection, data = leye.cs)
m55 <- polr(PecScore ~ Julian + ts.sunrise + AgIntensity + Detection, data = leye.cs)

### Five additive combinations ####

m56 <- polr(PecScore ~ Sex + Event + Julian + ts.sunrise + AgIntensity + Event * ts.sunrise, data = leye.cs)
m57 <- polr(PecScore ~ Sex + Event + Julian + ts.sunrise + Detection + Event * ts.sunrise, data = leye.cs)
m58 <- polr(PecScore ~ Sex + Event + Julian + AgIntensity + Detection, data = leye.cs)
m59 <- polr(PecScore ~ Sex + Event * ts.sunrise + AgIntensity + Detection, data = leye.cs)
m60 <- polr(PecScore ~ Sex + Julian + ts.sunrise + AgIntensity + Detection, data = leye.cs)
m61 <- polr(PecScore ~ Event + Julian + ts.sunrise + AgIntensity + Detection + Event*ts.sunrise, data = leye.cs)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:61)

models <- mget(model_names)

models$m.null <- m.null
models$m.global <- m.global

model_names <- c(model_names, "m.null", "m.global")

aictab(models, modnames = model_names)

# ---------------------------------------------------------------------------- #

# Top Model Summaries ####

summary(m16)
confint(m16)

# Julian date and agricultural intensity have a significant impact on PecScore.

summary(m39)
confint(m39)

# Julian date and agricultural intensity have a significant impact on PecScore.
