#----------------------------------#
#  Lesser Yellowlegs Fat Analysis  #
#          Created 11/26/2024      #          
#         Modified 12/02/2024      #
#----------------------------------#

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

leye$Fat <- factor(leye$Fat,
                   levels = c("0", "1", "2", "3", "4", "5"))

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

# filter birds that only contain fat information
leye <- leye %>% 
  filter(!is.na(Fat))

# combine fat categories to address small sample size
leye <- leye %>%
  mutate(Fat.G = case_when(
    Fat == 0 ~ "Low",       
    Fat == 1 | Fat == 2 | Fat == 3 ~ "Medium",  
    Fat == 4 | Fat == 5 ~ "High"       
  ))

leye$Fat.G <- factor(leye$Fat.G,
                   levels = c("Low", "Medium", "High"))

# Standardize continuous variables
leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

# Make Code for Tim
leye %>% 
  dplyr::select(Fat.G, AgIntensity, PercentAg, Julian, Detection, Individual,
                ts.sunrise, Species, DaysIntoSeason_S, Event) %>% 
  write.csv("leye_fat_analysis.csv")

# ---------------------------------------------------------------------------- #

# Data Visualization ####

my_theme <- theme_classic() +
  theme(text = element_text(size = 16),
        axis.title.y = element_text(margin = margin(r = 13)),
        axis.title.x = element_text(margin = margin(t = 13)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")

### ASG Presentation ####
ggplot(leye, aes(x = Fat.G, y = PercentAg)) + 
  geom_boxplot() +
  labs(x = "Fat Score Category", y = "Surrounding Agricultural Intensity (%)") +
  my_theme 


anova_result <- aov(PercentAg ~ Fat.G, data = leye)
summary(anova_result)

# leye$Fat <- as.numeric(leye$Fat)

# I don't think this is valid since fat is a factor?
ggplot(leye,
       aes(x = AgIntensity,
           y = Fat)) + geom_boxplot()

ggplot(leye, aes(x = Fat, fill = AgIntensity)) +
  geom_bar(position = "dodge") +
  labs(title = "Count of Fat by Agricultural Intensity",
       x = "Fat Category", y = "Count") +
  theme_light()


ggplot(leye, aes(x = Julian, y = Fat)) +
  geom_point()

ggplot(leye, aes(x = DaysIntoSeason_S, y = Fat)) +
  geom_point()


# ---------------------------------------------------------------------------- #

# Lesser Yellowlegs Modeling ####

## Univariate Analysis: Date ####

m1 <- polr(Fat.G ~ DaysIntoSeason_S, data = leye.cs)
m2 <- polr(Fat.G ~ Julian, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Julian performs better (wt = 0.82)
# But days into season is more valid (spring 2022 had fat scores almost completely of 0)

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Neonicotinoids ####

# Cannot use concentrations due to extremely skewed data

m1 <- polr(Fat.G ~ Detection, data = leye.cs)
m2 <- polr(Fat.G ~ LogNeonic, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Both perform similarly (log(neonic) performs best)
# Use detection for consistency

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Agricultural Intensity ####

m1 <- polr(Fat.G ~ PercentAg, data = leye.cs)
m2 <- polr(Fat.G ~ PercentAg + I(PercentAg^2), data = leye.cs)
m3 <- polr(Fat.G ~ AgIntensity, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Quadratic relationship performs the best

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Event and Capture Time ####

m1 <- polr(Fat.G ~ Event * ts.sunrise, data = leye.cs)
m2 <- polr(Fat.G ~ Event + ts.sunrise + I(ts.sunrise^2), data = leye.cs)
m3 <- polr(Fat.G ~ Event + ts.sunrise, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with interaction performs the worst.

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Detection & Sampling Event ####

m1 <- polr(Fat.G ~ Event * Detection, data = leye.cs)
m2 <- polr(Fat.G ~ Event + Detection, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model without interaction performs the best.

# ---------------------------------------------------------------------------- #

# INTERACTIONS/TRANSFORMATIONS INCLUDED: Event * Capture Time; PercentAg^2

# ---------------------------------------------------------------------------- #


## Single Models ####

m.global <- polr(Fat.G ~ Sex + Event*ts.sunrise + Julian + Detection + 
                   PercentAg + I(PercentAg^2), 
               data = leye.cs)

m.null <- polr(Fat.G ~ 1, data = leye.cs)

m1 <- polr(Fat.G ~ Sex, data = leye.cs)
m2 <- polr(Fat.G ~ Event, data = leye.cs)
m3 <- polr(Fat.G ~ Julian, data = leye.cs)
m4 <- polr(Fat.G ~ ts.sunrise, data = leye.cs)
m5 <- polr(Fat.G ~ PercentAg + I(PercentAg^2), data = leye.cs)
m6 <- polr(Fat.G ~ Detection, data = leye.cs)

## Additive Models ####

### Two additive combinations ####
m7 <- polr(Fat.G ~ Sex + Event, data = leye.cs)
m8 <- polr(Fat.G ~ Sex + Julian, data = leye.cs)
m9 <- polr(Fat.G ~ Sex + ts.sunrise, data = leye.cs)
m10 <- polr(Fat.G ~ Sex + PercentAg + I(PercentAg^2), data = leye.cs)
m11 <- polr(Fat.G ~ Sex + Detection, data = leye.cs)
m12 <- polr(Fat.G ~ Event + Julian, data = leye.cs)
m13 <- polr(Fat.G ~ Event * ts.sunrise, data = leye.cs)
m14 <- polr(Fat.G ~ Event + PercentAg + I(PercentAg^2), data = leye.cs)
m15 <- polr(Fat.G ~ Julian + ts.sunrise, data = leye.cs)
m16 <- polr(Fat.G ~ Julian + PercentAg + I(PercentAg^2), data = leye.cs)
m17 <- polr(Fat.G ~ Julian + Detection, data = leye.cs)
m18 <- polr(Fat.G ~ ts.sunrise + PercentAg + I(PercentAg^2), data = leye.cs)
m19 <- polr(Fat.G ~ ts.sunrise + Detection, data = leye.cs)
m20 <- polr(Fat.G ~ PercentAg + I(PercentAg^2) + Detection, data = leye.cs)

### Three additive combinations ####
m21 <- polr(Fat.G ~ Sex + Event + Julian, data = leye.cs)
m22 <- polr(Fat.G ~ Sex + Event * ts.sunrise, data = leye.cs)
m23 <- polr(Fat.G ~ Sex + Event + PercentAg + I(PercentAg^2), data = leye.cs)
m24 <- polr(Fat.G ~ Sex + Event + Detection, data = leye.cs)
m25 <- polr(Fat.G ~ Sex + Julian + ts.sunrise, data = leye.cs)
m26 <- polr(Fat.G ~ Sex + Julian + PercentAg + I(PercentAg^2), data = leye.cs)
m27 <- polr(Fat.G ~ Sex + Julian + Detection, data = leye.cs)
m28 <- polr(Fat.G ~ Sex + ts.sunrise + PercentAg + I(PercentAg^2), data = leye.cs)
m29 <- polr(Fat.G ~ Sex + ts.sunrise + Detection, data = leye.cs)
m30 <- polr(Fat.G ~ Sex + PercentAg + I(PercentAg^2) + Detection, data = leye.cs)
m31 <- polr(Fat.G ~ Event + Julian + ts.sunrise + Event * ts.sunrise, data = leye.cs)
m32 <- polr(Fat.G ~ Event + Julian + PercentAg + I(PercentAg^2), data = leye.cs)
m33 <- polr(Fat.G ~ Event + Julian + Detection, data = leye.cs)
m34 <- polr(Fat.G ~ Event * ts.sunrise + PercentAg + I(PercentAg^2), data = leye.cs)
m35 <- polr(Fat.G ~ Event * ts.sunrise + Detection, data = leye.cs)
m36 <- polr(Fat.G ~ Event + PercentAg + I(PercentAg^2) + Detection, data = leye.cs)
m37 <- polr(Fat.G ~ Julian + ts.sunrise + PercentAg + I(PercentAg^2), data = leye.cs)
m38 <- polr(Fat.G ~ Julian + ts.sunrise + Detection, data = leye.cs)
m39 <- polr(Fat.G ~ Julian + PercentAg + I(PercentAg^2) + Detection, data = leye.cs)
m40 <- polr(Fat.G ~ ts.sunrise + PercentAg + I(PercentAg^2) + Detection, data = leye.cs)

### Four additive combinations ####
m41 <- polr(Fat.G ~ Sex + Event + Julian + ts.sunrise + Event*ts.sunrise, data = leye.cs)
m42 <- polr(Fat.G ~ Sex + Event + Julian + PercentAg + I(PercentAg^2), data = leye.cs)
m43 <- polr(Fat.G ~ Sex + Event + Julian + Detection, data = leye.cs)
m44 <- polr(Fat.G ~ Sex + Event * ts.sunrise + PercentAg + I(PercentAg^2), data = leye.cs)
m45 <- polr(Fat.G ~ Sex + Event * ts.sunrise + Detection, data = leye.cs)
m46 <- polr(Fat.G ~ Sex + Event + PercentAg + I(PercentAg^2) + Detection, data = leye.cs)
m47 <- polr(Fat.G ~ Sex + Julian + ts.sunrise + PercentAg + I(PercentAg^2), data = leye.cs)
m48 <- polr(Fat.G ~ Sex + Julian + ts.sunrise + Detection, data = leye.cs)
m49 <- polr(Fat.G ~ Sex + Julian + PercentAg + I(PercentAg^2) + Detection, data = leye.cs)
m50 <- polr(Fat.G ~ Sex + ts.sunrise + PercentAg + I(PercentAg^2) + Detection, data = leye.cs)
m51 <- polr(Fat.G ~ Event + Julian + ts.sunrise + PercentAg + I(PercentAg^2) + Event*ts.sunrise, data = leye.cs)
m52 <- polr(Fat.G ~ Event + Julian + ts.sunrise + Detection + Event*ts.sunrise, data = leye.cs)
m53 <- polr(Fat.G ~ Event + Julian + PercentAg + I(PercentAg^2) + Detection, data = leye.cs)
m54 <- polr(Fat.G ~ Event * ts.sunrise + PercentAg + I(PercentAg^2) + Detection, data = leye.cs)
m55 <- polr(Fat.G ~ Julian + ts.sunrise + PercentAg + I(PercentAg^2) + Detection, data = leye.cs)

### Five additive combinations ####

m56 <- polr(Fat.G ~ Sex + Event + Julian + ts.sunrise + PercentAg + I(PercentAg^2) + Event * ts.sunrise, data = leye.cs)
m57 <- polr(Fat.G ~ Sex + Event + Julian + ts.sunrise + Detection + Event * ts.sunrise, data = leye.cs)
m58 <- polr(Fat.G ~ Sex + Event + Julian + PercentAg + I(PercentAg^2) + Detection, data = leye.cs)
m59 <- polr(Fat.G ~ Sex + Event * ts.sunrise + PercentAg + I(PercentAg^2) + Detection, data = leye.cs)
m60 <- polr(Fat.G ~ Sex + Julian + ts.sunrise + PercentAg + I(PercentAg^2) + Detection, data = leye.cs)
m61 <- polr(Fat.G ~ Event + Julian + ts.sunrise + PercentAg + I(PercentAg^2) + Detection + Event*ts.sunrise, data = leye.cs)

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

# Conclusion: Positive relationship between fat and Julian day
# Significant but non linear relationship between fat and agricultural intensity

# ---------------------------------------------------------------------------- #

# Plot Top Models ####


## Fat ~ Julian + PercentAg
m <- polr(as.factor(Fat) ~ Julian + PercentAg + I(PercentAg^2), data = leye)

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
    names_to = "Fat",           
    values_to = "Probability" 
  )

# Predicted Probabilities for Each Category of Fat
ggplot(d_long, aes(x = PercentAg, y = Probability, color = Fat)) +
  geom_line() +
  labs(x = "% Surrounding Agriculture", y = "Predicted Probability for Each Fat Score") +
  theme_light() +
  theme(legend.title = element_blank())



## Fat ~ Julian + AgIntensity ####

m <- polr(Fat ~ Julian + AgIntensity, data = leye, Hess = TRUE)

d <- expand.grid(AgIntensity = c("Low", "Moderate", "High"),
                 Julian = mean(leye$Julian))

d <- cbind(d, predict(m, newdata = d, type = "p"))

# transform data to long format
d_long <- d %>%
  pivot_longer(cols = `0`:`5`,   
               names_to = "Fat", 
               values_to = "Probability")

d_long <- d_long %>%
  mutate(Fat = factor(Fat, levels = 0:5))

# plot
# cannot figure out how to extract or add standard error bars to figure...
# must use delta method or bootstrapping
# something is also wrong because 0 probability for 0? not right
ggplot(d_long, aes(x = Fat, y = Probability)) +
  geom_point() +
  facet_wrap(~AgIntensity) +
  labs(title = "Distribution of Predicted Probabilities for Each Category of Fat",
       x = "Fat Score", y = "Predicted Probability") +
  theme_light() +
  theme(legend.title = element_blank())


# why zero?

table(leye$Fat)

# ---------------------------------------------------------------------------- #

## Single Models ####

m.global <- polr(Fat.G ~ Sex + Event*ts.sunrise + Julian + Detection + 
                   AgIntensity, 
                 data = leye.cs)

m.null <- polr(Fat.G ~ 1, data = leye.cs)

m1 <- polr(Fat.G ~ Sex, data = leye.cs)
m2 <- polr(Fat.G ~ Event, data = leye.cs)
m3 <- polr(Fat.G ~ Julian, data = leye.cs)
m4 <- polr(Fat.G ~ ts.sunrise, data = leye.cs)
m5 <- polr(Fat.G ~ AgIntensity, data = leye.cs)
m6 <- polr(Fat.G ~ Detection, data = leye.cs)

## Additive Models ####

### Two additive combinations ####
m7 <- polr(Fat.G ~ Sex + Event, data = leye.cs)
m8 <- polr(Fat.G ~ Sex + Julian, data = leye.cs)
m9 <- polr(Fat.G ~ Sex + ts.sunrise, data = leye.cs)
m10 <- polr(Fat.G ~ Sex + AgIntensity, data = leye.cs)
m11 <- polr(Fat.G ~ Sex + Detection, data = leye.cs)
m12 <- polr(Fat.G ~ Event + Julian, data = leye.cs)
m13 <- polr(Fat.G ~ Event * ts.sunrise, data = leye.cs)
m14 <- polr(Fat.G ~ Event + AgIntensity, data = leye.cs)
m15 <- polr(Fat.G ~ Julian + ts.sunrise, data = leye.cs)
m16 <- polr(Fat.G ~ Julian + AgIntensity, data = leye.cs)
m17 <- polr(Fat.G ~ Julian + Detection, data = leye.cs)
m18 <- polr(Fat.G ~ ts.sunrise + AgIntensity, data = leye.cs)
m19 <- polr(Fat.G ~ ts.sunrise + Detection, data = leye.cs)
m20 <- polr(Fat.G ~ AgIntensity + Detection, data = leye.cs)

### Three additive combinations ####
m21 <- polr(Fat.G ~ Sex + Event + Julian, data = leye.cs)
m22 <- polr(Fat.G ~ Sex + Event * ts.sunrise, data = leye.cs)
m23 <- polr(Fat.G ~ Sex + Event + AgIntensity, data = leye.cs)
m24 <- polr(Fat.G ~ Sex + Event + Detection, data = leye.cs)
m25 <- polr(Fat.G ~ Sex + Julian + ts.sunrise, data = leye.cs)
m26 <- polr(Fat.G ~ Sex + Julian + AgIntensity, data = leye.cs)
m27 <- polr(Fat.G ~ Sex + Julian + Detection, data = leye.cs)
m28 <- polr(Fat.G ~ Sex + ts.sunrise + AgIntensity, data = leye.cs)
m29 <- polr(Fat.G ~ Sex + ts.sunrise + Detection, data = leye.cs)
m30 <- polr(Fat.G ~ Sex + AgIntensity + Detection, data = leye.cs)
m31 <- polr(Fat.G ~ Event + Julian + ts.sunrise + Event * ts.sunrise, data = leye.cs)
m32 <- polr(Fat.G ~ Event + Julian + AgIntensity, data = leye.cs)
m33 <- polr(Fat.G ~ Event + Julian + Detection, data = leye.cs)
m34 <- polr(Fat.G ~ Event * ts.sunrise + AgIntensity, data = leye.cs)
m35 <- polr(Fat.G ~ Event * ts.sunrise + Detection, data = leye.cs)
m36 <- polr(Fat.G ~ Event + AgIntensity + Detection, data = leye.cs)
m37 <- polr(Fat.G ~ Julian + ts.sunrise + AgIntensity, data = leye.cs)
m38 <- polr(Fat.G ~ Julian + ts.sunrise + Detection, data = leye.cs)
m39 <- polr(Fat.G ~ Julian + AgIntensity + Detection, data = leye.cs)
m40 <- polr(Fat.G ~ ts.sunrise + AgIntensity + Detection, data = leye.cs)

### Four additive combinations ####
m41 <- polr(Fat.G ~ Sex + Event + Julian + ts.sunrise + Event*ts.sunrise, data = leye.cs)
m42 <- polr(Fat.G ~ Sex + Event + Julian + AgIntensity, data = leye.cs)
m43 <- polr(Fat.G ~ Sex + Event + Julian + Detection, data = leye.cs)
m44 <- polr(Fat.G ~ Sex + Event * ts.sunrise + AgIntensity, data = leye.cs)
m45 <- polr(Fat.G ~ Sex + Event * ts.sunrise + Detection, data = leye.cs)
m46 <- polr(Fat.G ~ Sex + Event + AgIntensity + Detection, data = leye.cs)
m47 <- polr(Fat.G ~ Sex + Julian + ts.sunrise + AgIntensity, data = leye.cs)
m48 <- polr(Fat.G ~ Sex + Julian + ts.sunrise + Detection, data = leye.cs)
m49 <- polr(Fat.G ~ Sex + Julian + AgIntensity + Detection, data = leye.cs)
m50 <- polr(Fat.G ~ Sex + ts.sunrise + AgIntensity + Detection, data = leye.cs)
m51 <- polr(Fat.G ~ Event + Julian + ts.sunrise + AgIntensity + Event*ts.sunrise, data = leye.cs)
m52 <- polr(Fat.G ~ Event + Julian + ts.sunrise + Detection + Event*ts.sunrise, data = leye.cs)
m53 <- polr(Fat.G ~ Event + Julian + AgIntensity + Detection, data = leye.cs)
m54 <- polr(Fat.G ~ Event * ts.sunrise + AgIntensity + Detection, data = leye.cs)
m55 <- polr(Fat.G ~ Julian + ts.sunrise + AgIntensity + Detection, data = leye.cs)

### Five additive combinations ####

m56 <- polr(Fat.G ~ Sex + Event + Julian + ts.sunrise + AgIntensity + Event * ts.sunrise, data = leye.cs)
m57 <- polr(Fat.G ~ Sex + Event + Julian + ts.sunrise + Detection + Event * ts.sunrise, data = leye.cs)
m58 <- polr(Fat.G ~ Sex + Event + Julian + AgIntensity + Detection, data = leye.cs)
m59 <- polr(Fat.G ~ Sex + Event * ts.sunrise + AgIntensity + Detection, data = leye.cs)
m60 <- polr(Fat.G ~ Sex + Julian + ts.sunrise + AgIntensity + Detection, data = leye.cs)
m61 <- polr(Fat.G ~ Event + Julian + ts.sunrise + AgIntensity + Detection + Event*ts.sunrise, data = leye.cs)

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

# Julian date and agricultural intensity have a significant impact on fat.

summary(m39)
confint(m39)

# Julian date and agricultural intensity have a significant impact on fat.


# ---------------------------------------------------------------------------- #

# Plot Top Model ####

## Fat ~ Julian + AgIntensity
m <- polr(Fat.G ~ Julian + AgIntensity, data = leye)

# hold julian date constant
d <- expand.grid(AgIntensity = c("Low", "Moderate", "High"),
                 Julian = mean(leye$Julian))

d <- cbind(d, predict(m, newdata = d, type = "p"))

# transform data to long format
d_long <- d %>%
  pivot_longer(cols = Low:High,   
               names_to = "Fat", 
               values_to = "Probability")

d_long <- d_long %>%
  mutate(Fat = factor(Fat, levels = c("Low", "Moderate", "High")))

# plot
# cannot figure out how to extract or add standard error bars to figure...
# must use delta method or bootstrapping
# something is also wrong because 0 probability for 0? not right
ggplot(d_long, aes(x = Fat, y = Probability)) +
  geom_point() +
  facet_wrap(~AgIntensity) +
  labs(title = "Distribution of Predicted Probabilities for Each Category of Fat",
       x = "Fat Score", y = "Predicted Probability") +
  theme_light() +
  theme(legend.title = element_blank())


## Fat ~ Julian + PercentAg
m <- polr(Fat.G ~ Julian + PercentAg + I(PercentAg^2), data = leye)

# hold julian date constant
d <- expand.grid(PercentAg = seq(min(leye$PercentAg),
                              max(leye$PercentAg), length = 100),
                 Julian = mean(leye$Julian))

d <- cbind(d, predict(m, newdata = d, type = "p"))

d_long <- d %>%
  pivot_longer(cols = Low:High,   
               names_to = "Fat", 
               values_to = "Probability")

d_long <- d_long %>%
  mutate(Fat = factor(Fat, levels = c("Low", "Moderate", "High")))


# Predicted Probabilities for Each Category of Fat
ggplot(d_long, aes(x = PercentAg, y = Probability, color = Fat)) +
  geom_line() +
  labs(x = "% Surrounding Agriculture", y = "Predicted Probability for Each Fat Score") +
  theme_light() +
  theme(legend.title = element_blank())



## Top Model ####
m <- polr(Fat.G ~ Julian + PercentAg + I(PercentAg^2), data = leye.cs)

# View results ####
summary(m)
confint(m)

#----------------------------------------------#

## Plot Top Model ####

### Fat ~ PercentAg ####

# Create new data
d <- expand.grid(PercentAg = seq(min(leye.cs$PercentAg),
                                 max(leye.cs$PercentAg), length = 100),
                 Julian = mean(leye.cs$Julian))

# Scale predictors
d$PercentAg <- scale(d$PercentAg, center = mean(leye.cs$PercentAg), scale = sd(leye.cs$PercentAg))
d$Julian <- scale(d$Julian, center = mean(leye.cs$Julian), scale = sd(leye.cs$Julian))

# Generate predictions
pred_probs <- predict(m, newdata = d, type = "probs", se.fit = TRUE)

d <- cbind(d, pred_probs)

# Bind the predictions and standard errors to the data frame
d <- cbind(d, pred_probs$fit)  # Predicted probabilities
d <- cbind(d, pred_probs$se.fit)  # Standard errors

# long format for easier plotting
# transform data to long format
d_long <- d %>%
  pivot_longer(cols = c("Low", "Medium", "High"),
               names_to = "Fat", 
               values_to = "Probability")

d_long <- d_long %>%
  mutate(Fat = factor(Fat, levels = c("Low", "Medium", "High")))

# Plot model
ggplot(d_long, aes(x = PercentAg, y = Probability, 
                   color = Fat)) +
  geom_line() + 
  labs(x = "Percent Agriculture", y = "Predicted Probability", color = "Fat Level") +
  theme_minimal() +
  theme(legend.position = "top")

