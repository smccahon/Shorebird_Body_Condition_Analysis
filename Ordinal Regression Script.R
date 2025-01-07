#----------------------------------#
#  Lesser Yellowlegs Fat Analysis  #
#         Created 12/02/2024       #          
#       Modified 12/02/2024        #
#----------------------------------#

# Load packages ####
library(tidyverse)
library(lme4)
library(MASS)
library(ordinal)

# MAIN Q: HOW DO I PLOT THIS MODEL?

# CUMULATIVE LINK MIXED MODELS (CLMM) ####

## Read data ####
birds <- read.csv("bird_fat_analysis.csv")

## Convert variables to a factor ####
birds$Fat.G <- factor(birds$Fat.G,
                      levels = c("Low", "Medium", "High"))

birds$Detection <- as.factor(birds$Detection)

birds$Species <- as.factor(birds$Species)

## Standardize continuous variables ####
birds.cs <- birds %>%
  mutate(across(where(is.numeric), scale))

## Top Model ####
m <- clmm(Fat.G ~ Event * ts.sunrise + DaysIntoSeason_S + Detection + (1|Species),
          data = birds.cs)

# View results ####
summary(m)
cbind(summary(m)$coefficients, confint(m))

## Plot Top Model ####

### Fat ~ Detection ####
d <- expand.grid(Event = c("Fall 2023"),
                 ts.sunrise = mean(birds.cs$ts.sunrise),
                 Detection = c("Detection", "Non-detection"),
                 Species = unique(birds.cs$Species))

# Scale predictors
d$ts.sunrise <- scale(d$ts.sunrise)

# Generate predictions
# Code doesn't work...
pred_probs <- predict(m, newdata = d, type = "probs", se.fit = TRUE)

d <- cbind(d, pred_probs)

#----------------------------------------------------#

# MAIN Q: HOW DO I ADD STANDARD ERRORS/CONFIDENCE INTERVALS?

# CUMULATIVE LINK MIXED MODELS (CLMM) ####

## Read data ####
leye <- read.csv("leye_fat_analysis.csv")

## Convert variables to a factor ####
leye$Fat.G <- factor(leye$Fat.G,
                      levels = c("Low", "Medium", "High"))

leye$Detection <- as.factor(leye$Detection)

## Standardize continuous variables ####
leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

## Fat ~ Julian + AgIntensity
m <- polr(Fat.G ~ Julian + AgIntensity, data = leye.cs)

# hold julian date constant
d <- expand.grid(AgIntensity = c("Low", "Moderate", "High"),
                 Julian = mean(leye.cs$Julian))

# Scale predictors
d$Julian <- scale(d$Julian, center = mean(leye.cs$Julian), scale = sd(leye.cs$Julian))

d <- cbind(d, predict(m, newdata = d, type = "p"))

# transform data to long format
d_long <- d %>%
  pivot_longer(cols = Low:High,   
               names_to = "Fat", 
               values_to = "Probability")

d_long <- d_long %>%
  mutate(Fat = factor(Fat, levels = c("Low", "Medium", "High")))

# Plot Model
# How do I add standard errors?
ggplot(d_long, aes(x = Fat, y = Probability)) +
  geom_point() +
  facet_wrap(~AgIntensity) +
  labs(title = "Distribution of Predicted Probabilities for Each Category of Fat",
       x = "Fat Score", y = "Predicted Probability") +
  theme_light() +
  theme(legend.title = element_blank())
