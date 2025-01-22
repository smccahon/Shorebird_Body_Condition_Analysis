#-----------------------------------------#
#  Lesser Yellowlegs Metabolite Analysis  #
#          Created 11/20/2024             #          
#         Modified 01/22/2025             #
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

# filter birds that only contain metabolite information
 leye <- leye %>% 
   filter(!is.na(Tri) & !is.na(Beta))
 
# Standardize continuous variables
leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

# ---------------------------------------------------------------------------- #

# Perform PCA ####

# Subset the dataset to only include 'Tri' and 'Beta'
leye_subset <- leye.cs[, c("Tri", "Beta")]

# Remove rows with NAs for PCA
leye_subset_clean <- leye_subset[complete.cases(leye_subset), ]

# Run PCA on the cleaned data
pca_result <- prcomp(leye_subset_clean, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# View the PCA scores (principal components)
pca_scores <- pca_result$x

# Merge PCA scores back to the original dataset
# Create a data frame to store the PCA scores for the rows with no missing data
leye.cs$PC1 <- NA  # Initialize with NA values
leye.cs$PC2 <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
leye.cs[complete.cases(leye_subset), "PC1"] <- pca_scores[, 1]
leye.cs[complete.cases(leye_subset), "PC2"] <- pca_scores[, 2]

# View the updated dataset with PCA scores
print(leye.cs)

# Merge PCA scores back to the original LEYE unstandardized dataset
# Create a data frame to store the PCA scores for the rows with no missing data
leye$PC1 <- NA  # Initialize with NA values
leye$PC2 <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
leye[complete.cases(leye_subset), "PC1"] <- pca_scores[, 1]
leye[complete.cases(leye_subset), "PC2"] <- pca_scores[, 2]

# View the updated dataset with PCA scores
print(leye)

# View the first few rows of the updated data frame
head(leye)

# ---------------------------------------------------------------------------- #

# Data Visualization ####

# Fattening Index ~ Detection (no relationship)
ggplot(leye.cs, 
       aes(x = Detection, y = PC1, fill = Detection)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(color = Detection),
              position = position_jitterdodge(jitter.width = 0.2),  
              size = 2, 
              alpha = 0.8) + 
  theme_light() +
  labs(x = "Neonicotinoid Detection", 
       y = "Fattening Index") +  
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

# Uric Acid ~ PercentAg (no clear relationship)
ggplot(leye, 
       aes(x = PercentAg, y = Uric)) + 
  geom_point(size = 2, alpha = 0.4, color = "black") +
  labs(x = "Surrounding Agricultural Intensity", 
       y = "LEYE Uric Acid Level (umol/L)") +  
  theme_classic() +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.position = "none",
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = scales::percent)


# Fattening index ~ Percent Ag ####
ggplot(leye, 
       aes(x = PercentAg, y = PC1)) + 
  geom_point(size = 2) + 
  labs(x = "Surrounding Agricultural Intensity", 
       y = "Lesser Yellowlegs Fattening Index") +  
  theme_classic() +
  theme(text = element_text(size = 20),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed",
             size = 1) +
  geom_smooth(method = "lm", color = "black", size = 1)  +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = scales::percent)

# ---------------------------------------------------------------------------- #

# Lesser Yellowlegs Modeling ####

## Univariate Analysis: Date ####

m1 <- lm(PC1 ~ DaysIntoSeason_S, data = leye.cs)
m2 <- lm(PC1 ~ Julian, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Julian performs better (wt = 0.54) but use date into trapping season for consistency

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Neonicotinoids ####

m1 <- lm(PC1 ~ Detection, data = leye.cs)
m2 <- lm(PC1 ~ OverallNeonic, data = leye.cs)
m3 <- lm(PC1 ~ LogNeonic, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# All perform similarly (detection performs best)
# Use detection for consistency

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Agricultural Intensity ####
# Can't use ag category for LEYE because there are only three birds in the low category

m1 <- lm(PC1 ~ PercentAg, data = leye.cs)
m2 <- lm(PC1 ~ PercentAg + I(PercentAg^2), data = leye.cs)
m3 <- lm(PC1 ~ LogAg, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Transformation does not significantly improve model fit. Use regular PercentAg

summary(m1)
confint(m1)
summary(m3)
confint(m3)

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Event and Capture Time ####

m1 <- lm(PC1 ~ Event * ts.sunrise, data = leye.cs)
m2 <- lm(PC1 ~ Event + ts.sunrise + I(ts.sunrise^2), data = leye.cs)
m3 <- lm(PC1 ~ Event + ts.sunrise, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model without interaction performs the best.
# Model with interaction performs the worst.

# ---------------------------------------------------------------------------- #

## Best Time Variable: Time or Time^2 ####

m1 <- lm(PC1 ~ ts.sunrise, data = leye.cs)
m2 <- lm(PC1 ~ ts.sunrise + I(ts.sunrise^2), data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with linear transformation performed best


# ---------------------------------------------------------------------------- #

# INTERACTIONS INCLUDED: NONE

# ---------------------------------------------------------------------------- #

## Single Models ####

m.global <- lm(PC1 ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + Detection + PercentAg + Age, 
               data = leye.cs)

m.null <- lm(PC1 ~ 1, data = leye.cs)

m1 <- lm(PC1 ~ Sex, data = leye.cs)
m2 <- lm(PC1 ~ Event, data = leye.cs)
m3 <- lm(PC1 ~ DaysIntoSeason_S, data = leye.cs)
m4 <- lm(PC1 ~ ts.sunrise, data = leye.cs)
m5 <- lm(PC1 ~ PercentAg, data = leye.cs)
m6 <- lm(PC1 ~ Detection, data = leye.cs)
m7 <- lm(PC1 ~ Age, data = leye.cs)

## Additive Models ####

### Two additive combinations ####
m8 <- lm(PC1 ~ Sex + Event, data = leye.cs)
m9 <- lm(PC1 ~ Sex + DaysIntoSeason_S, data = leye.cs)
m10 <- lm(PC1 ~ Sex + ts.sunrise, data = leye.cs)
m11 <- lm(PC1 ~ Sex + PercentAg, data = leye.cs)
m12 <- lm(PC1 ~ Sex + Detection, data = leye.cs)
m13 <- lm(PC1 ~ Event + DaysIntoSeason_S, data = leye.cs)
m14 <- lm(PC1 ~ Event + ts.sunrise, data = leye.cs)
m15 <- lm(PC1 ~ Event + PercentAg, data = leye.cs)
m16 <- lm(PC1 ~ DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m17 <- lm(PC1 ~ DaysIntoSeason_S + PercentAg, data = leye.cs)
m18 <- lm(PC1 ~ DaysIntoSeason_S + Detection, data = leye.cs)
m19 <- lm(PC1 ~ ts.sunrise + PercentAg, data = leye.cs)
m20 <- lm(PC1 ~ ts.sunrise + Detection, data = leye.cs)
m21 <- lm(PC1 ~ PercentAg + Detection, data = leye.cs)
m22 <- lm(PC1 ~ Age + Event, data = leye.cs)
m23 <- lm(PC1 ~ Age + DaysIntoSeason_S, data = leye.cs)
m24 <- lm(PC1 ~ Age + ts.sunrise, data = leye.cs)
m25 <- lm(PC1 ~ Age + PercentAg, data = leye.cs)
m26 <- lm(PC1 ~ Age + Detection, data = leye.cs)

### Three additive combinations ####
m27 <- lm(PC1 ~ Sex + Event + DaysIntoSeason_S, data = leye.cs)
m28 <- lm(PC1 ~ Sex + Event + ts.sunrise, data = leye.cs)
m29 <- lm(PC1 ~ Sex + Event + PercentAg, data = leye.cs)
m30 <- lm(PC1 ~ Sex + Event + Detection, data = leye.cs)
m31 <- lm(PC1 ~ Sex + DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m32 <- lm(PC1 ~ Sex + DaysIntoSeason_S + PercentAg, data = leye.cs)
m33 <- lm(PC1 ~ Sex + DaysIntoSeason_S + Detection, data = leye.cs)
m34 <- lm(PC1 ~ Sex + ts.sunrise + PercentAg, data = leye.cs)
m35 <- lm(PC1 ~ Sex + ts.sunrise + Detection, data = leye.cs)
m36 <- lm(PC1 ~ Sex + PercentAg + Detection, data = leye.cs)
m37 <- lm(PC1 ~ Event + DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m38 <- lm(PC1 ~ Event + DaysIntoSeason_S + PercentAg, data = leye.cs)
m39 <- lm(PC1 ~ Event + DaysIntoSeason_S + Detection, data = leye.cs)
m40 <- lm(PC1 ~ Event + ts.sunrise + PercentAg, data = leye.cs)
m41 <- lm(PC1 ~ Event + ts.sunrise + Detection, data = leye.cs)
m42 <- lm(PC1 ~ Event + PercentAg + Detection, data = leye.cs)
m43 <- lm(PC1 ~ DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m44 <- lm(PC1 ~ DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m45 <- lm(PC1 ~ DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m46 <- lm(PC1 ~ ts.sunrise + PercentAg + Detection, data = leye.cs)
m47 <- lm(PC1 ~ Age + Event + DaysIntoSeason_S, data = leye.cs)
m48 <- lm(PC1 ~ Age + Event + ts.sunrise, data = leye.cs)
m49 <- lm(PC1 ~ Age + Event + PercentAg, data = leye.cs)
m50 <- lm(PC1 ~ Age + Event + Detection, data = leye.cs)
m51 <- lm(PC1 ~ Age + DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m52 <- lm(PC1 ~ Age + DaysIntoSeason_S + PercentAg, data = leye.cs)
m53 <- lm(PC1 ~ Age + DaysIntoSeason_S + Detection, data = leye.cs)
m54 <- lm(PC1 ~ Age + ts.sunrise + PercentAg, data = leye.cs)
m55 <- lm(PC1 ~ Age + ts.sunrise + Detection, data = leye.cs)
m56 <- lm(PC1 ~ Age + PercentAg + Detection, data = leye.cs)

### Four additive combinations ####
m57 <- lm(PC1 ~ Sex + Event + DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m58 <- lm(PC1 ~ Sex + Event + DaysIntoSeason_S + PercentAg, data = leye.cs)
m59 <- lm(PC1 ~ Sex + Event + DaysIntoSeason_S + Detection, data = leye.cs)
m60 <- lm(PC1 ~ Sex + Event + ts.sunrise + PercentAg, data = leye.cs)
m61 <- lm(PC1 ~ Sex + Event + ts.sunrise + Detection, data = leye.cs)
m62 <- lm(PC1 ~ Sex + Event + PercentAg + Detection, data = leye.cs)
m63 <- lm(PC1 ~ Sex + DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m64 <- lm(PC1 ~ Sex + DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m65 <- lm(PC1 ~ Sex + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m66 <- lm(PC1 ~ Sex + ts.sunrise + PercentAg + Detection, data = leye.cs)
m67 <- lm(PC1 ~ Event + DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m68 <- lm(PC1 ~ Event + DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m69 <- lm(PC1 ~ Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m70 <- lm(PC1 ~ Event + ts.sunrise + PercentAg + Detection, data = leye.cs)
m71 <- lm(PC1 ~ DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)
m72 <- lm(PC1 ~ Age + Event + DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m73 <- lm(PC1 ~ Age + Event + DaysIntoSeason_S + PercentAg, data = leye.cs)
m74 <- lm(PC1 ~ Age + Event + DaysIntoSeason_S + Detection, data = leye.cs)
m75 <- lm(PC1 ~ Age + Event + ts.sunrise + PercentAg, data = leye.cs)
m76 <- lm(PC1 ~ Age + Event + ts.sunrise + Detection, data = leye.cs)
m77 <- lm(PC1 ~ Age + Event + PercentAg + Detection, data = leye.cs)
m78 <- lm(PC1 ~ Age + DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m79 <- lm(PC1 ~ Age + DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m80 <- lm(PC1 ~ Age + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m81 <- lm(PC1 ~ Age + ts.sunrise + PercentAg + Detection, data = leye.cs)

### Five additive combinations ####

m82 <- lm(PC1 ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m83 <- lm(PC1 ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m84 <- lm(PC1 ~ Sex + Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m85 <- lm(PC1 ~ Sex + Event + ts.sunrise + PercentAg + Detection, data = leye.cs)
m86 <- lm(PC1 ~ Sex + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)
m87 <- lm(PC1 ~ Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)
m88 <- lm(PC1 ~ Age + Event + DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m89 <- lm(PC1 ~ Age + Event + DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m90 <- lm(PC1 ~ Age + Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m91 <- lm(PC1 ~ Age + Event + ts.sunrise + PercentAg + Detection, data = leye.cs)
m92 <- lm(PC1 ~ Age + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:92)

models <- mget(model_names)

models$m.null <- m.null
models$m.global <- m.global

model_names <- c(model_names, "m.null", "m.global")

aictab(models, modnames = model_names)

# NEW SUMMARIES 2025-01-20 ####
confint(m18) # time and ag significant (negative)...
confint(m37) # time and ag significant (negative)...
confint(m28) # time and ag significant (negative)...

# NEW SUMMARIES 2025-01-22 ####
confint(m19) # time and ag significant (negative)...
confint(m43) # time and ag significant (negative)...
confint(m34) # time and ag significant (negative)...
confint(m46) # time and ag significant (negative)...

# ---------------------------------------------------------------------------- #

# Top Model Summaries ####

m4 <- lm(PC1 ~ ts.sunrise, data = leye.cs)
cbind(summary(m4)$coefficients, confint(m4))

# Capture time has a significant impact on fattening index. 

cbind(summary(m9)$coefficients, confint(m9))

# Capture time has a significant impact on fattening index. 
# Sex does not have significant impact on fattening index.

### Goodness of fit ####
summary(m4)
par(mfrow = c(2,2))
plot(m4)

### Are any other models good? ####
confint(m18)
confint(m19)
confint(m15)

# ---------------------------------------------------------------------------- #

# Graph Top Models ####

## Fattening Index ~ Time + Ag ####
m <- lm(PC1 ~ ts.sunrise, data = leye)

# holding ag intensity constant
d <- expand.grid(ts.sunrise = seq(min(leye$ts.sunrise), 
                                 max(leye$ts.sunrise),
                                 length.out = 1000))

d <- cbind(d, predict(m, newdata = d, interval = "confidence"))

ggplot(d, aes(x = ts.sunrise, y = fit)) + 
  geom_line(aes(color = "tan2"), size = 1, show.legend = FALSE) +  
  theme_light() +
  labs(x = "Time of Capture since Sunrise (min)",
       y = "Lesser Yellowlegs Fattening Index") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = LogNeonic),  
              alpha = 0.25, color = NA, show.legend = FALSE,
              fill = "tan1") +
  geom_point(data = leye, aes(x = ts.sunrise, y = PC1), 
             color = "tan2") +
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 0, y = max(d$fit), 
           label = "Sunrise", angle = 90, 
           vjust = -0.5, hjust = 0.5,
           size = 5) + 
  geom_hline(size = 1, color = "red", linetype = "dashed", yintercept = 0)

## Fattening Index ~ PercentAg ####
m <- lm(PC1 ~ PercentAg, data = leye)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg), 
                                  max(leye$PercentAg),
                                  length.out = 1000))

d <- cbind(d, predict(m, newdata = d, interval = "confidence"))

ggplot(d, aes(x = PercentAg, y = fit)) + 
  geom_line(aes(color = "tan2"), size = 1, show.legend = FALSE) +  
  theme_light() +
  labs(x = "% Surrounding Agriculture",
       y = "Lesser Yellowlegs Fattening Index") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = LogNeonic),  
              alpha = 0.25, color = NA, show.legend = FALSE,
              fill = "tan1") +
  geom_point(data = leye, aes(x = PercentAg, y = PC1), 
             color = "tan2") +
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none") + 
  geom_hline(size = 1, color = "red", linetype = "dashed", yintercept = 0)

# ---------------------------------------------------------------------------- #

# GAMs ####

library(splines)
library(mgcv)
library(MASS)

# class example
knots <- seq(2.4, 57.6, length = 27)
m <- gam(accel ~ s(times, bs = "bs", k = length(knots), sp = 0),
         data = mcycle, knots = list(x = knots))

## LEYE GAM (standardized) ####
knots <- seq(min(leye.cs$ts.sunrise), max(leye.cs$ts.sunrise), length = 6)

m <- gam(PC1 ~ s(ts.sunrise, bs = "bs", k = length(knots), sp = 0), 
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

ggplot(leye.cs,
       aes(x = ts.sunrise,
           y = PC1)) +
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
        legend.position = "none")
  
## LEYE GAM (unstandardized) ####
knots <- seq(min(leye$ts.sunrise), max(leye$ts.sunrise), length = 6)

m <- gam(PC1 ~ s(ts.sunrise, bs = "bs", k = length(knots), sp = 0), 
         data = leye, knots = list(x = knots))

# Plot model
d <- expand.grid(ts.sunrise = seq(min(leye$ts.sunrise),
                                  max(leye$ts.sunrise),
                                  length = 1000))

d$yhat <- predict(m, newdata = d)

d$se <- predict(m, newdata = d, se.fit = TRUE)$se.fit
d$lower <- d$yhat - 2*d$se
d$upper <- d$yhat + 2*d$se

ggplot(d,
       aes(x = ts.sunrise,
           y = PC1)) +
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

summary(m)

# AIC Comparison for GAM vs GLM ####

m1 <- gam(PC1 ~ s(ts.sunrise), data = leye.cs)

knots <- seq(min(leye.cs$ts.sunrise), max(leye.cs$ts.sunrise), length = 6)
m2 <- gam(PC1 ~ s(ts.sunrise, bs = "bs", k = length(knots), sp = 0), 
             data = leye.cs, knots = list(x = knots))

aic_m1 <- AIC(m1)
aic_m2 <- AIC(m2)

# Compare AIC values
aic_values <- c(m1 = aic_m1, m2 = aic_m2)

# GAM model slightly better fit than glm


# GAM for agricultural intensity ####
m <- gam(PC1 ~ s(PercentAg), data = leye)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg),
                                  max(leye$PercentAg),
                                  length = 1000))

d$yhat <- predict(m, newdata = d)

d$se <- predict(m, newdata = d, se.fit = TRUE)$se.fit
d$lower <- d$yhat - 2*d$se
d$upper <- d$yhat + 2*d$se

ggplot(d, aes(x = PercentAg, y = yhat)) +
  geom_line(color = "tan2") +  
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "tan1", alpha = 0.25) +  # Confidence intervals
  geom_point(data = leye, aes(x = PercentAg, y = PC1), color = "tan2") +  # Actual data points
  theme_light() +
  labs(x = "% Surrounding Agriculture", y = "Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 14, margin = margin(t = 13)),
        axis.title.y = element_text(size = 14, margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red",
             size = 1)

# ---------------------------------------------------------------------------- #

# Model Assumptions ####
m <- lm(PC1 ~ ts.sunrise, data = leye.cs)
summary(m)

library(splines)
m3 <- lm(PC1 ~ ns(ts.sunrise, df = 4), data = leye.cs)
summary(m3)

par(mfrow = c(2,2))
plot(m)

# Cyclical pattern in residuals...

# ASSESS WHETHER I SHOULD DO GAM OR ANOTHER NONLINEAR MODEL

m <- lm(PC1 ~ ts.sunrise + PercentAg, data = leye.cs)
summary(m)
par(mfrow = c(2,2))
plot(m)

# ---------------------------------------------------------------------------- #

# Plasma Metabolite Profiles at Different Capture Sites ####

# ---------------------------------------------------------------------------- #

# LAB MEETING DETECTION GRAPH ####

m <- lm(PC1 ~ Detection + ts.sunrise, data = leye)

d <- expand.grid(Detection = unique(leye$Detection),                    
                 ts.sunrise = mean(leye$ts.sunrise)) 

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
       y = "Predicted Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "red", size = 1)

summary(m)

# LAB MEETING AG INTENSITY GRAPH ####

m <- lm(PC1 ~ PercentAg + ts.sunrise, data = leye)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg),
                                 max(leye$PercentAg),
                                 length = 1000),                    
                 ts.sunrise = mean(leye$ts.sunrise)) 

predictions <- predict(m, newdata = d, se.fit = TRUE)

d$predicted_Mass <- predictions$fit

d$lwr <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upr <- d$predicted_Mass + 1.96 * predictions$se.fit

ggplot(d, aes(x = (PercentAg * 100), y = predicted_Mass)) +
  geom_line(size = 0.8, col = "black") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_light() +
  labs(x = "Surrounding Agricultural Intensity (%)", 
       y = "Predicted Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "red", size = 1)

summary(m)

# ---------------------------------------------------------------------------- #

# Age in different ag ####

leye$Age <- factor(leye$Age, levels = c("Juvenile", "Adult"))

ggplot(leye, aes(x = Age, y = PercentAg)) +
  geom_boxplot() +
  theme_light() +
  labs(y = "% Agriculture Surrounding Wetland", 
       x = "Age") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  theme(legend.position = "none")

# Refueling rates of juveniles versus adults ####
ggplot(leye, aes(x = Age, y = PC1)) +
  geom_boxplot() +
  theme_light() +
  labs(y = "Lesser Yellowlegs Fattening Index", 
       x = "Age") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "red", size = 1)

















