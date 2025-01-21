#---------------------------------#
# All Species Body Mass Analysis  #
#         Created 11/11/2024      #          
#         Modified 01/20/2025     #
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
library(nlme)
library(lattice)
library(MuMIn)

# ---------------------------------------------------------------------------- #

# Data Processing & Manipulation ####

# Read data
setwd("processed_data")
birds <- read.csv("Shorebird_Data_Cleaned_2024-12-9.csv")

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

# Logarithmic transformation of mass
birds <- birds %>% 
  mutate(LogMass = log10(Mass))

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


# Only include birds with mass 
birds <- birds %>%
  filter(!is.na(Mass))

# Only include species with at least three individuals
birds <- birds %>% 
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

birds$AgCategory <- factor(birds$AgCategory,
                           levels = c("Low", "Moderate", "High"))

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

# Standardize continuous variables and remove attributes
birds.cs <- birds %>%
  mutate(across(where(is.numeric), ~ scale(.) %>% `attributes<-`(NULL)))

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

# Mass ~ Agricultural Intensity
ggplot(data = birds,
       aes(x = AgCategory,
           y = Mass)) + 
  geom_boxplot() + theme_classic()

ggplot(data = birds,
       aes(x = PercentAg,
           y = Mass)) + 
  geom_point() + theme_classic()

# Mass ~ Migratory status
ggplot(data = birds,
       aes(x = MigStatus,
           y = Mass)) + 
  geom_boxplot() + theme_classic()


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

## Univariate Analysis: Agricultural Intensity ####

m1 <- lmer(Mass ~ AgCategory + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(Mass ~ PercentAg + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Agricultural intensity as a continuous variable performs significantly better.

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

# Model with interaction performs significantly better

# ---------------------------------------------------------------------------- #

# Address heteroscedasticity ####

m.global <- lmer(Mass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + 
                   MigStatus + PercentAg + Detection +
                   MigStatus * DaysIntoSeason_S +
                   Event * ts.sunrise + (1 | Species), REML = FALSE, data = birds.cs,
                   na.action = na.exclude)

fitted_vals <- fitted(m.global)  
residuals_vals <- resid(m.global)  
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# Plot residuals vs Event (CLEAR VIOLATION)
boxplot(residuals_vals ~ birds.cs$Event, xlab = "Event", ylab = "Residuals")

# Plot residuals vs MigStatus (CLEAR VIOLATION)
boxplot(residuals_vals ~ birds.cs$MigStatus, xlab = "MigStatus", ylab = "Residuals")

# Plot residuals vs Sex (slight violation)
boxplot(residuals_vals ~ birds.cs$Sex, xlab = "Sex", ylab = "Residuals")

# Plot residuals vs Detection (pretty good)
boxplot(residuals_vals ~ birds.cs$Detection, xlab = "Detection", ylab = "Residuals")

# Plot residuals vs Capture Time (pretty good)
plot(birds.cs$ts.sunrise, residuals_vals, xlab = "Capture Time", ylab = "Residuals")

# Plot residuals vs Date (good)
plot(birds.cs$DaysIntoSeason_S, residuals_vals, xlab = "Date", ylab = "Residuals")

# Plot residuals vs Agricultural Intensity (slight violation)
plot(birds.cs$PercentAg, residuals_vals, xlab = "Agricultural Intensity", ylab = "Residuals")

# Address heteroscedasticity: Log Transformation ####

m.global <- lmer(LogMass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + 
                   MigStatus + PercentAg + Detection +
                   MigStatus * DaysIntoSeason_S +
                   Event * ts.sunrise + (1 | Species), REML = FALSE, data = birds.cs,
                 na.action = na.exclude)

fitted_vals <- fitted(m.global)  
residuals_vals <- resid(m.global)  
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# Plot residuals vs Event (pretty good)
boxplot(residuals_vals ~ birds.cs$Event, xlab = "Event", ylab = "Residuals")

# Plot residuals vs MigStatus (pretty good)
boxplot(residuals_vals ~ birds.cs$MigStatus, xlab = "MigStatus", ylab = "Residuals")

# Plot residuals vs Sex (pretty good)
boxplot(residuals_vals ~ birds.cs$Sex, xlab = "Sex", ylab = "Residuals")

# Plot residuals vs Detection (pretty good)
boxplot(residuals_vals ~ birds.cs$Detection, xlab = "Detection", ylab = "Residuals")

# Plot residuals vs Capture Time (pretty good)
plot(birds.cs$ts.sunrise, residuals_vals, xlab = "Capture Time", ylab = "Residuals")

# Plot residuals vs Date (pretty good)
plot(birds.cs$DaysIntoSeason_S, residuals_vals, xlab = "Date", ylab = "Residuals")

# Plot residuals vs Agricultural Intensity (slight violation)
plot(birds.cs$PercentAg, residuals_vals, xlab = "Agricultural Intensity", ylab = "Residuals")

# Use a different variance structure and account for agricultural intensity issue
m.varFixed <- lme(LogMass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + 
                    MigStatus + PercentAg + Detection +
                    MigStatus * DaysIntoSeason_S +
                    Event * ts.sunrise,
                  random = ~ 1 | Species,
                  na.action = na.exclude,
                  method = "REML",
                  weights = varFixed(form = ~ PercentAg),
                  data = birds.cs)

fitted_vals <- fitted(m.varFixed)  
residuals_vals <- resid(m.varFixed)  
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# Plot residuals vs Agricultural Intensity (slight violation)
plot(birds.cs$PercentAg, residuals_vals, xlab = "Agricultural Intensity", ylab = "Residuals")

# Conclusion: none of these variance structures resolve issue with agricultural intensity
# It's only a slight violation so it's okay

#--------------#

## Identify optimal variance structure ####

m.lme <- lme(Mass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + 
                   MigStatus + PercentAg + Detection +
                   MigStatus * DaysIntoSeason_S +
                   Event * ts.sunrise,
             random = ~ 1 | Species, data = birds.cs,
             na.action = na.exclude)

fitted_vals <- fitted(m.lme)  
residuals_vals <- resid(m.lme)  
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

### Variance structures ####

# VarFixed (Fixed variance): Performs better than nothing
m.varFixed <- lme(Mass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + 
               MigStatus + PercentAg + Detection +
               MigStatus * DaysIntoSeason_S +
               Event * ts.sunrise,
               random = ~ 1 | Species,
               na.action = na.exclude,
               method = "REML",
               weights = varFixed(~ PercentAg),
               data = birds.cs)

anova(m.lme, m.varFixed)

aictab(list(m.lme, m.varFixed), modnames = c("m.lme", "m.varFixed"))

# Did it fix our issue? Nope
fitted_vals <- fitted(m.varFixed)  
residuals_vals <- resid(m.varFixed)  
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# VarIdent (Different variances per stratum)

m.varIdentEvent <- lme(Mass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + 
                    MigStatus + PercentAg + Detection +
                    MigStatus * DaysIntoSeason_S +
                    Event * ts.sunrise,
                  random = ~ 1 | Species,
                  na.action = na.exclude,
                  method = "REML",
                  weights = varIdent(form = ~ 1 | Event),
                  data = birds.cs)

m.varIdentMigStatus <- lme(Mass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + 
                    MigStatus + PercentAg + Detection +
                    MigStatus * DaysIntoSeason_S +
                    Event * ts.sunrise,
                  random = ~ 1 | Species,
                  na.action = na.exclude,
                  method = "REML",
                  weights = varIdent(form = ~ 1 | MigStatus),
                  data = birds.cs)

m.varIdentSpecies <- lme(Mass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + 
                             MigStatus + PercentAg + Detection +
                             MigStatus * DaysIntoSeason_S +
                             Event * ts.sunrise,
                           random = ~ 1 | Species,
                           na.action = na.exclude,
                           method = "REML",
                           weights = varIdent(form = ~ 1 | Species),
                           data = birds.cs)

# m.varIdentSpecies is best
anova(m.lme, m.varFixed, m.varIdentEvent, m.varIdentMigStatus, m.varIdentSpecies)

# Did it fix our issue? No
fitted_vals <- fitted(m.varIdentSpecies)  
residuals_vals <- resid(m.varIdentSpecies)  
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# VarPower (Power of the variance covariate)

# VarExp (Exponential of the variance covariate)

# VarConstPower (Constant plus power of the variance covariate)

# VarComb (A combination of variance functions)

# Log transformation of response variable
m.log <- lme(LogMass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + 
               MigStatus + PercentAg + Detection +
               MigStatus * DaysIntoSeason_S +
               Event * ts.sunrise,
             random = ~ 1 | Species, data = birds.cs,
             na.action = na.exclude)

# much better, solved all my issues
fitted_vals <- fitted(m.log)  
residuals_vals <- resid(m.log)  
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# Plot residuals vs Agricultural Intensity (slight violation but not bad)
plot(birds.cs$PercentAg, residuals_vals, xlab = "Agricultural Intensity", ylab = "Residuals")

# ---------------------------------------------------------------------------- #

m.global <- lmer(LogMass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + 
                     MigStatus + PercentAg + Detection +
                     MigStatus * DaysIntoSeason_S +
                     Event * ts.sunrise + (1 | Species), REML = FALSE, data = birds.cs)

m.null <- lmer(LogMass ~ 1 + (1|Species), data = birds.cs, REML = FALSE)

## Single Covariate Models ####

m1 <- lmer(LogMass ~ Sex + (1|Species), data = birds.cs, REML = FALSE)
m2 <- lmer(LogMass ~ Event + (1|Species), data = birds.cs, REML = FALSE)
m3 <- lmer(LogMass ~ ts.sunrise + (1|Species), data = birds.cs, REML = FALSE)
m4 <- lmer(LogMass ~ DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m5 <- lmer(LogMass ~ MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m6 <- lmer(LogMass ~ PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m7 <- lmer(LogMass ~ Detection + (1|Species), data = birds.cs, REML = FALSE)

##  Additive Models ####

### Two additive combinations ####

m8 <- lmer(LogMass ~ Sex + Event + (1|Species), data = birds.cs, REML = FALSE)
m9 <- lmer(LogMass ~ Sex + ts.sunrise + (1|Species), data = birds.cs, REML = FALSE)
m10 <- lmer(LogMass ~ Sex + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m11 <- lmer(LogMass ~ Sex + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m12 <- lmer(LogMass ~ Sex + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m13 <- lmer(LogMass ~ Sex + Detection + (1|Species), data = birds.cs, REML = FALSE)
m14 <- lmer(LogMass ~ Event + ts.sunrise + (1|Species) + Event*ts.sunrise, data = birds.cs, REML = FALSE)
m15 <- lmer(LogMass ~ Event + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m16 <- lmer(LogMass ~ Event + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m17 <- lmer(LogMass ~ Event + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m18 <- lmer(LogMass ~ Event + Detection + (1|Species), data = birds.cs, REML = FALSE)
m19 <- lmer(LogMass ~ ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m20 <- lmer(LogMass ~ ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m21 <- lmer(LogMass ~ ts.sunrise + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m22 <- lmer(LogMass ~ ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m23 <- lmer(LogMass ~ DaysIntoSeason_S + MigStatus + (1|Species) + DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m24 <- lmer(LogMass ~ DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m25 <- lmer(LogMass ~ DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m26 <- lmer(LogMass ~ MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m27 <- lmer(LogMass ~ MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m28 <- lmer(LogMass ~ PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)

### Three Additive Combinations ####

m29 <- lmer(LogMass ~ Sex + Event + ts.sunrise + (1|Species) + Event * ts.sunrise, data = birds.cs, REML = FALSE)
m30 <- lmer(LogMass ~ Sex + Event + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m31 <- lmer(LogMass ~ Sex + Event + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m32 <- lmer(LogMass ~ Sex + Event + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m33 <- lmer(LogMass ~ Sex + Event + Detection + (1|Species), data = birds.cs, REML = FALSE)
m34 <- lmer(LogMass ~ Sex + ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m35 <- lmer(LogMass ~ Sex + ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m36 <- lmer(LogMass ~ Sex + ts.sunrise + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m37 <- lmer(LogMass ~ Sex + ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m38 <- lmer(LogMass ~ Sex + DaysIntoSeason_S + MigStatus + (1|Species) + DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m39 <- lmer(LogMass ~ Sex + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m40 <- lmer(LogMass ~ Sex + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m41 <- lmer(LogMass ~ Sex + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m42 <- lmer(LogMass ~ Sex + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m43 <- lmer(LogMass ~ Sex + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m44 <- lmer(LogMass ~ Event * ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m45 <- lmer(LogMass ~ Event * ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m46 <- lmer(LogMass ~ Event * ts.sunrise + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m47 <- lmer(LogMass ~ Event * ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m48 <- lmer(LogMass ~ Event + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m49 <- lmer(LogMass ~ Event + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m50 <- lmer(LogMass ~ Event + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m51 <- lmer(LogMass ~ Event + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m52 <- lmer(LogMass ~ Event + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m53 <- lmer(LogMass ~ Event + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m54 <- lmer(LogMass ~ ts.sunrise + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m55 <- lmer(LogMass ~ ts.sunrise + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m56 <- lmer(LogMass ~ ts.sunrise + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m57 <- lmer(LogMass ~ ts.sunrise + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m58 <- lmer(LogMass ~ ts.sunrise + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m59 <- lmer(LogMass ~ ts.sunrise + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m60 <- lmer(LogMass ~ DaysIntoSeason_S * MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m61 <- lmer(LogMass ~ DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m62 <- lmer(LogMass ~ DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m63 <- lmer(LogMass ~ MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)

### Four Additive Combinations ####

m64 <- lmer(LogMass ~ Sex + Event * ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m65 <- lmer(LogMass ~ Sex + Event * ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m66 <- lmer(LogMass ~ Sex + Event * ts.sunrise + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m67 <- lmer(LogMass ~ Sex + Event * ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m68 <- lmer(LogMass ~ Sex + Event + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m69 <- lmer(LogMass ~ Sex + Event + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m70 <- lmer(LogMass ~ Sex + Event + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m71 <- lmer(LogMass ~ Sex + Event + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m72 <- lmer(LogMass ~ Sex + Event + MigStatus + Detection + (1|Species) + Event * MigStatus, data = birds.cs, REML = FALSE)
m73 <- lmer(LogMass ~ Sex + Event + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m74 <- lmer(LogMass ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m75 <- lmer(LogMass ~ Sex + ts.sunrise + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m76 <- lmer(LogMass ~ Sex + ts.sunrise + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m77 <- lmer(LogMass ~ Sex + ts.sunrise + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m78 <- lmer(LogMass ~ Sex + ts.sunrise + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m79 <- lmer(LogMass ~ Sex + ts.sunrise + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m80 <- lmer(LogMass ~ Sex + DaysIntoSeason_S * MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m81 <- lmer(LogMass ~ Sex + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m82 <- lmer(LogMass ~ Sex + DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m83 <- lmer(LogMass ~ Sex + MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m84 <- lmer(LogMass ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m85 <- lmer(LogMass ~ Event * ts.sunrise + DaysIntoSeason_S + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m86 <- lmer(LogMass ~ Event * ts.sunrise + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m87 <- lmer(LogMass ~ Event * ts.sunrise + MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m88 <- lmer(LogMass ~ Event * ts.sunrise + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m89 <- lmer(LogMass ~ Event * ts.sunrise + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m90 <- lmer(LogMass ~ Event + DaysIntoSeason_S * MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m91 <- lmer(LogMass ~ Event + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m92 <- lmer(LogMass ~ Event + DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m93 <- lmer(LogMass ~ Event + MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m94 <- lmer(LogMass ~ ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m95 <- lmer(LogMass ~ ts.sunrise + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m96 <- lmer(LogMass ~ ts.sunrise + DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m97 <- lmer(LogMass ~ ts.sunrise + MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m98 <- lmer(LogMass ~ DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)

### Five Additive Combinations ####

m99 <- lmer(LogMass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + MigStatus + (1|Species) +
               DaysIntoSeason_S * MigStatus + Event * ts.sunrise, data = birds.cs, REML = FALSE)
m100 <- lmer(LogMass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + PercentAg + (1|Species) +
               Event * ts.sunrise, data = birds.cs, REML = FALSE)
m101 <- lmer(LogMass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + Detection + (1|Species) +
               Event * ts.sunrise, data = birds.cs, REML = FALSE)
m102 <- lmer(LogMass ~ Sex + Event + ts.sunrise + MigStatus + PercentAg + (1|Species) +
               Event * ts.sunrise, data = birds.cs, REML = FALSE)
m103 <- lmer(LogMass ~ Sex + Event + ts.sunrise + MigStatus + Detection + (1|Species) +
               Event * ts.sunrise, data = birds.cs, REML = FALSE)
m104 <- lmer(LogMass ~ Sex + Event + ts.sunrise + PercentAg + Detection + (1|Species) +
               Event * ts.sunrise, data = birds.cs, REML = FALSE)
m105 <- lmer(LogMass ~ Sex + Event + DaysIntoSeason_S + MigStatus + PercentAg + (1|Species) +
               DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m106 <- lmer(LogMass ~ Sex + Event + DaysIntoSeason_S + MigStatus + Detection + (1|Species) +
                + DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m107 <- lmer(LogMass ~ Sex + Event + DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m108 <- lmer(LogMass ~ Sex + Event + MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m109 <- lmer(LogMass ~ Sex + ts.sunrise + DaysIntoSeason_S + MigStatus + PercentAg + (1|Species) +
               DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m110 <- lmer(LogMass ~ Sex + ts.sunrise + DaysIntoSeason_S + MigStatus + Detection + (1|Species) +
               DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m111 <- lmer(LogMass ~ Sex + ts.sunrise + DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m112 <- lmer(LogMass ~ Sex + ts.sunrise + MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m113 <- lmer(LogMass ~ Sex + DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m114 <- lmer(LogMass ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m115 <- lmer(LogMass ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m116 <- lmer(LogMass ~ Event * ts.sunrise + DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m117 <- lmer(LogMass ~ Event * ts.sunrise + MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m118 <- lmer(LogMass ~ Event + DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m119 <- lmer(LogMass ~ ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)

### Six Additive Combinations ####

m120 <- lmer(LogMass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + MigStatus + PercentAg + (1|Species) +
               Event * ts.sunrise + DaysIntoSeason_S * MigStatus, data = birds.cs, REML = FALSE)
m121 <- lmer(LogMass ~ Sex + Event * ts.sunrise + DaysIntoSeason_S * MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m122 <- lmer(LogMass ~ Sex + Event * ts.sunrise + DaysIntoSeason_S + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m123 <- lmer(LogMass ~ Sex + Event * ts.sunrise + MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m124 <- lmer(LogMass ~ Sex + Event + DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m125 <- lmer(LogMass ~ Sex + ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)
m126 <- lmer(LogMass ~ Event * ts.sunrise + DaysIntoSeason_S * MigStatus + PercentAg + Detection + (1|Species), data = birds.cs, REML = FALSE)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:126)

models <- mget(model_names)

models$m.null <- m.null
models$m.global <- m.global

model_names <- c(model_names, "m.null", "m.global")

aictab(models, modnames = model_names)

# NEW TOP RESULTS 2025-01-20 ####
summary(m84)
confint(m84) # event * time significant

summary(m114)
confint(m114) # percentag in model but not significant; event * time significant

summary(m44)
confint(m44) # date and event * time significant

summary(m45)
confint(m45) # migratory status and event * time significant


# ---------------------------------------------------------------------------- #

# Top Model Summaries ####

summary(m45)
confint(m45)
r.squaredGLMM(m45)
# M: 0.318, C: 0.964

summary(m87)
confint(m87)
r.squaredGLMM(m87)
# M: 0.317, C: 0.964

summary(m84)
confint(m84)
r.squaredGLMM(m84)
# M: 0.315, C: 0.966

summary(m14)
confint(m14)
r.squaredGLMM(m14)
# M: 0.00899, C: 0.966

summary(m114)
confint(m114)
r.squaredGLMM(m114)
# M: 0.319, C: 0.966

summary(m44)
confint(m44)
r.squaredGLMM(m44)
# M: 0.00944, C: 0.967

summary(m45)
confint(m45)
r.squaredGLMM(m45)
# M: 0.318, C: 0.964

# ---------------------------------------------------------------------------- #

# Graph top models ####

m103 <- lmer(Mass ~ Sex + Event + ts.sunrise + MigStatus + Detection + (1|Species) +
               Event * ts.sunrise + Event * Detection, data = birds, REML = FALSE)

# keeping everything constant but event, time, and species
d <- expand.grid(Event = c("Fall 2021", "Spring 2022", "Spring 2023", "Fall 2023"),
                 ts.sunrise = seq(min(birds$ts.sunrise), max(birds$ts.sunrise), 
                                  length = 1000),
                 MigStatus = c("Migratory"),
                 Detection = c("Detection"),
                 Sex = c("Male"),
                 Species = unique(birds$Species))

d$Mass <- predict(m103, newdata = d, re.form = NULL)

# points not being added correctly...likely because my expand.grid subsets data
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
           label = "Sunrise", angle = 90, vjust = -0.5, hjust = 1)


# keeping everything constant but migrant status
d <- expand.grid(Event = c("Fall 2023"),
                 ts.sunrise = mean(birds$ts.sunrise),
                 MigStatus = c("Migratory", "Resident"),
                 Detection = c("Detection"),
                 Sex = c("Male"),
                 Species = unique(birds$Species))

d$Mass <- predict(m103, newdata = d, re.form = NULL)


ggplot(d, aes(x = MigStatus, y = Mass, color = MigStatus)) +
  geom_boxplot() +
  theme_light() +
  labs(x = "Migrant Status", 
       y = "Body Mass (g)", 
       color = "Migratory Status") +
  theme(legend.position = "right") +
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10)) +
  scale_color_viridis_d(begin = 0.2, 
                        end = 0.8, 
                        option = "inferno")

## Other Graphing Code ####

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

## Mass ~ Detection * Event
# keep everything constant but detection and event
d <- expand.grid(Event = c("Fall 2021","Spring 2022","Spring 2023", "Fall 2023"),
                 ts.sunrise = mean(birds$ts.sunrise),
                 MigStatus = c("Migratory"),
                 Detection = c("Detection", "Non-detection"),
                 Sex = c("Male"),
                 Species = unique(birds$Species))

d$Mass <- predict(m103, newdata = d, re.form = NULL)


ggplot(d, aes(x = Detection, y = Mass, fill = Event)) +
  geom_boxplot() +
  theme_light() +
  labs(x = "Neonicotinoid Detection", 
       y = "Body Mass (g)", 
       fill = "Sampling Occasion") +
  theme(legend.position = "right") +
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10)) +
  scale_color_viridis_d(begin = 0.2, 
                        end = 0.8, 
                        option = "inferno")

## LAB MEETING DETECTION GRAPH SPRING 2025 ####
m <- lmer(LogMass ~ Detection + Event * ts.sunrise + MigStatus + Sex + (1 | Species), 
          data = birds, REML = FALSE)

d <- expand.grid(ts.sunrise = mean(birds$ts.sunrise),   
                 Event = c("Fall 2023"),                    
                 Detection = c("Detection", "Non-detection"),
                 Sex = c("Female"),
                 MigStatus = c("Resident"),
                 Species = unique(birds$Species)) 

d$Detection <- factor(d$Detection, levels = levels(birds$Detection))
d$Species <- factor(d$Species, levels = levels(birds$Species))

predictions <- predict(m, newdata = d, se.fit = TRUE, re.form = NULL)

d$predicted_Mass <- predictions$fit

d$lower_CI <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upper_CI <- d$predicted_Mass + 1.96 * predictions$se.fit

d$Species <- factor(d$Species, levels = c("Lesser Yellowlegs", 
                                          "Pectoral Sandpiper", 
                                          "Killdeer", 
                                          "American Avocet", 
                                          "Long-billed Dowitcher", 
                                          "Willet", 
                                          "Wilson's Phalarope", 
                                          "Semipalmated Sandpiper",
                                          "Least Sandpiper"))


ggplot(d, aes(x = Detection, y = predicted_Mass)) +
  geom_point(size = 3, col = "black") +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.1,
                col = "black",
                size = 1) +
  theme_light() +
  facet_wrap(~ Species) +
  labs(x = NULL, 
       y = "Log(Predicted Body Mass (g))") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none",
        strip.text = element_text(size = 18)) +
  theme(legend.position = "none")

library(lmerTest)
summary(m)
r.squaredGLMM(m)


## LAB MEETING AG GRAPH SPRING 2025 ####
m <- lmer(LogMass ~ PercentAg + Event * ts.sunrise + MigStatus + (1 | Species), 
        data = birds, REML = FALSE)

d <- expand.grid(PercentAg = seq(min(birds$PercentAg), 
                                 max(birds$PercentAg), 
                                 length = 1000),   
                 Event = c("Fall 2023"),                    
                 ts.sunrise = mean(birds$ts.sunrise),
                 MigStatus = c("Resident"),
                 Species = unique(birds$Species)) 

# scale continuous variables
# d$PercentAg <- scale(d$PercentAg, center = TRUE, scale = TRUE)
# d$ts.sunrise <- scale(d$ts.sunrise, center = TRUE, scale = TRUE)

predictions <- predict(m, newdata = d, se.fit = TRUE, re.form = NULL)

d$fit <- predictions$fit

d$lwr <- d$fit - 1.96 * predictions$se.fit
d$upr <- d$fit + 1.96 * predictions$se.fit

# SIZE
d$Species_group <- factor(
  ifelse(d$Species %in% c("Wilson's Phalarope", "Semipalmated Sandpiper", "Least Sandpiper"), "Small", 
         ifelse(d$Species %in% c("Lesser Yellowlegs", "Pectoral Sandpiper", "Killdeer"), "Medium", 
                ifelse(d$Species %in% c("American Avocet", "Long-billed Dowitcher", "Willet"), "Large", 
                       "Other Groups"))))

d$Species_group <- factor(d$Species_group, levels = c("Small", "Medium", "Large"))

custom_colors <- c(
  "Wilson's Phalarope" = "coral2", 
  "Semipalmated Sandpiper" = "skyblue2", 
  "Least Sandpiper" = "olivedrab3",
  "Lesser Yellowlegs" = "coral2",
  "Pectoral Sandpiper" = "skyblue2",
  "Killdeer" = "olivedrab3",
  "American Avocet" = "coral2",
  "Long-billed Dowitcher" = "skyblue2",
  "Willet" = "olivedrab3"  # 
)

# 3 GROUPS
ggplot(d, aes(x = (PercentAg * 100), y = fit, color = Species)) +
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Species), 
              alpha = 0.4, color = NA, show.legend = FALSE) +
  theme_light() +
  facet_wrap(~ Species_group, ncol = 3) + 
  labs(x = "Surrounding Agricultural Intensity (%)", 
       y = "Log(Predicted Body Mass (g))") +  
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18)) +
        # legend.position = "none") +
  # theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  scale_color_manual(values = custom_colors) + 
  scale_fill_manual(values = custom_colors)

# ALL: THIS ONE FOR LAB MEETING ####
d$Species <- factor(d$Species, levels = c("Lesser Yellowlegs", 
                                          "Pectoral Sandpiper", 
                                          "Killdeer", 
                                          "American Avocet", 
                                          "Long-billed Dowitcher", 
                                          "Willet", 
                                          "Wilson's Phalarope", 
                                          "Semipalmated Sandpiper",
                                          "Least Sandpiper"))

ggplot(d, aes(x = (PercentAg * 100), y = fit)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.2, show.legend = FALSE) +
  theme_light() +
  facet_wrap(~ Species) +
  labs(x = "Surrounding Agricultural Intensity (%)", 
       y = "Log(Predicted Body Mass (g))") + 
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none",
        strip.text = element_text(size = 18)) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 100, by = 20))

summary(m)

# derive p value
library(lmerTest)
summary(m)

# goodness of fit
r.squaredGLMM(m)




# ---------------------------------------------------------------------------- #

# Model Assumptions ####

# CLEAR VIOLATION OF HOMOSCEDASTICITY
M1 <- lmer(LogMass ~ Event * ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
M1 <- lmer(Mass ~ 1 + (1|Species), data= birds.cs, REML = FALSE)

# Set up the 2x2 plot layout
op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))

# Plot residuals vs fitted values (CLEAR VIOLATION)
fitted_vals <- fitted(M1)  # Get fitted values from the model
residuals_vals <- resid(M1)  # Get residuals from the model
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")  # Add a horizontal line at zero for reference

# Plot histogram of residuals (satisfied)
hist(residuals_vals, xlab = "Residuals", main = "")

# Plot residuals vs MigStatus (CLEAR VIOLATION)
boxplot(residuals_vals ~ birds.cs$Event, xlab = "Event", ylab = "Residuals")

# Plot residuals vs Capture Time
plot(birds.cs$ts.sunrise, residuals_vals, xlab = "Capture Time", ylab = "Residuals")

# Restore original plot settings
par(op)

# Plot residuals vs MigStatus (CLEAR VIOLATION)
boxplot(residuals_vals ~ birds.cs$Event, xlab = "Event", ylab = "Residuals")
