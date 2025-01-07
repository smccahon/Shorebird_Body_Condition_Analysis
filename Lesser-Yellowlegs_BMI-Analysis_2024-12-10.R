#----------------------------------#
#  Lesser Yellowlegs BMI Analysis  #
#          Created 12/10/2024      #          
#         Modified 01/06/2025      #
#----------------------------------#

# load packages
library(ggplot2)
library(AICcmodavg)
library(tidyverse)

# ---------------------------------------------------------------------------- #

# Data Processing & Manipulation ####

# Read data
setwd("processed_data")
birds <- read.csv("Shorebird_Data_Cleaned_2024-12-9.csv")

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

# Perform PCA: Females ####
# Peig and Green 2009 & Bajracharya 2022

# Subset the dataset to only include females
leye_female <- subset(leye.cs, Sex == "Female")

# Subset the dataset to only include tarsus, wing length, and bill length
leye_female <- leye_female[, c("Wing", "Culmen", "DiagTarsus")]

# Run PCA on the subsetted data
pca_result <- prcomp(leye_female, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# Extract residuals 
PC1_female <- pca_result$x[, 1]

# Regress PC1 against body mass
m_female <- lm(PC1_female ~ Mass, data = subset(leye.cs, Sex == "Female"))

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.f <- resid(m_female)

# Add size-corrected mass back into the full datasest
leye.cs$sc.mass <- NA  # Initialize with NA values
leye.cs[leye.cs$Sex == "Female", "sc.mass"] <- sc.mass.f


# Perform PCA: Males ####

# Subset the dataset to only include males
leye_male <- subset(leye.cs, Sex == "Male")

# Subset the dataset to only include tarsus, wing length, and bill length
leye_male <- leye_male[, c("Wing", "Culmen", "DiagTarsus")]

# Run PCA on the subsetted data
pca_result <- prcomp(leye_male, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# Extract residuals 
PC1_male <- pca_result$x[, 1]

# Regress PC1 against body mass
m_male <- lm(PC1_male ~ Mass, data = subset(leye.cs, Sex == "Male"))

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.m <- resid(m_male)

# Add size-corrected mass back into the full datasest
leye.cs[leye.cs$Sex == "Male", "sc.mass"] <- sc.mass.m

# Add size-corrected mass back into the full datasest
leye[leye.cs$Sex == "Male", "sc.mass"] <- sc.mass.m

# Add size-corrected mass back into the full datasest
leye[leye.cs$Sex == "Female", "sc.mass"] <- sc.mass.f

# ---------------------------------------------------------------------------- #

# Data Visualization ####

# Sex-specific body condition ~ Detection
ggplot(leye.cs, 
       aes(x = Detection, y = sc.mass, fill = Detection)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(color = Detection),
              position = position_jitterdodge(jitter.width = 0.2),  
              size = 2, 
              alpha = 0.8) + 
  theme_light() +
  labs(x = "Neonicotinoid Detection", 
       y = "Size-corrected Body Mass") +  
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

# Fattening Index ~ PercentAg ####
ggplot(leye, 
       aes(x = PercentAg, y = sc.mass)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed",
             size = 1) +
  geom_point(color = "black", size = 2) +  
  theme_light() +
  labs(x = "Surrounding Agricultural Intensity", 
       y = "Lesser Yellowlegs Body Condition") +  
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none") +
  xlim(0,1)

# for presentation
ggplot(leye.cs, 
       aes(x = AgCategory, y = sc.mass)) + 
  geom_boxplot(outlier.shape = NA, fill = "#F1CCB5") + 
  geom_point(size = 2, alpha = 0.4, color = "black") +
  labs(x = "Surrounding Agricultural Intensity", 
       y = "Lesser Yellowlegs Body Condition") +  
  theme_classic() +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.position = "none",
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed",
             size = 1) +
  scale_x_discrete(labels = c("Low <25%", "Moderate (25-50%)", "High (>75%)"))


# Fattening index ~ Percent Ag ####
ggplot(leye, 
       aes(x = PercentAg, y = sc.mass)) + 
  geom_point(size = 2) + 
  labs(x = "Surrounding Agricultural Intensity", 
       y = "Lesser Yellowlegs Body Condition") +  
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

m1 <- lm(sc.mass ~ DaysIntoSeason_S, data = leye.cs)
m2 <- lm(sc.mass ~ Julian, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Date into season performs better

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Agricultural Intensity ####
# Can't use ag category for LEYE because there are only three birds in the low category

m1 <- lm(sc.mass ~ PercentAg, data = leye.cs)
m2 <- lm(sc.mass ~ PercentAg + I(PercentAg^2), data = leye.cs)
m3 <- lm(sc.mass ~ LogAg, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Transformation does not significantly improve model fit. Use linear approximation of PercentAg

summary(m1)
confint(m1)
summary(m3)
confint(m3)

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Event and Capture Time ####

m1 <- lm(sc.mass ~ Event * ts.sunrise, data = leye.cs)
m2 <- lm(sc.mass ~ Event + I(ts.sunrise^2), data = leye.cs)
m3 <- lm(sc.mass ~ Event + ts.sunrise, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model without interaction performs the best.

# ---------------------------------------------------------------------------- #

# INTERACTIONS INCLUDED: NONE

# ---------------------------------------------------------------------------- #

## Single Models ####

m.global <- lm(sc.mass ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + Detection + PercentAg, 
               data = leye.cs)

m.null <- lm(sc.mass ~ 1, data = leye.cs)

m1 <- lm(sc.mass ~ Sex, data = leye.cs)
m2 <- lm(sc.mass ~ Event, data = leye.cs)
m3 <- lm(sc.mass ~ DaysIntoSeason_S, data = leye.cs)
m4 <- lm(sc.mass ~ ts.sunrise, data = leye.cs)
m5 <- lm(sc.mass ~ PercentAg, data = leye.cs)
m6 <- lm(sc.mass ~ Detection, data = leye.cs)

## Additive Models ####

### Two additive combinations ####
m7 <- lm(sc.mass ~ Sex + Event, data = leye.cs)
m8 <- lm(sc.mass ~ Sex + DaysIntoSeason_S, data = leye.cs)
m9 <- lm(sc.mass ~ Sex + ts.sunrise, data = leye.cs)
m10 <- lm(sc.mass ~ Sex + PercentAg, data = leye.cs)
m11 <- lm(sc.mass ~ Sex + Detection, data = leye.cs)
m12 <- lm(sc.mass ~ Event + DaysIntoSeason_S, data = leye.cs)
m13 <- lm(sc.mass ~ Event + ts.sunrise, data = leye.cs)
m14 <- lm(sc.mass ~ Event + PercentAg, data = leye.cs)
m15 <- lm(sc.mass ~ DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m16 <- lm(sc.mass ~ DaysIntoSeason_S + PercentAg, data = leye.cs)
m17 <- lm(sc.mass ~ DaysIntoSeason_S + Detection, data = leye.cs)
m18 <- lm(sc.mass ~ ts.sunrise + PercentAg, data = leye.cs)
m19 <- lm(sc.mass ~ ts.sunrise + Detection, data = leye.cs)
m20 <- lm(sc.mass ~ PercentAg + Detection, data = leye.cs)

### Three additive combinations ####
m21 <- lm(sc.mass ~ Sex + Event + DaysIntoSeason_S, data = leye.cs)
m22 <- lm(sc.mass ~ Sex + Event + ts.sunrise, data = leye.cs)
m23 <- lm(sc.mass ~ Sex + Event + PercentAg, data = leye.cs)
m24 <- lm(sc.mass ~ Sex + Event + Detection, data = leye.cs)
m25 <- lm(sc.mass ~ Sex + DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m26 <- lm(sc.mass ~ Sex + DaysIntoSeason_S + PercentAg, data = leye.cs)
m27 <- lm(sc.mass ~ Sex + DaysIntoSeason_S + Detection, data = leye.cs)
m28 <- lm(sc.mass ~ Sex + ts.sunrise + PercentAg, data = leye.cs)
m29 <- lm(sc.mass ~ Sex + ts.sunrise + Detection, data = leye.cs)
m30 <- lm(sc.mass ~ Sex + PercentAg + Detection, data = leye.cs)
m31 <- lm(sc.mass ~ Event + DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m32 <- lm(sc.mass ~ Event + DaysIntoSeason_S + PercentAg, data = leye.cs)
m33 <- lm(sc.mass ~ Event + DaysIntoSeason_S + Detection, data = leye.cs)
m34 <- lm(sc.mass ~ Event + ts.sunrise + PercentAg, data = leye.cs)
m35 <- lm(sc.mass ~ Event + ts.sunrise + Detection, data = leye.cs)
m36 <- lm(sc.mass ~ Event + PercentAg + Detection, data = leye.cs)
m37 <- lm(sc.mass ~ DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m38 <- lm(sc.mass ~ DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m39 <- lm(sc.mass ~ DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m40 <- lm(sc.mass ~ ts.sunrise + PercentAg + Detection, data = leye.cs)

### Four additive combinations ####
m41 <- lm(sc.mass ~ Sex + Event + DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m42 <- lm(sc.mass ~ Sex + Event + DaysIntoSeason_S + PercentAg, data = leye.cs)
m43 <- lm(sc.mass ~ Sex + Event + DaysIntoSeason_S + Detection, data = leye.cs)
m44 <- lm(sc.mass ~ Sex + Event + ts.sunrise + PercentAg, data = leye.cs)
m45 <- lm(sc.mass ~ Sex + Event + ts.sunrise + Detection, data = leye.cs)
m46 <- lm(sc.mass ~ Sex + Event + PercentAg + Detection, data = leye.cs)
m47 <- lm(sc.mass ~ Sex + DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m48 <- lm(sc.mass ~ Sex + DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m49 <- lm(sc.mass ~ Sex + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m50 <- lm(sc.mass ~ Sex + ts.sunrise + PercentAg + Detection, data = leye.cs)
m51 <- lm(sc.mass ~ Event + DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m52 <- lm(sc.mass ~ Event + DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m53 <- lm(sc.mass ~ Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m54 <- lm(sc.mass ~ Event + ts.sunrise + PercentAg + Detection, data = leye.cs)
m55 <- lm(sc.mass ~ DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)

### Five additive combinations ####

m56 <- lm(sc.mass ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m57 <- lm(sc.mass ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m58 <- lm(sc.mass ~ Sex + Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m59 <- lm(sc.mass ~ Sex + Event + ts.sunrise + PercentAg + Detection, data = leye.cs)
m60 <- lm(sc.mass ~ Sex + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)
m61 <- lm(sc.mass ~ Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)

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

confint(m4)
confint(m5)
confint(m3)

# No significant covariate effects on sex-specific, size-corrected body mass

# ---------------------------------------------------------------------------- #

# Top Model Figuers ####

## LAB MEETING DETECTION GRAPH ####
m <- lm(sc.mass ~ Detection, data = leye)

d <- expand.grid(Detection = unique(leye$Detection)) 

predictions <- predict(m, newdata = d, se.fit = TRUE)

d$predicted_Mass <- predictions$fit

d$lower_CI <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upper_CI <- d$predicted_Mass + 1.96 * predictions$se.fit

ggplot(d, aes(x = Detection, y = predicted_Mass)) +
  geom_point(size = 5, col = "black") + 
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.1,
                col = "black",
                size = 1) +
  theme_light() +
  labs(x = NULL, 
       y = "Predicted Body Condition Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1)







## LAB MEETING AG GRAPH ####
m <- lm(sc.mass ~ PercentAg, data = leye)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg), 
                                 max(leye$PercentAg), 
                                 length = 1000))

d <- cbind(d, predict(m, newdata = d, interval = "confidence"))

ggplot(d, aes(x = (PercentAg * 100), y = fit)) +
  geom_line(size = 0.8, col = "black") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_light() +
  labs(x = "Surrounding Agricultural Intensity (%)", 
       y = "Predicted Body Condition Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", color = "red")

summary(m)


