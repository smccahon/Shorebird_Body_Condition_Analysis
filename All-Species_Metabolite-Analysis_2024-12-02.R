#-----------------------------------#
#  All Species Metabolite Analysis  #
#          Created 12/02/2024       #          
#         Modified 01/07/2025       #
#-----------------------------------#

# load packages
library(dplyr)
library(ggplot2)
library(AICcmodavg)
library(tidyverse)
library(lme4)
library(car)

# ---------------------------------------------------------------------------- #

# Read data
setwd("processed_data")
birds <- read.csv("Shorebird_Data_Cleaned_2024-12-9.csv")

# Make neonicotinoid detection column (Detection/Non-detection)
birds$Detection <- ifelse(birds$OverallNeonic > 0, "Detection", "Non-detection")

# Reorder factor variables
birds$Sex <- factor(birds$Sex,
                    levels = c("M", "F"),
                    labels = c("Male", "Female"))

birds$Detection <- as.factor(birds$Detection)

birds$AgCategory <- factor(birds$AgCategory,
                           levels = c("Low", "Moderate", "High"))

birds$Species <- factor(birds$Species, 
                        levels = c("LesserYellowlegs", "Killdeer", "Willet", 
                                   "PectoralSandpiper", "WIPH", "LeastSandpiper", 
                                   "LongbilledDowitcher", "AmericanAvocet",
                                   "SemipalmatedSandpiper", "MarbledGodwit",
                                   "ShortbilledDowitcher", "GreaterYellowlegs"),
                        labels = c("Lesser Yellowlegs", "Killdeer", "Willet", 
                                   "Pectoral Sandpiper", "Wilson's Phalarope", 
                                   "Least Sandpiper","Long-billed Dowitcher", "American Avocet",
                                   "Semipalmated Sandpiper", "Marbled Godwit",
                                   "Short-billed Dowitcher", "Greater Yellowlegs"))

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

# Only include birds with metabolite data
birds <- birds %>% 
  filter(!is.na(Tri) & !is.na(Beta))

# Only include species with at least three individuals
birds <- birds %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

table(birds$Species)

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

# Standardize continuous variables
birds.cs <- birds %>%
  mutate(across(where(is.numeric), scale))

# ---------------------------------------------------------------------------- #

# Perform PCA ####

# Subset the dataset to only include 'Tri' and 'Beta'
birds_subset <- birds.cs[, c("Tri", "Beta")]

# Remove rows with NAs for PCA
birds_subset_clean <- birds_subset[complete.cases(birds_subset), ]

# Run PCA on the cleaned data
pca_result <- prcomp(birds_subset_clean, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# View the PCA scores (principal components)
pca_scores <- pca_result$x

# Merge PCA scores back to the original dataset
# Create a data frame to store the PCA scores for the rows with no missing data
birds.cs$PC1 <- NA  # Initialize with NA values
birds.cs$PC2 <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
birds.cs[complete.cases(birds_subset), "PC1"] <- pca_scores[, 1]
birds.cs[complete.cases(birds_subset), "PC2"] <- pca_scores[, 2]

# View the updated dataset with PCA scores
print(birds.cs)

# Merge PCA scores back to the original birds unstandardized dataset
# Create a data frame to store the PCA scores for the rows with no missing data
birds$PC1 <- NA  # Initialize with NA values
birds$PC2 <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
birds[complete.cases(birds_subset), "PC1"] <- pca_scores[, 1]
birds[complete.cases(birds_subset), "PC2"] <- pca_scores[, 2]

# View the updated dataset with PCA scores
print(birds)

# View the first few rows of the updated data frame
head(birds)

# ---------------------------------------------------------------------------- #

# Data Visualization ####

## Fattening Index ~ Detection ####
# No clear relationship
ggplot(birds, 
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
                        alpha = 0.8) +
  geom_hline(color = "red", linetype = "dashed", yintercept = 0, size = 1)

## Fattening Index ~ AgCategory ####
# low seems to be higher than high
ggplot(birds, 
       aes(x = AgCategory, y = PC1)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed",
             size = 1) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(color = ts.sunrise),
              position = position_jitterdodge(jitter.width = 0.2),  
              size = 2, 
              alpha = 0.8) + 
  theme_light() +
  labs(x = "Agricultural Intensity", 
       y = "All Species Fattening Index") +  
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none") +
  scale_color_viridis_c(option = "inferno", 
                        alpha = 0.8)

# significant different in agricultural intensities
kruskal.test(PC1 ~ AgCategory, data = birds)

# no significant difference in low and high agricultural intensities (only low and medium)
AgCategory.lh <- subset(birds, AgCategory %in% c("Low", "High"))
kruskal.test(PC1 ~ AgCategory, data = AgCategory.lh)

# check if there's a different result using a t.test
t.test(PC1 ~ AgCategory, data = AgCategory.lh) # nope

# ASG Presentation
# low seems to be higher than high
ggplot(birds, 
       aes(x = AgCategory, y = PC1)) + 
  geom_boxplot(outlier.shape = NA, fill = "#F1CCB5") + 
  geom_point(size = 2, alpha = 0.4, color = "black") +
  labs(x = "Surrounding Agricultural Intensity", 
       y = "All Species Fattening Index") +  
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

ggplot(birds, 
       aes(x = AgCategory, y = PC1)) + 
  geom_boxplot(outlier.shape = NA, fill = "#F1CCB5") +
  labs(x = "Surrounding Agricultural Intensity", 
       y = "All Species Fattening Index") +  
  theme_light() +
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


ggplot(birds, 
       aes(x = PercentAg, y = PC1)) + 
  geom_point(size = 2) + 
  labs(x = "Surrounding Agricultural Intensity", 
       y = "All Species Fattening Index") +  
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


# Fattening index ~ Percent Ag ####
ggplot(birds, 
       aes(x = PercentAg, y = PC1)) + 
  geom_point(size = 2) + 
  labs(x = "Surrounding Agricultural Intensity", 
       y = "All Species Fattening Index") +  
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

## Univariate Analysis: Date ####

m1 <- lmer(PC1 ~ DaysIntoSeason_S + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(PC1 ~ Julian + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Both perform similarily; use date into season for consistency

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Neonicotinoid ####

m1 <- lmer(PC1 ~ OverallNeonic + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(PC1 ~ LogNeonic + (1 | Species), data = birds.cs, REML = FALSE)
m3 <- lmer(PC1 ~ Detection + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Detection and LogNeonic perform similarily
# Log transformation is necessary to address issues of homoscedasticity and normality
# Detection makes sense biologically
# Just use detection for consistency with other models

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Event and Capture Time ####

m1 <- lmer(PC1 ~ Event * ts.sunrise + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(PC1 ~ Event + I(ts.sunrise^2) + (1 | Species), data = birds.cs, REML = FALSE)
m3 <- lmer(PC1 ~ Event + ts.sunrise + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Simple model without interaction by far performs the best.

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Date * Migratory Status ####

m1 <- lmer(PC1 ~ DaysIntoSeason_S + MigStatus + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(PC1 ~ DaysIntoSeason_S * MigStatus + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with interaction performs slightly better (wt = 52%)

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Agricultural Intensity ####

m1 <- lmer(PC1 ~ PercentAg + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(PC1 ~ PercentAg + I(PercentAg^2) + (1 | Species), data = birds.cs, REML = FALSE)
m3 <- lmer(PC1 ~ AgCategory + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# AgCategory as a category performs best

# ---------------------------------------------------------------------------- #

# INTERACTIONS INCLUDED: NONE; AgCategory (categorical)

# ---------------------------------------------------------------------------- #

m.global <- lmer(PC1 ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + 
                   MigStatus + AgCategory + Detection + (1 | Species), 
                   REML = FALSE, data = birds.cs)

m.null <- lmer(PC1 ~ 1 + (1|Species), data = birds.cs, REML = FALSE)

## Single Covariate Models ####

m1 <- lmer(PC1 ~ Sex + (1|Species), data = birds.cs, REML = FALSE)
m2 <- lmer(PC1 ~ Event + (1|Species), data = birds.cs, REML = FALSE)
m3 <- lmer(PC1 ~ ts.sunrise + (1|Species), data = birds.cs, REML = FALSE)
m4 <- lmer(PC1 ~ DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m5 <- lmer(PC1 ~ MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m6 <- lmer(PC1 ~ AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m7 <- lmer(PC1 ~ Detection + (1|Species), data = birds.cs, REML = FALSE)

##  Additive Models ####

### Two additive combinations ####

m8 <- lmer(PC1 ~ Sex + Event + (1|Species), data = birds.cs, REML = FALSE)
m9 <- lmer(PC1 ~ Sex + ts.sunrise + (1|Species), data = birds.cs, REML = FALSE)
m10 <- lmer(PC1 ~ Sex + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m11 <- lmer(PC1 ~ Sex + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m12 <- lmer(PC1 ~ Sex + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m13 <- lmer(PC1 ~ Sex + Detection + (1|Species), data = birds.cs, REML = FALSE)
m14 <- lmer(PC1 ~ Event + ts.sunrise + (1|Species), data = birds.cs, REML = FALSE)
m15 <- lmer(PC1 ~ Event + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m16 <- lmer(PC1 ~ Event + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m17 <- lmer(PC1 ~ Event + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m18 <- lmer(PC1 ~ Event + Detection + (1|Species), data = birds.cs, REML = FALSE)
m19 <- lmer(PC1 ~ ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m20 <- lmer(PC1 ~ ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m21 <- lmer(PC1 ~ ts.sunrise + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m22 <- lmer(PC1 ~ ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m23 <- lmer(PC1 ~ DaysIntoSeason_S + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m24 <- lmer(PC1 ~ DaysIntoSeason_S + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m25 <- lmer(PC1 ~ DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m26 <- lmer(PC1 ~ MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m27 <- lmer(PC1 ~ MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m28 <- lmer(PC1 ~ AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)

### Three Additive Combinations ####

m29 <- lmer(PC1 ~ Sex + Event + ts.sunrise + (1|Species), data = birds.cs, REML = FALSE)
m30 <- lmer(PC1 ~ Sex + Event + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m31 <- lmer(PC1 ~ Sex + Event + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m32 <- lmer(PC1 ~ Sex + Event + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m33 <- lmer(PC1 ~ Sex + Event + Detection + (1|Species), data = birds.cs, REML = FALSE)
m34 <- lmer(PC1 ~ Sex + ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m35 <- lmer(PC1 ~ Sex + ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m36 <- lmer(PC1 ~ Sex + ts.sunrise + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m37 <- lmer(PC1 ~ Sex + ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m38 <- lmer(PC1 ~ Sex + DaysIntoSeason_S + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m39 <- lmer(PC1 ~ Sex + DaysIntoSeason_S + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m40 <- lmer(PC1 ~ Sex + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m41 <- lmer(PC1 ~ Sex + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m42 <- lmer(PC1 ~ Sex + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m43 <- lmer(PC1 ~ Sex + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m44 <- lmer(PC1 ~ Event + ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m45 <- lmer(PC1 ~ Event + ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m46 <- lmer(PC1 ~ Event + ts.sunrise + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m47 <- lmer(PC1 ~ Event + ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m48 <- lmer(PC1 ~ Event + DaysIntoSeason_S + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m49 <- lmer(PC1 ~ Event + DaysIntoSeason_S + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m50 <- lmer(PC1 ~ Event + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m51 <- lmer(PC1 ~ Event + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m52 <- lmer(PC1 ~ Event + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m53 <- lmer(PC1 ~ Event + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m54 <- lmer(PC1 ~ ts.sunrise + DaysIntoSeason_S + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m55 <- lmer(PC1 ~ ts.sunrise + DaysIntoSeason_S + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m56 <- lmer(PC1 ~ ts.sunrise + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m57 <- lmer(PC1 ~ ts.sunrise + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m58 <- lmer(PC1 ~ ts.sunrise + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m59 <- lmer(PC1 ~ ts.sunrise + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m60 <- lmer(PC1 ~ DaysIntoSeason_S + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m61 <- lmer(PC1 ~ DaysIntoSeason_S + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m62 <- lmer(PC1 ~ DaysIntoSeason_S + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m63 <- lmer(PC1 ~ MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)

### Four Additive Combinations ####

m64 <- lmer(PC1 ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + (1|Species), data = birds.cs, REML = FALSE)
m65 <- lmer(PC1 ~ Sex + Event + ts.sunrise + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m66 <- lmer(PC1 ~ Sex + Event + ts.sunrise + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m67 <- lmer(PC1 ~ Sex + Event + ts.sunrise + Detection + (1|Species), data = birds.cs, REML = FALSE)
m68 <- lmer(PC1 ~ Sex + Event + DaysIntoSeason_S + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m69 <- lmer(PC1 ~ Sex + Event + DaysIntoSeason_S + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m70 <- lmer(PC1 ~ Sex + Event + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m71 <- lmer(PC1 ~ Sex + Event + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m72 <- lmer(PC1 ~ Sex + Event + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m73 <- lmer(PC1 ~ Sex + Event + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m74 <- lmer(PC1 ~ Sex + ts.sunrise + DaysIntoSeason_S + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m75 <- lmer(PC1 ~ Sex + ts.sunrise + DaysIntoSeason_S + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m76 <- lmer(PC1 ~ Sex + ts.sunrise + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m77 <- lmer(PC1 ~ Sex + ts.sunrise + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m78 <- lmer(PC1 ~ Sex + ts.sunrise + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m79 <- lmer(PC1 ~ Sex + ts.sunrise + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m80 <- lmer(PC1 ~ Sex + DaysIntoSeason_S + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m81 <- lmer(PC1 ~ Sex + DaysIntoSeason_S + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m82 <- lmer(PC1 ~ Sex + DaysIntoSeason_S + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m83 <- lmer(PC1 ~ Sex + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m84 <- lmer(PC1 ~ Event + ts.sunrise + DaysIntoSeason_S + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m85 <- lmer(PC1 ~ Event + ts.sunrise + DaysIntoSeason_S + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m86 <- lmer(PC1 ~ Event + ts.sunrise + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m87 <- lmer(PC1 ~ Event + ts.sunrise + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m88 <- lmer(PC1 ~ Event + ts.sunrise + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m89 <- lmer(PC1 ~ Event + ts.sunrise + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m90 <- lmer(PC1 ~ Event + DaysIntoSeason_S + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m91 <- lmer(PC1 ~ Event + DaysIntoSeason_S + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m92 <- lmer(PC1 ~ Event + DaysIntoSeason_S + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m93 <- lmer(PC1 ~ Event + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m94 <- lmer(PC1 ~ ts.sunrise + DaysIntoSeason_S + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m95 <- lmer(PC1 ~ ts.sunrise + DaysIntoSeason_S + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m96 <- lmer(PC1 ~ ts.sunrise + DaysIntoSeason_S + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m97 <- lmer(PC1 ~ ts.sunrise + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m98 <- lmer(PC1 ~ DaysIntoSeason_S + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)

### Five Additive Combinations ####

m99 <- lmer(PC1 ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m100 <- lmer(PC1 ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m101 <- lmer(PC1 ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + Detection + (1|Species), data = birds.cs, REML = FALSE)
m102 <- lmer(PC1 ~ Sex + Event + ts.sunrise + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m103 <- lmer(PC1 ~ Sex + Event + ts.sunrise + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m104 <- lmer(PC1 ~ Sex + Event + ts.sunrise + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m105 <- lmer(PC1 ~ Sex + Event + DaysIntoSeason_S + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m106 <- lmer(PC1 ~ Sex + Event + DaysIntoSeason_S + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m107 <- lmer(PC1 ~ Sex + Event + DaysIntoSeason_S + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m108 <- lmer(PC1 ~ Sex + Event + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m109 <- lmer(PC1 ~ Sex + ts.sunrise + DaysIntoSeason_S + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m110 <- lmer(PC1 ~ Sex + ts.sunrise + DaysIntoSeason_S + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m111 <- lmer(PC1 ~ Sex + ts.sunrise + DaysIntoSeason_S + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m112 <- lmer(PC1 ~ Sex + ts.sunrise + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m113 <- lmer(PC1 ~ Sex + DaysIntoSeason_S + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m114 <- lmer(PC1 ~ Event + ts.sunrise + DaysIntoSeason_S + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m115 <- lmer(PC1 ~ Event + ts.sunrise + DaysIntoSeason_S + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m116 <- lmer(PC1 ~ Event + ts.sunrise + DaysIntoSeason_S + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m117 <- lmer(PC1 ~ Event + ts.sunrise + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m118 <- lmer(PC1 ~ Event + DaysIntoSeason_S + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m119 <- lmer(PC1 ~ ts.sunrise + DaysIntoSeason_S + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)

### Six Additive Combinations ####

m120 <- lmer(PC1 ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + MigStatus + AgCategory + (1|Species), data = birds.cs, REML = FALSE)
m121 <- lmer(PC1 ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + MigStatus + Detection + (1|Species), data = birds.cs, REML = FALSE)
m122 <- lmer(PC1 ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m123 <- lmer(PC1 ~ Sex + Event + ts.sunrise + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m124 <- lmer(PC1 ~ Sex + Event + DaysIntoSeason_S + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m125 <- lmer(PC1 ~ Sex + ts.sunrise + DaysIntoSeason_S + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)
m126 <- lmer(PC1 ~ Event + ts.sunrise + DaysIntoSeason_S + MigStatus + AgCategory + Detection + (1|Species), data = birds.cs, REML = FALSE)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:126)

models <- mget(model_names)

models$m.null <- m.null
models$m.global <- m.global

model_names <- c(model_names, "m.null", "m.global")

aictab(models, modnames = model_names)

# ---------------------------------------------------------------------------- #

# Top Model Summaries ####

# Significant difference between low and moderate ag intensity.
summary(m6)
confint(m6)
r.squaredGLMM(m6)
# M: 0.0898, C: 0.264

m6 <- lmer(PC1 ~ AgCategory + (1|Species), data = birds, REML = FALSE)

d <- expand.grid(
  AgCategory = unique(birds$AgCategory),
  Species = unique(birds$Species))

predictions <- predict(m6, newdata = d, type = "response", se.fit = TRUE) 

d$fit <- predictions$fit

d$lower_CI <- predictions$fit - 1.96 * predictions$se.fit  # Lower CI
d$upper_CI <- predictions$fit + 1.96 * predictions$se.fit  # Upper CI

ggplot(d, aes(x = AgCategory, y = fit)) +
  geom_point(size = 5) +  # Points showing predicted values
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.1,
                col = "black",
                size = 1) +  # Add confidence intervals
  theme_classic() +
  facet_wrap(~ Species, scales = "free_y") +
  labs(x = "Agricultural Intensity Category", 
       y = "Predicted Mean Fattening Index") +
  theme(axis.title.x = element_text(size = 16,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 16,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none") +
  theme(legend.position = "none") +
  geom_hline(size = 1, linetype = "dashed", col = "red", yintercept = 0)

# no heterogeneity
fitted_vals <- fitted(m6)  
residuals_vals <- resid(m6)  
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# Significant difference between low and moderate ag intensity. No effect of capture time.
summary(m21)
confint(m21)
r.squaredGLMM(m21)
# M: 0.121, C: 0.25

m21 <- lmer(PC1 ~ ts.sunrise + AgCategory + (1|Species), data = birds, REML = FALSE)

d <- expand.grid(
  AgCategory = unique(birds$AgCategory),
  Species = unique(birds$Species),
  ts.sunrise = mean(birds$ts.sunrise))

predictions <- predict(m6, newdata = d, type = "response", se.fit = TRUE) 

d$fit <- predictions$fit

d$lower_CI <- predictions$fit - 1.96 * predictions$se.fit  # Lower CI
d$upper_CI <- predictions$fit + 1.96 * predictions$se.fit  # Upper CI

ggplot(d, aes(x = AgCategory, y = fit)) +
  geom_point(size = 3) +  # Points showing predicted values
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.1,
                col = "black",
                size = 0.8) +  # Add confidence intervals
  theme_light() +
  facet_wrap(~ Species, scales = "free_y") +
  labs(x = "Agricultural Intensity Category", 
       y = "Predicted Mean Fattening Index") +
  theme(axis.title.x = element_text(size = 16,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 16,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none") +
  theme(legend.position = "none") +
  geom_hline(size = 1, linetype = "dashed", col = "red", yintercept = 0)

# no heterogeneity
fitted_vals <- fitted(m21)  
residuals_vals <- resid(m21)  
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# No significant effects on fattening index.
summary(m24)
confint(m24)
r.squaredGLMM(m24)
# M: 0.109, C: 0.264

# No significant effects on fattening index.
summary(m55)
confint(m55)
r.squaredGLMM(m55)
# M: 0.146, C: 0.257

# Significant difference between low and moderate ag intensity. No effect of detection on fattening index.
summary(m28)
confint(m28)
r.squaredGLMM(m28)
# M: 0.0971, C: 0.278

# No significant effects on fattening index.
summary(m96)
confint(m96)
r.squaredGLMM(m96)
# M: 0.168, C: 0.265

# no heterogeneity
fitted_vals <- fitted(m96)  
residuals_vals <- resid(m96)  
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# Significant difference between low and moderate ag intensity. 
# No significant effect of capture time or detection on fattening index.
summary(m59)
confint(m59)
r.squaredGLMM(m59)
# M: 0.132, C: 0.263

# no heterogeneity
fitted_vals <- fitted(m59)  
residuals_vals <- resid(m59)  
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# No significant effects on fattening index.
summary(m62)
confint(m62)
r.squaredGLMM(m62)
# M: 0.122, C: 0.274

# ---------------------------------------------------------------------------- #

# Model assumptions ####

m.global <- lmer(PC1 ~ Sex + Event + ts.sunrise + DaysIntoSeason_S + 
                   MigStatus + PercentAg + Detection +
                   MigStatus * DaysIntoSeason_S +
                   Event * ts.sunrise + (1 | Species), REML = FALSE, data = birds.cs,
                 na.action = na.exclude)

# no severe violations of heteroscedasticity
fitted_vals <- fitted(m.global)  
residuals_vals <- resid(m.global)  
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

### Note, the overall residuals do not exhibit heteroscedasticity
### Although individual predictors do, the overall residuals do not
### Heteroscedasticity does not significantly affect the overall model's predictive performance

# CLEAR VIOLATION OF HETEROSCEDASTICITY
boxplot(residuals_vals ~ birds.cs$Event, xlab = "Event", ylab = "Residuals")

# Plot residuals vs MigStatus (pretty good)
boxplot(residuals_vals ~ birds.cs$MigStatus, xlab = "MigStatus", ylab = "Residuals")

# Plot residuals vs Sex (slight violation)
boxplot(residuals_vals ~ birds.cs$Sex, xlab = "Sex", ylab = "Residuals")

# Plot residuals vs Detection (slight violation)
boxplot(residuals_vals ~ birds.cs$Detection, xlab = "Detection", ylab = "Residuals")

# Plot residuals vs Capture Time (pretty good)
plot(birds.cs$ts.sunrise, residuals_vals, xlab = "Capture Time", ylab = "Residuals")

# Plot residuals vs Date (pretty good)
plot(birds.cs$DaysIntoSeason_S, residuals_vals, xlab = "Date", ylab = "Residuals")

# Plot residuals vs Agricultural Intensity (clear violation)
plot(birds.cs$PercentAg, residuals_vals, xlab = "Agricultural Intensity", ylab = "Residuals")




# ---------------------------------------------------------------------------- #

# Removing two outliers in high category does not affect inference

birds.nooutliers <- birds.cs %>% 
  filter(PC1 < 3)

m <- lmer(PC1 ~ AgCategory + (1 | Species), REML = FALSE, data = birds.nooutliers)
m1 <- lmer(PC1 ~ AgCategory + (1 | Species), REML = FALSE, data = birds.cs)

summary(m)
confint(m)
summary(m1)
confint(m1)

ggplot(birds.nooutliers,
       aes(x = AgCategory, y = PC1)) +
  geom_boxplot()

# ---------------------------------------------------------------------------- #

# Plasma Metabolite Profiles ####

# correlation between triglyceride and B-hydroxy levels

# overall
cor.test(birds$Tri, birds$Beta)


cor_results <- birds %>%
  group_by(Species) %>%
  summarise(
    cor_test = list(cor.test(Tri, Beta)),  # Perform correlation test
    .groups = "drop"
  ) %>%
  mutate(
    cor_value = sapply(cor_test, function(x) x$estimate),   # Extract correlation coefficient
    p_value = sapply(cor_test, function(x) x$p.value)       # Extract p-value
  ) %>%
  select(Species, cor_value, p_value)

print(cor_results)

# code above works
leye <- subset(birds, Species == "Lesser Yellowlegs")
cor.test(leye$Tri, leye$Beta)

# PCA for each species
leye <- subset(birds, Species == "Lesser Yellowlegs")
leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

lesa <- subset(birds, Species == "Least Sandpiper")
lesa.cs <- lesa %>%
  mutate(across(where(is.numeric), scale))

lbdo <- subset(birds, Species == "Long-billed Dowitcher")
lbdo.cs <- lbdo %>%
  mutate(across(where(is.numeric), scale))

birds_subset <- lbdo.cs[, c("Tri", "Beta")]

# Remove rows with NAs for PCA
birds_subset_clean <- birds_subset[complete.cases(birds_subset), ]

# Run PCA on the cleaned data
pca_result <- prcomp(birds_subset_clean, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# View the PCA scores (principal components)
pca_scores <- pca_result$x

# Extract the factor score coefficients (loadings)
pca_result$rotation




# correlation between body mass and TRIG & BHOB ####
cor.test(birds$Mass, birds$Beta) # significantly negatively correlated

cor.test(birds$Mass, birds$Tri) # significantly positively correlated

# ---------------------------------------------------------------------------- #

# LAB MEETING DETECTION GRAPH ####

m <- lmer(PC1 ~ Detection + AgCategory + (1 | Species), data = birds, REML = FALSE)

d <- expand.grid(
  AgCategory = c("Moderate"),
  Detection = unique(birds$Detection),
  Species = unique(birds$Species))

d$Species <- factor(d$Species, levels = c("Lesser Yellowlegs", 
                                          "Pectoral Sandpiper", 
                                          "Killdeer", 
                                          "American Avocet", 
                                          "Long-billed Dowitcher", 
                                          "Willet", 
                                          "Wilson's Phalarope", 
                                          "Semipalmated Sandpiper",
                                          "Least Sandpiper"))

predictions <- predict(m, newdata = d, type = "response", se.fit = TRUE) 

d$fit <- predictions$fit

d$lower_CI <- predictions$fit - 1.96 * predictions$se.fit  # Lower CI
d$upper_CI <- predictions$fit + 1.96 * predictions$se.fit  # Upper CI

ggplot(d, aes(x = Detection, y = fit)) +
  geom_point(size = 3) +  # Points showing predicted values
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.1,
                col = "black",
                size = 1) +  # Add confidence intervals
  theme_light() +
  facet_wrap(~ Species) +
  labs(x = NULL, 
       y = "Predicted Fattening Index") +
  theme(axis.title.x = element_text(size = 16,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 16,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        strip.text = element_text(size = 18)) +
  theme(legend.position = "none") +
  geom_hline(size = 1, linetype = "dashed", col = "red", yintercept = 0)

library(lmerTest)
library(MuMIn)

summary(m)
r.squaredGLMM(m)

# LAB MEETING AG GRAPH ####

m <- lmer(PC1 ~ AgCategory + (1 | Species), data = birds, REML = FALSE)

d <- expand.grid(
  AgCategory = unique(birds$AgCategory),
  Species = unique(birds$Species))

d$Species <- factor(d$Species, levels = c("Lesser Yellowlegs", 
                                          "Pectoral Sandpiper", 
                                          "Killdeer", 
                                          "American Avocet", 
                                          "Long-billed Dowitcher", 
                                          "Willet", 
                                          "Wilson's Phalarope", 
                                          "Semipalmated Sandpiper",
                                          "Least Sandpiper"))

predictions <- predict(m, newdata = d, type = "response", se.fit = TRUE) 

d$fit <- predictions$fit

d$lower_CI <- predictions$fit - 1.96 * predictions$se.fit  # Lower CI
d$upper_CI <- predictions$fit + 1.96 * predictions$se.fit  # Upper CI

ggplot(d, aes(x = AgCategory, y = fit)) +
  geom_point(size = 3) +  # Points showing predicted values
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.1,
                col = "black",
                size = 1) +  # Add confidence intervals
  theme_light() +
  facet_wrap(~ Species) +
  labs(x = "Agricultural Intensity", 
       y = "Predicted Fattening Index") +
  theme(axis.title.x = element_text(size = 16,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 16,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        strip.text = element_text(size = 18)) +
  theme(legend.position = "none") +
  geom_hline(size = 1, linetype = "dashed", col = "red", yintercept = 0)

summary(m)
r.squaredGLMM(m)
