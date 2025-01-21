#----------------------------#
#  All Species BMI Analysis  #
#     Created 12/11/2024     #          
#    Modified 01/20/2024     #
#----------------------------#

# load packages
library(ggplot2)
library(AICcmodavg)
library(tidyverse)

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

birds$Species <- factor(birds$Species, 
                        levels = c("LesserYellowlegs", "Killdeer", "Willet", 
                                   "PectoralSandpiper", "WIPH", "LeastSandpiper", 
                                   "AmericanAvocet", "LongbilledDowitcher", 
                                   "SemipalmatedSandpiper", 
                                   "MarbledGodwit", "ShortbilledDowitcher",
                                   "GreaterYellowlegs"),
                        labels = c("Lesser Yellowlegs", "Killdeer", "Willet", 
                                   "Pectoral Sandpiper", "Wilson's Phalarope", 
                                   "Least Sandpiper", "American Avocet", 
                                   "Long-billed Dowitcher", "Semipalmated Sandpiper", 
                                   "Marbled Godwit",  "Short-billed Dowitcher",
                                   "Greater Yellowlegs"))

birds$AgCategory <- factor(birds$AgCategory,
                           levels = c("Low", "Moderate", "High"))

# Only include species with at least three individuals that also have mass, tarsus, wing length, and bill length

birds <- birds %>%
  filter(!is.na(Mass))

birds <- birds %>%
  filter(!is.na(Culmen))

birds <- birds %>%
  filter(!is.na(DiagTarsus))

birds <- birds %>%
  filter(!is.na(Wing))

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

birds.cs <- birds %>%
  mutate(across(where(is.numeric), scale))

# ---------------------------------------------------------------------------- #

# Subset data for Lesser Yellowlegs
leye <- subset(birds, Species %in% c("Lesser Yellowlegs"))

leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

leye_female <- subset(leye.cs, Sex == "Female")
leye_male <- subset(leye.cs, Sex == "Male")

leye_female <- leye_female[, c("Wing", "Culmen", "DiagTarsus")]
leye_male <- leye_male[, c("Wing", "Culmen", "DiagTarsus")]

# Subset data for Least Sandpiper
lesa <- subset(birds, Species %in% c("Least Sandpiper"))

lesa.cs <- lesa %>%
  mutate(across(where(is.numeric), scale))

lesa.cs.subset <- lesa.cs[, c("Wing", "Culmen", "DiagTarsus")]

# Subset data for Semipalmated Sandpiper
sesa <- subset(birds, Species %in% c("Semipalmated Sandpiper"))

sesa.cs <- sesa %>%
  mutate(across(where(is.numeric), scale))

sesa.cs.subset <- sesa.cs[, c("Wing", "Culmen", "DiagTarsus")]

# Subset data for Pectoral Sandpiper
pesa <- subset(birds, Species %in% c("Pectoral Sandpiper"))

pesa.cs <- pesa %>%
  mutate(across(where(is.numeric), scale))

pesa.cs.subset <- pesa.cs[, c("Wing", "Culmen", "DiagTarsus")]

# Subset data for Wilson's Phalarope
wiph <- subset(birds, Species %in% c("Wilson's Phalarope"))

wiph.cs <- wiph %>%
  mutate(across(where(is.numeric), scale))

wiph_female <- subset(wiph.cs, Sex == "Female")
wiph_male <- subset(wiph.cs, Sex == "Male")

wiph_female <- wiph_female[, c("Wing", "Culmen", "DiagTarsus")]
wiph_male <- wiph_male[, c("Wing", "Culmen", "DiagTarsus")]

# Subset data for American Avocet
amav <- subset(birds, Species %in% c("American Avocet"))

amav.cs <- amav %>%
  mutate(across(where(is.numeric), scale))

amav.cs.subset <- amav.cs[, c("Wing", "Culmen", "DiagTarsus")]

# Subset data for Killdeer
kill <- subset(birds, Species %in% c("Killdeer"))

kill.cs <- kill %>%
  mutate(across(where(is.numeric), scale))

kill.cs.subset <- kill.cs[, c("Wing", "Culmen", "DiagTarsus")]

# Subset data for Willet
will <- subset(birds, Species %in% c("Willet"))

will.cs <- will %>%
  mutate(across(where(is.numeric), scale))

will_female <- subset(will.cs, Sex == "Female")
will_male <- subset(will.cs, Sex == "Male")

will_female <- will_female[, c("Wing", "Culmen", "DiagTarsus")]
will_male <- will_male[, c("Wing", "Culmen", "DiagTarsus")]

# Subset data for Long-billed Dowitcher
lbdo <- subset(birds, Species %in% c("Long-billed Dowitcher"))

lbdo.cs <- lbdo %>%
  mutate(across(where(is.numeric), scale))

lbdo.cs.subset <- lbdo.cs[, c("Wing", "Culmen", "DiagTarsus")]

# ---------------------------------------------------------------------------- #

# Perform PCA ####
## Lesser Yellowlegs: Females ####

# Run PCA on the subsetted data
pca_result_f <- prcomp(leye_female, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result_f)

# Extract residuals 
PC1_leye_female <- pca_result_f$x[, 1]

# Regress PC1 against body mass
m_leye_female <- lm(PC1_leye_female ~ Mass, data = subset(leye.cs, Sex == "Female"))

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.leye.f <- resid(m_leye_female)

# Add size-corrected mass back into the full dataset
leye.cs$sc.mass <- NA
leye$sc.mass <- NA
leye.cs[leye.cs$Sex == "Female", "sc.mass"] <- sc.mass.leye.f
leye[leye$Sex == "Female", "sc.mass"] <- sc.mass.leye.f

## Lesser Yellowlegs: Males ####

# Run PCA on the subsetted data
pca_result_m <- prcomp(leye_male, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result_m)

# Extract residuals 
PC1_leye_male <- pca_result_m$x[, 1]

# Regress PC1 against body mass
m_leye_male <- lm(PC1_leye_male ~ Mass, data = subset(leye.cs, Sex == "Male"))

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.leye.m <- resid(m_leye_male)

# Add size-corrected mass back into the full datasest
leye.cs[leye.cs$Sex == "Male", "sc.mass"] <- sc.mass.leye.m
leye[leye$Sex == "Male", "sc.mass"] <- sc.mass.leye.m

# Add leye size corrected mass to full dataset

#---------#

## Wilson's Phalaropes: Females ####

# Run PCA on the subsetted data
pca_result_f <- prcomp(wiph_female, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result_f)

# Extract residuals 
PC1_wiph_female <- pca_result_f$x[, 1]

# Regress PC1 against body mass
m_wiph_female <- lm(PC1_wiph_female ~ Mass, data = subset(wiph.cs, Sex == "Female"))

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.wiph.f <- resid(m_wiph_female)

# Add size-corrected mass back into the full dataset
wiph.cs$sc.mass <- NA
wiph$sc.mass <- NA
wiph.cs[wiph.cs$Sex == "Female", "sc.mass"] <- sc.mass.wiph.f
wiph[wiph$Sex == "Female", "sc.mass"] <- sc.mass.wiph.f

## Wilson's Phalaropes: Males ####

# Run PCA on the subsetted data
pca_result_m <- prcomp(wiph_male, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result_m)

# Extract residuals 
PC1_wiph_male <- pca_result_m$x[, 1]

# Regress PC1 against body mass
m_wiph_male <- lm(PC1_wiph_male ~ Mass, data = subset(wiph.cs, Sex == "Male"))

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.wiph.m <- resid(m_wiph_male)

# Add size-corrected mass back into the full datasest
wiph.cs[wiph.cs$Sex == "Male", "sc.mass"] <- sc.mass.wiph.m
wiph[wiph$Sex == "Male", "sc.mass"] <- sc.mass.wiph.m

#---------#

## Willet: Females ####

# Run PCA on the subsetted data
pca_result_f <- prcomp(will_female, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result_f)

# Extract residuals 
PC1_will_female <- pca_result_f$x[, 1]

# Regress PC1 against body mass
m_will_female <- lm(PC1_will_female ~ Mass, data = subset(will.cs, Sex == "Female"))

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.will.f <- resid(m_will_female)

# Add size-corrected mass back into the full dataset
will.cs$sc.mass <- NA
will$sc.mass <- NA
will.cs[will.cs$Sex == "Female", "sc.mass"] <- sc.mass.will.f
will[will$Sex == "Female", "sc.mass"] <- sc.mass.will.f

## Willet: Males ####

# Run PCA on the subsetted data
pca_result_m <- prcomp(will_male, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result_m)

# Extract residuals 
PC1_will_male <- pca_result_m$x[, 1]

# Regress PC1 against body mass
m_will_male <- lm(PC1_will_male ~ Mass, data = subset(will.cs, Sex == "Male"))

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.will.m <- resid(m_will_male)

# Add size-corrected mass back into the full datasest
will.cs[will.cs$Sex == "Male", "sc.mass"] <- sc.mass.will.m
will[will$Sex == "Male", "sc.mass"] <- sc.mass.will.m

#---------#

## Least Sandpiper ####

# Run PCA on the subsetted data
pca_result <- prcomp(lesa.cs.subset, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# Extract residuals 
PC1_lesa <- pca_result$x[, 1]

# Regress PC1 against body mass
m_lesa <- lm(PC1_lesa ~ Mass, data = lesa.cs)

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.lesa <- resid(m_lesa)

# Add size-corrected mass back into the full datasest
lesa.cs$sc.mass <- NA
lesa$sc.mass <- NA
lesa.cs["sc.mass"] <- sc.mass.lesa
lesa["sc.mass"] <- sc.mass.lesa

#---------#

## Semipalmated Sandpiper ####

# Run PCA on the subsetted data
pca_result <- prcomp(sesa.cs.subset, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# Extract residuals 
PC1_sesa <- pca_result$x[, 1]

# Regress PC1 against body mass
m_sesa <- lm(PC1_sesa ~ Mass, data = sesa.cs)

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.sesa <- resid(m_sesa)

# Add size-corrected mass back into the full datasest
sesa.cs$sc.mass <- NA
sesa$sc.mass <- NA
sesa.cs["sc.mass"] <- sc.mass.sesa
sesa["sc.mass"] <- sc.mass.sesa

#---------#

## Pectoral Sandpiper ####

# Run PCA on the subsetted data
pca_result <- prcomp(pesa.cs.subset, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# Extract residuals 
PC1_pesa <- pca_result$x[, 1]

# Regress PC1 against body mass
m_pesa <- lm(PC1_pesa ~ Mass, data = pesa.cs)

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.pesa <- resid(m_pesa)

# Add size-corrected mass back into the full datasest
pesa.cs$sc.mass <- NA
pesa$sc.mass <- NA
pesa.cs["sc.mass"] <- sc.mass.pesa
pesa["sc.mass"] <- sc.mass.pesa

#---------#

## Killdeer ####

# Run PCA on the subsetted data
pca_result <- prcomp(kill.cs.subset, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# Extract residuals 
PC1_kill <- pca_result$x[, 1]

# Regress PC1 against body mass
m_kill <- lm(PC1_kill ~ Mass, data = kill.cs)

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.kill <- resid(m_kill)

# Add size-corrected mass back into the full datasest
kill.cs$sc.mass <- NA
kill$sc.mass <- NA
kill.cs["sc.mass"] <- sc.mass.kill
kill["sc.mass"] <- sc.mass.kill

#---------#

## American Avocet ####

# Run PCA on the subsetted data
pca_result <- prcomp(amav.cs.subset, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# Extract residuals 
PC1_amav <- pca_result$x[, 1]

# Regress PC1 against body mass
m_amav <- lm(PC1_amav ~ Mass, data = amav.cs)

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.amav <- resid(m_amav)

# Add size-corrected mass back into the full datasest
amav.cs$sc.mass <- NA
amav$sc.mass <- NA
amav.cs["sc.mass"] <- sc.mass.amav
amav["sc.mass"] <- sc.mass.amav

#---------#

## Long-billed Dowitcher ####

# Run PCA on the subsetted data
pca_result <- prcomp(lbdo.cs.subset, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# Extract residuals 
PC1_lbdo <- pca_result$x[, 1]

# Regress PC1 against body mass
m_lbdo <- lm(PC1_lbdo ~ Mass, data = lbdo.cs)

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.lbdo <- resid(m_lbdo)

# Add size-corrected mass back into the full datasest
lbdo.cs$sc.mass <- NA
lbdo$sc.mass <- NA
lbdo.cs["sc.mass"] <- sc.mass.lbdo
lbdo["sc.mass"] <- sc.mass.lbdo

# ---------------------------------------------------------------------------- #

# Combine all species together

birds.cs <- rbind(leye.cs, pesa.cs, lbdo.cs, amav.cs, kill.cs, lesa.cs, 
                     will.cs, sesa.cs, wiph.cs)

birds <- rbind(leye, pesa, lbdo, amav, kill, lesa, 
                  will, sesa, wiph)

# ---------------------------------------------------------------------------- #

# Data Visualization ####

# Sex-specific body condition ~ Detection
ggplot(birds.cs %>% filter(!is.na(Detection)),
       aes(x = Detection, y = sc.mass, fill = Detection)) +
  geom_boxplot(outlier.shape = NA) + 
  theme_light() +
  labs(x = "", 
       y = "All Species Body Condition Index") +  
  theme(axis.title.x = element_text(size = 14),
                                    # margin = margin(t = 13)),
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
  scale_x_discrete(labels = c("Detection" = "Neonicotinoid Detection", 
                              "Non-detection" = "Neonicotinoid Non-detection")) +
  geom_hline(yintercept = 0, linetype = "dashed", color ="red", size = 1)

# Body Condition ~ PercentAg ####
ggplot(birds, 
       aes(x = PercentAg, y = sc.mass, color = Species)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed",
             size = 1) +
  geom_point(size = 2) +  
  theme_light() +
  labs(x = "Surrounding Agricultural Intensity", 
       y = "Body Condition Index") +  
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = scales::percent) +
  scale_color_viridis_d(begin = 0.2, 
                       end = 0.8, 
                       option = "inferno",
                       alpha = 0.8)

# for presentation
ggplot(birds.cs, 
       aes(x = AgCategory, y = sc.mass)) + 
  geom_boxplot(outlier.shape = NA, fill = "#F1CCB5") + 
  labs(x = "Surrounding Agricultural Intensity", 
       y = "All Species Body Condition Index") +  
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


# Body Condition ~ Percent Ag ####
ggplot(birds, 
       aes(x = PercentAg, y = sc.mass)) + 
  geom_point(size = 2) + 
  labs(x = "Surrounding Agricultural Intensity", 
       y = "All Species Body Condition Index") +  
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
    breaks = seq(0, 1, by = 0.25),
    labels = scales::percent)

# Body Condition ~ Capture Time ####
ggplot(birds, 
       aes(x = ts.sunrise, y = sc.mass)) + 
  geom_point(size = 2) + 
  labs(x = "Capture Time Since Sunrise (min)", 
       y = "All Species Body Condition Index") +  
  theme_classic() +
  theme(text = element_text(size = 20),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed",
             size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1)

  
  
  
  
  
  
  
  
  
  
m51 <- lm(Mass ~ Event * ts.sunrise + LogNeonic, data = leye)

d <- expand.grid(Event = c("Fall 2021", "Spring 2022", "Fall 2023"),
                 ts.sunrise = seq(min(leye$ts.sunrise), 
                                  max(leye$ts.sunrise), 
                                  length = 1000),
                 LogNeonic = mean(leye$LogNeonic))

d <- cbind(d, predict(m51, newdata = d, interval = "confidence"))

ggplot(d, aes(x = ts.sunrise, y = fit, col = Event)) + 
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Event), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_light() +
  labs(x = "Time of Capture since Sunrise (min)",
       y = "Lesser Yellowlegs Body Mass (g)",
       color = "Sampling Occasion") + 
  theme(legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        legend.position = c(0.98, 0.05),
        legend.justification = c(0.98, 0.05)) +
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_color_viridis_d(begin = 0.2, 
                        end = 0.8, 
                        option = "inferno") +
  scale_fill_viridis_d(begin = 0.2, 
                       end = 0.8, 
                       option = "inferno") +
  geom_point(data = leye, aes(x = ts.sunrise, y = Mass)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 0, y = max(d$fit), 
           label = "Sunrise", angle = 90, 
           vjust = -0.5, hjust = 0.5,
           size = 5)




# ---------------------------------------------------------------------------- #

## RANDOM EFFECTS NOT NEEDED BECAUSE MASS IS ALREADY SCALED TO SPECIES

## Interaction Analysis: Event and Capture Time ####

m1 <- lm(sc.mass ~ Event * ts.sunrise, data = birds.cs)
m2 <- lm(sc.mass ~ Event + ts.sunrise, data = birds.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model without interaction performs the best

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Time

m1 <- lm(sc.mass ~ ts.sunrise, data = birds.cs)
m2 <- lm(sc.mass ~ ts.sunrise + I(ts.sunrise^2), data = birds.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# quadratic relationship performs the best

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Date * Migratory Status ####

m1 <- lm(sc.mass ~ DaysIntoSeason_S + MigStatus, data = birds.cs)
m2 <- lm(sc.mass ~ DaysIntoSeason_S * MigStatus, data = birds.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model without interaction performs significantly better

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Agricultural Intensity ####

m1 <- lm(sc.mass ~ PercentAg, data = birds.cs)
m2 <- lm(sc.mass ~ PercentAg + I(PercentAg^2), data = birds.cs)
m3 <- lm(sc.mass ~ AgCategory, data = birds.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Agricultural intensity as a continuous variable performs significantly better

# ---------------------------------------------------------------------------- #

# INTERACTIONS INCLUDED: None; Agricultural Intensity (continuous); Quadratic

# ---------------------------------------------------------------------------- #

m.global <- lm(sc.mass ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + 
                   MigStatus + PercentAg + Detection, 
                  data = birds.cs)

m.null <- lm(sc.mass ~ 1, data = birds.cs)

## Single Covariate Models ####

m1 <- lm(sc.mass ~ Sex, data = birds.cs)
m2 <- lm(sc.mass ~ Event, data = birds.cs)
m3 <- lm(sc.mass ~ ts.sunrise + I(ts.sunrise^2), data = birds.cs)
m4 <- lm(sc.mass ~ DaysIntoSeason_S, data = birds.cs)
m5 <- lm(sc.mass ~ MigStatus, data = birds.cs)
m6 <- lm(sc.mass ~ PercentAg, data = birds.cs)
m7 <- lm(sc.mass ~ Detection, data = birds.cs)

##  Additive Models ####

### Two additive combinations ####

m8 <- lm(sc.mass ~ Sex + Event, data = birds.cs)
m9 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2), data = birds.cs)
m10 <- lm(sc.mass ~ Sex + DaysIntoSeason_S, data = birds.cs)
m11 <- lm(sc.mass ~ Sex + MigStatus, data = birds.cs)
m12 <- lm(sc.mass ~ Sex + PercentAg, data = birds.cs)
m13 <- lm(sc.mass ~ Sex + Detection, data = birds.cs)
m14 <- lm(sc.mass ~ Event + ts.sunrise + I(ts.sunrise^2), data = birds.cs)
m15 <- lm(sc.mass ~ Event + DaysIntoSeason_S, data = birds.cs)
m16 <- lm(sc.mass ~ Event + MigStatus, data = birds.cs)
m17 <- lm(sc.mass ~ Event + PercentAg, data = birds.cs)
m18 <- lm(sc.mass ~ Event + Detection, data = birds.cs)
m19 <- lm(sc.mass ~ ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S, data = birds.cs)
m20 <- lm(sc.mass ~ ts.sunrise + I(ts.sunrise^2) + MigStatus, data = birds.cs)
m21 <- lm(sc.mass ~ ts.sunrise + I(ts.sunrise^2) + PercentAg, data = birds.cs)
m22 <- lm(sc.mass ~ ts.sunrise + I(ts.sunrise^2) + Detection, data = birds.cs)
m23 <- lm(sc.mass ~ DaysIntoSeason_S + MigStatus, data = birds.cs)
m24 <- lm(sc.mass ~ DaysIntoSeason_S + PercentAg, data = birds.cs)
m25 <- lm(sc.mass ~ DaysIntoSeason_S + Detection, data = birds.cs)
m26 <- lm(sc.mass ~ MigStatus + PercentAg, data = birds.cs)
m27 <- lm(sc.mass ~ MigStatus + Detection, data = birds.cs)
m28 <- lm(sc.mass ~ PercentAg + Detection, data = birds.cs)

### Three Additive Combinations ####

m29 <- lm(sc.mass ~ Sex + Event + ts.sunrise + I(ts.sunrise^2), data = birds.cs)
m30 <- lm(sc.mass ~ Sex + Event + DaysIntoSeason_S, data = birds.cs)
m31 <- lm(sc.mass ~ Sex + Event + MigStatus, data = birds.cs)
m32 <- lm(sc.mass ~ Sex + Event + PercentAg, data = birds.cs)
m33 <- lm(sc.mass ~ Sex + Event + Detection, data = birds.cs)
m34 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S, data = birds.cs)
m35 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2) + MigStatus, data = birds.cs)
m36 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2) + PercentAg, data = birds.cs)
m37 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2) + Detection, data = birds.cs)
m38 <- lm(sc.mass ~ Sex + DaysIntoSeason_S + MigStatus, data = birds.cs)
m39 <- lm(sc.mass ~ Sex + DaysIntoSeason_S + PercentAg, data = birds.cs)
m40 <- lm(sc.mass ~ Sex + DaysIntoSeason_S + Detection, data = birds.cs)
m41 <- lm(sc.mass ~ Sex + MigStatus + PercentAg, data = birds.cs)
m42 <- lm(sc.mass ~ Sex + MigStatus + Detection, data = birds.cs)
m43 <- lm(sc.mass ~ Sex + PercentAg + Detection, data = birds.cs)
m44 <- lm(sc.mass ~ Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S, data = birds.cs)
m45 <- lm(sc.mass ~ Event + ts.sunrise + I(ts.sunrise^2) + MigStatus, data = birds.cs)
m46 <- lm(sc.mass ~ Event + ts.sunrise + I(ts.sunrise^2) + PercentAg, data = birds.cs)
m47 <- lm(sc.mass ~ Event + ts.sunrise + I(ts.sunrise^2) + Detection, data = birds.cs)
m48 <- lm(sc.mass ~ Event + DaysIntoSeason_S + MigStatus, data = birds.cs)
m49 <- lm(sc.mass ~ Event + DaysIntoSeason_S + PercentAg, data = birds.cs)
m50 <- lm(sc.mass ~ Event + DaysIntoSeason_S + Detection, data = birds.cs)
m51 <- lm(sc.mass ~ Event + MigStatus + PercentAg, data = birds.cs)
m52 <- lm(sc.mass ~ Event + MigStatus + Detection, data = birds.cs)
m53 <- lm(sc.mass ~ Event + PercentAg + Detection, data = birds.cs)
m54 <- lm(sc.mass ~ ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus, data = birds.cs)
m55 <- lm(sc.mass ~ ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + PercentAg, data = birds.cs)
m56 <- lm(sc.mass ~ ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + Detection, data = birds.cs)
m57 <- lm(sc.mass ~ ts.sunrise + I(ts.sunrise^2) + MigStatus + PercentAg, data = birds.cs)
m58 <- lm(sc.mass ~ ts.sunrise + I(ts.sunrise^2) + MigStatus + Detection, data = birds.cs)
m59 <- lm(sc.mass ~ ts.sunrise + I(ts.sunrise^2) + PercentAg + Detection, data = birds.cs)
m60 <- lm(sc.mass ~ DaysIntoSeason_S + MigStatus + PercentAg, data = birds.cs)
m61 <- lm(sc.mass ~ DaysIntoSeason_S + MigStatus + Detection, data = birds.cs)
m62 <- lm(sc.mass ~ DaysIntoSeason_S + PercentAg + Detection, data = birds.cs)
m63 <- lm(sc.mass ~ MigStatus + PercentAg + Detection, data = birds.cs)

### Four Additive Combinations ####

m64 <- lm(sc.mass ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S, data = birds.cs)
m65 <- lm(sc.mass ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + MigStatus, data = birds.cs)
m66 <- lm(sc.mass ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + PercentAg, data = birds.cs)
m67 <- lm(sc.mass ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + Detection, data = birds.cs)
m68 <- lm(sc.mass ~ Sex + Event + DaysIntoSeason_S + MigStatus, data = birds.cs)
m69 <- lm(sc.mass ~ Sex + Event + DaysIntoSeason_S + PercentAg, data = birds.cs)
m70 <- lm(sc.mass ~ Sex + Event + DaysIntoSeason_S + Detection, data = birds.cs)
m71 <- lm(sc.mass ~ Sex + Event + MigStatus + PercentAg, data = birds.cs)
m72 <- lm(sc.mass ~ Sex + Event + MigStatus + Detection, data = birds.cs)
m73 <- lm(sc.mass ~ Sex + Event + PercentAg + Detection, data = birds.cs)
m74 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus, data = birds.cs)
m75 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + PercentAg, data = birds.cs)
m76 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + Detection, data = birds.cs)
m77 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2) + MigStatus + PercentAg, data = birds.cs)
m78 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2) + MigStatus + Detection, data = birds.cs)
m79 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2) + PercentAg + Detection, data = birds.cs)
m80 <- lm(sc.mass ~ Sex + DaysIntoSeason_S + MigStatus + PercentAg, data = birds.cs)
m81 <- lm(sc.mass ~ Sex + DaysIntoSeason_S + MigStatus + Detection, data = birds.cs)
m82 <- lm(sc.mass ~ Sex + DaysIntoSeason_S + PercentAg + Detection, data = birds.cs)
m83 <- lm(sc.mass ~ Sex + MigStatus + PercentAg + Detection, data = birds.cs)
m84 <- lm(sc.mass ~ Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus, data = birds.cs)
m85 <- lm(sc.mass ~ Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + PercentAg, data = birds.cs)
m86 <- lm(sc.mass ~ Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + Detection, data = birds.cs)
m87 <- lm(sc.mass ~ Event + ts.sunrise + I(ts.sunrise^2) + MigStatus + PercentAg, data = birds.cs)
m88 <- lm(sc.mass ~ Event + ts.sunrise + I(ts.sunrise^2) + MigStatus + Detection, data = birds.cs)
m89 <- lm(sc.mass ~ Event + ts.sunrise + I(ts.sunrise^2) + PercentAg + Detection, data = birds.cs)
m90 <- lm(sc.mass ~ Event + DaysIntoSeason_S + MigStatus + PercentAg, data = birds.cs)
m91 <- lm(sc.mass ~ Event + DaysIntoSeason_S + MigStatus + Detection, data = birds.cs)
m92 <- lm(sc.mass ~ Event + DaysIntoSeason_S + PercentAg + Detection, data = birds.cs)
m93 <- lm(sc.mass ~ Event + MigStatus + PercentAg + Detection, data = birds.cs)
m94 <- lm(sc.mass ~ ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus + PercentAg, data = birds.cs)
m95 <- lm(sc.mass ~ ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus + Detection, data = birds.cs)
m96 <- lm(sc.mass ~ ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + PercentAg + Detection, data = birds.cs)
m97 <- lm(sc.mass ~ ts.sunrise + I(ts.sunrise^2) + MigStatus + PercentAg + Detection, data = birds.cs)
m98 <- lm(sc.mass ~ DaysIntoSeason_S + MigStatus + PercentAg + Detection, data = birds.cs)

### Five Additive Combinations ####

m99 <- lm(sc.mass ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus, data = birds.cs)
m100 <- lm(sc.mass ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + PercentAg, data = birds.cs)
m101 <- lm(sc.mass ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + Detection, data = birds.cs)
m102 <- lm(sc.mass ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + MigStatus + PercentAg, data = birds.cs)
m103 <- lm(sc.mass ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + MigStatus + Detection, data = birds.cs)
m104 <- lm(sc.mass ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + PercentAg + Detection, data = birds.cs)
m105 <- lm(sc.mass ~ Sex + Event + DaysIntoSeason_S + MigStatus + PercentAg, data = birds.cs)
m106 <- lm(sc.mass ~ Sex + Event + DaysIntoSeason_S + MigStatus + Detection, data = birds.cs)
m107 <- lm(sc.mass ~ Sex + Event + DaysIntoSeason_S + PercentAg + Detection, data = birds.cs)
m108 <- lm(sc.mass ~ Sex + Event + MigStatus + PercentAg + Detection, data = birds.cs)
m109 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus + PercentAg, data = birds.cs)
m110 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus + Detection, data = birds.cs)
m111 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + PercentAg + Detection, data = birds.cs)
m112 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2) + MigStatus + PercentAg + Detection, data = birds.cs)
m113 <- lm(sc.mass ~ Sex + DaysIntoSeason_S + MigStatus + PercentAg + Detection, data = birds.cs)
m114 <- lm(sc.mass ~ Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus + PercentAg, data = birds.cs)
m115 <- lm(sc.mass ~ Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus + Detection, data = birds.cs)
m116 <- lm(sc.mass ~ Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + PercentAg + Detection, data = birds.cs)
m117 <- lm(sc.mass ~ Event + ts.sunrise + I(ts.sunrise^2) + MigStatus + PercentAg + Detection, data = birds.cs)
m118 <- lm(sc.mass ~ Event + DaysIntoSeason_S + MigStatus + PercentAg + Detection, data = birds.cs)
m119 <- lm(sc.mass ~ ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus + PercentAg + Detection, data = birds.cs)

### Six Additive Combinations ####

m120 <- lm(sc.mass ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus + PercentAg, data = birds.cs)
m121 <- lm(sc.mass ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus + Detection, data = birds.cs)
m122 <- lm(sc.mass ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + PercentAg + Detection, data = birds.cs)
m123 <- lm(sc.mass ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + MigStatus + PercentAg + Detection, data = birds.cs)
m124 <- lm(sc.mass ~ Sex + Event + DaysIntoSeason_S + MigStatus + PercentAg + Detection, data = birds.cs)
m125 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus + PercentAg + Detection, data = birds.cs)
m126 <- lm(sc.mass ~ Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + MigStatus + PercentAg + Detection, data = birds.cs)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:126)

models <- mget(model_names)

models$m.null <- m.null
models$m.global <- m.global

model_names <- c(model_names, "m.null", "m.global")

aictab(models, modnames = model_names)

# NEW RESULTS 2025-01-20 ####



# ---------------------------------------------------------------------------- #

# Top Model Summaries & Assumptions ####

# capture time significant
summary(m37)
confint(m37)

plot(predict(m37), rstudent(m37))

# severe violation of homoscedasticity
# DID NOT APPROPRIATELY ADDRESS BUT NOT A MAIN QUESTION OF INTEREST
# HETEROSCEDASTICITY MAY SIMPLY BE A PROBLEM OF NO RELATIONSHIPS
fitted_vals <- fitted(m37)  
residuals_vals <- resid(m37)  
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# Plot residuals vs Sex (pretty good)
boxplot(residuals_vals ~ birds.c$Sex, xlab = "Sex", ylab = "Residuals")

# Plot residuals vs Detection (pretty good)
boxplot(residuals_vals ~ birds.c$Detection, xlab = "Detection", ylab = "Residuals")

# Plot residuals vs Capture Time (pretty good)
plot(birds.c$ts.sunrise, residuals_vals, xlab = "Capture Time", ylab = "Residuals")

#other assumptions
par(mfrow = c(2,2))
plot(m37)

# Weighted least squares
birds.c <- birds.cs[complete.cases(birds.cs$Sex, birds.cs$Event, birds.cs$ts.sunrise, 
                                   birds.cs$DaysIntoSeason_S, birds.cs$MigStatus, 
                                   birds.cs$PercentAg, birds.cs$Detection), ]

m37 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2) + Detection, data = birds.c)

fitted_vals <- fitted(m37)
weights <- 1 / fitted_vals^2

m37 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2) + Detection, data = birds.c,
          weights = weights)

par(mfrow = c(2,2))
plot(m37)

#----------------------#

# capture time significant
summary(m22)
confint(m22)

par(mfrow = c(2,2))
plot(m22)

# capture time significant
summary(m76)
confint(m76)

par(mfrow = c(2,2))
plot(m76)

# nothing significant
summary(m13)
confint(m13)

# capture time significant
summary(m56)
confint(m56)

# capture time significant
summary(m78)
confint(m78)

# detection not significant
summary(m7)
confint(m7)

# ---------------------------------------------------------------------------- #

# Graph Top Models
m37 <- lm(sc.mass ~ Sex + ts.sunrise + I(ts.sunrise^2) + Detection, data = birds)

d <- expand.grid(Detection = c("Detection", "Non-detection"),
                 Sex = c("Female"),
                 ts.sunrise = seq(min(birds$ts.sunrise), 
                                  max(birds$ts.sunrise), 
                                  length = 2000))

d <- cbind(d, predict(m37, newdata = d, interval = "confidence"))

ggplot(d, aes(x = ts.sunrise, y = fit)) + 
  geom_line(color = "black", size = 0) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_light() +
  labs(x = "Time of Capture since Sunrise (min)",
       y = "Body Condition Index") + 
  theme(legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        legend.position = c(0.98, 0.05),
        legend.justification = c(0.98, 0.05)) +
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  geom_point(data = birds, aes(x = ts.sunrise, y = sc.mass)) +
  geom_vline(xintercept = 0, linetype = "dashed", 
             color = "black", size = 1) +
  annotate("text", x = 0, y = 3, 
           label = "Sunrise", angle = 90, 
           vjust = -0.5, hjust = 0.5,
           size = 5) +
  geom_hline(linetype = "dashed", color = "red", size = 1, yintercept = 0)

# ---------------------------------------------------------------------------- #

# Model assumptions ####

m.global <- lm(sc.mass ~ Sex + Event + ts.sunrise + I(ts.sunrise^2) + DaysIntoSeason_S + 
                 MigStatus + PercentAg + Detection, 
               data = birds.cs)

# slight violation of homoscedasticity
fitted_vals <- fitted(m.global)  
residuals_vals <- resid(m.global)  
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# satisfied
hist(residuals_vals, xlab = "Residuals", main = "")

birds.c <- birds.cs[complete.cases(birds.cs$Sex, birds.cs$Event, birds.cs$ts.sunrise, 
                                             birds.cs$DaysIntoSeason_S, birds.cs$MigStatus, 
                                             birds.cs$PercentAg, birds.cs$Detection), ]

# Plot residuals vs Event (pretty good)
boxplot(residuals_vals ~ birds.c$Event, xlab = "Event", ylab = "Residuals")

# Plot residuals vs MigStatus (slight violation)
boxplot(residuals_vals ~ birds.c$MigStatus, xlab = "MigStatus", ylab = "Residuals")

# Plot residuals vs Sex (pretty good)
boxplot(residuals_vals ~ birds.c$Sex, xlab = "Sex", ylab = "Residuals")

# Plot residuals vs Detection (pretty good)
boxplot(residuals_vals ~ birds.c$Detection, xlab = "Detection", ylab = "Residuals")

# Plot residuals vs Capture Time (pretty good)
plot(birds.c$ts.sunrise, residuals_vals, xlab = "Capture Time", ylab = "Residuals")

# Plot residuals vs Date (pretty good)
plot(birds.cs$DaysIntoSeason_S, residuals_vals, xlab = "Date", ylab = "Residuals")

# Plot residuals vs Agricultural Intensity (violation)
plot(birds.c$PercentAg, residuals_vals, xlab = "Agricultural Intensity", ylab = "Residuals")

# ---------------------------------------------------------------------------- #

# TOP MODEL PLOTS FOR LAB MEETING

## LAB MEETING PLOT DETECTION ####
m <- lm(sc.mass ~ Detection + ts.sunrise + I(ts.sunrise^2) + Sex, data = birds)

d <- expand.grid(Detection = c("Detection", "Non-detection"),                    
                 ts.sunrise = mean(leye$ts.sunrise),
                 Sex = "Female") 

predictions <- predict(m, newdata = d, se.fit = TRUE)

d$predicted_Mass <- predictions$fit

d$lower_CI <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upper_CI <- d$predicted_Mass + 1.96 * predictions$se.fit

ggplot(d, aes(x = Detection, y = predicted_Mass, color = Detection)) +
  geom_point(size = 5, col = "black") +  # Points showing predicted values
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.1,
                col = "black",
                size = 1) +  # Add confidence intervals
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
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "red", size = 1)

## LAB MEETING PLOT AG ####
m <- lm(sc.mass ~ PercentAg + Detection + ts.sunrise + I(ts.sunrise^2) + Sex, data = birds)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg), 
                                 max(leye$PercentAg), 
                                 length = 1000),
                 Detection = c("Detection"),
                 Sex = "Female",                    
                 ts.sunrise = mean(leye$ts.sunrise)) 

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
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "red", size = 1)
