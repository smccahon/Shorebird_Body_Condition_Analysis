#---------------------------------------------#
#  Lesser Yellowlegs Pectoral Muscle Analysis #
#           Created 11/11/2024                #          
#          Modified 11/22/2024                #
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
birds <- read.csv("LM_ShorebirdsALLNeg.csv")

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
ggplot(data = leye, aes(x = Detection, y = PecSizeBest)) + 
  theme_classic() + 
  geom_boxplot() + 
  geom_jitter() +
  labs(x = "Neonicotinoid Detection",
       y = "Lesser Yellowlegs Body PecSizeBest")

# PecSizeBest ~ Migration Date
leye$Event <- factor(leye$Event, levels = c("Fall_2021", "Spring_2022", "Fall_2023"))
ggplot(data = leye, aes(x = MigDate, y = PecSizeBest, color = Event)) +
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

# ---------------------------------------------------------------------------- #

# Lesser Yellowlegs Modeling ####

## Univariate Analysis: Date ####

m1 <- lm(PecSizeBest ~ MigDate, data = leye.cs)
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
# Conclusion: try both

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

## Interaction Analysis: LogNeonic & Sampling Event ####

m1 <- lm(PecSizeBest ~ Event * LogNeonic, data = leye.cs)
m2 <- lm(PecSizeBest ~ Event + LogNeonic, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model without interaction performs significantly better

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Detection & Sampling Event ####

m1 <- lm(PecSizeBest ~ Event * Detection, data = leye.cs)
m2 <- lm(PecSizeBest ~ Event + Detection, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Cannot consider interaction because one season has 100% detection

# ---------------------------------------------------------------------------- #

## Single Models ####

m.global.d <- lm(PecSizeBest ~ Sex + Event*ts.sunrise + MigDate + Detection + PercentAg, 
                 data = leye.cs)

m.global.l <- lm(PecSizeBest ~ Sex + Event*ts.sunrise + MigDate + LogNeonic + PercentAg, 
                 data = leye.cs)

m.null <- lm(PecSizeBest ~ 1, data = leye.cs)

m1 <- lm(PecSizeBest ~ Sex, data = leye.cs)
m2 <- lm(PecSizeBest ~ Event, data = leye.cs)
m3 <- lm(PecSizeBest ~ MigDate, data = leye.cs)
m4 <- lm(PecSizeBest ~ ts.sunrise, data = leye.cs)
m5 <- lm(PecSizeBest ~ PercentAg, data = leye.cs)
m6 <- lm(PecSizeBest ~ Detection, data = leye.cs)
m7 <- lm(PecSizeBest ~ LogNeonic, data = leye.cs)

## Additive Models ####

### Two additive combinations ####
m8 <- lm(PecSizeBest ~ Sex + Event, data = leye.cs)
m9 <- lm(PecSizeBest ~ Sex + MigDate, data = leye.cs)
m10 <- lm(PecSizeBest ~ Sex + ts.sunrise, data = leye.cs)
m11 <- lm(PecSizeBest ~ Sex + PercentAg, data = leye.cs)
m12 <- lm(PecSizeBest ~ Sex + Detection, data = leye.cs)
m13 <- lm(PecSizeBest ~ Sex + LogNeonic, data = leye.cs)
m14 <- lm(PecSizeBest ~ Event + MigDate, data = leye.cs)
m15 <- lm(PecSizeBest ~ Event * ts.sunrise, data = leye.cs)
m16 <- lm(PecSizeBest ~ Event + PercentAg, data = leye.cs)
m17 <- lm(PecSizeBest ~ Event + Detection, data = leye.cs)
m18 <- lm(PecSizeBest ~ Event + LogNeonic, data = leye.cs)
m19 <- lm(PecSizeBest ~ MigDate + ts.sunrise, data = leye.cs)
m20 <- lm(PecSizeBest ~ MigDate + PercentAg, data = leye.cs)
m21 <- lm(PecSizeBest ~ MigDate + Detection, data = leye.cs)
m22 <- lm(PecSizeBest ~ MigDate + LogNeonic, data = leye.cs)
m23 <- lm(PecSizeBest ~ ts.sunrise + PercentAg, data = leye.cs)
m24 <- lm(PecSizeBest ~ ts.sunrise + Detection, data = leye.cs)
m25 <- lm(PecSizeBest ~ ts.sunrise + LogNeonic, data = leye.cs)
m26 <- lm(PecSizeBest ~ PercentAg + Detection, data = leye.cs)
m27 <- lm(PecSizeBest ~ PercentAg + LogNeonic, data = leye.cs)

### Three additive combinations ####
m28 <- lm(PecSizeBest ~ Sex + Event + MigDate, data = leye.cs)
m29 <- lm(PecSizeBest ~ Sex + Event * ts.sunrise, data = leye.cs)
m30 <- lm(PecSizeBest ~ Sex + Event + PercentAg, data = leye.cs)
m31 <- lm(PecSizeBest ~ Sex + Event + Detection, data = leye.cs)
m32 <- lm(PecSizeBest ~ Sex + MigDate + ts.sunrise, data = leye.cs)
m33 <- lm(PecSizeBest ~ Sex + MigDate + PercentAg, data = leye.cs)
m34 <- lm(PecSizeBest ~ Sex + MigDate + Detection, data = leye.cs)
m35 <- lm(PecSizeBest ~ Sex + ts.sunrise + PercentAg, data = leye.cs)
m36 <- lm(PecSizeBest ~ Sex + ts.sunrise + Detection, data = leye.cs)
m37 <- lm(PecSizeBest ~ Sex + PercentAg + Detection, data = leye.cs)
m38 <- lm(PecSizeBest ~ Event + MigDate + ts.sunrise + Event * ts.sunrise, data = leye.cs)
m39 <- lm(PecSizeBest ~ Event + MigDate + PercentAg, data = leye.cs)
m40 <- lm(PecSizeBest ~ Event + MigDate + Detection, data = leye.cs)
m41 <- lm(PecSizeBest ~ Event * ts.sunrise + PercentAg, data = leye.cs)
m42 <- lm(PecSizeBest ~ Event * ts.sunrise + Detection, data = leye.cs)
m43 <- lm(PecSizeBest ~ Event + PercentAg + Detection, data = leye.cs)
m44 <- lm(PecSizeBest ~ MigDate + ts.sunrise + PercentAg, data = leye.cs)
m45 <- lm(PecSizeBest ~ MigDate + ts.sunrise + Detection, data = leye.cs)
m46 <- lm(PecSizeBest ~ MigDate + PercentAg + Detection, data = leye.cs)
m47 <- lm(PecSizeBest ~ ts.sunrise + PercentAg + Detection, data = leye.cs)
m48 <- lm(PecSizeBest ~ MigDate + ts.sunrise + LogNeonic, data = leye.cs)
m49 <- lm(PecSizeBest ~ MigDate + PercentAg + LogNeonic, data = leye.cs)
m50 <- lm(PecSizeBest ~ ts.sunrise + PercentAg + LogNeonic, data = leye.cs)
m51 <- lm(PecSizeBest ~ Event * ts.sunrise + LogNeonic, data = leye.cs)
m52 <- lm(PecSizeBest ~ Event + PercentAg + LogNeonic, data = leye.cs)
m53 <- lm(PecSizeBest ~ Event + MigDate + LogNeonic, data = leye.cs)
m54 <- lm(PecSizeBest ~ Sex + ts.sunrise + LogNeonic, data = leye.cs)
m55 <- lm(PecSizeBest ~ Sex + PercentAg + LogNeonic, data = leye.cs)
m56 <- lm(PecSizeBest ~ Sex + MigDate + LogNeonic, data = leye.cs)
m57 <- lm(PecSizeBest ~ Sex + Event + LogNeonic, data = leye.cs)

### Four additive combinations ####
m58 <- lm(PecSizeBest ~ Sex + Event + MigDate + ts.sunrise + Event*ts.sunrise, data = leye.cs)
m59 <- lm(PecSizeBest ~ Sex + Event + MigDate + PercentAg, data = leye.cs)
m60 <- lm(PecSizeBest ~ Sex + Event + MigDate + Detection, data = leye.cs)
m61 <- lm(PecSizeBest ~ Sex + Event * ts.sunrise + PercentAg, data = leye.cs)
m62 <- lm(PecSizeBest ~ Sex + Event * ts.sunrise + Detection, data = leye.cs)
m63 <- lm(PecSizeBest ~ Sex + Event + PercentAg + Detection, data = leye.cs)
m64 <- lm(PecSizeBest ~ Sex + MigDate + ts.sunrise + PercentAg, data = leye.cs)
m65 <- lm(PecSizeBest ~ Sex + MigDate + ts.sunrise + Detection, data = leye.cs)
m66 <- lm(PecSizeBest ~ Sex + MigDate + PercentAg + Detection, data = leye.cs)
m67 <- lm(PecSizeBest ~ Sex + ts.sunrise + PercentAg + Detection, data = leye.cs)
m68 <- lm(PecSizeBest ~ Event + MigDate + ts.sunrise + PercentAg + Event*ts.sunrise, data = leye.cs)
m69 <- lm(PecSizeBest ~ Event + MigDate + ts.sunrise + Detection + Event*ts.sunrise, data = leye.cs)
m70 <- lm(PecSizeBest ~ Event + MigDate + PercentAg + Detection, data = leye.cs)
m71 <- lm(PecSizeBest ~ Event * ts.sunrise + PercentAg + Detection, data = leye.cs)
m72 <- lm(PecSizeBest ~ MigDate + ts.sunrise + PercentAg + Detection, data = leye.cs)
m73 <- lm(PecSizeBest ~ Sex + Event + MigDate + LogNeonic, data = leye.cs)
m74 <- lm(PecSizeBest ~ Sex + Event * ts.sunrise + LogNeonic, data = leye.cs)
m75 <- lm(PecSizeBest ~ Sex + Event + PercentAg + LogNeonic, data = leye.cs)
m76 <- lm(PecSizeBest ~ Sex + MigDate + ts.sunrise + LogNeonic, data = leye.cs)
m77 <- lm(PecSizeBest ~ Sex + MigDate + PercentAg + LogNeonic, data = leye.cs)
m78 <- lm(PecSizeBest ~ Sex + ts.sunrise + PercentAg + LogNeonic, data = leye.cs)
m79 <- lm(PecSizeBest ~ Event + MigDate + ts.sunrise + LogNeonic + Event*ts.sunrise, data = leye.cs)
m80 <- lm(PecSizeBest ~ Event + MigDate + PercentAg + LogNeonic, data = leye.cs)
m81 <- lm(PecSizeBest ~ Event * ts.sunrise + PercentAg + LogNeonic, data = leye.cs)
m82 <- lm(PecSizeBest ~ MigDate + ts.sunrise + PercentAg + LogNeonic, data = leye.cs)

### Five additive combinations ####

m83 <- lm(PecSizeBest ~ Sex + Event + MigDate + ts.sunrise + PercentAg + Event * ts.sunrise, data = leye.cs)
m84 <- lm(PecSizeBest ~ Sex + Event + MigDate + ts.sunrise + Detection + Event * ts.sunrise, data = leye.cs)
m85 <- lm(PecSizeBest ~ Sex + Event + MigDate + PercentAg + Detection, data = leye.cs)
m86 <- lm(PecSizeBest ~ Sex + Event * ts.sunrise + PercentAg + Detection, data = leye.cs)
m87 <- lm(PecSizeBest ~ Sex + MigDate + ts.sunrise + PercentAg + Detection, data = leye.cs)
m88 <- lm(PecSizeBest ~ Event + MigDate + ts.sunrise + PercentAg + Detection + Event*ts.sunrise, data = leye.cs)
m89 <- lm(PecSizeBest ~ Sex + Event + MigDate + ts.sunrise + LogNeonic + Event * ts.sunrise, data = leye.cs)
m90 <- lm(PecSizeBest ~ Sex + Event + MigDate + PercentAg + LogNeonic, data = leye.cs)
m91 <- lm(PecSizeBest ~ Sex + Event * ts.sunrise + PercentAg + LogNeonic, data = leye.cs)
m92 <- lm(PecSizeBest ~ Sex + MigDate + ts.sunrise + PercentAg + LogNeonic, data = leye.cs)
m93 <- lm(PecSizeBest ~ Event + MigDate + ts.sunrise + PercentAg + LogNeonic + Event*ts.sunrise, data = leye.cs)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:93)

models <- mget(model_names)

models$m.null <- m.null
models$m.global.d <- m.global.d
models$m.global.l <- m.global.l

model_names <- c(model_names, "m.null", "m.global.d", "m.global.l")

aictab(models, modnames = model_names)

# ---------------------------------------------------------------------------- #

# Top Model Summaries: All Models ####
m7 <- lm(PecSizeBest ~ LogNeonic, data = leye.cs)
cbind(summary(m7)$coefficients, confint(m7))

# Significant positive relationship between log(concentration) and pectoral muscle size
# Opposite of prediction...

m27 <- lm(PecSizeBest ~ PercentAg + LogNeonic, data = leye.cs)
cbind(summary(m27)$coefficients, confint(m27))

# Significant positive relationship between log(concentration) and pectoral muscle size
# No significant impact of agricultural intensity on pectoral size.

# ---------------------------------------------------------------------------- #

# Graph Top Models ####

# Pectoral ~ LogNeonic
m <- lm(PecSizeBest ~ LogNeonic, data = leye)

d <- expand.grid(LogNeonic = seq(min(leye$LogNeonic), 
                                 max(leye$LogNeonic),
                                 length.out = 1000))

d <- cbind(d, predict(m, newdata = d, interval = "confidence"))

ggplot(d, aes(x = LogNeonic, y = fit)) + 
  geom_line(aes(color = "tan2"), size = 1, show.legend = FALSE) +  
  theme_light() +
  labs(x = "Log(Neonicotinoid Concentration [ug/L])",
       y = expression("Lesser Yellowlegs Breast Muscle Size" ~~~ (mm[score]))) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = LogNeonic),  
              alpha = 0.25, color = NA, show.legend = FALSE,
              fill = "tan1") +
  geom_point(data = leye, aes(x = LogNeonic, y = PecSizeBest), color = "tan2")


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

# Pectoral ~ LogNeonic + Agricultural Intensity

m <- lm(PecSizeBest ~ LogNeonic + PercentAg, data = leye)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg), 
                                 max(leye$PercentAg),
                                 length.out = 1000),
                 LogNeonic = mean(leye$LogNeonic))

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

# ---------------------------------------------------------------------------- #

# Model assumptions ####
m <- lm(PecSizeBest ~ LogNeonic, data = leye.cs)

par(mfrow = c(2,2))
plot(m)

# Some deviation of assumptions...but not super severe? Worth investigating further...
