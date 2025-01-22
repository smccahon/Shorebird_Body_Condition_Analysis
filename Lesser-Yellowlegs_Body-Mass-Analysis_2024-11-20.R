#---------------------------------------#
# Lesser Yellowlegs Body Mass Analysis  #
#         Created 11/11/2024            #          
#        Modified 01/22/2025            #
#---------------------------------------#

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


# Standardize continuous variables
leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

# ---------------------------------------------------------------------------- #

# Data Visualization ####

# Mass ~ Concentration
ggplot(data = leye, aes(x = OverallNeonic, y = Mass)) + 
  theme_classic() + 
  geom_point() + 
  labs(x = "Neonicotinoid Concentration (ug/L)", 
       y = "Lesser Yellowlegs Body Mass")

# Mass ~ Log(Concentration)
ggplot(data = leye, aes(x = LogNeonic, y = Mass)) +
  theme_classic() +
  geom_point() + 
  labs(x = "Log(Concentration)",
       y = "Lesser Yellowlegs Body Mass")

# Mass ~ Detection
ggplot(data = leye, aes(x = Detection, y = Mass)) + 
  theme_classic() + 
  geom_boxplot() + 
  geom_jitter() +
  labs(x = "Neonicotinoid Detection",
       y = "Lesser Yellowlegs Body Mass")

# Mass ~ Migration Date
leye$Event <- factor(leye$Event, levels = c("Fall_2021", "Spring_2022", "Fall_2023"))
ggplot(data = leye, aes(x = MigDate, y = Mass, color = Event)) +
  theme_classic() +
  geom_point() +
  labs(x = "Date into Migration Season",
       y = "Lesser Yellowlegs Body Mass",
       color = "Sampling Event") +
  scale_color_manual(
    values = c("orange", "purple", "green"),
    labels = c("Fall 2021", "Spring 2022", "Fall 2023")
  )

# Mass ~ Julian
ggplot(data = leye, aes(x = Julian, y = Mass, color = Event)) +
  theme_classic() +
  geom_point() +
  labs(x = "Julian Day",
       y = "Lesser Yellowlegs Body Mass",
       color = "Sampling Event") +
  scale_color_manual(
    values = c("orange", "purple", "green"),
    labels = c("Fall 2021", "Spring 2022", "Fall 2023")
  )

# Mass ~ Sex
ggplot(data = leye, aes(x = Sex, y = Mass)) + 
  theme_classic() + 
  geom_boxplot() + 
  geom_jitter() +
  labs(x = "Sex",
       y = "Lesser Yellowlegs Body Mass")

# Mass ~ Sampling Event
ggplot(data = leye, aes(x = Event, y = Mass, color = Event)) +
  theme_classic() +
  geom_boxplot() +
  geom_jitter() + 
  labs(x = "Sampling Event",
       y = "Lesser Yellowlegs Body Mass",
       color = "Sampling Event") +
  scale_color_manual(
    values = c("orange", "purple", "green"),
    labels = c("Fall 2021", "Spring 2022", "Fall 2023")
  ) + 
  scale_x_discrete(
    labels = c("Fall_2021" = "Fall 2021", "Spring_2022" = "Spring 2022", "Fall_2023" = "Fall 2023")
  )

# Mass ~ Capture Time
ggplot(data = leye, aes(x = ts.sunrise, y = Mass)) +
  theme_classic() +
  geom_smooth(method = "lm",  color = "blue") +
  geom_point() +
  labs(x = "Capture Time (min.)",
       y = "Lesser Yellowlegs Body Mass (g)")

# Mass ~ Capture Time * Event
ggplot(data = leye, aes(x = ts.sunrise, y = Mass, col = Event)) +
  theme_classic() +
  geom_smooth(method = "lm",  color = "blue") +
  geom_point() +
  labs(x = "Capture Time (min.)",
       y = "Lesser Yellowlegs Body Mass (g)")



# ---------------------------------------------------------------------------- #

# Lesser Yellowlegs Modeling ####

## Univariate Analysis: Date ####

m1 <- lm(Mass ~ DaysIntoSeason_S, data = leye.cs)
m2 <- lm(Mass ~ Julian, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Date into trapping season performs significantly better than Julian date

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Neonicotinoids ####

m1 <- lm(Mass ~ Detection, data = leye.cs)
m2 <- lm(Mass ~ OverallNeonic, data = leye.cs)
m3 <- lm(Mass ~ LogNeonic, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Detection and Log transformation of concentrations perform significantly better
# Conclusion: include detection only.

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Agricultural Intensity ####

m1 <- lm(Mass ~ AgCategory, data = leye.cs)
m2 <- lm(Mass ~ PercentAg, data = leye.cs)
m3 <- lm(Mass ~ LogAg, data = leye.cs)
m4 <- lm(Mass ~ PercentAg + I(PercentAg^2), data = leye.cs)

model_names <- paste0("m", 1:4)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Percent Ag as a linear, continuous variable performs the best (no transformation)

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Event and Capture Time ####

m1 <- lm(Mass ~ Event * ts.sunrise, data = leye.cs)
m2 <- lm(Mass ~ Event + I(ts.sunrise^2), data = leye.cs)
m3 <- lm(Mass ~ Event + ts.sunrise, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with interaction by far performs the best.

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Detection & Sampling Event ####

m1 <- lm(Mass ~ Event * Detection, data = leye.cs)
m2 <- lm(Mass ~ Event + Detection, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model without interaction performs significantly better

# ---------------------------------------------------------------------------- #

# INTERACTIONS INCLUDED: Event * Time, PercentAg

# ---------------------------------------------------------------------------- #

## Single Models ####

m.global <- lm(Mass ~ Sex + Event * ts.sunrise + DaysIntoSeason_S + Detection + PercentAg + Age, 
               data = leye.cs)

m.null <- lm(Mass ~ 1, data = leye.cs)

m1 <- lm(Mass ~ Sex, data = leye.cs)
m2 <- lm(Mass ~ Event, data = leye.cs)
m3 <- lm(Mass ~ DaysIntoSeason_S, data = leye.cs)
m4 <- lm(Mass ~ ts.sunrise, data = leye.cs)
m5 <- lm(Mass ~ PercentAg, data = leye.cs)
m6 <- lm(Mass ~ Detection, data = leye.cs)
m7 <- lm(Mass ~ Age, data = leye.cs)

## Additive Models ####

### Two additive combinations ####
m8 <- lm(Mass ~ Sex + Event, data = leye.cs)
m9 <- lm(Mass ~ Sex + DaysIntoSeason_S, data = leye.cs)
m10 <- lm(Mass ~ Sex + ts.sunrise, data = leye.cs)
m11 <- lm(Mass ~ Sex + PercentAg, data = leye.cs)
m12 <- lm(Mass ~ Sex + Detection, data = leye.cs)
m13 <- lm(Mass ~ Event + DaysIntoSeason_S, data = leye.cs)
m14 <- lm(Mass ~ Event * ts.sunrise, data = leye.cs)
m15 <- lm(Mass ~ Event + PercentAg, data = leye.cs)
m16 <- lm(Mass ~ DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m17 <- lm(Mass ~ DaysIntoSeason_S + PercentAg, data = leye.cs)
m18 <- lm(Mass ~ DaysIntoSeason_S + Detection, data = leye.cs)
m19 <- lm(Mass ~ ts.sunrise + PercentAg, data = leye.cs)
m20 <- lm(Mass ~ ts.sunrise + Detection, data = leye.cs)
m21 <- lm(Mass ~ PercentAg + Detection, data = leye.cs)
m22 <- lm(Mass ~ Age + Event, data = leye.cs)
m23 <- lm(Mass ~ Age + DaysIntoSeason_S, data = leye.cs)
m24 <- lm(Mass ~ Age + ts.sunrise, data = leye.cs)
m25 <- lm(Mass ~ Age + PercentAg, data = leye.cs)
m26 <- lm(Mass ~ Age + Detection, data = leye.cs)

### Three additive combinations ####
m27 <- lm(Mass ~ Sex + Event + DaysIntoSeason_S, data = leye.cs)
m28 <- lm(Mass ~ Sex + Event * ts.sunrise, data = leye.cs)
m29 <- lm(Mass ~ Sex + Event + PercentAg, data = leye.cs)
m30 <- lm(Mass ~ Sex + Event + Detection, data = leye.cs)
m31 <- lm(Mass ~ Sex + DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m32 <- lm(Mass ~ Sex + DaysIntoSeason_S + PercentAg, data = leye.cs)
m33 <- lm(Mass ~ Sex + DaysIntoSeason_S + Detection, data = leye.cs)
m34 <- lm(Mass ~ Sex + ts.sunrise + PercentAg, data = leye.cs)
m35 <- lm(Mass ~ Sex + ts.sunrise + Detection, data = leye.cs)
m36 <- lm(Mass ~ Sex + PercentAg + Detection, data = leye.cs)
m37 <- lm(Mass ~ Event + DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m38 <- lm(Mass ~ Event + DaysIntoSeason_S + PercentAg, data = leye.cs)
m39 <- lm(Mass ~ Event + DaysIntoSeason_S + Detection, data = leye.cs)
m40 <- lm(Mass ~ Event * ts.sunrise + PercentAg, data = leye.cs)
m41 <- lm(Mass ~ Event * ts.sunrise + Detection, data = leye.cs)
m42 <- lm(Mass ~ Event + PercentAg + Detection, data = leye.cs)
m43 <- lm(Mass ~ DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m44 <- lm(Mass ~ DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m45 <- lm(Mass ~ DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m46 <- lm(Mass ~ ts.sunrise + PercentAg + Detection, data = leye.cs)
m47 <- lm(Mass ~ Age + Event + DaysIntoSeason_S, data = leye.cs)
m48 <- lm(Mass ~ Age + Event * ts.sunrise, data = leye.cs)
m49 <- lm(Mass ~ Age + Event + PercentAg, data = leye.cs)
m50 <- lm(Mass ~ Age + Event + Detection, data = leye.cs)
m51 <- lm(Mass ~ Age + DaysIntoSeason_S + ts.sunrise, data = leye.cs)
m52 <- lm(Mass ~ Age + DaysIntoSeason_S + PercentAg, data = leye.cs)
m53 <- lm(Mass ~ Age + DaysIntoSeason_S + Detection, data = leye.cs)
m54 <- lm(Mass ~ Age + ts.sunrise + PercentAg, data = leye.cs)
m55 <- lm(Mass ~ Age + ts.sunrise + Detection, data = leye.cs)
m56 <- lm(Mass ~ Age + PercentAg + Detection, data = leye.cs)

### Four additive combinations ####
m57 <- lm(Mass ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + Event * ts.sunrise, data = leye.cs)
m58 <- lm(Mass ~ Sex + Event + DaysIntoSeason_S + PercentAg, data = leye.cs)
m59 <- lm(Mass ~ Sex + Event + DaysIntoSeason_S + Detection, data = leye.cs)
m60 <- lm(Mass ~ Sex + Event * ts.sunrise + PercentAg, data = leye.cs)
m61 <- lm(Mass ~ Sex + Event * ts.sunrise + Detection, data = leye.cs)
m62 <- lm(Mass ~ Sex + Event + PercentAg + Detection, data = leye.cs)
m63 <- lm(Mass ~ Sex + DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m64 <- lm(Mass ~ Sex + DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m65 <- lm(Mass ~ Sex + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m66 <- lm(Mass ~ Sex + ts.sunrise + PercentAg + Detection, data = leye.cs)
m67 <- lm(Mass ~ Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Event * ts.sunrise, data = leye.cs)
m68 <- lm(Mass ~ Event + DaysIntoSeason_S + ts.sunrise + Detection + Event * ts.sunrise, data = leye.cs)
m69 <- lm(Mass ~ Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m70 <- lm(Mass ~ Event * ts.sunrise + PercentAg + Detection, data = leye.cs)
m71 <- lm(Mass ~ DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)
m72 <- lm(Mass ~ Age + Event + DaysIntoSeason_S + ts.sunrise + Event * ts.sunrise, data = leye.cs)
m73 <- lm(Mass ~ Age + Event + DaysIntoSeason_S + PercentAg, data = leye.cs)
m74 <- lm(Mass ~ Age + Event + DaysIntoSeason_S + Detection, data = leye.cs)
m75 <- lm(Mass ~ Age + Event * ts.sunrise + PercentAg, data = leye.cs)
m76 <- lm(Mass ~ Age + Event * ts.sunrise + Detection, data = leye.cs)
m77 <- lm(Mass ~ Age + Event + PercentAg + Detection, data = leye.cs)
m78 <- lm(Mass ~ Age + DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs)
m79 <- lm(Mass ~ Age + DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs)
m80 <- lm(Mass ~ Age + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m81 <- lm(Mass ~ Age + ts.sunrise + PercentAg + Detection, data = leye.cs)

### Five additive combinations ####

m82 <- lm(Mass ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Event * ts.sunrise, data = leye.cs)
m83 <- lm(Mass ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + Detection + Event * ts.sunrise, data = leye.cs)
m84 <- lm(Mass ~ Sex + Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m85 <- lm(Mass ~ Sex + Event * ts.sunrise + PercentAg + Detection, data = leye.cs)
m86 <- lm(Mass ~ Sex + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)
m87 <- lm(Mass ~ Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection + Event * ts.sunrise, data = leye.cs)
m88 <- lm(Mass ~ Age + Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Event * ts.sunrise, data = leye.cs)
m89 <- lm(Mass ~ Age + Event + DaysIntoSeason_S + ts.sunrise + Detection + Event * ts.sunrise, data = leye.cs)
m90 <- lm(Mass ~ Age + Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs)
m91 <- lm(Mass ~ Age + Event * ts.sunrise + PercentAg + Detection, data = leye.cs)
m92 <- lm(Mass ~ Age + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:92)

models <- mget(model_names)

models$m.null <- m.null
models$m.global <- m.global

model_names <- c(model_names, "m.null", "m.global")

aictab(models, modnames = model_names)

# NEW RESULTS AFTER INCORPORATING AGE 2025-01-22 ####
confint(m48) # event * time significant
confint(m14) # event * time significant
confint(m41) # event * time significant
confint(m76) # event * time significant
confint(m72) # event * time significant

# ---------------------------------------------------------------------------- #

### Top Model Summaries ####

cbind(summary(m13)$coefficients, confint(m13))
cbind(summary(m36)$coefficients, confint(m36))
cbind(summary(m35)$coefficients, confint(m35))

### Goodness of Fit ####

summary(m13)
summary(m36)
summary(m35)

# ---------------------------------------------------------------------------- #

# Graph Top Models ####

## Mass ~ Event * Time + LogNeonic ####

# Hold LogNeonic constant
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


# Hold Event * Time constant
m51 <- lm(Mass ~ Event * ts.sunrise + LogNeonic, data = leye)

d <- expand.grid(Event = c("Fall 2023"),
                 LogNeonic = seq(min(leye$LogNeonic), 
                                  max(leye$LogNeonic), 
                                  length = 1000),
                 ts.sunrise = mean(leye$ts.sunrise))

d <- cbind(d, predict(m51, newdata = d, interval = "confidence"))

ggplot(d, aes(x = LogNeonic, y = fit)) + 
  geom_line(size = 0.8,
            col = "darkgoldenrod4") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, color = NA, show.legend = FALSE,
              fill = "darkgoldenrod1") +
  theme_light() +
  labs(x = "Log(Neonicotinoid Concentration [ug/L])",
       y = "Lesser Yellowlegs Body Mass (g)") +
  theme(axis.title.x = element_text(size = 14,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 14,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  geom_point(data = leye, aes(x = LogNeonic, y = Mass), col = "darkgoldenrod3")


## Mass ~ Event * Time: THIS ONE FOR THESIS ####
m15 <- lm(Mass ~ Event * ts.sunrise, data = leye)

d_list <- list()

# Loop through each event and create a prediction grid for each event
for (event in unique(leye$Event)) {
  event_data <- leye[leye$Event == event, ]
  d_event <- expand.grid(Event = event,
                         ts.sunrise = seq(min(event_data$ts.sunrise), 
                                          max(event_data$ts.sunrise), 
                                          length = 1000))
  d_event <- cbind(d_event, predict(m15, newdata = d_event, interval = "confidence"))
  d_list[[event]] <- d_event
}

d <- expand.grid(Event = event,
                 ts.sunrise = seq(min(event_data$ts.sunrise), 
                                  max(event_data$ts.sunrise), 
                                  length = 1000))

d <- cbind(d_event, predict(m15, newdata = d_event, interval = "confidence"))

d <- do.call(rbind, d_list)

ggplot(d, aes(x = ts.sunrise, y = fit, col = Event)) + 
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Event), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_classic() +
  labs(x = "Time of Capture since Sunrise (min)",
       y = "Lesser Yellowlegs Body Mass (g)",
       color = "Sampling Occasion") + 
  theme(legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1)) +
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

## Mass ~ Event * Time + Detection: THIS ONE FOR THESIS ####
## LAB MEETING DETECTION GRAPH SPRING 2025 ####
m <- lm(Mass ~ Event * ts.sunrise + Detection, data = leye)

d <- expand.grid(Detection = unique(leye$Detection),   
                 Event = "Fall 2023",                    
                 ts.sunrise = mean(leye$ts.sunrise)) 

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
       y = "Predicted Lesser Yellowlegs Body Mass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  theme(legend.position = "none") +
  ylim(78, 100)
  # scale_x_discrete(labels = c("Detection" = "Neonicotinoid Detection", 
                              # "Non-detection" = "Neonicotinoid Non-detection"))

## LAB MEETING AG GRAPH SPRING 2025 ####
m <- lm(Mass ~ Event * ts.sunrise + PercentAg, data = leye)
m1 <- lm(Mass ~ Event * ts.sunrise + PercentAg, data = leye.cs)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg), 
                                 max(leye$PercentAg), 
                                 length = 1000),   
                 Event = "Fall 2023",                    
                 ts.sunrise = mean(leye$ts.sunrise)) 

d <- cbind(d, predict(m, newdata = d, interval = "confidence"))

ggplot(d, aes(x = (PercentAg * 100), y = fit)) +
  geom_line(size = 0.8, col = "black") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_light() +
  labs(x = "Surrounding Agricultural Intensity (%)", 
       y = "Predicted Lesser Yellowlegs Body Mass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 100, by = 20))
  
  

summary(m)



# BOXPLOT WITH NO MODELING
ggplot(leye, aes(x = Detection, y = Mass, fill = Detection)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(color = Detection),
               position = position_jitterdodge(jitter.width = 0.2),  
              size = 2, 
              alpha = 0.8) + 
  theme_light() +
  labs(x = "Neonicotinoid Detection", 
       y = "Lesser Yellowlegs Body Mass (g)") +  
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


# ---------------------------------------------------------------------------- #

# Model Assumptions #### REVIEW AND REDO

## Mass ~ Event * Time ####
m <- lm(Mass ~ Event * ts.sunrise, data = leye.cs)
par(mfrow = c(2, 2))
plot(m)

# Plot 1: Residuals vs Fitted Values (Check for linearity and homoscedasticity)
plot(fitted(m), resid(m), 
     main = "Residuals vs Fitted", 
     xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")  # Add horizontal line at 0

boxplot(residuals(m) ~ leye$Event, 
        main = "Residuals by Sampling Event",
        xlab = "Sampling Event",
        ylab = "Residuals")

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

# Conclusion: Few influential points

# Calculate Cook's Distance
cooks_dist <- cooks.distance(m)

# Set the threshold for influential points
threshold <- 4 / (length(resid(m)) - length(coef(m)))

# Identify the indices of influential points (Cook's Distance > threshold)
influential_points <- which(cooks_dist > threshold)

# Print the indices of influential points
influential_points

# If you want to see the corresponding data points in your dataset:
leye[influential_points, ]


# Plot residuals by categorical variable to assess homoscedasticity
m <- lm(Mass ~ Event * ts.sunrise, data = leye.cs)
boxplot(residuals(m) ~ leye$Event, 
        main = "Residuals by Sampling Event",
        xlab = "Sampling Event",
        ylab = "Residuals")

# Plot residuals by fitted values to assess homoscedasticity
residuals_m <- resid(m)
fitted_values_m <- fitted(m)

ggplot(data = data.frame(fitted_values_m, residuals_m), aes(x = fitted_values_m, y = residuals_m)) +
  geom_point(color = "blue", size = 2) +  # Plot points
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line at 0
  theme_classic() +
  labs(x = "Fitted Values", 
       y = "Residuals",
       title = "Residuals vs Fitted Values") +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

# Conclusion: No severe violation of homoscedasticity

# Plot QQPlot to assess normality of residuals
qqnorm(residuals_m)
qqline(residuals_m, col = "red")
shapiro.test(residuals_m)

# Conclusion: No significant deviation from normality

# ---------------------------------------------------------------------------- #

## Mass ~ Event * Time + Detection ####

# Plot residuals by categorical variable to assess homoscedasticity
m <- lm(Mass ~ Event * ts.sunrise + Detection, data = leye.cs)

par(mfrow = c(1, 1))

boxplot(residuals(m) ~ leye$Event, 
        main = "Residuals by Sampling Event",
        xlab = "Sampling Event",
        ylab = "Residuals")

boxplot(residuals(m) ~ leye$Detection, 
        main = "Residuals by Detection",
        xlab = "Detection",
        ylab = "Residuals")

# Conclusion: No severe violations of homoscedasticity

# Plot 1: Residuals vs Fitted Values (Check for linearity and homoscedasticity)
plot(fitted(m), resid(m), 
     main = "Residuals vs Fitted", 
     xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")  # Add horizontal line at 0

boxplot(residuals(m) ~ leye$Event, 
        main = "Residuals by Sampling Event",
        xlab = "Sampling Event",
        ylab = "Residuals")

boxplot(residuals(m) ~ leye$Detection, 
        main = "Residuals by Detection",
        xlab = "Detection",
        ylab = "Residuals")

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

# Conclusion: Four influential points but I'm not sure if they're worth omitting..

# Calculate Cook's Distance
cooks_dist <- cooks.distance(m)

# Set the threshold for influential points
threshold <- 4 / (length(resid(m)) - length(coef(m)))

# Identify the indices of influential points (Cook's Distance > threshold)
influential_points <- which(cooks_dist > threshold)

# Print the indices of influential points
influential_points

# If you want to see the corresponding data points in your dataset:
leye.cs[influential_points, ]

# Plot residuals by categorical variable to assess homoscedasticity
m <- lm(Mass ~ Event * ts.sunrise, data = leye.cs)
boxplot(residuals(m) ~ leye$Event, 
        main = "Residuals by Sampling Event",
        xlab = "Sampling Event",
        ylab = "Residuals")

# Plot residuals by fitted values to assess homoscedasticity
residuals_m <- resid(m)
fitted_values_m <- fitted(m)

ggplot(data = data.frame(fitted_values_m, residuals_m), aes(x = fitted_values_m, y = residuals_m)) +
  geom_point(color = "blue", size = 2) +  # Plot points
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line at 0
  theme_classic() +
  labs(x = "Fitted Values", 
       y = "Residuals",
       title = "Residuals vs Fitted Values") +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

# Conclusion: No severe violation of homoscedasticity

# Plot QQPlot to assess normality of residuals
qqnorm(residuals_m)
qqline(residuals_m, col = "red")
shapiro.test(residuals_m)

# Conclusion: No significant deviation from normality

# ---------------------------------------------------------------------------- #

# Omit influential points and re-run analysis: Mass ~ Event x Time ####

m <- lm(Mass ~ Event * ts.sunrise, data = leye)

# Calculate Cook's Distance
cooks_dist <- cooks.distance(m)

# Set the threshold for influential points
threshold <- 4 / (length(resid(m)) - length(coef(m)))

# Identify the indices of influential points (Cook's Distance > threshold)
influential_points <- which(cooks_dist > threshold)

# Print the indices of influential points
influential_points

# If you want to see the corresponding data points in your dataset:
leye[influential_points, ]

# Remove the influential points from the dataset
leye_without_influential <- leye[-influential_points, ]

# Rerun the model without influential points
m_without_influential <- lm(Mass ~ Event * ts.sunrise, data = leye_without_influential)

# View the summary of the new model
cbind(summary(m_without_influential)$coefficients, confint(m_without_influential))
cbind(summary(m)$coefficients, confint(m))

# Interaction no longer significant...graph it and then check assumptions again!

d_list <- list()

# Loop through each event and create a prediction grid for each event
for (event in unique(leye_without_influential$Event)) {
  event_data <- leye_without_influential[leye_without_influential$Event == event, ]
  d_event <- expand.grid(Event = event,
                         ts.sunrise = seq(min(event_data$ts.sunrise), 
                                          max(event_data$ts.sunrise), 
                                          length = 1000))
  d_event <- cbind(d_event, predict(m_without_influential, newdata = d_event, interval = "confidence"))
  d_list[[event]] <- d_event
}

d <- expand.grid(Event = event,
                 ts.sunrise = seq(min(event_data$ts.sunrise), 
                                  max(event_data$ts.sunrise), 
                                  length = 1000))

d <- cbind(d_event, predict(m_without_influential, newdata = d_event, interval = "confidence"))

d <- do.call(rbind, d_list)

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
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1)) +
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
  geom_point(data = leye_without_influential, aes(x = ts.sunrise, y = Mass)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 0, y = max(d$fit), 
           label = "Sunrise", angle = 90, 
           vjust = -0.5, hjust = 0.5,
           size = 5)

# Plot model without influential points ####
par(mfrow = c(2, 2))
plot(m_without_influential)

# 

