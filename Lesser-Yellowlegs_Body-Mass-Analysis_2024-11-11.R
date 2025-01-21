#---------------------------------------#
# Lesser Yellowlegs Body Mass Analysis  #
#         Created 11/11/2024            #          
#        Modified 11/20/2024            #
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

## Interaction Analysis: Event and Capture Time ####

m1 <- lm(Mass ~ Event * ts.sunrise, data = leye.cs)
m2 <- lm(Mass ~ Event + I(ts.sunrise^2), data = leye.cs)
m3 <- lm(Mass ~ Event + ts.sunrise, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with interaction by far performs the best.

# ---------------------------------------------------------------------------- #

## Interaction Analysis: LogNeonic & Sampling Event ####

m1 <- lm(Mass ~ Event * LogNeonic, data = leye.cs)
m2 <- lm(Mass ~ Event + LogNeonic, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model without interaction performs significantly better

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Detection & Sampling Event ####

m1 <- lm(Mass ~ Event * Detection, data = leye.cs)
m2 <- lm(Mass ~ Event + Detection, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model without interaction performs significantly better

# ---------------------------------------------------------------------------- #

## Single Models ####

m.global <- lm(Mass ~ Sex + Event*ts.sunrise + MigDate + Detection + PercentAg, 
                 data = leye.cs)


m.null <- lm(Mass ~ 1, data = leye.cs)

m1 <- lm(Mass ~ Sex, data = leye.cs)
m2 <- lm(Mass ~ Event, data = leye.cs)
m3 <- lm(Mass ~ MigDate, data = leye.cs)
m4 <- lm(Mass ~ ts.sunrise, data = leye.cs)
m5 <- lm(Mass ~ PercentAg, data = leye.cs)
m6 <- lm(Mass ~ Detection, data = leye.cs)
m7 <- lm(Mass ~ LogNeonic, data = leye.cs)

## Additive Models ####

### Two additive combinations ####
m8 <- lm(Mass ~ Sex + Event, data = leye.cs)
m9 <- lm(Mass ~ Sex + MigDate, data = leye.cs)
m10 <- lm(Mass ~ Sex + ts.sunrise, data = leye.cs)
m11 <- lm(Mass ~ Sex + PercentAg, data = leye.cs)
m12 <- lm(Mass ~ Sex + Detection, data = leye.cs)
m13 <- lm(Mass ~ Sex + LogNeonic, data = leye.cs)
m14 <- lm(Mass ~ Event + MigDate, data = leye.cs)
m15 <- lm(Mass ~ Event * ts.sunrise, data = leye.cs)
m16 <- lm(Mass ~ Event + PercentAg, data = leye.cs)
m17 <- lm(Mass ~ Event + Detection, data = leye.cs)
m18 <- lm(Mass ~ Event + LogNeonic, data = leye.cs)
m19 <- lm(Mass ~ MigDate + ts.sunrise, data = leye.cs)
m20 <- lm(Mass ~ MigDate + PercentAg, data = leye.cs)
m21 <- lm(Mass ~ MigDate + Detection, data = leye.cs)
m22 <- lm(Mass ~ MigDate + LogNeonic, data = leye.cs)
m23 <- lm(Mass ~ ts.sunrise + PercentAg, data = leye.cs)
m24 <- lm(Mass ~ ts.sunrise + Detection, data = leye.cs)
m25 <- lm(Mass ~ ts.sunrise + LogNeonic, data = leye.cs)
m26 <- lm(Mass ~ PercentAg + Detection, data = leye.cs)
m27 <- lm(Mass ~ PercentAg + LogNeonic, data = leye.cs)

### Three additive combinations ####
m28 <- lm(Mass ~ Sex + Event + MigDate, data = leye.cs)
m29 <- lm(Mass ~ Sex + Event * ts.sunrise, data = leye.cs)
m30 <- lm(Mass ~ Sex + Event + PercentAg, data = leye.cs)
m31 <- lm(Mass ~ Sex + Event + Detection, data = leye.cs)
m32 <- lm(Mass ~ Sex + MigDate + ts.sunrise, data = leye.cs)
m33 <- lm(Mass ~ Sex + MigDate + PercentAg, data = leye.cs)
m34 <- lm(Mass ~ Sex + MigDate + Detection, data = leye.cs)
m35 <- lm(Mass ~ Sex + ts.sunrise + PercentAg, data = leye.cs)
m36 <- lm(Mass ~ Sex + ts.sunrise + Detection, data = leye.cs)
m37 <- lm(Mass ~ Sex + PercentAg + Detection, data = leye.cs)
m38 <- lm(Mass ~ Event + MigDate + ts.sunrise + Event * ts.sunrise, data = leye.cs)
m39 <- lm(Mass ~ Event + MigDate + PercentAg, data = leye.cs)
m40 <- lm(Mass ~ Event + MigDate + Detection, data = leye.cs)
m41 <- lm(Mass ~ Event * ts.sunrise + PercentAg, data = leye.cs)
m42 <- lm(Mass ~ Event * ts.sunrise + Detection, data = leye.cs)
m43 <- lm(Mass ~ Event + PercentAg + Detection, data = leye.cs)
m44 <- lm(Mass ~ MigDate + ts.sunrise + PercentAg, data = leye.cs)
m45 <- lm(Mass ~ MigDate + ts.sunrise + Detection, data = leye.cs)
m46 <- lm(Mass ~ MigDate + PercentAg + Detection, data = leye.cs)
m47 <- lm(Mass ~ ts.sunrise + PercentAg + Detection, data = leye.cs)
m48 <- lm(Mass ~ MigDate + ts.sunrise + LogNeonic, data = leye.cs)
m49 <- lm(Mass ~ MigDate + PercentAg + LogNeonic, data = leye.cs)
m50 <- lm(Mass ~ ts.sunrise + PercentAg + LogNeonic, data = leye.cs)
m51 <- lm(Mass ~ Event * ts.sunrise + LogNeonic, data = leye.cs)
m52 <- lm(Mass ~ Event + PercentAg + LogNeonic, data = leye.cs)
m53 <- lm(Mass ~ Event + MigDate + LogNeonic, data = leye.cs)
m54 <- lm(Mass ~ Sex + ts.sunrise + LogNeonic, data = leye.cs)
m55 <- lm(Mass ~ Sex + PercentAg + LogNeonic, data = leye.cs)
m56 <- lm(Mass ~ Sex + MigDate + LogNeonic, data = leye.cs)
m57 <- lm(Mass ~ Sex + Event + LogNeonic, data = leye.cs)

### Four additive combinations ####
m58 <- lm(Mass ~ Sex + Event + MigDate + ts.sunrise + Event*ts.sunrise, data = leye.cs)
m59 <- lm(Mass ~ Sex + Event + MigDate + PercentAg, data = leye.cs)
m60 <- lm(Mass ~ Sex + Event + MigDate + Detection, data = leye.cs)
m61 <- lm(Mass ~ Sex + Event * ts.sunrise + PercentAg, data = leye.cs)
m62 <- lm(Mass ~ Sex + Event * ts.sunrise + Detection, data = leye.cs)
m63 <- lm(Mass ~ Sex + Event + PercentAg + Detection, data = leye.cs)
m64 <- lm(Mass ~ Sex + MigDate + ts.sunrise + PercentAg, data = leye.cs)
m65 <- lm(Mass ~ Sex + MigDate + ts.sunrise + Detection, data = leye.cs)
m66 <- lm(Mass ~ Sex + MigDate + PercentAg + Detection, data = leye.cs)
m67 <- lm(Mass ~ Sex + ts.sunrise + PercentAg + Detection, data = leye.cs)
m68 <- lm(Mass ~ Event + MigDate + ts.sunrise + PercentAg + Event*ts.sunrise, data = leye.cs)
m69 <- lm(Mass ~ Event + MigDate + ts.sunrise + Detection + Event*ts.sunrise, data = leye.cs)
m70 <- lm(Mass ~ Event + MigDate + PercentAg + Detection, data = leye.cs)
m71 <- lm(Mass ~ Event * ts.sunrise + PercentAg + Detection, data = leye.cs)
m72 <- lm(Mass ~ MigDate + ts.sunrise + PercentAg + Detection, data = leye.cs)
m73 <- lm(Mass ~ Sex + Event + MigDate + LogNeonic, data = leye.cs)
m74 <- lm(Mass ~ Sex + Event * ts.sunrise + LogNeonic, data = leye.cs)
m75 <- lm(Mass ~ Sex + Event + PercentAg + LogNeonic, data = leye.cs)
m76 <- lm(Mass ~ Sex + MigDate + ts.sunrise + LogNeonic, data = leye.cs)
m77 <- lm(Mass ~ Sex + MigDate + PercentAg + LogNeonic, data = leye.cs)
m78 <- lm(Mass ~ Sex + ts.sunrise + PercentAg + LogNeonic, data = leye.cs)
m79 <- lm(Mass ~ Event + MigDate + ts.sunrise + LogNeonic + Event*ts.sunrise, data = leye.cs)
m80 <- lm(Mass ~ Event + MigDate + PercentAg + LogNeonic, data = leye.cs)
m81 <- lm(Mass ~ Event * ts.sunrise + PercentAg + LogNeonic, data = leye.cs)
m82 <- lm(Mass ~ MigDate + ts.sunrise + PercentAg + LogNeonic, data = leye.cs)
 
### Five additive combinations ####

m83 <- lm(Mass ~ Sex + Event + MigDate + ts.sunrise + PercentAg + Event * ts.sunrise, data = leye.cs)
m84 <- lm(Mass ~ Sex + Event + MigDate + ts.sunrise + Detection + Event * ts.sunrise, data = leye.cs)
m85 <- lm(Mass ~ Sex + Event + MigDate + PercentAg + Detection, data = leye.cs)
m86 <- lm(Mass ~ Sex + Event * ts.sunrise + PercentAg + Detection, data = leye.cs)
m87 <- lm(Mass ~ Sex + MigDate + ts.sunrise + PercentAg + Detection, data = leye.cs)
m88 <- lm(Mass ~ Event + MigDate + ts.sunrise + PercentAg + Detection + Event*ts.sunrise, data = leye.cs)
m89 <- lm(Mass ~ Sex + Event + MigDate + ts.sunrise + LogNeonic + Event * ts.sunrise, data = leye.cs)
m90 <- lm(Mass ~ Sex + Event + MigDate + PercentAg + LogNeonic, data = leye.cs)
m91 <- lm(Mass ~ Sex + Event * ts.sunrise + PercentAg + LogNeonic, data = leye.cs)
m92 <- lm(Mass ~ Sex + MigDate + ts.sunrise + PercentAg + LogNeonic, data = leye.cs)
m93 <- lm(Mass ~ Event + MigDate + ts.sunrise + PercentAg + LogNeonic + Event*ts.sunrise, data = leye.cs)

# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:93)

models <- mget(model_names)

models$m.null <- m.null
models$m.global <- m.global

model_names <- c(model_names, "m.null", "m.global")

aictab(models, modnames = model_names)

# NEW RESULTS 2025-01-20 ####


# ---------------------------------------------------------------------------- #

### Top Model Summaries ####

summary(m42)$coefficients
confint(m42)

cbind(summary(m51)$coefficients, confint(m51))

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


## Mass ~ Event * Time ####
m15 <- lm(Mass ~ Event * ts.sunrise, data = leye)

d <- expand.grid(Event = c("Fall 2021", "Spring 2022", "Fall 2023"),
                 ts.sunrise = seq(min(leye$ts.sunrise), 
                                  max(leye$ts.sunrise), 
                                  length = 1000))

d <- cbind(d, predict(m15, newdata = d, interval = "confidence"))

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

## Mass ~ Event * Time + Detection ####

m42 <- lm(Mass ~ Event * ts.sunrise + Detection, data = leye)

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

## Mass ~ Event * Time + LogNeonic (m51) ####

# Plot residuals by categorical variable to assess homoscedasticity
m <- lm(Mass ~ Event * ts.sunrise + LogNeonic, data = leye.cs)
residuals_m <- resid(m)
fitted_values_m <- fitted(m)

boxplot(residuals(m) ~ leye$Event, 
        main = "Residuals by Sampling Event",
        xlab = "Sampling Event",
        ylab = "Residuals")

# Conclusion: No violation of homoscedasticity

# Normality of residuals
qqnorm(residuals_m)
qqline(residuals_m, col = "red")
shapiro.test(residuals_m)

# Conclusion: Some significant deviation from normality.


## Mass ~ Event + Detection (m5)

# Plot residuals by categorical variable to assess homoscedasticity
m <- lm(Mass ~ Event + Detection, data = leye)
residuals_m <- resid(m)
fitted_values_m <- fitted(m)

boxplot(residuals(m) ~ leye$Event, 
        main = "Residuals by Sampling Event",
        xlab = "Sampling Event",
        ylab = "Residuals")

boxplot(residuals(m) ~ leye$Detection, 
        main = "Residuals by Sampling Event",
        xlab = "Detection",
        ylab = "Residuals")

# Conclusion: No violation of homoscedasticity

# Normality of residuals
qqnorm(residuals_m)
qqline(residuals_m, col = "red")
shapiro.test(residuals_m)

# Conclusion: Some significant deviation of normality


## Mass ~ Event + LogNeonic (m18)
m <- lm(Mass ~ Event + LogNeonic, data = leye)
residuals_m <- resid(m)
fitted_values_m <- fitted(m)

par(mfrow = c(2,2))
plot(m)

boxplot(residuals(m) ~ leye$Event, 
        main = "Residuals by Sampling Event",
        xlab = "Sampling Event",
        ylab = "Residuals")

# Greater residual variance in fall 2021 but not extreme

# Normality of residuals
qqnorm(residuals_m)
qqline(residuals_m, col = "red")
shapiro.test(residuals_m)

# Conclusion: Some significant deviation of normality

