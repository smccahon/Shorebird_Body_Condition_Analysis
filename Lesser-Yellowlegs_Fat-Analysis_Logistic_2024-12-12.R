#----------------------------------#
#  Lesser Yellowlegs Fat Analysis  #
#          Created 12/11/2024      #          
#         Modified 01/20/2025      #
#----------------------------------#

# MODIFED TO USE FAT AS ORDINAL (LOGSTIC REGRESSION, 2 GROUPS: Low and High)

# load packages
library(ggplot2)
library(AICcmodavg)
library(tidyverse)
library(pROC)
library(splines)

# ---------------------------------------------------------------------------- #

# Data Processing & Manipulation ####

# Read data
setwd("processed_data")
birds <- read.csv("Shorebird_Data_Cleaned_2025-01-20.csv")

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

leye$Fat <- factor(leye$Fat,
                   levels = c("0", "1", "2", "3", "4", "5"))

# Convert Fat to two levels (low and high)
leye <- leye %>%
  mutate(Fat.G = case_when(
    Fat == 0 | Fat == 1 | Fat == 2 ~ "Low",       
    Fat == 3 | Fat == 4 | Fat == 5 ~ "High"))

leye$Fat.G <- factor(leye$Fat.G,
                     levels = c("Low", "High"))

leye$Fat.G_binary <- ifelse(leye$Fat.G == "Low", 0, 1)

# Subset data to only include fall 2021 and fall 2023
leye_fall <- leye[leye$Event %in% c("Fall 2021", "Fall 2023"), ]

# Standardize continuous variables
leye.cs <- leye_fall %>%
  mutate(across(where(is.numeric), scale))

# ---------------------------------------------------------------------------- #

# Data Visualization

my_theme <- theme_classic() +
  theme(text = element_text(size = 16),
        axis.title.y = element_text(margin = margin(r = 13)),
        axis.title.x = element_text(margin = margin(t = 13)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")

# Fat.G ~ Detection
ggplot(leye_fall, aes(x = Fat.G, fill = Detection)) +
  geom_bar(position = "dodge") +  
  labs(x = "Lesser Yellowlegs Fat Score", 
       y = "Frequency") +
  theme_classic() +
  theme(text = element_text(size = 16),
        axis.title.y = element_text(margin = margin(r = 13)),
        axis.title.x = element_text(margin = margin(t = 13)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

m <- glm(Fat.G ~ Detection, data = leye, family = "binomial")
summary(m)
confint(m)

# Fat.G ~ Agricultural Intensity
ggplot(leye_fall, aes(x = Fat.G, y = PercentAg)) +
  geom_boxplot() +  
  labs(y = "Surrounding Agricultural Intensity", 
       x = "Lesser Yellowlegs Fat Score") +
  my_theme +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = scales::percent)

# Fat.G ~ log(Agricultural Intensity)
ggplot(leye_fall, aes(x = Fat.G, y = LogAg)) +
  geom_boxplot() +  
  labs(y = "Surrounding Agricultural Intensity", 
       x = "Lesser Yellowlegs Fat Score") +
   my_theme


# Fat.G ~ Season
ggplot(leye_fall, aes(x = Fat.G, fill = Event)) +
  geom_bar(position = "dodge") +  
  labs(x = "Lesser Yellowlegs Fat Score", 
       y = "Frequency") +
  theme_classic() +
  theme(text = element_text(size = 16),
        axis.title.y = element_text(margin = margin(r = 13)),
        axis.title.x = element_text(margin = margin(t = 13)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

# Fat ~ Time
ggplot(leye_fall,
       aes(x = Fat.G, y = ts.sunrise)) +
  geom_boxplot() +  
  labs(x = "Lesser Yellowlegs Fat Score", 
       y = "Capture Time Since Sunrise (min)") +
  theme_classic() +
  theme(text = element_text(size = 16),
        axis.title.y = element_text(margin = margin(r = 13)),
        axis.title.x = element_text(margin = margin(t = 13)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

  

# ---------------------------------------------------------------------------- #

## Interaction Analysis: Event and Capture Time ####

m1 <- glm(Fat.G ~ Event * ts.sunrise, data = leye.cs, family = "binomial")
m2 <- glm(Fat.G ~ Event + ts.sunrise, data = leye.cs, family = "binomial")

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with interaction performs significantly better.

# ---------------------------------------------------------------------------- #

## Univariate Analysis: Time

m1 <- glm(Fat.G ~ ts.sunrise, data = leye.cs, family = "binomial")
m2 <- glm(Fat.G ~ ts.sunrise + I(ts.sunrise^2), data = leye.cs, family = "binomial")
m3 <- glm(Fat.G ~ ts.sunrise + I(ts.sunrise^2) + I(ts.sunrise^3), data = leye.cs, family = "binomial")
m4 <- glm(Fat.G ~ ns(ts.sunrise, df = 3) + Sex + Event * ts.sunrise, family = "binomial", data = leye.cs)

model_names <- paste0("m", 1:4)

models <- mget(model_names)

aictab(models, modnames = model_names)

# model with the splines for time performs best

## Splines and interaction? ####
m1 <- glm(Fat.G ~ Event * ts.sunrise, data = leye.cs, family = "binomial")
m2 <- glm(Fat.G ~ Event + ts.sunrise, data = leye.cs, family = "binomial")
m3 <- glm(Fat.G ~ Event + ns(ts.sunrise, df = 3), data = leye.cs, family = "binomial")

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model with interaction performs significantly better
# ---------------------------------------------------------------------------- #

# INTERACTIONS/TRANSFORMATIONS INCLUDED: Event * ts.sunrise

# ---------------------------------------------------------------------------- #

## Single Models ####
## Single Models ####

m.global <- glm(Fat.G ~ Sex + Event * ts.sunrise + DaysIntoSeason_S + Detection + PercentAg + Age, 
               data = leye.cs, family = "binomial")

m.null <- glm(Fat.G ~ 1, data = leye.cs, family = "binomial")

m1 <- glm(Fat.G ~ Sex, data = leye.cs, family = "binomial")
m2 <- glm(Fat.G ~ Event, data = leye.cs, family = "binomial")
m3 <- glm(Fat.G ~ DaysIntoSeason_S, data = leye.cs, family = "binomial")
m4 <- glm(Fat.G ~ ts.sunrise, data = leye.cs, family = "binomial")
m5 <- glm(Fat.G ~ PercentAg, data = leye.cs, family = "binomial")
m6 <- glm(Fat.G ~ Detection, data = leye.cs, family = "binomial")
m7 <- glm(Fat.G ~ Age, data = leye.cs, family = "binomial")

## Additive Models ####

### Two additive combinations ####
m8 <- glm(Fat.G ~ Sex + Event, data = leye.cs, family = "binomial")
m9 <- glm(Fat.G ~ Sex + DaysIntoSeason_S, data = leye.cs, family = "binomial")
m10 <- glm(Fat.G ~ Sex + ts.sunrise, data = leye.cs, family = "binomial")
m11 <- glm(Fat.G ~ Sex + PercentAg, data = leye.cs, family = "binomial")
m12 <- glm(Fat.G ~ Sex + Detection, data = leye.cs, family = "binomial")
m13 <- glm(Fat.G ~ Event + DaysIntoSeason_S, data = leye.cs, family = "binomial")
m14 <- glm(Fat.G ~ Event * ts.sunrise, data = leye.cs, family = "binomial")
m15 <- glm(Fat.G ~ Event + PercentAg, data = leye.cs, family = "binomial")
m16 <- glm(Fat.G ~ DaysIntoSeason_S + ts.sunrise, data = leye.cs, family = "binomial")
m17 <- glm(Fat.G ~ DaysIntoSeason_S + PercentAg, data = leye.cs, family = "binomial")
m18 <- glm(Fat.G ~ DaysIntoSeason_S + Detection, data = leye.cs, family = "binomial")
m19 <- glm(Fat.G ~ ts.sunrise + PercentAg, data = leye.cs, family = "binomial")
m20 <- glm(Fat.G ~ ts.sunrise + Detection, data = leye.cs, family = "binomial")
m21 <- glm(Fat.G ~ PercentAg + Detection, data = leye.cs, family = "binomial")
m22 <- glm(Fat.G ~ Age + Event, data = leye.cs, family = "binomial")
m23 <- glm(Fat.G ~ Age + DaysIntoSeason_S, data = leye.cs, family = "binomial")
m24 <- glm(Fat.G ~ Age + ts.sunrise, data = leye.cs, family = "binomial")
m25 <- glm(Fat.G ~ Age + PercentAg, data = leye.cs, family = "binomial")
m26 <- glm(Fat.G ~ Age + Detection, data = leye.cs, family = "binomial")

### Three additive combinations ####
m27 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S, data = leye.cs, family = "binomial")
m28 <- glm(Fat.G ~ Sex + Event * ts.sunrise, data = leye.cs, family = "binomial")
m29 <- glm(Fat.G ~ Sex + Event + PercentAg, data = leye.cs, family = "binomial")
m30 <- glm(Fat.G ~ Sex + Event + Detection, data = leye.cs, family = "binomial")
m31 <- glm(Fat.G ~ Sex + DaysIntoSeason_S + ts.sunrise, data = leye.cs, family = "binomial")
m32 <- glm(Fat.G ~ Sex + DaysIntoSeason_S + PercentAg, data = leye.cs, family = "binomial")
m33 <- glm(Fat.G ~ Sex + DaysIntoSeason_S + Detection, data = leye.cs, family = "binomial")
m34 <- glm(Fat.G ~ Sex + ts.sunrise + PercentAg, data = leye.cs, family = "binomial")
m35 <- glm(Fat.G ~ Sex + ts.sunrise + Detection, data = leye.cs, family = "binomial")
m36 <- glm(Fat.G ~ Sex + PercentAg + Detection, data = leye.cs, family = "binomial")
m37 <- glm(Fat.G ~ Event + DaysIntoSeason_S + ts.sunrise, data = leye.cs, family = "binomial")
m38 <- glm(Fat.G ~ Event + DaysIntoSeason_S + PercentAg, data = leye.cs, family = "binomial")
m39 <- glm(Fat.G ~ Event + DaysIntoSeason_S + Detection, data = leye.cs, family = "binomial")
m40 <- glm(Fat.G ~ Event * ts.sunrise + PercentAg, data = leye.cs, family = "binomial")
m41 <- glm(Fat.G ~ Event * ts.sunrise + Detection, data = leye.cs, family = "binomial")
m42 <- glm(Fat.G ~ Event + PercentAg + Detection, data = leye.cs, family = "binomial")
m43 <- glm(Fat.G ~ DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs, family = "binomial")
m44 <- glm(Fat.G ~ DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs, family = "binomial")
m45 <- glm(Fat.G ~ DaysIntoSeason_S + PercentAg + Detection, data = leye.cs, family = "binomial")
m46 <- glm(Fat.G ~ ts.sunrise + PercentAg + Detection, data = leye.cs, family = "binomial")
m47 <- glm(Fat.G ~ Age + Event + DaysIntoSeason_S, data = leye.cs, family = "binomial")
m48 <- glm(Fat.G ~ Age + Event * ts.sunrise, data = leye.cs, family = "binomial")
m49 <- glm(Fat.G ~ Age + Event + PercentAg, data = leye.cs, family = "binomial")
m50 <- glm(Fat.G ~ Age + Event + Detection, data = leye.cs, family = "binomial")
m51 <- glm(Fat.G ~ Age + DaysIntoSeason_S + ts.sunrise, data = leye.cs, family = "binomial")
m52 <- glm(Fat.G ~ Age + DaysIntoSeason_S + PercentAg, data = leye.cs, family = "binomial")
m53 <- glm(Fat.G ~ Age + DaysIntoSeason_S + Detection, data = leye.cs, family = "binomial")
m54 <- glm(Fat.G ~ Age + ts.sunrise + PercentAg, data = leye.cs, family = "binomial")
m55 <- glm(Fat.G ~ Age + ts.sunrise + Detection, data = leye.cs, family = "binomial")
m56 <- glm(Fat.G ~ Age + PercentAg + Detection, data = leye.cs, family = "binomial")

### Four additive combinations ####
m57 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + Event * ts.sunrise, data = leye.cs, family = "binomial")
m58 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S + PercentAg, data = leye.cs, family = "binomial")
m59 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S + Detection, data = leye.cs, family = "binomial")
m60 <- glm(Fat.G ~ Sex + Event * ts.sunrise + PercentAg, data = leye.cs, family = "binomial")
m61 <- glm(Fat.G ~ Sex + Event * ts.sunrise + Detection, data = leye.cs, family = "binomial")
m62 <- glm(Fat.G ~ Sex + Event + PercentAg + Detection, data = leye.cs, family = "binomial")
m63 <- glm(Fat.G ~ Sex + DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs, family = "binomial")
m64 <- glm(Fat.G ~ Sex + DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs, family = "binomial")
m65 <- glm(Fat.G ~ Sex + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs, family = "binomial")
m66 <- glm(Fat.G ~ Sex + ts.sunrise + PercentAg + Detection, data = leye.cs, family = "binomial")
m67 <- glm(Fat.G ~ Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Event * ts.sunrise, data = leye.cs, family = "binomial")
m68 <- glm(Fat.G ~ Event + DaysIntoSeason_S + ts.sunrise + Detection + Event * ts.sunrise, data = leye.cs, family = "binomial")
m69 <- glm(Fat.G ~ Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs, family = "binomial")
m70 <- glm(Fat.G ~ Event * ts.sunrise + PercentAg + Detection, data = leye.cs, family = "binomial")
m71 <- glm(Fat.G ~ DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs, family = "binomial")
m72 <- glm(Fat.G ~ Age + Event + DaysIntoSeason_S + ts.sunrise + Event * ts.sunrise, data = leye.cs, family = "binomial")
m73 <- glm(Fat.G ~ Age + Event + DaysIntoSeason_S + PercentAg, data = leye.cs, family = "binomial")
m74 <- glm(Fat.G ~ Age + Event + DaysIntoSeason_S + Detection, data = leye.cs, family = "binomial")
m75 <- glm(Fat.G ~ Age + Event * ts.sunrise + PercentAg, data = leye.cs, family = "binomial")
m76 <- glm(Fat.G ~ Age + Event * ts.sunrise + Detection, data = leye.cs, family = "binomial")
m77 <- glm(Fat.G ~ Age + Event + PercentAg + Detection, data = leye.cs, family = "binomial")
m78 <- glm(Fat.G ~ Age + DaysIntoSeason_S + ts.sunrise + PercentAg, data = leye.cs, family = "binomial")
m79 <- glm(Fat.G ~ Age + DaysIntoSeason_S + ts.sunrise + Detection, data = leye.cs, family = "binomial")
m80 <- glm(Fat.G ~ Age + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs, family = "binomial")
m81 <- glm(Fat.G ~ Age + ts.sunrise + PercentAg + Detection, data = leye.cs, family = "binomial")

### Five additive combinations ####

m82 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Event * ts.sunrise, data = leye.cs, family = "binomial")
m83 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S + ts.sunrise + Detection + Event * ts.sunrise, data = leye.cs, family = "binomial")
m84 <- glm(Fat.G ~ Sex + Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs, family = "binomial")
m85 <- glm(Fat.G ~ Sex + Event * ts.sunrise + PercentAg + Detection, data = leye.cs, family = "binomial")
m86 <- glm(Fat.G ~ Sex + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs, family = "binomial")
m87 <- glm(Fat.G ~ Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection + Event * ts.sunrise, data = leye.cs, family = "binomial")
m88 <- glm(Fat.G ~ Age + Event + DaysIntoSeason_S + ts.sunrise + PercentAg + Event * ts.sunrise, data = leye.cs, family = "binomial")
m89 <- glm(Fat.G ~ Age + Event + DaysIntoSeason_S + ts.sunrise + Detection + Event * ts.sunrise, data = leye.cs, family = "binomial")
m90 <- glm(Fat.G ~ Age + Event + DaysIntoSeason_S + PercentAg + Detection, data = leye.cs, family = "binomial")
m91 <- glm(Fat.G ~ Age + Event * ts.sunrise + PercentAg + Detection, data = leye.cs, family = "binomial")
m92 <- glm(Fat.G ~ Age + DaysIntoSeason_S + ts.sunrise + PercentAg + Detection, data = leye.cs, family = "binomial")


# ---------------------------------------------------------------------------- #

### AIC Model Selection ####

model_names <- paste0("m", 1:92)

models <- mget(model_names)

models$m.null <- m.null
models$m.global <- m.global

model_names <- c(model_names, "m.null", "m.global")

aictab(models, modnames = model_names)

#NEW RESULTS 2025-01-20 ####
confint(m13) # event*time signficant
confint(m23) # event*time significant
confint(m3) # date significant
confint(m2) # event significant
confint(m8) # date significant
confint(m32) # event*time signficant
confint(m7) # nothing significant
confint(m42) # event*time significant
confint(m16) # nothing significant

# NEW RESULTS AFTER INCLUDING AGE 2025-01-22 ####
confint(m23) # juveniles have a higher probability of having less fat; date and fat are positively correlated
summary(m23)

ggplot(leye, aes(x = Age, fill = Fat.G)) + geom_bar(position = "dodge")

confint(m72) # juveniles have a higher probability of having less fat; event * time significant
summary(m72)

confint(m47) # juveniles have a higher probability of having less fat;
summary(m47)

confint(m53) # juveniles have a higher probability of having less fat; date and fat are positively correlated
summary(m53)

confint(m14) # event *time significant
confint(m28) # event * time significant
confint(m48) # event * time significant
confint(m51) # juveniles have a higher probability of having less fat; date and fat are positively correlated
confint(m3) # date and fat are positively correlated
confint(m22) # event significant
confint(m2) 

 # ---------------------------------------------------------------------------- #

# Top Model Summaries ####

summary(m45)
confint(m45)
# Event * time is significant; birds in high ag have significantly higher fat scores

summary(m35)
confint(m35)
# Event * time is significant; birds in high ag have significantly higher fat scores

summary(m14)
confint(m14)
# Birds in fall 2023 had significantly higher fat scores compared to fall 2021

summary(m24)
confint(m24)
# Birds in high ag have significantly higher fat scores

summary(m13)
confint(m13)
# Event * time is significant

summary(m23)
confint(m23)
# Event * time is significant

## Goodness of Fit ####
predicted_probs <- predict(m45, type = "response")
roc_curve <- roc(leye.cs$Fat.G, predicted_probs)
plot(roc_curve)
auc(roc_curve)

# ---------------------------------------------------------------------------- #

# Model assumptions ####

## Linear relationship between explanatory variables and logit of response variable ####
# satistifed: agricultural intensity has a linear relationship with logit of Fat
# not satisfied: capture time has a non-linear relationship with logit of Fat 
# (splines of time improve fit but model with interaction and no splines is still a better fit)

# Fit logistic regression model
m45 <- glm(Fat.G ~ PercentAg + Sex + Event * ts.sunrise, family = "binomial", data = leye_fall)
summary(m45)
confint(m45)

# Box-Tidwell transformation for PercentAg
leye_fall$log_PercentAg <- log(leye_fall$PercentAg + 1)  # Adding a small constant to avoid log(0)

# Fit model with interaction terms (PercentAg and its log transformation)
m_bt <- glm(Fat.G ~ PercentAg * log_PercentAg + Sex + Event * ts.sunrise, family = "binomial", data = leye_fall)

# View the summary of the Box-Tidwell model
summary(m_bt)

# ---------------------------------------------------------------------------- #

# Plot Top Models

m35 <- glm(Fat.G ~ Event * ts.sunrise + PercentAg, data = leye, family = "binomial")

d <- data.frame(
  Event = c("Fall 2023"),
  PercentAg = seq(min(leye$PercentAg), max(leye$PercentAg), length = 1000),
  ts.sunrise = mean(leye$ts.sunrise)
)

predictions <- predict(m35, newdata = d, type = "response", se.fit = TRUE) 

d$fit <- predictions$fit

d$lower_CI <- predictions$fit - 1.96 * predictions$se.fit  # Lower CI
d$upper_CI <- predictions$fit + 1.96 * predictions$se.fit  # Upper CI

ggplot(d, aes(x = PercentAg, y = fit)) +
  geom_line(aes(color = "tan2"), size = 1, show.legend = FALSE) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), 
              alpha = 0.25, color = NA, show.legend = FALSE,
              fill = "tan1") +  # Confidence intervals as ribbon
  theme_classic() +
  labs(x = "Surrounding Agricultural Intensity (%)", 
       y = "Probability of Lesser Yellowlegs Having High Fat") +
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 13)),
        axis.title.y = element_text(size = 16, margin = margin(r = 13)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none") +
  theme(legend.position = "none")  +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = scales::percent)


# LAB MEETING AG GRAPH ####
m <- glm(Fat.G ~ PercentAg + Event * ts.sunrise + Sex, data = leye_fall, family = "binomial")

d <- expand.grid(PercentAg = seq(min(leye_fall$PercentAg),
                                 max(leye_fall$PercentAg),
                                 length = 1000),
                 Event = c("Fall 2023"),
                 ts.sunrise = mean(leye_fall$ts.sunrise),
                 Sex = c("Female"))

predictions <- predict(m, newdata = d, type = "response", se.fit = TRUE) 

d$fit <- predictions$fit

d$lwr <- predictions$fit - 1.96 * predictions$se.fit  # Lower CI
d$upr <- predictions$fit + 1.96 * predictions$se.fit  # Upper CI

ggplot(d, aes(x = PercentAg, y = fit)) +
  geom_line(size = 0.8, col = "black") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_light() +
  labs(x ="Surrounding Agricultural Intensity (%)", 
       y = "P(High Fat)") +
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

# Get predicted probabilities
predicted_probs <- predict(m, type = "response")

# Create a ROC curve
roc_curve <- roc(leye_fall$Fat.G, predicted_probs)

# AUC value
auc_value <- auc(roc_curve)
auc_value



# LAB MEETING DETECTION GRAPH ####

m <- glm(Fat.G ~ Detection + PercentAg + Event, data = leye_fall, family = "binomial")

d <- expand.grid(PercentAg = mean(leye_fall$PercentAg),
                 Event = c("Fall 2023"),
                 Detection = c("Detection", "Non-detection"))

predictions <- predict(m, newdata = d, type = "response", se.fit = TRUE) 

d$fit <- predictions$fit

d$lwr <- predictions$fit - 1.96 * predictions$se.fit  # Lower CI
d$upr <- predictions$fit + 1.96 * predictions$se.fit  # Upper CI

ggplot(d, aes(x = Detection, y = fit)) +
  geom_point(size = 5, col = "black") +  # Points showing predicted values
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1,
                col = "black",
                size = 1) +  # Add confidence intervals
  theme_light() +
  labs(x = NULL, 
       y = "P(High Fat)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  theme(legend.position = "none")

summary(m)

# Get predicted probabilities
predicted_probs <- predict(m, type = "response")

# Create a ROC curve
roc_curve <- roc(leye_fall$Fat.G, predicted_probs)

# AUC value
auc_value <- auc(roc_curve)
auc_value
