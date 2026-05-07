################# Packages #################
library(tidyr)
library(dplyr)
library(car)
library(lme4)
library(glmmTMB)
library(tidyverse)
library(corrplot)
library(dotwhisker)
library(car)
library(gllvm)
library(MCMCglmm)
library(brms)
library(MasterBayes)
library(ggpedigree)
library(kinship2)
library(corrplot)
library(lubridate)
library(lme4)
library(DHARMa)
library(mgcv)
library(gridExtra)

install.packages("gridExtra")

########### Merging data ##################
# Dataframe 1
all_juveniles<- read.csv("All juveniles.csv")
View(all_juveniles)
# Change Anolis.ID to AnolisID
colnames(all_juveniles)[which(names(all_juveniles) == "Anolis.ID")] <- "AnolisID"
# Dataframe 2
xray_dat <- read.csv("xrays2.csv")
View(xray_dat)
# Create new dataframe which merges both data
new_dat <- merge(xray_dat, all_juveniles, by="AnolisID", all = TRUE)
View(new_dat)

################## Cleaning and filtering combined dataset ########################
# Remove data with NAs for Date column              
new_dat_clean <- new_dat %>% filter(!is.na(Date))
View(new_dat_clean) 
# Remove data without NAs for Notes.x
df_filtered <- new_dat_clean[is.na(new_dat_clean$Notes.x), ]
View(df_filtered)
# Select specific columns for analyses and filter for Andros (A) and San Salvador (S) populations
df_filtered_sub <- df_filtered %>% 
  dplyr::select(c("AnolisID","Order","Date","X.ray..","Sex.x","Mass","Meas",
                  "SVL","Head.Length","Head.Width","Pectoral","Pelvic","Humerus",
                  "Radius","Ulna","Femur","Tibia","Date.of.Hatch","Mother","Father",
                  "Age.at.Death","Egg.Weight","Date.Egg.Found","Mother.population",
                  "Father.Population","Population")) %>%
  dplyr::filter(Population == c("A","S"))
View(df_filtered_sub)
# Create object of morphological variables to select
morph_vars <- c("Head.Length","Head.Width","Pectoral","Pelvic","Humerus",
                "Radius","Ulna","Femur","Tibia")
morph_vars2 <- c("SVL","Head.Length","Head.Width","Pectoral","Pelvic","Humerus",
                "Radius","Ulna","Femur","Tibia")

# Take log of each morphological variables
df_filtered_sub_log <- df_filtered_sub %>% dplyr::mutate(across(c(SVL,Head.Length,Head.Width,Pectoral,
                                                                  Pelvic,Humerus,Radius,Ulna,Femur,Tibia), log))
View(df_filtered_sub_log)

################ Is growth rate different between male and female juveniles for each island population? ############
# Create new dataframe with growth rates for each individual
class(df_filtered_sub_log$Date)
class(df_filtered_sub_log$Date.of.Hatch)
# Convert Date and Date.of.Hatch to ISO format
df_filtered_sub_log$Date <- mdy(df_filtered_sub_log$Date)
df_filtered_sub_log$Date.of.Hatch <- mdy(df_filtered_sub_log$Date.of.Hatch)
# Order the data based on AnolisID and Date
df_filtered_sub_log_ordered <- df_filtered_sub_log %>% arrange(AnolisID, Date)
View(df_filtered_sub_log_ordered)
# Create "Age" column which is Date (date of measurement) - Date.of.Hatch
df_filtered_sub_log_new <- df_filtered_sub_log_ordered %>% mutate(Age = as.numeric(Date - Date.of.Hatch))
View(df_filtered_sub_log_new)
# Create mixed effects model comparing growth rates of different morphological variables between sex for each population

###### Filter for San Salvador population -------
df_model_san <- df_filtered_sub_log_new %>% dplyr::filter(Population == "S")
View(df_model_san)

# Plots the data showing growth trajectories of either males or females for each trait
# SVL and Age (body size seems to be different in males and females at final age)
ggplot(df_model_san, aes(Age, SVL, fill=Sex.x)) +
  geom_point() + 
  geom_smooth() +
  theme_bw()

# Linear mixed effects model - SVL
# Clean data for each variable
df_clean_san <- df_model_san %>% dplyr::filter(!is.na(SVL),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
View(df_clean_san)
# Model 1 - SVL
model1 <- lmer(SVL ~ Age * Sex.x + (1 | AnolisID),
               data = df_clean_san)
summary(model1) # Anolis ID (random effect) variance is 0
# boundary (singular) fit: see help('isSingular')
# I can rop random effect
lm1 <- lm(SVL ~ Age*Sex.x, data=df_clean_san)
shapiro.test(lm1$residuals) # p-value < 2.2e-16
qqnorm(lm1$residuals)
par(mfrow = c(2,2))
plot(lm1) # Residuals vs fitted are non-linear


# Clean data for each variable
df_clean_san <- df_model_san %>% dplyr::filter(!is.na(SVL),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
View(df_clean_san)
# Try using a gam (generalized additive model)
# Convert sex to a factor instead of character
df_clean_san$Sex.x <- as.factor(df_clean_san$Sex.x)
df_clean_san$AnolisID <- as.factor(df_clean_san$AnolisID)
# GAM - SVL -----
model1_gam <- gam(SVL ~ Sex.x + s(Age, by = Sex.x), data = df_clean_san, method = "ML")
summary(model1_gam)
model1.2_gam <- gam(SVL ~ s(Age, by=Sex.x), data = df_clean_san, method = "ML")
summary(model1.2_gam)
model1.3_gam <- gam(SVL ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_san, method="ML")
summary(model1.3_gam)
# Model comparison
AIC(model1_gam, model1.2_gam, model1.3_gam) # first model is the best with lowest AIC
# Create new date object to store model predicted values
newdat <- expand.grid(Age = seq(min(df_clean_san$Age, na.rm = TRUE),max(df_clean_san$Age, na.rm = TRUE), length.out = 200),
                      Sex.x = unique(df_clean_san$Sex.x,))
# Store predicted values of SVL from the model
newdat$SVL_pred <- predict(model1_gam, newdata = newdat, type="response")
# Display model and measured values using ggplot
salvador_svl <- ggplot() +
  geom_point(
    data = df_clean_san,
    aes(x = Age, y = SVL, color = Sex.x),
    alpha = 0.3,
    position = position_jitter(width = 0.1)) +
  geom_line(data = newdat,aes(x = Age, y = SVL_pred, color = Sex.x, group = Sex.x),
    linewidth = 1) +
  theme_classic() +
  labs(x = "Age", y = "SVL", color = "Sex", fill = "Sex")


# GAM - Head Width -----
# Clean data for each variable
df_model_san <- df_filtered_sub_log_new %>% dplyr::filter(Population == "S")
df_clean_san <- df_model_san %>% dplyr::filter(!is.na(Head.Width),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_san$Sex.x <- as.factor(df_clean_san$Sex.x)
df_clean_san$AnolisID <- as.factor(df_clean_san$AnolisID)
# Gam models
model2_gam <- gam(Head.Width ~ Sex.x + s(Age, by = Sex.x), data = df_clean_san, method = "ML")
summary(model2_gam)
model2.2_gam <- gam(Head.Width ~ s(Age, by=Sex.x), data = df_clean_san, method = "ML")
summary(model2.2_gam)
model2.3_gam <- gam(Head.Width ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_san, method="ML")
summary(model2.3_gam)
# Model comparison
AIC(model2_gam, model2.2_gam, model2.3_gam) # Third model is best supported by lowest AIC
# Create new date object to store model predicted values
newdat <- expand.grid(Age = seq(min(df_clean_san$Age, na.rm = TRUE),max(df_clean_san$Age, na.rm = TRUE), length.out = 200),
                      Sex.x = unique(df_clean_san$Sex.x,))
# Store predicted values of SVL from the model
newdat$HW_pred <- predict(model2.3_gam, newdata = newdat, type="response")
# Display model and measured values using ggplot
ggplot() +
  geom_point(
    data = df_clean_san,
    aes(x = Age, y = Head.Width, color = Sex.x),
    alpha = 0.3,
    position = position_jitter(width = 0.1)) +
  geom_line(data = newdat, aes(x = Age, y = HW_pred, color = Sex.x, group = Sex.x),
    linewidth = 1) +
  theme_classic() +
  labs(x = "Age", y = "Head Width", color = "Sex", fill = "Sex")


# GAM - Head Length -----
# Clean data for each variable
df_model_san <- df_filtered_sub_log_new %>% dplyr::filter(Population == "S")
df_clean_san <- df_model_san %>% dplyr::filter(!is.na(Head.Length),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_san$Sex.x <- as.factor(df_clean_san$Sex.x)
df_clean_san$AnolisID <- as.factor(df_clean_san$AnolisID)
# Gam models
model3_gam <- gam(Head.Length ~ Sex.x + s(Age, by = Sex.x), data = df_clean_san, method = "ML")
summary(model3_gam)
model3.2_gam <- gam(Head.Length ~ s(Age, by=Sex.x), data = df_clean_san, method = "ML")
summary(model3.2_gam)
model3.3_gam <- gam(Head.Length ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_san, method="ML")
summary(model3.3_gam)
# Model comparison
AIC(model3_gam, model3.2_gam, model3.3_gam) # Third model is best supported by lowest AIC

# GAM - Pectoral Width -----
# Clean data for each variable
df_model_san <- df_filtered_sub_log_new %>% dplyr::filter(Population == "S")
df_clean_san <- df_model_san %>% dplyr::filter(!is.na(Pectoral),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_san$Sex.x <- as.factor(df_clean_san$Sex.x)
df_clean_san$AnolisID <- as.factor(df_clean_san$AnolisID)
# Gam models
model4_gam <- gam(Pectoral ~ Sex.x + s(Age, by = Sex.x), data = df_clean_san, method = "ML")
summary(model4_gam)
model4.2_gam <- gam(Pectoral ~ s(Age, by=Sex.x), data = df_clean_san, method = "ML")
summary(model4.3_gam)
model4.3_gam <- gam(Pectoral ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_san, method="ML")
summary(model4.3_gam)
# Model comparison
AIC(model4_gam, model4.2_gam, model4.3_gam) # Third model is best supported by lowest AIC

# GAM - Pelvic Width -----
# Clean data for each variable
df_model_san <- df_filtered_sub_log_new %>% dplyr::filter(Population == "S")
df_clean_san <- df_model_san %>% dplyr::filter(!is.na(Pelvic),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_san$Sex.x <- as.factor(df_clean_san$Sex.x)
df_clean_san$AnolisID <- as.factor(df_clean_san$AnolisID)
# Gam models
model5_gam <- gam(Pelvic ~ Sex.x + s(Age, by = Sex.x), data = df_clean_san, method = "ML")
summary(model5_gam)
model5.2_gam <- gam(Pelvic ~ s(Age, by=Sex.x), data = df_clean_san, method = "ML")
summary(model5.2_gam)
model5.3_gam <- gam(Pelvic ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_san, method="ML")
summary(model5.3_gam)
# Model comparison
AIC(model5_gam, model5.2_gam, model5.3_gam) # Third model is best supported by lowest AIC

# GAM - Humerus Length -----
# Clean data for each variable
df_model_san <- df_filtered_sub_log_new %>% dplyr::filter(Population == "S")
df_clean_san <- df_model_san %>% dplyr::filter(!is.na(Humerus),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_san$Sex.x <- as.factor(df_clean_san$Sex.x)
df_clean_san$AnolisID <- as.factor(df_clean_san$AnolisID)
# Gam models
model6_gam <- gam(Humerus ~ Sex.x + s(Age, by = Sex.x), data = df_clean_san, method = "ML")
summary(model6_gam)
model6.2_gam <- gam(Humerus ~ s(Age, by=Sex.x), data = df_clean_san, method = "ML")
summary(model6.2_gam)
model6.3_gam <- gam(Humerus ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_san, method="ML")
summary(model6.3_gam)
# Model comparison
AIC(model6_gam, model6.2_gam, model6.3_gam) # Third model is best supported by lowest AIC

# GAM - Radius Length -----
# Clean data for each variable
df_model_san <- df_filtered_sub_log_new %>% dplyr::filter(Population == "S")
df_clean_san <- df_model_san %>% dplyr::filter(!is.na(Radius),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_san$Sex.x <- as.factor(df_clean_san$Sex.x)
df_clean_san$AnolisID <- as.factor(df_clean_san$AnolisID)
# Gam models
model7_gam <- gam(Radius ~ Sex.x + s(Age, by = Sex.x), data = df_clean_san, method = "ML")
summary(model7_gam)
model7.2_gam <- gam(Radius ~ s(Age, by=Sex.x), data = df_clean_san, method = "ML")
summary(model7.2_gam)
model7.3_gam <- gam(Radius ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_san, method="ML")
summary(model7.3_gam)
# Model comparison
AIC(model7_gam, model7.2_gam, model7.3_gam) # Third model is best supported by lowest AIC

# GAM - Ulna Length -----
# Clean data for each variable
df_model_san <- df_filtered_sub_log_new %>% dplyr::filter(Population == "S")
df_clean_san <- df_model_san %>% dplyr::filter(!is.na(Ulna),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_san$Sex.x <- as.factor(df_clean_san$Sex.x)
df_clean_san$AnolisID <- as.factor(df_clean_san$AnolisID)
# Gam models
model8_gam <- gam(Ulna ~ Sex.x + s(Age, by = Sex.x), data = df_clean_san, method = "ML")
summary(model8_gam)
model8.2_gam <- gam(Ulna ~ s(Age, by=Sex.x), data = df_clean_san, method = "ML")
summary(model8.2_gam)
model8.3_gam <- gam(Ulna ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_san, method="ML")
summary(model8.3_gam)
# Model comparison
AIC(model8_gam, model8.2_gam, model8.3_gam) # Third model is best supported by lowest AIC

# GAM - Femur Length -----
# Clean data for each variable
df_model_san <- df_filtered_sub_log_new %>% dplyr::filter(Population == "S")
df_clean_san <- df_model_san %>% dplyr::filter(!is.na(Femur),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_san$Sex.x <- as.factor(df_clean_san$Sex.x)
df_clean_san$AnolisID <- as.factor(df_clean_san$AnolisID)
# Gam models
model9_gam <- gam(Femur ~ Sex.x + s(Age, by = Sex.x), data = df_clean_san, method = "ML")
summary(model9_gam)
model9.2_gam <- gam(Femur ~ s(Age, by=Sex.x), data = df_clean_san, method = "ML")
summary(model9.2_gam)
model9.3_gam <- gam(Femur ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_san, method="ML")
summary(model9.3_gam)
# Model comparison
AIC(model9_gam, model9.2_gam, model9.3_gam) # Third model is best supported by lowest AIC

# GAM - Tibia Length -----
# Clean data for each variable
df_model_san <- df_filtered_sub_log_new %>% dplyr::filter(Population == "S")
df_clean_san <- df_model_san %>% dplyr::filter(!is.na(Tibia),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_san$Sex.x <- as.factor(df_clean_san$Sex.x)
df_clean_san$AnolisID <- as.factor(df_clean_san$AnolisID)
# Gam models
model10_gam <- gam(Tibia ~ Sex.x + s(Age, by = Sex.x), data = df_clean_san, method = "ML")
summary(model9_gam)
model10.2_gam <- gam(Tibia ~ s(Age, by=Sex.x), data = df_clean_san, method = "ML")
summary(model9.2_gam)
model10.3_gam <- gam(Tibia ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_san, method="ML")
summary(model9.3_gam)
# Model comparison
AIC(model10_gam, model10.2_gam, model10.3_gam) # Third model is best supported by lowest AIC

##### Filter for Andros population-----
df_model_and <- df_filtered_sub_log_new %>% dplyr::filter(Population == "A")
View(df_model_and)

# Plots the data
# SVL and Age (body size seems to be closer between sexes at final age)
ggplot(df_model_and, aes(Age, SVL, fill=Sex.x)) +
  geom_point() + 
  geom_smooth() +
  theme_bw()

# GAM - SVL -----
# Clean data for each variable
df_model_and <- df_filtered_sub_log_new %>% dplyr::filter(Population == "A")
df_clean_and <- df_model_and %>% dplyr::filter(!is.na(SVL),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_and$Sex.x <- as.factor(df_clean_and$Sex.x)
df_clean_and$AnolisID <- as.factor(df_clean_and$AnolisID)
# Gam models
model1_and_gam <- gam(SVL ~ Sex.x + s(Age, by = Sex.x), data = df_clean_and, method = "ML")
summary(model1_and_gam)
model1.2_and_gam <- gam(SVL ~ s(Age, by=Sex.x), data = df_clean_and, method = "ML")
summary(model1.2_and_gam)
model1.3_and_gam <- gam(SVL ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_and, method="ML")
summary(model1.3_and_gam)
# Model comparison
AIC(model1_and_gam, model1.2_and_gam, model1.3_and_gam) # Third model is best supported by lowest AIC
# Create new date object to store model predicted values
newdat <- expand.grid(Age = seq(min(df_clean_and$Age, na.rm = TRUE),max(df_clean_and$Age, na.rm = TRUE), length.out = 200),
                      Sex.x = unique(df_clean_and$Sex.x,))
# Store predicted values of SVL from the model
newdat$SVL_pred <- predict(model1_and_gam, newdata = newdat, type="response")
# Display model and measured values using ggplot
andros_svl <- ggplot() +
  geom_point(
    data = df_clean_and,
    aes(x = Age, y = SVL, color = Sex.x),
    alpha = 0.3,
    position = position_jitter(width = 0.1)) +
  geom_line(data = newdat,aes(x = Age, y = SVL_pred, color = Sex.x, group = Sex.x),
            linewidth = 1) +
  theme_classic() +
  labs(x = "Age", y = "SVL", color = "Sex", fill = "Sex")



# GAM - Head Width -----
# Clean data for each variable
df_model_and <- df_filtered_sub_log_new %>% dplyr::filter(Population == "A")
df_clean_and <- df_model_and %>% dplyr::filter(!is.na(Head.Length),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_and$Sex.x <- as.factor(df_clean_and$Sex.x)
df_clean_and$AnolisID <- as.factor(df_clean_and$AnolisID)
View(df_clean_and)
# Gam models
model2_and_gam <- gam(Head.Width ~ Sex.x + s(Age, by = Sex.x), data = df_clean_and, method = "ML")
summary(model2_and_gam)
model2.2_and_gam <- gam(Head.Width ~ s(Age, by=Sex.x), data = df_clean_and, method = "ML")
summary(model2.2_and_gam)
model2.3_and_gam <- gam(Head.Width ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_and, method="ML")
summary(model2.3_and_gam)
# Model comparison
AIC(model2_and_gam, model2.2_and_gam, model2.3_and_gam) # Third model is best supported by lowest AIC

# GAM - Head Length -----
# Clean data for each variable
df_model_and <- df_filtered_sub_log_new %>% dplyr::filter(Population == "A")
df_clean_and <- df_model_and %>% dplyr::filter(!is.na(Head.Length),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_and$Sex.x <- as.factor(df_clean_and$Sex.x)
df_clean_and$AnolisID <- as.factor(df_clean_and$AnolisID)
# Gam models
model3_and_gam <- gam(Head.Length ~ Sex.x + s(Age, by = Sex.x), data = df_clean_and, method = "ML")
summary(model3_and_gam)
model3.2_and_gam <- gam(Head.Length ~ s(Age, by=Sex.x), data = df_clean_and, method = "ML")
summary(model3.2_and_gam)
model3.3_and_gam <- gam(Head.Length ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_and, method="ML")
summary(model3.3_and_gam)
# Model comparison
AIC(model3_and_gam, model3.2_and_gam, model3.3_and_gam) # third model is best supported by lowest AIC

# GAM - Pectoral Width -----
# Clean data for each variable
df_model_and <- df_filtered_sub_log_new %>% dplyr::filter(Population == "A")
df_clean_and <- df_model_and %>% dplyr::filter(!is.na(Pectoral),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_and$Sex.x <- as.factor(df_clean_and$Sex.x)
df_clean_and$AnolisID <- as.factor(df_clean_and$AnolisID)
# Gam models
model4_and_gam <- gam(Pectoral ~ Sex.x + s(Age, by = Sex.x), data = df_clean_and, method = "ML")
summary(model4_and_gam)
model4.2_and_gam <- gam(Pectoral ~ s(Age, by=Sex.x), data = df_clean_and, method = "ML")
summary(model4.2_and_gam)
model4.3_and_gam <- gam(Pectoral ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_and, method="ML")
summary(model4.3_and_gam)
# Model comparison
AIC(model4_and_gam, model4.2_and_gam, model4.3_and_gam) # first model is best supported by lowest AIC

# GAM - Pelvic Width -----
# Clean data for each variable
df_model_and <- df_filtered_sub_log_new %>% dplyr::filter(Population == "A")
df_clean_and <- df_model_and %>% dplyr::filter(!is.na(Pelvic),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_and$Sex.x <- as.factor(df_clean_and$Sex.x)
df_clean_and$AnolisID <- as.factor(df_clean_and$AnolisID)
# Gam models
model5_and_gam <- gam(Pelvic ~ Sex.x + s(Age, by = Sex.x), data = df_clean_and, method = "ML")
summary(model5_and_gam)
model5.2_and_gam <- gam(Pelvic ~ s(Age, by=Sex.x), data = df_clean_and, method = "ML")
summary(model5.2_and_gam)
model5.3_and_gam <- gam(Pelvic ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_and, method="ML")
summary(model5.3_and_gam)
# Model comparison
AIC(model5_and_gam, model5.2_and_gam, model5.3_and_gam) # third model is best supported by lowest AIC

# GAM - Humerus Length -----
# Clean data for each variable
df_model_and <- df_filtered_sub_log_new %>% dplyr::filter(Population == "A")
df_clean_and <- df_model_and %>% dplyr::filter(!is.na(Humerus),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_and$Sex.x <- as.factor(df_clean_and$Sex.x)
df_clean_and$AnolisID <- as.factor(df_clean_and$AnolisID)
# Gam models
model6_and_gam <- gam(Humerus ~ Sex.x + s(Age, by = Sex.x), data = df_clean_and, method = "ML")
summary(model6_and_gam)
model6.2_and_gam <- gam(Humerus ~ s(Age, by=Sex.x), data = df_clean_and, method = "ML")
summary(model6.2_and_gam)
model6.3_and_gam <- gam(Humerus ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_and, method="ML")
summary(model6.3_and_gam)
# Model comparison
AIC(model6_and_gam, model6.2_and_gam, model6.3_and_gam) # third model is best supported by lowest AIC

# GAM - Radius Length -----
# Clean data for each variable
df_model_and <- df_filtered_sub_log_new %>% dplyr::filter(Population == "A")
df_clean_and <- df_model_and %>% dplyr::filter(!is.na(Radius),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_and$Sex.x <- as.factor(df_clean_and$Sex.x)
df_clean_and$AnolisID <- as.factor(df_clean_and$AnolisID)
# Gam models
model7_and_gam <- gam(Radius ~ Sex.x + s(Age, by = Sex.x), data = df_clean_and, method = "ML")
summary(model7_and_gam)
model7.2_and_gam <- gam(Radius ~ s(Age, by=Sex.x), data = df_clean_and, method = "ML")
summary(model7.2_and_gam)
model7.3_and_gam <- gam(Radius ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_and, method="ML")
summary(model7.3_and_gam)
# Model comparison
AIC(model7_and_gam, model7.2_and_gam, model7.3_and_gam) # third model is best supported by lowest AIC


# GAM - Ulna Length -----
# Clean data for each variable
df_model_and <- df_filtered_sub_log_new %>% dplyr::filter(Population == "A")
df_clean_and <- df_model_and %>% dplyr::filter(!is.na(Ulna),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_and$Sex.x <- as.factor(df_clean_and$Sex.x)
df_clean_and$AnolisID <- as.factor(df_clean_and$AnolisID)
# Gam models
model8_and_gam <- gam(Ulna ~ Sex.x + s(Age, by = Sex.x), data = df_clean_and, method = "ML")
summary(model8_and_gam)
model8.2_and_gam <- gam(Ulna ~ s(Age, by=Sex.x), data = df_clean_and, method = "ML")
summary(model8.2_and_gam)
model8.3_and_gam <- gam(Ulna ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_and, method="ML")
summary(model8.3_and_gam)
# Model comparison
AIC(model8_and_gam, model8.2_and_gam, model8.3_and_gam) # third model is best supported by lowest AIC

# GAM - Femur Length -----
# Clean data for each variable
df_model_and <- df_filtered_sub_log_new %>% dplyr::filter(Population == "A")
df_clean_and <- df_model_and %>% dplyr::filter(!is.na(Femur),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_and$Sex.x <- as.factor(df_clean_and$Sex.x)
df_clean_and$AnolisID <- as.factor(df_clean_and$AnolisID)
# Gam models
model9_and_gam <- gam(Femur ~ Sex.x + s(Age, by = Sex.x), data = df_clean_and, method = "ML")
summary(model9_and_gam)
model9.2_and_gam <- gam(Femur ~ s(Age, by=Sex.x), data = df_clean_and, method = "ML")
summary(model9.2_and_gam)
model9.3_and_gam <- gam(Femur ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_and, method="ML")
summary(model9.3_and_gam)
# Model comparison
AIC(model9_and_gam, model9.2_and_gam, model9.3_and_gam) # third model is best supported by lowest AIC

# GAM - Tibia Length -----
# Clean data for each variable
df_model_and <- df_filtered_sub_log_new %>% dplyr::filter(Population == "A")
df_clean_and <- df_model_and %>% dplyr::filter(!is.na(Tibia),!is.na(Age),!is.na(Sex.x),!is.na(AnolisID))
# Convert sex to a factor instead of character
df_clean_and$Sex.x <- as.factor(df_clean_and$Sex.x)
df_clean_and$AnolisID <- as.factor(df_clean_and$AnolisID)
# Gam models
model10_and_gam <- gam(Tibia ~ Sex.x + s(Age, by = Sex.x), data = df_clean_and, method = "ML")
summary(model10_and_gam)
model10.2_and_gam <- gam(Tibia ~ s(Age, by=Sex.x), data = df_clean_and, method = "ML")
summary(model10.2_and_gam)
model10.3_and_gam <- gam(Tibia ~ Sex.x + s(Age, Sex.x, bs="fs"), data=df_clean_and, method="ML")
summary(model10.3_and_gam)
# Model comparison
AIC(model10_and_gam, model10.2_and_gam, model10.3_and_gam) # third model is best supported by lowest AIC

# Combine ggplots for SVL
grid.arrange(salvador_svl_new, andros_svl_new, ncol = 2)

salvador_svl_new <- salvador_svl + 
  theme(legend.position = "none") +
  annotate("text", x=120, y=4, label="San Salvador", size=5) +
  ylab("Snout-vent length") +
  xlab("Days since hatch")
andros_svl_new <- andros_svl +
  theme(axis.title.y = element_blank(), labs(y="Days since hatch")) +
  annotate("text", x=60, y=4, label="Andros", size=5) +
  ylab("Snout-vent length") +
  xlab("Days since hatch")

################## Are there morphological differences between sex at the end of development ? ##############
View(df_filtered_sub_log)
###### Filter for San Salvador population -------
df_model_san <- df_filtered_sub_log_new %>% dplyr::filter(Population == "S")
View(df_model_san)

### Filter for SVL
df_model_san_svl <- df_model_san %>% dplyr::select("AnolisID", "SVL", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_san_svl2 <- df_model_san_svl %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
View(df_model_san_svl2)
# Only keep rows with ages ranging from 600 to 900 days
df_model_san_svl3 <- df_model_san_svl2 %>%
  filter(Age >= 600 & Age <= 900)
View(df_model_san_svl3)
# Does svl significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - SVL ~ Sex 
lm_svl_san <- lm(SVL ~ Sex.x, data=df_model_san_svl3)
summary(lm_svl_san)
# Test for normality
shapiro.test(lm_svl_san$residuals) # data is not normally distributed.
par(mfrow = c(2,2))
plot(lm_svl_san) # Residuals vs fitted is linear
# Try running generalized linear model
glm_svl_san <- glm(SVL ~ Sex.x, data=df_model_san_svl3)
summary(glm_svl_san)
glm_svl_san2 <- glm(SVL ~ Sex.x + Age, data=df_model_san_svl3)
summary(glm_svl_san2) # no significant effect of age


### Filter for Head Width
df_model_san_hw <- df_model_san %>% dplyr::select("AnolisID", "Head.Width", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_san_hw2 <- df_model_san_hw %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
# Only keep rows with ages ranging from 600 to 900 days
df_model_san_hw3 <- df_model_san_hw2 %>%
  filter(Age >= 600 & Age <= 900)
# Does head width significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - Head width ~ Sex 
# Try running generalized linear model
glm_hw_san <- glm(Head.Width ~ Sex.x, data=df_model_san_hw3)
summary(glm_hw_san)
glm_hw_san2 <- glm(Head.Width ~ Sex.x + Age, data=df_model_san_hw3)
summary(glm_hw_san2) # significant effect of age

### Filter for Head Length
df_model_san_hl <- df_model_san %>% dplyr::select("AnolisID", "Head.Length", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_san_hl2 <- df_model_san_hl %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
# Only keep rows with ages ranging from 600 to 900 days
df_model_san_hl3 <- df_model_san_hl2 %>%
  filter(Age >= 600 & Age <= 900)
# Does head width significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - Head width ~ Sex 
# Try running generalized linear model
glm_hl_san <- glm(Head.Length ~ Sex.x, data=df_model_san_hl3)
summary(glm_hl_san)
glm_hl_san2 <- glm(Head.Length ~ Sex.x + Age, data=df_model_san_hl3)
summary(glm_hl_san2) # mo significant effect of age

### Filter for Pectoral width
df_model_san_pw <- df_model_san %>% dplyr::select("AnolisID", "Pectoral", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_san_pw2 <- df_model_san_pw %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
# Only keep rows with ages ranging from 600 to 900 days
df_model_san_pw3 <- df_model_san_pw2 %>%
  filter(Age >= 600 & Age <= 900)
# Does head width significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - Head width ~ Sex 
# Try running generalized linear model
glm_pw_san <- glm(Pectoral ~ Sex.x, data=df_model_san_pw3)
summary(glm_pw_san)
glm_pw_san2 <- glm(Pectoral ~ Sex.x + Age, data=df_model_san_pw3)
summary(glm_pw_san2) # no significant effect of age

### Filter for Pelvic width
df_model_san_pvw <- df_model_san %>% dplyr::select("AnolisID", "Pelvic", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_san_pvw2 <- df_model_san_pvw %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
# Only keep rows with ages ranging from 600 to 900 days
df_model_san_pvw3 <- df_model_san_pvw2 %>%
  filter(Age >= 600 & Age <= 900)
# Does head width significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - Head width ~ Sex 
# Try running generalized linear model
glm_pvw_san <- glm(Pelvic ~ Sex.x, data=df_model_san_pvw3)
summary(glm_pvw_san)
glm_pvw_san2 <- glm(Pelvic ~ Sex.x + Age, data=df_model_san_pvw3)
summary(glm_pvw_san2) # no significant effect of age

### Filter for Humerus length
df_model_san_hl <- df_model_san %>% dplyr::select("AnolisID", "Humerus", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_san_hl2 <- df_model_san_hl %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
# Only keep rows with ages ranging from 600 to 900 days
df_model_san_hl3 <- df_model_san_hl2 %>%
  filter(Age >= 600 & Age <= 900)
# Does head width significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - Head width ~ Sex 
# Try running generalized linear model
glm_hl_san <- glm(Humerus ~ Sex.x, data=df_model_san_hl3)
summary(glm_hl_san)
glm_hl_san2 <- glm(Humerus ~ Sex.x + Age, data=df_model_san_hl3)
summary(glm_hl_san2) # no significant effect of age


### Filter for Radius length
df_model_san_rl <- df_model_san %>% dplyr::select("AnolisID", "Radius", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_san_rl2 <- df_model_san_rl %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
# Only keep rows with ages ranging from 600 to 900 days
df_model_san_rl3 <- df_model_san_rl2 %>%
  filter(Age >= 600 & Age <= 900)
# Does head width significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - Head width ~ Sex 
# Try running generalized linear model
glm_rl_san <- glm(Radius ~ Sex.x, data=df_model_san_rl3)
summary(glm_rl_san)
glm_rl_san2 <- glm(Radius ~ Sex.x + Age, data=df_model_san_rl3)
summary(glm_rl_san2) # no significant effect of age

### Filter for Ulna length
df_model_san_ul <- df_model_san %>% dplyr::select("AnolisID", "Ulna", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_san_ul2 <- df_model_san_ul %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
# Only keep rows with ages ranging from 600 to 900 days
df_model_san_ul3 <- df_model_san_ul2 %>%
  filter(Age >= 600 & Age <= 900)
# Does head width significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - Head width ~ Sex 
# Try running generalized linear model
glm_ul_san <- glm(Ulna ~ Sex.x, data=df_model_san_ul3)
summary(glm_ul_san)
glm_ul_san2 <- glm(Ulna ~ Sex.x + Age, data=df_model_san_ul3)
summary(glm_ul_san2) # no significant effect of age

### Filter for Femur length
df_model_san_fl <- df_model_san %>% dplyr::select("AnolisID", "Femur", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_san_fl2 <- df_model_san_fl %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
# Only keep rows with ages ranging from 600 to 900 days
df_model_san_fl3 <- df_model_san_fl2 %>%
  filter(Age >= 600 & Age <= 900)
# Does head width significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - Head width ~ Sex 
# Try running generalized linear model
glm_fl_san <- glm(Femur ~ Sex.x, data=df_model_san_fl3)
summary(glm_fl_san)
glm_fl_san2 <- glm(Femur ~ Sex.x + Age, data=df_model_san_fl3)
summary(glm_fl_san2) # no significant effect of age

### Filter for Tibia length
df_model_san_tl <- df_model_san %>% dplyr::select("AnolisID", "Tibia", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_san_tl2 <- df_model_san_tl %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
# Only keep rows with ages ranging from 600 to 900 days
df_model_san_tl3 <- df_model_san_tl2 %>%
  filter(Age >= 600 & Age <= 900)
# Does head width significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - Head width ~ Sex 
# Try running generalized linear model
glm_tl_san <- glm(Tibia ~ Sex.x, data=df_model_san_tl3)
summary(glm_tl_san)
glm_tl_san2 <- glm(Tibia ~ Sex.x + Age, data=df_model_san_tl3)
summary(glm_tl_san2) # no significant effect of age

###### Filter for Andros population -------
df_model_and <- df_filtered_sub_log_new %>% dplyr::filter(Population == "A")
View(df_model_and)

### Filter for SVL
df_model_and_svl <- df_model_and %>% dplyr::select("AnolisID", "SVL", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_and_svl2 <- df_model_and_svl %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
View(df_model_and_svl2)
# Only keep rows with ages ranging from 600 to 900 days
df_model_and_svl3 <- df_model_and_svl2 %>%
  filter(Age >= 600 & Age <= 900)
View(df_model_and_svl3)
# Does svl significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - SVL ~ Sex 
lm_svl_and <- lm(SVL ~ Sex.x, data=df_model_and_svl3)
summary(lm_svl_and)
# Test for normality
shapiro.test(lm_svl_and$residuals) # data is not normally distributed.
par(mfrow = c(2,2))
plot(lm_svl_and) # Residuals vs fitted is linear
# Try running generalized linear model
glm_svl_and <- glm(SVL ~ Sex.x, data=df_model_and_svl3)
summary(glm_svl_and)
glm_svl_and2 <- glm(SVL ~ Sex.x + Age, data=df_model_and_svl3)
summary(glm_svl_and2) # significant effect of age

### Filter for Head Width
df_model_and_hw <- df_model_and %>% dplyr::select("AnolisID", "Head.Width", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_and_hw2 <- df_model_and_hw %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
# Only keep rows with ages ranging from 600 to 900 days
df_model_and_hw3 <- df_model_and_hw2 %>%
  filter(Age >= 600 & Age <= 900)
# Does head width significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - Head width ~ Sex 
# Try running generalized linear model
glm_hw_and <- glm(Head.Width ~ Sex.x, data=df_model_and_hw3)
summary(glm_hw_and)
glm_hw_and2 <- glm(Head.Width ~ Sex.x + Age, data=df_model_and_hw3)
summary(glm_hw_and2) # significant effect of age

### Filter for Head Length
df_model_and_hl <- df_model_and %>% dplyr::select("AnolisID", "Head.Length", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_and_hl2 <- df_model_and_hl %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
# Only keep rows with ages ranging from 600 to 900 days
df_model_and_hl3 <- df_model_and_hl2 %>%
  filter(Age >= 600 & Age <= 900)
# Does head width significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - Head width ~ Sex 
# Try running generalized linear model
glm_hl_and <- glm(Head.Length ~ Sex.x, data=df_model_and_hl3)
summary(glm_hl_and)
glm_hl_and2 <- glm(Head.Length ~ Sex.x + Age, data=df_model_and_hl3)
summary(glm_hl_and2) # no significant effect of age

### Filter for Pectoral width
df_model_and_pw <- df_model_and %>% dplyr::select("AnolisID", "Pectoral", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_and_pw2 <- df_model_and_pw %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
# Only keep rows with ages ranging from 600 to 900 days
df_model_and_pw3 <- df_model_and_pw2 %>%
  filter(Age >= 600 & Age <= 900)
# Does head width significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - Head width ~ Sex 
# Try running generalized linear model
glm_pw_and <- glm(Pectoral ~ Sex.x, data=df_model_and_pw3)
summary(glm_pw_and)
glm_pw_and2 <- glm(Pectoral ~ Sex.x + Age, data=df_model_and_pw3)
summary(glm_pw_san2) # no significant effect of age

### Filter for Pelvic width
df_model_and_pvw <- df_model_and %>% dplyr::select("AnolisID", "Pelvic", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_and_pvw2 <- df_model_and_pvw %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
# Only keep rows with ages ranging from 600 to 900 days
df_model_and_pvw3 <- df_model_and_pvw2 %>%
  filter(Age >= 600 & Age <= 900)
# Does head width significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - Head width ~ Sex 
# Try running generalized linear model
glm_pvw_and <- glm(Pelvic ~ Sex.x, data=df_model_and_pvw3)
summary(glm_pvw_and)
glm_pvw_and2 <- glm(Pelvic ~ Sex.x + Age, data=df_model_and_pvw3)
summary(glm_pvw_and2) # significant effect of age

### Filter for Humerus length
df_model_and_hl <- df_model_and %>% dplyr::select("AnolisID", "Humerus", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_and_hl2 <- df_model_and_hl %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
# Only keep rows with ages ranging from 600 to 900 days
df_model_and_hl3 <- df_model_and_hl2 %>%
  filter(Age >= 600 & Age <= 900)
# Does head width significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - Head width ~ Sex 
# Try running generalized linear model
glm_hl_and <- glm(Humerus ~ Sex.x, data=df_model_and_hl3)
summary(glm_hl_and)
glm_hl_and2 <- glm(Humerus ~ Sex.x + Age, data=df_model_and_hl3)
summary(glm_hl_and2) # no significant effect of age


### Filter for Radius length
df_model_and_rl <- df_model_and %>% dplyr::select("AnolisID", "Radius", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_and_rl2 <- df_model_and_rl %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
# Only keep rows with ages ranging from 600 to 900 days
df_model_and_rl3 <- df_model_and_rl2 %>%
  filter(Age >= 600 & Age <= 900)
# Does head width significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - Head width ~ Sex 
# Try running generalized linear model
glm_rl_and <- glm(Radius ~ Sex.x, data=df_model_and_rl3)
summary(glm_rl_and)
glm_rl_and2 <- glm(Radius ~ Sex.x + Age, data=df_model_and_rl3)
summary(glm_rl_and2) # no significant effect of age

### Filter for Ulna length
df_model_and_ul <- df_model_and %>% dplyr::select("AnolisID", "Ulna", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_and_ul2 <- df_model_and_ul %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
# Only keep rows with ages ranging from 600 to 900 days
df_model_and_ul3 <- df_model_and_ul2 %>%
  filter(Age >= 600 & Age <= 900)
# Does head width significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - Head width ~ Sex 
# Try running generalized linear model
glm_ul_and <- glm(Ulna ~ Sex.x, data=df_model_and_ul3)
summary(glm_ul_and)
glm_ul_and2 <- glm(Ulna ~ Sex.x + Age, data=df_model_and_ul3)
summary(glm_ul_and2) # no significant effect of age

### Filter for Femur length
df_model_and_fl <- df_model_and %>% dplyr::select("AnolisID", "Femur", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_and_fl2 <- df_model_and_fl %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
# Only keep rows with ages ranging from 600 to 900 days
df_model_and_fl3 <- df_model_and_fl2 %>%
  filter(Age >= 600 & Age <= 900)
# Does head width significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - Head width ~ Sex 
# Try running generalized linear model
glm_fl_and <- glm(Femur ~ Sex.x, data=df_model_and_fl3)
summary(glm_fl_and)
glm_fl_and2 <- glm(Femur ~ Sex.x + Age, data=df_model_and_fl3)
summary(glm_fl_and2) # significant effect of age

### Filter for Tibia length
df_model_and_tl <- df_model_and %>% dplyr::select("AnolisID", "Tibia", "Sex.x", "Age")
### Group by individual and select largest age 
df_model_and_tl2 <- df_model_and_tl %>%
  group_by(AnolisID) %>%
  slice_max(Age, n = 1) %>%
  ungroup()
# Only keep rows with ages ranging from 600 to 900 days
df_model_and_tl3 <- df_model_and_tl2 %>%
  filter(Age >= 600 & Age <= 900)
# Does head width significantly differ between males and females late in development (Age 600 - 900 days)?
# Model - Head width ~ Sex 
# Try running generalized linear model
glm_tl_and <- glm(Tibia ~ Sex.x, data=df_model_and_tl3)
summary(glm_tl_and)
glm_tl_and2 <- glm(Tibia ~ Sex.x + Age, data=df_model_and_tl3)
summary(glm_tl_and2) # no significant effect of age

############################## Multivariate Generalised Linear Mixed Models using MCMCglmm ###########################
###### Andros population -----
# Filter for Population = A  and Select AnolisID, Mother, Father, Sex, SVL and other morphological traits
df_filtered_sub_log_A <- df_filtered_sub_log %>% dplyr::filter(Population == "A") %>%
  dplyr::select(c("AnolisID","Sex.x","SVL","Head.Length","Head.Width","Pectoral","Pelvic","Humerus",
                  "Radius","Ulna","Femur","Tibia","Mother","Father"))
View(df_filtered_sub_log_A)

# Change column names of AnolisID, Father and Mother 
colnames(df_filtered_sub_log_A)[which(names(df_filtered_sub_log_A) == "AnolisID")] <- "id"
colnames(df_filtered_sub_log_A)[which(names(df_filtered_sub_log_A) == "Father")] <- "sire"
colnames(df_filtered_sub_log_A)[which(names(df_filtered_sub_log_A) == "Mother")] <- "dam"
# make id, dam, sire factors
df_filtered_sub_log_A$id <- as.factor(df_filtered_sub_log_A$id)
df_filtered_sub_log_A$dam <- as.factor(df_filtered_sub_log_A$dam)
df_filtered_sub_log_A$sire <- as.factor(df_filtered_sub_log_A$sire)
# Need separate columns later to estimate addititive and maternal genetic effects
df_filtered_sub_log_A$ANIMAL <- as.factor(df_filtered_sub_log_A$id)
df_filtered_sub_log_A$MOTHER <- as.factor(df_filtered_sub_log_A$dam)
# Repeat for connection to the extra matrices in the models
df_filtered_sub_log_A$ANIMAL2 <- as.factor(df_filtered_sub_log_A$id)
df_filtered_sub_log_A$MOTHER2 <- as.factor(df_filtered_sub_log_A$dam)
# Need to ensure other values are correctly read in as numeric/factors
df_filtered_sub_log_A$SVL <- as.numeric(df_filtered_sub_log_A$SVL)
df_filtered_sub_log_A$Head.Length <- as.numeric(df_filtered_sub_log_A$Head.Length)
df_filtered_sub_log_A$Head.Width <- as.numeric(df_filtered_sub_log_A$Head.Width)
df_filtered_sub_log_A$Pectoral <- as.numeric(df_filtered_sub_log_A$Pectoral)
df_filtered_sub_log_A$Pelvic <- as.numeric(df_filtered_sub_log_A$Pelvic)
df_filtered_sub_log_A$Humerus <- as.numeric(df_filtered_sub_log_A$Humerus)
df_filtered_sub_log_A$Radius <- as.numeric(df_filtered_sub_log_A$Radius)
df_filtered_sub_log_A$Ulna <- as.numeric(df_filtered_sub_log_A$Ulna)
df_filtered_sub_log_A$Femur <- as.numeric(df_filtered_sub_log_A$Femur)
df_filtered_sub_log_A$Tibia <- as.numeric(df_filtered_sub_log_A$Tibia)

df_filtered_sub_log_A$Sex.x <- as.factor(df_filtered_sub_log_A$Sex.x) # Can be done for all the factors in the data
View(df_filtered_sub_log_A)

# Pedigree for the data
pedigree <- df_filtered_sub_log_A[,c("id","dam","sire")]
View(pedigree) # 885 rows
# Remove duplicate rows
pedigree_unique <- unique(pedigree)
View(pedigree_unique) # 392 rows

# Gets the inverse of the relatedness matrix (how is the relatedness matrix being created, G?)
# inverseA requires every parent (sire or dam) to also exist as a rwo in the pedigree
Ainv <- inverseA(pedigree_unique)$Ainv # Error - individuals appearing as dams but not in pedigree
# Check for missing dams in id
missing_dams <- setdiff(pedigree_unique$dam, pedigree_unique$id)
missing_dams # 90 dams not in id column
missing_sires <- setdiff(pedigree_unique$sire, pedigree_unique$id)
missing_sires # 60 sires not in id column
# Create new rows with missing sires and dams
missing_parents <- unique(c(missing_dams, missing_sires))
missing_parents
new_rows <- data.frame(id=missing_parents, sire=NA, dam=NA)
pedigree_unique_new <- rbind(pedigree_unique, new_rows)
View(pedigree_unique_new) # 542 rows
# Try running inverseA again
Ainv <- inverseA(pedigree_unique_new)$Ainv # Error - dams appearing before their offspring, try orderPed from MasterBayes
# Orders the pedigree so that parents come before their offspring
pedigree_unique_new_ord <- orderPed(pedigree_unique_new)
View(pedigree_unique_new_ord)
# Try running inverseA again
Ainv_obj <- inverseA(pedigree_unique_new_ord)
View(Ainv_obj)
Ainv <- Ainv_obj$Ainv # Use IDs returned by inverseA(); this fixes error below (line 174)
Ainv
# Basic model - includes only additive genetic variance and residual variance on the trait
# What is nu? What is V?
prior1 <- list(G = list(G1 = list(V = 1, nu = 0.002)),
               R = list(V = 1, nu = 0.002))

MCMC_basic <- MCMCglmm(SVL ~ 1, random = ~ ANIMAL, data = df_filtered_sub_log_A, ginverse = list(ANIMAL=Ainv), prior = prior1)
# Error - ANIMAL ginvere must have non-null rownames
# The inverse relationship matrix (Ainv) does not have row names
rownames(Ainv) # NULL; matrix inside of ginverse list should have rownames that exactly match the levels of your ANIMAL column
View(Ainv)
# It runs now
# Basic model for SVL
MCMC_basic <- MCMCglmm(SVL ~ 1, random = ~ ANIMAL, data = df_filtered_sub_log_A, ginverse = list(ANIMAL=Ainv), prior = prior1)
summary(MCMC_basic)
MCMC_basic$ginverse
# The values given under the G-structure in the model summary are the values for the random effects 
# In this case the additive genetic variance (the ANIMAL term) 
# R-structure gives the residual variance estimates (units)
# The effective sample sizes should be over 1000. This can be improved by running the model for a longer
# period of time and increasing the thinning interval to reduce autocorrelation. This can be done by
# altering the nitt (number of iteractions) argument, burnin (the burn-in period)
# and thin (the thinning interval - how often iterations are sampled)

View(df_filtered_sub_log_A)
# Model with SVL and Head length 
MCMC_model1 <- MCMCglmm(Head.Length ~ SVL, random = ~ ANIMAL, data = df_filtered_sub_log_A, ginverse = list(ANIMAL=Ainv), prior = prior1,
                        nitt=13000*10, thin=10*10, burnin=3000*10)
summary(MCMC_model1)
# Model with SVL and Head length 
MCMC_model2 <- MCMCglmm(Head.Width ~ SVL, random = ~ ANIMAL, data = df_filtered_sub_log_A, ginverse = list(ANIMAL=Ainv), prior = prior1,
                        nitt=13000*10, thin=10*10, burnin=3000*10)
summary(MCMC_model2)
# Model with SVL and Pectoral
MCMC_model3 <- MCMCglmm(Pectoral ~ SVL, random = ~ ANIMAL, data = df_filtered_sub_log_A, ginverse = list(ANIMAL=Ainv), prior = prior1,
                        nitt=13000*10, thin=10*10, burnin=3000*10)
summary(MCMC_model3)
# Model with SVL and Pelvic
MCMC_model4 <- MCMCglmm(Pelvic ~ SVL, random = ~ ANIMAL, data = df_filtered_sub_log_A, ginverse = list(ANIMAL=Ainv), prior = prior1,
                        nitt=13000*10, thin=10*10, burnin=3000*10)
summary(MCMC_model4)
# Model with SVL and Humerus
MCMC_model5 <- MCMCglmm(Humerus ~ SVL, random = ~ ANIMAL, data = df_filtered_sub_log_A, ginverse = list(ANIMAL=Ainv), prior = prior1,
                        nitt=13000*10, thin=10*10, burnin=3000*10)
summary(MCMC_model5)
# Model with SVL and Radius
MCMC_model6 <- MCMCglmm(Radius ~ SVL, random = ~ ANIMAL, data = df_filtered_sub_log_A, ginverse = list(ANIMAL=Ainv), prior = prior1,
                        nitt=13000*10, thin=10*10, burnin=3000*10)
summary(MCMC_model6)
# Model with SVL and Ulna
MCMC_model7 <- MCMCglmm(Ulna ~ SVL, random = ~ ANIMAL, data = df_filtered_sub_log_A, ginverse = list(ANIMAL=Ainv), prior = prior1,
                        nitt=13000*10, thin=10*10, burnin=3000*10)
summary(MCMC_model7)
# Model with SVL and Femur
MCMC_model8 <- MCMCglmm(Femur ~ SVL, random = ~ ANIMAL, data = df_filtered_sub_log_A, ginverse = list(ANIMAL=Ainv), prior = prior1,
                        nitt=13000*10, thin=10*10, burnin=3000*10)
summary(MCMC_model8)
# Model with SVL and 
MCMC_model9 <- MCMCglmm(Tibia ~ SVL, random = ~ ANIMAL, data = df_filtered_sub_log_A, ginverse = list(ANIMAL=Ainv), prior = prior1,
                        nitt=13000*10, thin=10*10, burnin=3000*10)
summary(MCMC_model9)

# Filter out NAs





# Multiple response variable using cbind?
MCMC_multi_model <- MCMCglmm(cbind(Head.Length, Head.Width, Pectoral)~SVL, random=~ANIMAL,data = df_filtered_sub_log_A, ginverse = list(ANIMAL=Ainv), prior = prior1)

# Bad ggplot2
ggplot(df_filtered_sub_log_A, aes(SVL, Femur)) +
  geom_point(shape = 1, size = 1) +
  geom_smooth(linewidth=1) +
  geom_point(shape=5, size=5) +
  geom_smooth(linewidth=5) +
  geom_point(shape=10) +
  geom_smooth(shape=10, size=10) +
  geom_point(shape=12, size=20) +
  geom_smooth(linewidth=20)+
  geom_area() +
  geom_bin_2d() +
  geom_bin2d() +
  geom_blank() +
  geom_col() +
  geom_density_2d() +
  geom_function() +
  geom_hex() +
  theme_dark() +
  theme(axis.text.x = element_text(family="mono", size = 75, color="orange"),
        axis.text.y = element_text(family="Times New Roman", size = 150, color="green"),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0),
        legend.key.size = unit(4.5, "cm"),
        legend.text = element_text(family="Arial", size = 50,color="red"),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0)) 

###### San Salvador population -----
# Filter for Population = S and Select AnolisID, Mother, Father, Sex, SVL and other morphological traits
df_filtered_sub_log_S <- df_filtered_sub_log %>% dplyr::filter(Population == "S") %>%
  dplyr::select(c("AnolisID","Sex.x","SVL","Head.Length","Head.Width","Pectoral","Pelvic","Humerus",
                  "Radius","Ulna","Femur","Tibia","Mother","Father"))
View(df_filtered_sub_log_S)

# Change column names of AnolisID, Father and Mother 
colnames(df_filtered_sub_log_S)[which(names(df_filtered_sub_log_S) == "AnolisID")] <- "id"
colnames(df_filtered_sub_log_S)[which(names(df_filtered_sub_log_S) == "Father")] <- "sire"
colnames(df_filtered_sub_log_S)[which(names(df_filtered_sub_log_S) == "Mother")] <- "dam"
# make id, dam, sire factors
df_filtered_sub_log_S$id <- as.factor(df_filtered_sub_log_S$id)
df_filtered_sub_log_S$dam <- as.factor(df_filtered_sub_log_S$dam)
df_filtered_sub_log_S$sire <- as.factor(df_filtered_sub_log_S$sire)
# Need separate columns later to estimate addititive and maternal genetic effects
df_filtered_sub_log_S$ANIMAL <- as.factor(df_filtered_sub_log_S$id)
df_filtered_sub_log_S$MOTHER <- as.factor(df_filtered_sub_log_S$dam)
# Repeat for connection to the extra matrices in the models
df_filtered_sub_log_S$ANIMAL2 <- as.factor(df_filtered_sub_log_S$id)
df_filtered_sub_log_S$MOTHER2 <- as.factor(df_filtered_sub_log_S$dam)
# Need to ensure other values are correctly read in as numeric/factors
df_filtered_sub_log_S$SVL <- as.numeric(df_filtered_sub_log_S$SVL)
df_filtered_sub_log_S$Head.Length <- as.numeric(df_filtered_sub_log_S$Head.Length)
df_filtered_sub_log_S$Head.Width <- as.numeric(df_filtered_sub_log_S$Head.Width)
df_filtered_sub_log_S$Pectoral <- as.numeric(df_filtered_sub_log_S$Pectoral)
df_filtered_sub_log_S$Pelvic <- as.numeric(df_filtered_sub_log_S$Pelvic)
df_filtered_sub_log_S$Humerus <- as.numeric(df_filtered_sub_log_S$Humerus)
df_filtered_sub_log_S$Radius <- as.numeric(df_filtered_sub_log_S$Radius)
df_filtered_sub_log_S$Ulna <- as.numeric(df_filtered_sub_log_S$Ulna)
df_filtered_sub_log_S$Femur <- as.numeric(df_filtered_sub_log_S$Femur)
df_filtered_sub_log_S$Tibia <- as.numeric(df_filtered_sub_log_S$Tibia)

df_filtered_sub_log_S$Sex.x <- as.factor(df_filtered_sub_log_S$Sex.x) # Can be done for all the factors in the data
View(df_filtered_sub_log_S)

# Pedigree for the data
pedigreeS <- df_filtered_sub_log_S[,c("id","dam","sire")]
View(pedigreeS) # 1,378 rows
# Remove duplicate rows
pedigreeS_unique <- unique(pedigreeS)
View(pedigreeS_unique) # 504 rows

# Gets the inverse of the relatedness matrix (how is the relatedness matrix being created, G?)
# inverseA requires every parent (sire or dam) to also exist as a rwo in the pedigree
AinvS <- inverseA(pedigreeS_unique)$Ainv # Error - individuals appearing as dams but not in pedigree
# Check for missing dams in id
missing_dams <- setdiff(pedigreeS_unique$dam, pedigreeS_unique$id)
missing_dams # 110 dams not in id column
missing_sires <- setdiff(pedigreeS_unique$sire, pedigreeS_unique$id)
missing_sires # 52 sires not in id column
# Create new rows with missing sires and dams
missing_parents <- unique(c(missing_dams, missing_sires))
missing_parents
new_rows <- data.frame(id=missing_parents, sire=NA, dam=NA)
pedigreeS_unique_new <- rbind(pedigreeS_unique, new_rows)
View(pedigreeS_unique_new) # 666 rows
# Try running inverseA again
Ainv_S <- inverseA(pedigreeS_unique_new)$Ainv # Error - dams appearing before their offspring, try orderPed from MasterBayes
# Orders the pedigree so that parents come before their offspring
pedigreeS_unique_new_ord <- orderPed(pedigreeS_unique_new)
View(pedigreeS_unique_new_ord)
# Try running inverseA again
Ainv_obj_S <- inverseA(pedigreeS_unique_new_ord)
View(Ainv_obj)
Ainv_S <- Ainv_obj_S$Ainv # Use IDs returned by inverseA(); this fixes error below (line 174)
Ainv_S
# Basic model - includes only additive genetic variance and residual variance on the trait
# What is nu? What is V?
prior1 <- list(G = list(G1 = list(V = 1, nu = 0.002)),
               R = list(V = 1, nu = 0.002))

# Model with SVL and Head length 
MCMC_model1_S <- MCMCglmm(Head.Length ~ SVL, random = ~ ANIMAL, data = df_filtered_sub_log_S, ginverse = list(ANIMAL=Ainv_S), prior = prior1,
                          nitt=13000*10, thin=10*10, burnin=3000*10)
summary(MCMC_model1_S)
# Model with SVL and Head length 
MCMC_model2_S <- MCMCglmm(Head.Width ~ SVL, random = ~ ANIMAL, data = df_filtered_sub_log_S, ginverse = list(ANIMAL=Ainv_S), prior = prior1,
                          nitt=13000*10, thin=10*10, burnin=3000*10)
summary(MCMC_model2_S)
# Model with SVL and Pectoral
MCMC_model3_S <- MCMCglmm(Pectoral ~ SVL, random = ~ ANIMAL, data = df_filtered_sub_log_S, ginverse = list(ANIMAL=Ainv_S), prior = prior1,
                          nitt=13000*10, thin=10*10, burnin=3000*10)
summary(MCMC_model3_S)
# Model with SVL and Pelvic
MCMC_model4_S <- MCMCglmm(Pelvic ~ SVL, random = ~ ANIMAL, data = df_filtered_sub_log_S, ginverse = list(ANIMAL=Ainv_S), prior = prior1,
                          nitt=13000*10, thin=10*10, burnin=3000*10)
summary(MCMC_model4_S)
# Model with SVL and Humerus
MCMC_model5_S <- MCMCglmm(Humerus ~ SVL, random = ~ ANIMAL, data = df_filtered_sub_log_S, ginverse = list(ANIMAL=Ainv_S), prior = prior1,
                          nitt=13000*10, thin=10*10, burnin=3000*10)
summary(MCMC_model5_S)
# Model with SVL and Radius
MCMC_model6_S <- MCMCglmm(Radius ~ SVL, random = ~ ANIMAL, data = df_filtered_sub_log_S, ginverse = list(ANIMAL=Ainv_S), prior = prior1,
                          nitt=13000*10, thin=10*10, burnin=3000*10)
summary(MCMC_model6_S)
# Model with SVL and Ulna
MCMC_model7_S <- MCMCglmm(Ulna ~ SVL, random = ~ ANIMAL, data = df_filtered_sub_log_S, ginverse = list(ANIMAL=Ainv_S), prior = prior1,
                          nitt=13000*10, thin=10*10, burnin=3000*10)
summary(MCMC_model7_S)
# Model with SVL and Femur
MCMC_model8_S <- MCMCglmm(Femur ~ SVL, random = ~ ANIMAL, data = df_filtered_sub_log_S, ginverse = list(ANIMAL=Ainv_S), prior = prior1,
                          nitt=13000*10, thin=10*10, burnin=3000*10)
summary(MCMC_model8_S)
# Model with SVL and 
MCMC_model9_S <- MCMCglmm(Tibia ~ SVL, random = ~ ANIMAL, data = df_filtered_sub_log_S, ginverse = list(ANIMAL=Ainv_S), prior = prior1,
                          nitt=13000*10, thin=10*10, burnin=3000*10)
summary(MCMC_model9_S)













