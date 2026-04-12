################# Packages ##############
library(dplyr)
library(ggplot2)
library(tidyr)
library(performance)
library(see)
library(effects)
library(emmeans)
library(DHARMa)
library(AICcmodavg)
######## Norris 2024 dataset
norris_dat <- read.csv("Norris2024_master_data.csv")
### Filter out individuals that died before or during experiments
norris_dat_filtered <- norris_dat %>% dplyr::filter(!(egg_id %in% c(32, 65, 107, 144, 162, 171, 225, 234,
                                                                    258, 297, 298, 307, 345, 359, 361, 365, 370,
                                                                    390, 402, 412, 416, 420, 429)))
# Remove the comments and notes columns
norris_dat_filtered_new <- norris_dat_filtered[-c(12,13)]
# Remove NAs from the entire dataset
norris_dat_filtered_new <- na.omit(norris_dat_filtered_new)

# What is the effect of incubation treatment, sex, and egg mass on hatchling svl (Body size)?
# Filter for ppen island source individuals -
norris_dat_filtered_new_openparent <- norris_dat_filtered_new %>% dplyr::filter(parental_island == "Open")
# Run Gaussian GLM
glm1 <- glm(hatch_svl ~ treatment+sex+egg_mass, norris_dat_filtered_new_openparent, family = "gaussian")
summary(glm1) # Treatment and Egg mass significantly effect hatchling body size. There is no statistically
# clear effect of sex on SVL
plot(allEffects(glm1)) # Individuals incubated in the cool, wet (shade) treatment had larger body size compared to 
# individuals incubated in the hot, dry treatment (open). Individuals with greater egg mass tended to have 
# larger body size across sex.
simulationOutput1 <- simulateResiduals(fittedModel = glm1, plot = T)
# Residuals closely aligned along the diagnal and seem to fit a normal distribution well. No significant dispersion. No outliers. 
# ggplot
dat.new <- expand.grid(egg_mass=seq(from = min(norris_dat_filtered_new_openparent$egg_mass),
                             to = max(norris_dat_filtered_new_openparent$egg_mass),length.out = 180),
                    treatment = unique(norris_dat_filtered_new_openparent$treatment),
                    sex = unique(norris_dat_filtered_new_openparent$sex))

dat.new$yhat <- predict(glm1,type="response",newdata = dat.new)
norris_dat_filtered_new_openparent$yhat2 <- predict(glm1,type="response")
plot1 <- ggplot(data=norris_dat_filtered_new_openparent,aes(x=treatment,y=hatch_svl,color=sex))+
  geom_point(size=2,shape =1) +
  geom_line(data=dat.new, aes(x=treatment,y=yhat,col = sex))
#geom_line(aes(x=date,y=yhat2,col = species))
plot1

# Run Gamma GLM
glm2 <- glm(hatch_svl ~ treatment+sex+egg_mass, norris_dat_filtered_new_openparent, family = "Gamma")
summary(glm2)
plot(allEffects(glm2)) # Similar results as the gaussian model. Incubation treatment and egg mass
# significantly affect hatchling body size but no statistically clear effect of sex on body size.
simulationOutput2 <- simulateResiduals(fittedModel = glm2, plot = T)
# Plotted residuals are very close along the diagnal. Only 1 outlier.
# ggplot
dat.new2 <- expand.grid(egg_mass=seq(from = min(norris_dat_filtered_new_openparent$egg_mass),
                                    to = max(norris_dat_filtered_new_openparent$egg_mass),length.out = 180),
                       treatment = unique(norris_dat_filtered_new_openparent$treatment),
                       sex = unique(norris_dat_filtered_new_openparent$sex))

dat.new2$yhat <- predict(glm2,type="response",newdata = dat.new2)
norris_dat_filtered_new_openparent$yhat2 <- predict(glm2,type="response")
plot2 <- ggplot(data=norris_dat_filtered_new_openparent,aes(x=treatment,y=hatch_svl,color=sex))+
  geom_point(size=2,shape =1) +
  geom_line(data=dat.new2, aes(x=treatment,y=yhat,col = sex))
plot2


#### Week 11 - Model Comparison
# Open island source individuals -
#### Model comparison - LIkelihood ratio test
lm1 <- lm(hatch_svl ~ 1, norris_dat_filtered_new_openparent)
lm2 <- lm(hatch_svl ~ treatment, norris_dat_filtered_new_openparent)
lm3 <- lm(hatch_svl ~ treatment + sex, norris_dat_filtered_new_openparent)
lm4 <- lm(hatch_svl ~ treatment + sex + egg_mass, norris_dat_filtered_new_openparent)

anova(lm1, lm2, lm3, lm4) # the best model is lm4?

#### Model comparison - AIC
AIC(lm1, lm2, lm3, lm4) # the best model is lm4 (lowest AIC = 491)
aictab(cand.set=list(lm1,lm2,lm3,lm4),modnames=c("lm1","lm2","lm3","lm4"))#AIC table









