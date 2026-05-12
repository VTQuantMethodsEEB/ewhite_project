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
library(lme4)
library(glmmTMB)

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
View(norris_dat_filtered_new)
# Filter for ppen island source individuals -
norris_dat_filtered_new_openparent <- norris_dat_filtered_new %>% dplyr::filter(parental_island == "Open")
# General linear mixed effects model (with lmer)
lm1 <- lmer(hatchling_max_speed~treatment + hatch_svl + (1|cage), data=norris_dat_filtered_new_openparent) 
summary(lm1)
#KL - singularity warning here - better to try glmmTMB
# Calculate p-values
# Number of individuals + number of cages + (number of treatments + continuous variable)
180 - 33 - (2+1) # 144
# Treatment
t.value <- 0.176 # absolute value
p.value <- 2*pt(t.value, df = 144, lower = FALSE)
p.value # 0.860; no statistically clear effect
# Body size
t.value <- 3.055
p.value <- 2*pt(t.value, df = 144, lower = FALSE)
p.value # 0.002 positive effect of hatchling svl on hatchling max speed.
















