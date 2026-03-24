################# Packages ##############
library(dplyr)
library(ggplot2)
library(tidyr)
library(performance)
library(see)
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

# Add log_water_loss column
#norris_dat_filtered_new_log <- norris_dat_filtered_new %>% mutate(log_water_loss = log(water_loss_g))
#View(norris_dat_filtered_new_log)

# Remove NAs
#norris_dat_filtered_new_log <- na.omit(norris_dat_filtered_new_log)

######### Effect of hatchling mass on maximum sprint speed for open parental source island
norris_dat_filtered_new_sunparent <- norris_dat_filtered_new %>% dplyr::filter(parental_island == "Open")

l1 <- lm(hatchling_max_speed~hatch_svl, data = norris_dat_filtered_new_sunparent)
lm(hatchling_max_speed~hatch_svl, data = norris_dat_filtered_new_sunparent)
summary(l1)
plot(l1)# Residuals and standardized residuals are very close to 0. 
check_model(l1) 


l1$residuals
hist(l1$residuals) # probably normal
shapiro.test(l1$residuals) # nope p-value = 0.0001
hist(norris_dat_filtered_new_sunparent$hatchling_max_speed)

##ggplot2##
library(ggplot2)
r=ggplot(data=norris_dat_filtered_new_sunparent, aes(x=hatch_svl, y=hatchling_max_speed))+ 
  geom_point()+
  stat_smooth(method = "lm")+
  theme_bw() + 
  theme(axis.title=element_text(size=20),axis.text=element_text(size=10),panel.grid = element_blank(), axis.line=element_line(),legend.position="top",legend.title=element_blank())
print(r)








