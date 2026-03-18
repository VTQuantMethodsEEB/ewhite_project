##################### Packages #####################
library(dplyr)
library(ggplot2)
library(egg)
library(tidyr)
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
norris_dat_filtered_new_log <- norris_dat_filtered_new %>% mutate(log_water_loss = log(water_loss_g))
View(norris_dat_filtered_new_log)

# Remove NAs
norris_dat_filtered_new_log <- na.omit(norris_dat_filtered_new_log)


################### Permutation test ######################
# Question: Do hatchlings from each environment (open or closed) show significant differences in water loss between incubation treatments?
# Select Open source island
norris_dat_filtered_open <- norris_dat_filtered_new_log %>% dplyr::filter(parental_island == "Open")
View(norris_dat_filtered_open)
# Select treatment and log water loss column
norris_dat_filtered_open_sub <- norris_dat_filtered_open %>% dplyr::select(c("treatment", "log_water_loss"))
View(norris_dat_filtered_open_sub)

set.seed(15)
res <- NA
 
for (i in 1:10000) {
  colonyboot <- sample(c(norris_dat_filtered_open_sub$treatment == "Open", norris_dat_filtered_open_sub$treatment == "Shade")) ## scramble
  ## pick out open & shade samples
  openboot <- colonyboot[1:length(norris_dat_filtered_open_sub$treatment[norris_dat_filtered_open_sub$treatment=="Open"])] #this says assign the first boot samples to open treatment
  shadeboot <- colonyboot[(length(norris_dat_filtered_open_sub$treatment[norris_dat_filtered_open_sub$treatment=="Open"])+1):length(norris_dat_filtered_open_sub$treatment)] #this says assign the rest of the observations to shade treatment
  
  ## compute & store difference in means
  res[i] <- mean(openboot)-mean(shadeboot) #calculate the difference in the open means and the shade treatment means
  #[i] says "where i", and i is a counter, after running this loop, i should be 1000
}
res

#what is our observed mean difference?
obs <- mean(norris_dat_filtered_open_sub$treatment=="Open")-mean(norris_dat_filtered_open_sub$treatment=="Shade")
obs

hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")

##so how do we get our p-value?
res[res>=obs]
length(res[res>=obs])
6682/10000 # p-value = 0.6682 
mean(res>=obs)
# No statistically clear difference in water loss between treatments for open island hatchlings


#####Shapiro-Wilk Test#######
#Are our data normally distributed?##
#The null hypothesis is that the data are normally distributed
#P<0.05 indicates NOT normal

# Filter for open treatment
norris_dat_filtered_open_sub_open <- norris_dat_filtered_open_sub %>% dplyr::filter(treatment == "Open")
View(norris_dat_filtered_open_sub_open)

shapiro.test.open.open.waterloss <- shapiro.test(norris_dat_filtered_open_sub_open$log_water_loss)
shapiro.test.open.open.waterloss # p-value = 0.160. normal


# Filter for shade treatment
norris_dat_filtered_open_sub_shade <- norris_dat_filtered_open_sub %>% dplyr::filter(treatment == "Shade")
View(norris_dat_filtered_open_sub_shade)

shapiro.test.open.shade.waterloss <- shapiro.test(norris_dat_filtered_open_sub_shade$log_water_loss)
shapiro.test.open.shade.waterloss # p-value = 1.28e-5. not normal

####### welch's t-test ######
tt.openisland.waterloss <- t.test(log_water_loss~treatment,data=norris_dat_filtered_open_sub)
tt.openisland.waterloss # p-value = 0.93. no statistically clear difference in log water loss between treatments for open island hatchlings

### Hypothesis 2: Does maximum sprint speed differ between treatments for shade island source hatchlings?
# Subset for shade parental island
norris_dat_filtered_new_shade <- norris_dat_filtered_new %>% dplyr::filter(parental_island == "Shade")
# Filter for shade treatment
norris_dat_filtered_new_shade_shade <- norris_dat_filtered_new %>% dplyr::filter(treatment == "Shade")
# Shapiro test
shapiro.test.shade.shade.maxspeed <- shapiro.test(norris_dat_filtered_new_shade_shade$hatchling_max_speed)
shapiro.test.shade.shade.maxspeed # data is not normally distributed 

# Welch's t-test 
tt.shadeisland.maxspeed <- t.test(hatchling_max_speed~treatment, data=norris_dat_filtered_new_shade)
tt.shadeisland.maxspeed # No statistically clear difference in max speed among hatchlings incubated across treatments
