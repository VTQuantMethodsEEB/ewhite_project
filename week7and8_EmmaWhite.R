################# Packages ##############
library(dplyr)
library(ggplot2)
library(tidyr)
library(performance)
library(see)
library(effects)
library(emmeans)
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

######### H1: Effect of hatchling mass on maximum sprint speed for open parental source island
norris_dat_filtered_new_openparent <- norris_dat_filtered_new %>% dplyr::filter(parental_island == "Open")

l1 <- lm(hatchling_max_speed~hatch_svl, data = norris_dat_filtered_new_openparent)
lm(hatchling_max_speed~hatch_svl, data = norris_dat_filtered_new_openparent)
summary(l1) # Significant positive relationship between hatchling svl and hatchling max speed (p-value = 0.002)

par(mfrow = c(2,2)) 
plot(l1)# Residuals vs fitted plot show values are close to the 0 line (maybe not much heteroscedasticity?). Residuals at both tails seem to not be close the line (maybe not normal)
# Leverage plot doesn't show any datapoints that could be really influential.
check_model(l1) 

l1$residuals
hist(l1$residuals) # looks close to normal
shapiro.test(l1$residuals) # nope p-value = 0.0001
hist(norris_dat_filtered_new_openparent$hatchling_max_speed)

##ggplot2##
library(ggplot2)
r=ggplot(data=norris_dat_filtered_new_openparent, aes(x=hatch_svl, y=hatchling_max_speed))+ 
  geom_point()+
  stat_smooth(method = "lm")+
  theme_bw() + 
  theme(axis.title=element_text(size=20),axis.text=element_text(size=10),panel.grid = element_blank(), axis.line=element_line(),legend.position="top",legend.title=element_blank())
print(r)


######## H2: What is the effect of treatment and parent island on hatchling max speed 
norris_dat_filtered_new_sub <- norris_dat_filtered_new %>% dplyr::select("hatchling_max_speed", "treatment", "parental_island")

# Create additive model
l2 <- lm(hatchling_max_speed~treatment+parental_island, data=norris_dat_filtered_new_sub)
summary(l2) # no statistically clear effect of treatment or parental island on hatchling speed (p-values > 0.05)
plot(allEffects(l2))
emmeans(l2, specs = ~treatment+parental_island)
# all contrasts
all_comps <- emmeans(l2, pairwise ~ treatment + parental_island)
all_comps$contrasts

# plotting the additive model
norris_dat_filtered_new_sub$yhat <- predict(l2)
pp <- with(norris_dat_filtered_new_sub, expand.grid(treatment = unique(treatment),
                                                    parental_island = unique(parental_island)))
pp$hatchling_max_speed <- predict(l2,newdata = pp)

ggplot(pp, aes(x=treatment,y=hatchling_max_speed,colour=parental_island))+
  geom_point(color="green")+
  geom_line(aes(group=parental_island))+
  geom_point(data=norris_dat_filtered_new, aes(x=treatment,y=hatchling_max_speed,colour=parental_island))



# Create interactive model
norris_dat_filtered_new_sub <- norris_dat_filtered_new %>% dplyr::select("hatchling_max_speed", "treatment", "parental_island")

l3 <- lm(hatchling_max_speed~treatment*parental_island, data=norris_dat_filtered_new_sub)
summary(l3) # no statistically clear effects of treatment, parental island, or the interaction
plot(allEffects(l3))
emmeans(l3, specs = ~treatment*parental_island)
# all contrasts
all_comps3 <- emmeans(l3, pairwise ~ treatment*parental_island)
all_comps3$contrasts

# plotting the interactive model
norris_dat_filtered_new_sub$yhat <- predict(l3)
pp <- with(norris_dat_filtered_new_sub, expand.grid(treatment = unique(treatment),
                                            parental_island = unique(parental_island)))
pp$hatchling_max_speed <- predict(l3,newdata = pp)
pp
ggplot(pp, aes(x=treatment,y=hatchling_max_speed,colour=parental_island))+
  geom_point(colour="green")+
  geom_line(aes(group=parental_island))+
  geom_point(data=norris_dat_filtered_new, aes(x=treatment,y=hatchling_max_speed,colour=parental_island))

#what does all this tell you?

######## H3: What is the effect of svl and mass on hatchling max speed 
norris_dat_filtered_new_openparent <- norris_dat_filtered_new %>% dplyr::filter(parental_island == "Open")

# Create additive model
l4 <- lm(hatchling_max_speed~hatch_svl+hatch_mass, data=norris_dat_filtered_new_openparent)
summary(l4) # statistically clear effect of mass on on hatchling speed (p-values < 0.05)
plot(allEffects(l4))

new.dat.combos <- with(norris_dat_filtered_new_openparent,
                       expand.grid(hatch_mass = seq(min(hatch_mass, na.rm=T), max(hatch_mass, na.rm=T), by=1),
                                   hatch_svl = seq(min(hatch_svl, na.rm=T), max(hatch_svl, na.rm=T), by=1)))

new.dat.combos$hatchling_max_speed <- predict(l4,newdata=new.dat.combos)

ggplot(new.dat.combos, aes(x=hatch_svl, y=hatchling_max_speed,colour=hatch_mass))+
  geom_line(aes(group=hatch_mass))+
  geom_point(data=norris_dat_filtered_new_openparent, aes(x=hatch_svl, y=hatchling_max_speed,colour=hatch_mass))

# Create interactive model
l5 <- lm(hatchling_max_speed~hatch_svl*hatch_mass, data=norris_dat_filtered_new_openparent)
summary(l5) # no statistically clear effect of mass, svl or interaction on on hatchling speed (p-values > 0.05)
plot(allEffects(l5))

new.dat.combos <- with(norris_dat_filtered_new_openparent,
                       expand.grid(hatch_mass = seq(min(hatch_mass, na.rm=T), max(hatch_mass, na.rm=T), by=1),
                                   hatch_svl = seq(min(hatch_svl, na.rm=T), max(hatch_svl, na.rm=T), by=1)))

new.dat.combos$hatchling_max_speed <- predict(l5,newdata=new.dat.combos)

ggplot(new.dat.combos, aes(x=hatch_svl, y=hatchling_max_speed,colour=hatch_mass))+
  geom_line(aes(group=hatch_mass))+
  geom_point(data=norris_dat_filtered_new_openparent, aes(x=hatch_svl, y=hatchling_max_speed,colour=hatch_mass))



