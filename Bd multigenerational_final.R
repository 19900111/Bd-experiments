library(tidyverse)
library(lme4)
library(jagsUI)
library(dplyr)
library(gridExtra)
library(drc)

############################# First transgenerational test (carry-over effects)
##### Scenario: exposing the first generation only and returning the second and third generations to the metal free environment (culture)

setwd("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions")

# transgenerational data

tg_dat <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Bd transgeneration test1.csv")

# make a subset of data for "Cu"

cu_dat <- filter(tg_dat, Metal == "Cu")

# calculate the colony fecundity (Cf = Zt/C)

cu_dat <- cu_dat %>%
  mutate(Colony_fecundity = Total_zoospores / Colony)

# calculate the sporangium fecundity (Sf = Zt/St)  

cu_dat <- cu_dat %>%
  mutate(Sporangium_fecundity = Total_zoospores / Sporangium)

# calculate the total number of zoospores (Zt = St * Sf) 

cu_dat <- cu_dat %>%
  mutate(Total_zoospores = Sporangium * Sporangium_fecundity)

# calculate density of sporangium (S = C * Sf) 

cu_dat <- cu_dat %>%
  mutate(Sporangium  = Colony * Sporangium_fecundity)


Cu_log10_3g <- write.csv(cu_dat,"C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Cu_log10_3g.csv", row.names = F)

Cu_log10_3g <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Cu_log10_3g.csv")

# Convert the values to the logartimic scale


Cu_log10_3g

View(Cu_log10_3g)

head(Cu_log10_3g)

## Plot the data with fitted regression lines 

x_breaks <- c(0, 0.2, 0.4, 1, 2, 4, 6)

# colony
colony <- ggplot(Cu_log10_3g, aes(x = Concentration, y = Colony, color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+
  labs(x = "Toxic unit", y = "Colony (scale log10)", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

colony


# colony fecundity

colony_f <- ggplot(Cu_log10_3g, aes(x = Concentration, y = Colony_fecundity , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+ 
  labs(x = "Toxic unit", y = "Colony fecundity (scale log10)", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

colony_f

# sporangium fecundity

sporangium_f <- ggplot(Cu_log10_3g, aes(x = Concentration, y = Sporangium_fecundity , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) +  scale_y_log10()+
  labs(x = "Toxic unit", y = "Sporangium fecundity (scale log10)", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

sporangium_f


# total zoospores 

total_zoospores <- ggplot(Cu_log10_3g, aes(x = Concentration, y = Total_zoospores , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+
  labs(x = "Toxic unit", y = "Total zoospores/mL(scale log10)", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

total_zoospores 

# sporangium

sporangium <- ggplot(Cu_log10_3g, aes(x = Concentration, y = Sporangium , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+
  labs(x = "Toxic unit", y = "Sporangium/mL(scale log10)", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

sporangium 


### Aggregate the plots

cu_all_plots <- grid.arrange(colony, colony_f, sporangium_f, total_zoospores, sporangium,   ncol = 2, nrow=3) 

ggsave("cu_all_plots.png", cu_all_plots, width = 10, height = 6, units = "in")


# make a subset of data for "Zn"

zn_dat <- filter(tg_dat, Metal == "Zn")

# calculate the colony fecundity (Cf = Zt/C)

zn_dat <- zn_dat %>%
  mutate(Colony_fecundity = Total_zoospores / Colony)

# calculate the sporangium fecundity (Sf = Zt/St)  

zn_dat <- zn_dat %>%
  mutate(Sporangium_fecundity = Total_zoospores / Sporangium)

# calculate the total number of zoospores (Zt = St * Sf) 

zn_dat <- zn_dat %>%
  mutate(Total_zoospores = Sporangium * Sporangium_fecundity)

# calculate density of sporangium (S = C * Sf) 

zn_dat <- zn_dat %>%
  mutate(Sporangium  = Colony * Sporangium_fecundity)


Zn_log10_3g <- write.csv(zn_dat,"C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Zn_log10_3g.csv", row.names = F)

Zn_log10_3g <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Zn_log10_3g.csv")

Zn_log10_3g

View(Zn_log10_3g)

head(Zn_log10_3g)

## Plot the data with fitted regression lines 

x_breaks <- c(0, 0.02, 0.2, 2, 4, 6, 8)

# colony
colony <- ggplot(Zn_log10_3g, aes(x = Concentration, y = Colony, color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+ 
  labs(x = "Toxic unit", y = "Colony (scale log10)", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

colony


# colony fecundity

colony_f <- ggplot(Zn_log10_3g, aes(x = Concentration, y = Colony_fecundity , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) +  scale_y_log10()+ 
  labs(x = "Toxic unit", y = "Colony fecundity (scale log10)", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

colony_f

# sporangium fecundity

sporangium_f <- ggplot(Zn_log10_3g, aes(x = Concentration, y = Sporangium_fecundity , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) +  scale_y_log10()+
  labs(x = "Toxic unit", y = "Sporangium fecundity (scale log10)", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

sporangium_f


# total zoospores 

total_zoospores <- ggplot(Zn_log10_3g, aes(x = Concentration, y = Total_zoospores , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+
  labs(x = "Toxic unit", y = "Total zoospores/mL(scale log10)", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

total_zoospores 

# sporangium

sporangium <- ggplot(Zn_log10_3g, aes(x = Concentration, y = Sporangium , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+
  labs(x = "Toxic unit", y = "Sporangium/mL(scale log10)", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

sporangium 


### Aggregate the plots

zn_all_plots <- grid.arrange(colony, colony_f, sporangium_f, total_zoospores, sporangium,   ncol = 2, nrow=3) 

ggsave("zn_all_plots.png", zn_all_plots, width = 10, height = 6, units = "in")


# make a subset of data for "Mixed"

mixed_dat <- filter(tg_dat, Metal == "Mixed")

# calculate the colony fecundity (Cf = Zt/C)

mixed_dat <- mixed_dat %>%
  mutate(Colony_fecundity = Total_zoospores / Colony)

# calculate the sporangium fecundity (Sf = Zt/St)  

mixed_dat <- mixed_dat %>%
  mutate(Sporangium_fecundity = Total_zoospores / Sporangium)

# calculate the total number of zoospores (Zt = St * Sf) 

mixed_dat <- mixed_dat %>%
  mutate(Total_zoospores = Sporangium * Sporangium_fecundity)

# calculate density of sporangium (S = C * Sf) 

mixed_dat <- mixed_dat %>%
  mutate(Sporangium  = Colony * Sporangium_fecundity)


Mixed_log10_3g <- write.csv(mixed_dat,"C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Mixed_log10_3g.csv", row.names = F)

Mixed_log10_3g <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Mixed_log10_3g.csv")

Mixed_log10_3g

View(Mixed_log10_3g)

head(Mixed_log10_3g)

## Plot the data with fitted regression lines 

x_breaks <- c(0, 0.22, 0.6, 3, 6, 10, 14)

# colony
colony <- ggplot(Mixed_log10_3g, aes(x = Concentration, y = Colony, color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+ 
  labs(x = "Total toxic unit", y = "Colony (scale log10)", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

colony


# colony fecundity

colony_f <- ggplot(Mixed_log10_3g, aes(x = Concentration, y = Colony_fecundity , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+ 
  labs(x = "Total toxic unit", y = "Colony fecundity (scale log10)", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

colony_f

# sporangium fecundity

sporangium_f <- ggplot(Mixed_log10_3g, aes(x = Concentration, y = Sporangium_fecundity , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) +  scale_y_log10()+ 
  labs(x = "Total toxic unit", y = "Sporangium fecundity (scale log10)", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

sporangium_f


# total zoospores 

total_zoospores <- ggplot(Mixed_log10_3g, aes(x = Concentration, y = Total_zoospores , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+
  labs(x = "Total toxic unit", y = "Total zoospores/mL(scale log10)", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

total_zoospores 

# sporangium

sporangium <- ggplot(Mixed_log10_3g, aes(x = Concentration, y = Sporangium , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+
  labs(x = "Total toxic unit", y = "Sporangium/mL(scale log10)", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

sporangium 


### Aggregate the plots

mixed_all_plots <- grid.arrange(colony, colony_f, sporangium_f, total_zoospores, sporangium,   ncol = 2, nrow=3) 

ggsave("mixed_all_plots.png", mixed_all_plots, width = 10, height = 6, units = "in")



############################# Second transgenerational test
##### Scenario 1: exposing the first generation to C,L,M, and H concentrations and cross exposing to the same concentrations in the second generation 

setwd("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions")

# transgenerational data

tg2_dat <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Bd transgeneration test2.csv")

head(tg2_dat)

################################ Effects of Cu

# make a subset of data for "Cu"

cu2_dat <- filter(tg2_dat, F1_metal == "Cu")

# calculate the colony fecundity (Cf = Zt/C)

view(cu2_dat)

cu2_dat <- cu2_dat %>%
  mutate(Colony_fecundity = Total_zoospores / Colony)

view(cu2_dat)

# calculate the sporangium fecundity (Sf = Zt/St)  

cu2_dat <- cu2_dat %>%
  mutate(Sporangium_fecundity = Total_zoospores / Sporangium)


Cu2_log10_3g <- write.csv(cu2_dat,"C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Cu2_log10_3g.csv", row.names = F)

Cu2_log10_3g <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Cu2_log10_3g.csv")


Cu2_log10_3g

View(Cu2_log10_3g)

head(Cu2_log10_3g)

## Plot the data with fitted regression lines 

x_breaks <- c(0, 0.2, 2, 6)

# colony
colony <- ggplot(Cu2_log10_3g, aes(x = F2_conc, y = Colony, color = factor(F1_conc))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+
  labs(x = "Toxic unit (F1)", y = "Colony (scale log10)", color = "Toxic unit (F2)") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

colony

# colony fecundity

colony_f <- ggplot(Cu2_log10_3g, aes(x = F2_conc, y = Colony_fecundity , color = factor(F1_conc))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) +  scale_y_log10()+
  labs(x = "Toxic unit (F1)", y = "Colony fecundity (scale log10)", color = "Toxic unit (F2)") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

colony_f

# sporangium fecundity

sporangium_f <- ggplot(Cu2_log10_3g, aes(x = F2_conc, y = Sporangium_fecundity , color = factor(F1_conc))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) +  scale_y_log10()+
  labs(x = "Toxic unit (F1)", y = "Sporangium fecundity (scale log10)", color = "Toxic unit (F2)") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

sporangium_f

view(Cu2_log10_3g)
# survived zoospores 

survived_zoospores <- ggplot(Cu2_log10_3g, aes(x = F2_conc, y = Survived_zoospores , color = factor(F1_conc))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+
  labs(x = "Toxic unit (F1)", y = "Survived zoospores/mL (scale log10)", color = "Toxic unit (F2)") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

survived_zoospores 

# sporangium

sporangium <- ggplot(Cu2_log10_3g, aes(x = F2_conc, y = Sporangium , color = factor(F1_conc))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+
  labs(x = "Toxic unit (F1)", y = "Sporangium/mL (scale log10)", color = "Toxic unit (F2)") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

sporangium 

view(Cu2_log10_3g)

### Aggregate the plots

cu2_all_plots <- grid.arrange(colony, colony_f, sporangium, sporangium_f, survived_zoospores,   ncol = 2, nrow=3) 

ggsave("cu2_all_plots.png", cu2_all_plots, width = 10, height = 6, units = "in")


################################ Effects of Zn

# make a subset of data for "Zn"

zn2_dat <- filter(tg2_dat, F1_metal == "Zn")

# calculate the colony fecundity (Cf = Zt/C)

view(zn2_dat)

zn2_dat <- zn2_dat %>%
  mutate(Colony_fecundity = Total_zoospores / Colony)

view(zn2_dat)

# calculate the sporangium fecundity (Sf = Zt/St)  

zn2_dat <- zn2_dat %>%
  mutate(Sporangium_fecundity = Total_zoospores / Sporangium)


Zn2_log10_3g <- write.csv(zn2_dat,"C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Zn2_log10_3g.csv", row.names = F)

Zn2_log10_3g <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Zn2_log10_3g.csv")


Zn2_log10_3g

View(Zn2_log10_3g)

head(Zn2_log10_3g)

## Plot the data with fitted regression lines 

x_breaks <- c(0, 0.02, 4, 8)

# colony
colony <- ggplot(Zn2_log10_3g, aes(x = F2_conc, y = Colony, color = factor(F1_conc))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10() +
  labs(x = "Toxic unit (F1)", y = "Colony (scale log10)", color = "Toxic unit (F2)") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

colony

view(Zn2_log10_3g)

# colony fecundity

colony_f <- ggplot(Zn2_log10_3g, aes(x = F2_conc, y = Colony_fecundity , color = factor(F1_conc))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+
  labs(x = "Toxic unit (F1)", y = "Colony fecundity (scale log10)", color = "Toxic unit (F2)") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

colony_f

# sporangium fecundity

sporangium_f <- ggplot(Zn2_log10_3g, aes(x = F2_conc, y = Sporangium_fecundity , color = factor(F1_conc))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+ 
  labs(x = "Toxic unit (F1)", y = "Sporangium fecundity (scale log10)", color = "Toxic unit (F2)") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

sporangium_f


# survived zoospores 

survived_zoospores <- ggplot(Zn2_log10_3g, aes(x = F2_conc, y = Survived_zoospores , color = factor(F1_conc))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+
  labs(x = "Toxic unit (F1)", y = "Survived zoospores/mL (scale log10)", color = "Toxic unit (F2)") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

survived_zoospores 

# sporangium

sporangium <- ggplot(Zn2_log10_3g, aes(x = F2_conc, y = Sporangium , color = factor(F1_conc))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+
  labs(x = "Toxic unit (F1)", y = "Sporangium/mL (scale log10)", color = "Toxic unit (F2)") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

sporangium 

view(Zn2_log10_3g)

### Aggregate the plots

Zn2_all_plots <- grid.arrange(colony, colony_f, sporangium, sporangium_f, survived_zoospores,   ncol = 2, nrow=3) 

ggsave("Zn2_all_plots.png", Zn2_all_plots, width = 10, height = 6, units = "in")



################################ Combination of Cu and Zn


# make a subset of data for "Mixed"

mixed2_dat <- filter(tg2_dat, F1_metal == "Mixed")

# calculate the colony fecundity (Cf = Zt/C)

view(mixed2_dat)

mixed2_dat <- mixed2_dat %>%
  mutate(Colony_fecundity = Total_zoospores / Colony)

view(mixed2_dat)

# calculate the sporangium fecundity (Sf = Zt/St)  

mixed2_dat <- mixed2_dat %>%
  mutate(Sporangium_fecundity = Total_zoospores / Sporangium)


Mixed2_log10_3g <- write.csv(mixed2_dat,"C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Mixed2_log10_3g.csv", row.names = F)

Mixed2_log10_3g <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Mixed2_log10_3g.csv")


Mixed2_log10_3g

View(Mixed2_log10_3g)

head(Mixed2_log10_3g)

## Plot the data with fitted regression lines 

x_breaks <- c(0, 0.22, 6, 14)

# colony
colony <- ggplot(Mixed2_log10_3g, aes(x = F2_conc, y = Colony, color = factor(F1_conc))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+  
  labs(x = "Total toxic unit (F1)", y = "Colony (scale log10)", color = "Toxic unit (F2)") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

colony

view(Mixed2_log10_3g)

# colony fecundity

colony_f <- ggplot(Mixed2_log10_3g, aes(x = F2_conc, y = Colony_fecundity , color = factor(F1_conc))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+ 
  labs(x = "Total toxic unit (F1)", y = "Colony fecundity (scale log10)", color = "Toxic unit (F2)") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

colony_f

# sporangium fecundity

sporangium_f <- ggplot(Mixed2_log10_3g, aes(x = F2_conc, y = Sporangium_fecundity , color = factor(F1_conc))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+ 
  labs(x = "Total toxic unit (F1)", y = "Sporangium fecundity (scale log10)", color = "Toxic unit (F2)") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

sporangium_f


# survived zoospores 

survived_zoospores <- ggplot(Mixed2_log10_3g, aes(x = F2_conc, y = Survived_zoospores , color = factor(F1_conc))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+ 
  labs(x = "Total toxic unit (F1)", y = "Survived zoospores/mL (scale log10)", color = "Toxic unit (F2)") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

survived_zoospores 

# sporangium

sporangium <- ggplot(Mixed2_log10_3g, aes(x = F2_conc, y = Sporangium , color = factor(F1_conc))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + scale_y_log10()+ 
  labs(x = "Total toxic unit (F1)", y = "Sporangium/mL (scale log10)", color = "Toxic unit (F2)") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

sporangium 

view(Mixed2_log10_3g)

### Aggregate the plots

Mixed2_all_plots <- grid.arrange(colony, colony_f, sporangium, sporangium_f, survived_zoospores,   ncol = 2, nrow=3) 

ggsave("Mixed2_all_plots.png", Mixed2_all_plots, width = 10, height = 6, units = "in")


############################# Third transgenerational test
##### Scenario 2: continue exposing Bd to C,L,M, and H concentrations in three generations 

setwd("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions")

# transgenerational data

tg3_dat <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Bd transgeneration test3.csv")

head(tg3_dat)


################################ Effects of Cu


# make a subset of data for "Cu"

cu3_dat <- filter(tg3_dat, Metal == "Cu")

# calculate the colony fecundity (Cf = Zt/C)

view(cu3_dat)

cu3_dat <- cu3_dat %>%
  mutate(Colony_fecundity = Total_zoospores / Colony)

view(cu3_dat)

# calculate the sporangium fecundity (Sf = Zt/St)  

cu3_dat <- cu3_dat %>%
  mutate(Sporangium_fecundity = Total_zoospores / Sporangium)

view(cu3_dat)

Cu3_log10_3g <- write.csv(cu3_dat,"C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Cu3_log10_3g.csv", row.names = F)

Cu3_log10_3g <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Cu3_log10_3g.csv")


Cu3_log10_3g

View(Cu3_log10_3g)

head(Cu3_log10_3g)


## Plot the data with fitted regression lines 

x_breaks <- c(0, 0.2, 2, 6)

# colony
colony <- ggplot(Cu3_log10_3g, aes(x = Concentration , y = Colony, color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + 
  labs(x = "Toxic unit", y = "Colony (count)", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

colony

view(Cu3_log10_3g)

# colony fecundity

colony_f <- ggplot(Cu3_log10_3g, aes(x = Concentration, y = Colony_fecundity , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) +  
  labs(x = "Toxic unit", y = "Colony fecundity", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

colony_f

view(Cu3_log10_3g)


# sporangium fecundity

sporangium_f <- ggplot(Cu3_log10_3g, aes(x = Concentration, y = Sporangium_fecundity , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) +  
  labs(x = "Toxic unit", y = "Sporangium fecundity", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

sporangium_f


# survived zoospores 

survived_zoospores <- ggplot(Cu3_log10_3g, aes(x = Concentration, y = Survived_zoospores , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + 
  labs(x = "Toxic unit", y = "Survived zoospores/mL", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

survived_zoospores 

# sporangium

sporangium <- ggplot(Cu3_log10_3g, aes(x = Concentration, y = Sporangium , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + 
  labs(x = "Toxic unit", y = "Sporangium/mL", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

sporangium 

view(Cu3_log10_3g)

### Aggregate the plots

Cu3_all_plots <- grid.arrange(colony, colony_f, sporangium, sporangium_f, survived_zoospores,   ncol = 2, nrow=3) 

ggsave("Cu3_all_plots.png", Cu3_all_plots, width = 10, height = 6, units = "in")


################################ Effects of Zn


# make a subset of data for "Zn"

zn3_dat <- filter(tg3_dat, Metal == "Zn")

# calculate the colony fecundity (Cf = Zt/C)

view(zn3_dat)

zn3_dat <- zn3_dat %>%
  mutate(Colony_fecundity = Total_zoospores / Colony)

view(zn3_dat)

# calculate the sporangium fecundity (Sf = Zt/St)  

zn3_dat <- zn3_dat %>%
  mutate(Sporangium_fecundity = Total_zoospores / Sporangium)

view(zn3_dat)

Zn3_log10_3g <- write.csv(zn3_dat,"C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Zn3_log10_3g.csv", row.names = F)

Zn3_log10_3g <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Zn3_log10_3g.csv")


Zn3_log10_3g

View(Zn3_log10_3g)

head(Zn3_log10_3g)

## Plot the data with fitted regression lines 

x_breaks <- c(0, 0.02, 4, 8)

# colony
colony <- ggplot(Zn3_log10_3g, aes(x = Concentration , y = Colony, color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "glm", se = TRUE, alpha = 0.3, size = 1) +  
  labs(x = "Toxic unit", y = "Colony (count)", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

colony

view(Zn3_log10_3g)

# colony fecundity

colony_f <- ggplot(Zn3_log10_3g, aes(x = Concentration, y = Colony_fecundity , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) +  
  labs(x = "Toxic unit", y = "Colony fecundity", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

colony_f

view(Zn3_log10_3g)


# sporangium fecundity

sporangium_f <- ggplot(Zn3_log10_3g, aes(x = Concentration, y = Sporangium_fecundity , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) +  
  labs(x = "Toxic unit", y = "Sporangium fecundity", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

sporangium_f


# survived zoospores 

survived_zoospores <- ggplot(Zn3_log10_3g, aes(x = Concentration, y = Survived_zoospores , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + 
  labs(x = "Toxic unit", y = "Survived zoospores/mL", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

survived_zoospores 

# sporangium

sporangium <- ggplot(Zn3_log10_3g, aes(x = Concentration, y = Sporangium , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + 
  labs(x = "Toxic unit", y = "Sporangium/mL", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

sporangium 

view(Zn3_log10_3g)

### Aggregate the plots

Zn3_all_plots <- grid.arrange(colony, colony_f, sporangium, sporangium_f, survived_zoospores,   ncol = 2, nrow=3) 

ggsave("Zn3_all_plots.png", Zn3_all_plots, width = 10, height = 6, units = "in")



################################ Combination of Cu and Zn


# make a subset of data for "Mixed"

mixed3_dat <- filter(tg3_dat, Metal == "Mixed")

# calculate the colony fecundity (Cf = Zt/C)

view(mixed3_dat)

mixed3_dat <- mixed3_dat %>%
  mutate(Colony_fecundity = Total_zoospores / Colony)

view(mixed3_dat)

# calculate the sporangium fecundity (Sf = Zt/St)  

mixed3_dat <- mixed3_dat %>%
  mutate(Sporangium_fecundity = Total_zoospores / Sporangium)

view(mixed3_dat)

Mixed3_log10_3g <- write.csv(mixed3_dat,"C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Mixed3_log10_3g.csv", row.names = F)

Mixed3_log10_3g <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Mixed3_log10_3g.csv")


Mixed3_log10_3g

View(Mixed3_log10_3g)

head(Mixed3_log10_3g)

## Plot the data with fitted regression lines 

x_breaks <- c(0, 0.22, 6, 14)

# colony
colony <- ggplot(Mixed3_log10_3g, aes(x = Concentration , y = Colony, color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + 
  labs(x = "Total toxic unit", y = "Colony (count)", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

colony

view(Mixed3_log10_3g)

# colony fecundity

colony_f <- ggplot(Mixed3_log10_3g, aes(x = Concentration, y = Colony_fecundity , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) +  
  labs(x = "Total toxic unit", y = "Colony fecundity", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

colony_f

view(Mixed3_log10_3g)


# sporangium fecundity

sporangium_f <- ggplot(Mixed3_log10_3g, aes(x = Concentration, y = Sporangium_fecundity , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) +  
  labs(x = "Total toxic unit", y = "Sporangium fecundity", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

sporangium_f


# survived zoospores 

survived_zoospores <- ggplot(Mixed3_log10_3g, aes(x = Concentration, y = Survived_zoospores , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + 
  labs(x = "Total toxic unit", y = "Survived zoospores/mL", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

survived_zoospores 

# sporangium

sporangium <- ggplot(Mixed3_log10_3g, aes(x = Concentration, y = Sporangium , color = factor(Generation))) +
  geom_point(size=0.75) + scale_x_continuous(breaks = x_breaks)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1) + 
  labs(x = "Total toxic unit", y = "Sporangium/mL", color = "Generation") +
  theme_bw()+ ggeasy::easy_all_text_size(6) + ggeasy::easy_all_text_color("black")

sporangium 

view(Mixed3_log10_3g)

### Aggregate the plots

Mixed3_all_plots <- grid.arrange(colony, colony_f, sporangium, sporangium_f, survived_zoospores,   ncol = 2, nrow=3) 

ggsave("Mixed3_all_plots.png", Mixed3_all_plots, width = 10, height = 6, units = "in")



###########################################################

library(ssdtools)
library(ggplot2)
library(readr)

# read dataset
# the file argument of read_csv() assumes the file is in your working directory. You may need to change the file path to correctly read your dataset.

dat_zn <- read_csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/SSD_96h_LC50_Zn_frogs.csv")

view(dat_zn)

# fix unacceptable column names
colnames(dat_zn) <- make.names(colnames(dat_zn))

# fit distributions
dist <- ssd_fit_dists(dat_zn, left = 'Concentration', dists = c('lgumbel'), silent = TRUE, reweight = FALSE, min_pmix = 0, nrow = 6L, computable = TRUE, at_boundary_ok = FALSE, rescale = TRUE)
# plot distributions
ssd_plot_cdf(dist, delta = Inf)
# goodness of fit table
ssd_gof(dist)

# save plot
# width and height are in inches, dpi (dots per inch) sets resolution
ggsave('fit_dist_plot.png', width = 8 , height = 6 , dpi = 300)

# plot model average
# to add confidence intervals set ci = TRUE in predict and ssd_plot
# we recommend using nboot = 10000 in predict, although this may take several minutes to run
pred <- predict(dist, nboot = 10L, ci = TRUE)
zn_ssdplot <- ssd_plot(dat_zn, pred, left = 'Concentration', label = 'Species', color = NULL, shape = NULL, hc = NULL, ci = TRUE,
                       shift_x = 2.5, xlab = 'Concentration (µg/L)', ylab = 'Species Affected (%)') + ggeasy::easy_all_text_size(14)
ggtitle('') +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, colour='black'),
        axis.text = element_text(color = 'black'),
        legend.key = element_rect(fill = NA, colour = NA)) +
  expand_limits(x = 28379.99999) +
  scale_color_brewer(palette = 'Dark2', name = '-none-') +
  scale_shape(name = NULL)


zn_ssdplot

### Aggregate the plots

final_ssd_plot <- grid.arrange(zn_ssdplot,   ncol = 1, nrow=1) 

ggsave("final_ssd_plot.png", final_ssd_plot, width = 10, height = 6, units = "in")



# save plot
# width and height are in inches, dpi (dots per inch) sets resolution
ggsave('model_average_plot.png', width = 8 , height = 6 , dpi = 600)

# get confidence limits
# use the nboot argument in ssd_hp to set the number of bootstrap samples
rbind(ssd_hp(dist, conc = 1850L, ci = TRUE, nboot = 1000L), ssd_hp(dist, conc = 1850L, ci = TRUE, average = FALSE, nboot = 1000L))


######################################## Dose-response curves 
### Cu plots


Cu_curves <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd_experiments_LC50/Zoospore_survival_96h_Cu.csv")

x_breaks <- c(0, 0.5, 1, 2.5, 5, 10, 15)


# Create the plot using ggplot2
plot_cu <- ggplot(Cu_curves, aes(x = Concentration, y = Survival_probability, color = Test, shape = Test)) +
  geom_point(size = 3) + scale_x_continuous(breaks = x_breaks) + scale_x_log10() + # Set the size of the data points
  geom_line() + # Connect the data points with lines
  scale_shape_manual(values = c(1, 2)) + # Set custom point shapes
  labs(
    x = "Cu concentration (µg/L)",
    y = "Survival probability",
    color = "Test", # Set the legend title
    shape = "Test"  # Set the legend title for point shapes
  ) +
  theme_minimal() + # Use a minimal theme (you can customize this)
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.title = element_blank(), 
        panel.grid.major = element_line(color = "gray"), 
        panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0, color = "black") + 
  geom_vline(xintercept = 0, color = "black") + ggeasy::easy_all_text_size(10) +ggeasy::easy_all_text_color("black")

plot_cu


#### Zn plot

Zn_curves <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd_experiments_LC50/Zoospore_survival_96h_Zn.csv")

x_breaks <- c(0, 5, 50, 500, 1000, 1500, 2000)


# Create the plot using ggplot2
plot_zn <- ggplot(Zn_curves, aes(x = Concentration, y = Survival_probability, color = Test, shape = Test)) +
  geom_point(size = 3) + scale_x_continuous(breaks = x_breaks) + scale_x_log10() + # Set the size of the data points
  geom_line() + # Connect the data points with lines
  scale_shape_manual(values = c(1, 2)) + # Set custom point shapes
  labs(
    x = "Zn concentration (µg/L)",
    y = "Survival probability",
    color = "Test", # Set the legend title
    shape = "Test"  # Set the legend title for point shapes
  ) +
  theme_minimal() + # Use a minimal theme (you can customize this)
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.title = element_blank(), 
        panel.grid.major = element_line(color = "gray"), 
        panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0, color = "black") + 
  geom_vline(xintercept = 0, color = "black") + ggeasy::easy_all_text_size(10) +ggeasy::easy_all_text_color("black")

plot_zn


#### Mixed plot

Mixed_curves <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd_experiments_LC50/Zoospore_survival_96h_Mixed.csv")

x_breaks <- c(0, 1, 2, 3, 4, 5, 6)
x_labels <- c("0", "0.5+5", "1+50", "2.5+500", "5+1000", "10+1500", "15+2000")

# Create the plot using ggplot2
plot_mixed <- ggplot(Mixed_curves, aes(x = Concentration, y = Survival_probability, color = Test, shape = Test)) +
  geom_point(size = 3) + scale_x_continuous(breaks = x_breaks, labels = x_labels) + scale_x_log10() + # Set the size of the data points
  geom_line() + # Connect the data points with lines
  scale_shape_manual(values = c(1, 2)) + # Set custom point shapes
  labs(
    x = "Cu+Zn concentration (µg/L)",
    y = "Survival probability",
    color = "Test", # Set the legend title
    shape = "Test"  # Set the legend title for point shapes
  ) +
  theme_minimal() + # Use a minimal theme (you can customize this)
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.title = element_blank(), 
        panel.grid.major = element_line(color = "gray"), 
        panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0, color = "black") + 
  geom_vline(xintercept = 0, color = "black") + ggeasy::easy_all_text_size(10) +ggeasy::easy_all_text_color("black")

plot_mixed


######## Save all dose-response curves in one page 

final_drc_plot <- grid.arrange(plot_cu,plot_zn, plot_mixed,   ncol = 1, nrow=3) 

ggsave("final_drc_plot.png", final_drc_plot, width = 6, height = 10, units = "in")



R.version.string 






