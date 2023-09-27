

library(tidyverse)
library(lme4)
library(jagsUI)

#--------------------------------------------------------------------------------
# function to extract variables from JAGS summary
ext <- function(in.sum, x) {
  a <- in.sum[substr(rownames(in.sum),1,nchar(x)) == x, ]
  if(is.vector(a)) a <- t(matrix(a))
  return(data.frame(mean=a[,1], lcl=a[,3], ucl=a[,7], lcl50=a[,4], ucl50=a[,6]))
}

#--------------------------------------------------------------------------------

setwd("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions")

# survival data
dat <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/Bd survival_raw data.csv")

glimpse(dat)

# number of live and dead cells per tube
dat <- dat %>%
  rowwise() %>% 
  mutate(live = sum(c_across(L1:L5), na.rm = T),
         dead = sum(c_across(D1:D5), na.rm = T))

dat <- dat %>%
  mutate(prop.live = live / (live + dead),
         id = paste(Experiment, Metal, Tube),
         dur.days = Duration / 24)

glimpse(dat)
table(dat$id)

# plot data on raw scale
ggplot(dat, aes(y = prop.live, x = dur.days, colour = factor(Concentration))) +
  geom_point() +
  stat_smooth(n = 10, se = F) +
  facet_grid(Experiment ~ Metal) +
  theme_bw() +
  theme(legend.position = "top")


# plot data on log scale
ggplot(dat, aes(y = prop.live, x = dur.days, colour = factor(Concentration))) +
  geom_point() +
  stat_smooth(n = 10, se = F) +
  facet_grid(Experiment ~ Metal) +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position = "top")

# in experiment 3, proportion live at time = 0
hist(dat$prop.live[dat$Experiment == 3 & dat$Duration == 0])

################################################################################
# fit exponential survival model using lme4

###################################################################################################################################################################################################################################################
### Cu test 3

# choose subset of data to fit model to
sub.dat <- filter(dat, Experiment == 3 & Metal == "Cu")

# plot
ggplot(sub.dat, aes())

# allow for repeated measures by id
m1 <- glmer(cbind(live, dead) ~ dur.days:factor(Concentration) + 
              (dur.days|id), family = binomial(link = "log"), 
            data = sub.dat)
anova(m1)
summary(m1)

# make predictions
# first get combinations of duration and concentration to predict to
dur <- seq(0, 4, 0.1)
conc <- unique(sub.dat$Concentration)
newd <- expand.grid(dur, conc) 
names(newd) <- c("dur.days", "Concentration")
newd$id <- "new"

newd$pred <- predict(m1, newdata = newd, type = "response", allow.new.levels = T)
head(newd)

ggplot(sub.dat, aes(y = prop.live, x = dur.days)) +
  geom_line(aes(group = id), colour = "grey") +
  geom_point(aes(colour = factor(Concentration)), size = 3) +
  geom_line(data = newd, aes(y = pred, colour = factor(Concentration)), lwd = 1.1) +
  theme_bw() +
  theme(legend.position = "top")


######################################################################################
#-------------------------------------------------------------------------------------
# fit exponential model using JAGS

# data for JAGS model
n.live <- sub.dat$live
n.total <- sub.dat$live + sub.dat$dead
dur <- sub.dat$dur.days
id <- as.numeric(factor(sub.dat$id))

# numbers of observations / groups
N <- length(n.live)
N.id <- max(id)

# the concentration associated with each tube (id)
a <- table(sub.dat$id, sub.dat$Concentration)
a
conc <- apply(a, 1, function(x) which(x > 0))
N.conc <- max(conc)

# JAGS model
mod <- "model
  {
  # Likelihood
  for(i in 1:N) {
    n.live[i] ~ dbinom(mu[i], n.total[i])
    mu[i] <- int * exp(-lambda[id[i]] * dur[i])
  }
  
  for(i in 1:N.id) {
    lambda[i] ~ dnorm(mu.lambda[conc[i]], tau)
  }

  int ~ dunif(0, 1)
  for(i in 1:N.conc) {
    mu.lambda[i] ~ dnorm(0, 0.01) T(0, )
  }
  
  tau <- 1/(sigma * sigma)
  sigma ~ dunif(0, 10)


}"

# write model
write(mod, "survival.txt")

mod <- jags(model = "survival.txt",
            data = list(n.live = n.live, n.total = n.total, dur = dur, id = id,
                        conc = conc, N = N, N.id = N.id, N.conc = N.conc),
            param = c( "int", "mu.lambda", "sigma"),
            n.chains = 3,
            n.iter =11000,
            n.burnin = 1000,
            parallel = T)

all.sum <- mod$summary
all.sum 

#------------------------------------------------------------------------
# extract lambda coefficients and plot
lam <- ext(all.sum, "mu.lambda") %>%
  rename(lambda = mean)

# add concentrations
lam$Concentration <- as.numeric(levels(factor(sub.dat$Concentration)))

# plot results
ggplot(lam, aes(y = lambda, x = Concentration)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0) +
  theme_bw()


#------------------------------------------------------------------------
# fit to data
dur <- seq(0, 4, 0.01)
conc <- lam$Concentration
newd <- expand.grid(dur, conc)
names(newd) <- c("dur.days", "Concentration")

# add the mean lambda value for each concentration
newd <- left_join(newd, dplyr::select(lam, lambda, Concentration)) 

# add the intercept term
newd$int <- ext(all.sum, "int")[1, 1]

# predict values and calculate hazard
newd <- newd %>%
  mutate(pred = int * exp(-lambda * dur.days),
         haz = lambda)

head(newd)

ggplot(sub.dat, aes(y = prop.live, x = dur.days)) +
  geom_line(aes(group = id), colour = "grey") +
  geom_point(aes(colour = factor(Concentration)), size = 3) +
  geom_line(data = newd, aes(y = pred, colour = factor(Concentration)), lwd = 1.1) +
  theme_bw() +
  theme(legend.position = "top")

#--------------------------------------------------------------------------
# plot the hazard functions

ggplot(newd, aes(y = haz, x = dur.days)) +
  geom_line(aes(colour = factor(Concentration)), lwd = 1.1) +
  theme_classic() +
  theme(legend.position = "top")

#---------------------------------------------------------------------------
# plot observed versus predicted
# add the int and lambda values for each concentraion to the data
pred.dat <- left_join(sub.dat, dplyr::select(lam, lambda, Concentration))
pred.dat$int <- ext(all.sum, "int")[1, 1]

pred.dat <- pred.dat %>%
  mutate(pred = int*(exp(-lambda*dur.days)))

glimpse(pred.dat)

ggplot(pred.dat, aes(y = prop.live, x = pred)) +
  geom_point(aes(colour = factor(Concentration)), size = 4) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic()


##############################################################################
#-----------------------------------------------------------------------------
# fit Weibull model using JAGS

# data for JAGS model
n.live <- sub.dat$live
n.total <- sub.dat$live + sub.dat$dead
dur <- sub.dat$dur.days
id <- as.numeric(factor(sub.dat$id))

# numbers of observations / groups
N <- length(n.live)
N.id <- max(id)

# the concentration associated with each tube (id)
a <- table(sub.dat$id, sub.dat$Concentration)
a
conc <- apply(a, 1, function(x) which(x > 0))

N.conc <- max(conc)

# JAGS model
mod <- "model
  {
  # Likelihood
  for(i in 1:N) {
    n.live[i] ~ dbinom(mu[i], n.total[i])
    mu[i] <- int * exp(-(lambda[id[i]] * dur[i])^k[id[i]])
  }
  
  for(i in 1:N.id) {
    lambda[i] ~ dnorm(mu.lambda[conc[i]], tau)
    k[i] ~ dnorm(mu.k[conc[i]], tau)
  }

  int ~ dunif(0, 1)
  for(i in 1:N.conc) {
    mu.lambda[i] ~ dnorm(0, 0.1) T(0, )
    mu.k[i] ~ dnorm(0, 0.1) T(0, )
  }
  
  tau <- 1/(sigma * sigma)
  sigma ~ dunif(0, 10)


}"

# write model
write(mod, "survival.txt")

mod <- jags(model = "survival.txt",
            data = list(n.live = n.live, n.total = n.total, dur = dur, id = id,
                        conc = conc, N = N, N.id = N.id, N.conc = N.conc),
            param = c( "int", "mu.k", "mu.lambda","sigma"),
            n.chains = 3,
            n.iter =11000,
            n.burnin = 1000,
            parallel = T)

all.sum <- mod$summary
all.sum 

#------------------------------------------------------------------------
# extract lambda coefficients and plot
lam <- ext(all.sum, "mu.lambda") %>%
  dplyr::rename(lambda = mean)

# add concentrations
lam$Concentration <- as.numeric(levels(factor(sub.dat$Concentration)))

x_breaks <- c(0, 0.2, 0.4, 1, 2, 4, 6)


# plot results
lambda_cu3 <- ggplot(lam, aes(y = lambda, x = Concentration)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0) +scale_x_continuous(breaks = x_breaks, labels = x_breaks) + theme_bw()+ xlab("Toxic unit") + ylab("Lambda") + ggeasy::easy_all_text_size(8) + ggeasy::easy_all_text_color("black")

lambda_cu3

view(lam)
table(lam$lambda)

#------------------------------------------------------------------------
# extract k coefficients and plot
k <- ext(all.sum, "mu.k") %>%
  dplyr::rename(k = mean)

# add concentrations
k$Concentration <- as.numeric(levels(factor(sub.dat$Concentration)))

# Define breaks and labels for the x-axis
x_breaks <- c(0.2, 0.4, 1, 2, 4, 6)
x_labels <- c("0.2", "0.4", "1", "2", "4", "6")


# plot results
k_cu3 <- ggplot(k, aes(y = k, x = Concentration)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = c(0.2, 6)) +
  theme_bw() +
  xlab("Toxic unit") +
  ylab("K") +
  ggeasy::easy_all_text_size(8) +
  ggeasy::easy_all_text_color("black")

k_cu3


#------------------------------------------------------------------------
# fit to data
dur <- seq(0, 4, 0.01)
conc <- lam$Concentration
newd <- expand.grid(dur, conc)
names(newd) <- c("dur.days", "Concentration")

# add the mean lambda value for each concentration
newd <- left_join(newd, dplyr::select(lam, lambda, Concentration)) 

# add the k value for each concentration
newd <- left_join(newd, dplyr::select(k, k, Concentration)) 

# add the intercept term
newd$int <- ext(all.sum, "int")[1, 1]

# predict values and calculate the hazard
newd <- newd %>%
  mutate(pred = int * exp(-(lambda * dur.days)^k),
         haz = (lambda*k) * (lambda*dur.days)^k)

head(newd)

# Define breaks and labels for the y-axis

y_labels <- c("0", "0.25", "0.5", "0.75", "1")


exp_3_cu <- ggplot(sub.dat, aes(y = prop.live, x = dur.days)) +
  geom_line(aes(group = id), colour = "grey") +
  geom_point(aes(colour = factor(Concentration)), size = 2) +
  geom_line(data = newd, aes(y = pred, colour = factor(Concentration)), lwd = 1.1) +
  theme_classic() + scale_y_continuous(labels = y_labels, limits = c(0, 1))+
  theme(legend.position = "top")+
  theme_classic()+labs(x= "Exposure duration (day)",  y= "Survival probability") + ggeasy::easy_add_legend_title("Toxic unit") + ggeasy::easy_all_text_size(10) + ggeasy::easy_all_text_color("black") 

exp_3_cu

#------------------------------------------------------------------------------
# plot the hazard functions

ggplot(newd, aes(y = haz, x = dur.days)) +
  geom_line(aes(colour = factor(Concentration)), lwd = 1.1) +
  theme_classic() +
  theme(legend.position = "top")

#---------------------------------------------------------------------------
# plot observed versus predicted
# add the int, k and lambda values for each concentration to the data
pred.dat <- left_join(sub.dat, dplyr::select(lam, lambda, Concentration))
pred.dat <- left_join(pred.dat, dplyr::select(k, k, Concentration))
pred.dat$int <- ext(all.sum, "int")[1, 1]

pred.dat <- pred.dat %>%
  mutate(pred = int*(exp(-(lambda*dur.days)^k)))


x_breaks <- c(0, 0.25, 0.50, 0.75, 1.00)
y_breaks <- c(0, 0.25, 0.50, 0.75, 1.00)


exp3_cu_pred <- ggplot(pred.dat, aes(y = prop.live, x = pred)) +
  geom_point(aes(colour = factor(Concentration)), size = 2) +
  geom_abline(intercept = 0, slope = 1) + scale_x_continuous(breaks = x_breaks) + scale_y_continuous(breaks = y_breaks) +
  theme_classic()+labs(x= "Predicted probability",  y= "Observed probability") + ggeasy::easy_add_legend_title("Toxic unit") + ggeasy::easy_all_text_size(10) + ggeasy::easy_all_text_color("black") 

exp3_cu_pred





###################################################################################################################################################################################################################################################
### Cu test 2

# choose subset of data to fit model to
sub.dat <- filter(dat, Experiment == 2 & Metal == "Cu")

# plot
ggplot(sub.dat, aes())

# allow for repeated measures by id
m1 <- glmer(cbind(live, dead) ~ dur.days:factor(Concentration) + 
              (dur.days|id), family = binomial(link = "log"), 
            data = sub.dat)
anova(m1)
summary(m1)

# make predictions
# first get combinations of duration and concentration to predict to
dur <- seq(0, 4, 0.1)
conc <- unique(sub.dat$Concentration)
newd <- expand.grid(dur, conc) 
names(newd) <- c("dur.days", "Concentration")
newd$id <- "new"

newd$pred <- predict(m1, newdata = newd, type = "response", allow.new.levels = T)
head(newd)

ggplot(sub.dat, aes(y = prop.live, x = dur.days)) +
  geom_line(aes(group = id), colour = "grey") +
  geom_point(aes(colour = factor(Concentration)), size = 3) +
  geom_line(data = newd, aes(y = pred, colour = factor(Concentration)), lwd = 1.1) +
  theme_bw() +
  theme(legend.position = "top")

######################################################################################
#-------------------------------------------------------------------------------------
# fit exponential model using JAGS

# data for JAGS model
n.live <- sub.dat$live
n.total <- sub.dat$live + sub.dat$dead
dur <- sub.dat$dur.days
id <- as.numeric(factor(sub.dat$id))

# numbers of observations / groups
N <- length(n.live)
N.id <- max(id)

# the concentration associated with each tube (id)
a <- table(sub.dat$id, sub.dat$Concentration)
a
conc <- apply(a, 1, function(x) which(x > 0))
N.conc <- max(conc)

# JAGS model
mod <- "model
  {
  # Likelihood
  for(i in 1:N) {
    n.live[i] ~ dbinom(mu[i], n.total[i])
    mu[i] <- int * exp(-lambda[id[i]] * dur[i])
  }
  
  for(i in 1:N.id) {
    lambda[i] ~ dnorm(mu.lambda[conc[i]], tau)
  }

  int ~ dunif(0, 1)
  for(i in 1:N.conc) {
    mu.lambda[i] ~ dnorm(0, 0.01) T(0, )
  }
  
  tau <- 1/(sigma * sigma)
  sigma ~ dunif(0, 10)


}"

# write model
write(mod, "survival.txt")

mod <- jags(model = "survival.txt",
            data = list(n.live = n.live, n.total = n.total, dur = dur, id = id,
                        conc = conc, N = N, N.id = N.id, N.conc = N.conc),
            param = c( "int", "mu.lambda", "sigma"),
            n.chains = 3,
            n.iter =11000,
            n.burnin = 1000,
            parallel = T)

all.sum <- mod$summary
all.sum 

#------------------------------------------------------------------------
# extract lambda coefficients and plot
lam <- ext(all.sum, "mu.lambda") %>%
  dplyr::rename(lambda = mean)

# add concentrations
lam$Concentration <- as.numeric(levels(factor(sub.dat$Concentration)))

# plot results
ggplot(lam, aes(y = lambda, x = Concentration)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0) +
  theme_bw()

#------------------------------------------------------------------------
# fit to data
dur <- seq(0, 4, 0.01)
conc <- lam$Concentration
newd <- expand.grid(dur, conc)
names(newd) <- c("dur.days", "Concentration")

# add the mean lambda value for each concentration
newd <- left_join(newd, dplyr::select(lam, lambda, Concentration)) 

# add the intercept term
newd$int <- ext(all.sum, "int")[1, 1]

# predict values and calculate hazard
newd <- newd %>%
  mutate(pred = int * exp(-lambda * dur.days),
         haz = lambda)

head(newd)

ggplot(sub.dat, aes(y = prop.live, x = dur.days)) +
  geom_line(aes(group = id), colour = "grey") +
  geom_point(aes(colour = factor(Concentration)), size = 3) +
  geom_line(data = newd, aes(y = pred, colour = factor(Concentration)), lwd = 1.1) +
  theme_bw() +
  theme(legend.position = "top")

#--------------------------------------------------------------------------
# plot the hazard functions

ggplot(newd, aes(y = haz, x = dur.days)) +
  geom_line(aes(colour = factor(Concentration)), lwd = 1.1) +
  theme_classic() +
  theme(legend.position = "top")

#---------------------------------------------------------------------------
# plot observed versus predicted
# add the int and lambda values for each concentraion to the data
pred.dat <- left_join(sub.dat, dplyr::select(lam, lambda, Concentration))
pred.dat$int <- ext(all.sum, "int")[1, 1]

pred.dat <- pred.dat %>%
  mutate(pred = int*(exp(-lambda*dur.days)))

glimpse(pred.dat)

ggplot(pred.dat, aes(y = prop.live, x = pred)) +
  geom_point(aes(colour = factor(Concentration)), size = 4) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic()


##############################################################################
#-----------------------------------------------------------------------------
# fit Weibull model using JAGS

# data for JAGS model
n.live <- sub.dat$live
n.total <- sub.dat$live + sub.dat$dead
dur <- sub.dat$dur.days
id <- as.numeric(factor(sub.dat$id))

# numbers of observations / groups
N <- length(n.live)
N.id <- max(id)

# the concentration associated with each tube (id)
a <- table(sub.dat$id, sub.dat$Concentration)
a
conc <- apply(a, 1, function(x) which(x > 0))

N.conc <- max(conc)

# JAGS model
mod <- "model
  {
  # Likelihood
  for(i in 1:N) {
    n.live[i] ~ dbinom(mu[i], n.total[i])
    mu[i] <- int * exp(-(lambda[id[i]] * dur[i])^k[id[i]])
  }
  
  for(i in 1:N.id) {
    lambda[i] ~ dnorm(mu.lambda[conc[i]], tau)
    k[i] ~ dnorm(mu.k[conc[i]], tau)
  }

  int ~ dunif(0, 1)
  for(i in 1:N.conc) {
    mu.lambda[i] ~ dnorm(0, 0.1) T(0, )
    mu.k[i] ~ dnorm(0, 0.1) T(0, )
  }
  
  tau <- 1/(sigma * sigma)
  sigma ~ dunif(0, 10)


}"

# write model
write(mod, "survival.txt")

mod <- jags(model = "survival.txt",
            data = list(n.live = n.live, n.total = n.total, dur = dur, id = id,
                        conc = conc, N = N, N.id = N.id, N.conc = N.conc),
            param = c( "int", "mu.k", "mu.lambda","sigma"),
            n.chains = 3,
            n.iter =11000,
            n.burnin = 1000,
            parallel = T)

all.sum <- mod$summary
all.sum 

#------------------------------------------------------------------------
# extract lambda coefficients and plot
lam <- ext(all.sum, "mu.lambda") %>%
  dplyr::rename(lambda = mean)

# add concentrations
lam$Concentration <- as.numeric(levels(factor(sub.dat$Concentration)))

x_breaks <- c(0, 0.2, 0.4, 1, 2, 4, 6)


# plot results
lambda_cu2 <- ggplot(lam, aes(y = lambda, x = Concentration)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0) +scale_x_continuous(breaks = x_breaks, labels = x_breaks) + theme_bw()+ xlab("Toxic unit") + ylab("Lambda") + ggeasy::easy_all_text_size(8) + ggeasy::easy_all_text_color("black")

lambda_cu2

view(lam)
#------------------------------------------------------------------------
# extract k coefficients and plot
k <- ext(all.sum, "mu.k") %>%
  dplyr::rename(k = mean)

# add concentrations
k$Concentration <- as.numeric(levels(factor(sub.dat$Concentration)))

# Define breaks and labels for the x-axis
x_breaks <- c(0.2, 0.4, 1, 2, 4, 6)
x_labels <- c("0.2", "0.4", "1", "2", "4", "6")

View(k)

# plot results
k_cu2 <- ggplot(k, aes(y = k, x = Concentration)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = c(0.2, 6)) +
  theme_bw() +
  xlab("Toxic unit") +
  ylab("K") +
  ggeasy::easy_all_text_size(8) +
  ggeasy::easy_all_text_color("black")

k_cu2


#------------------------------------------------------------------------
# fit to data
dur <- seq(0, 4, 0.01)
conc <- lam$Concentration
newd <- expand.grid(dur, conc)
names(newd) <- c("dur.days", "Concentration")

# add the mean lambda value for each concentration
newd <- left_join(newd, dplyr::select(lam, lambda, Concentration)) 

# add the k value for each concentration
newd <- left_join(newd, dplyr::select(k, k, Concentration)) 

# add the intercept term
newd$int <- ext(all.sum, "int")[1, 1]

# predict values and calculate the hazard
newd <- newd %>%
  mutate(pred = int * exp(-(lambda * dur.days)^k),
         haz = (lambda*k) * (lambda*dur.days)^k)

head(newd)

# Define breaks and labels for the y-axis

y_labels <- c("0", "0.25", "0.5", "0.75", "1")


exp_2_cu <- ggplot(sub.dat, aes(y = prop.live, x = dur.days)) +
  geom_line(aes(group = id), colour = "grey") +
  geom_point(aes(colour = factor(Concentration)), size = 2) +
  geom_line(data = newd, aes(y = pred, colour = factor(Concentration)), lwd = 1.1) +
  theme_classic() + scale_y_continuous(labels = y_labels, limits = c(0, 1))+
  theme(legend.position = "top")+
  theme_classic()+labs(x= "Exposure duration (day)",  y= "Survival probability") + ggeasy::easy_add_legend_title("Toxic unit") + ggeasy::easy_all_text_size(10) + ggeasy::easy_all_text_color("black") 

exp_2_cu

#------------------------------------------------------------------------------
# plot the hazard functions

ggplot(newd, aes(y = haz, x = dur.days)) +
  geom_line(aes(colour = factor(Concentration)), lwd = 1.1) +
  theme_classic() +
  theme(legend.position = "top")

#---------------------------------------------------------------------------
# plot observed versus predicted
# add the int, k and lambda values for each concentration to the data
pred.dat <- left_join(sub.dat, dplyr::select(lam, lambda, Concentration))
pred.dat <- left_join(pred.dat, dplyr::select(k, k, Concentration))
pred.dat$int <- ext(all.sum, "int")[1, 1]

pred.dat <- pred.dat %>%
  mutate(pred = int*(exp(-(lambda*dur.days)^k)))


x_breaks <- c(0, 0.25, 0.50, 0.75, 1.00)
y_breaks <- c(0, 0.25, 0.50, 0.75, 1.00)


exp2_cu_pred <- ggplot(pred.dat, aes(y = prop.live, x = pred)) +
  geom_point(aes(colour = factor(Concentration)), size = 2) +
  geom_abline(intercept = 0, slope = 1) + scale_x_continuous(breaks = x_breaks) + scale_y_continuous(breaks = y_breaks) +
  theme_classic()+labs(x= "Predicted probability",  y= "Observed probability") + ggeasy::easy_add_legend_title("Toxic unit") + ggeasy::easy_all_text_size(10) + ggeasy::easy_all_text_color("black") 

exp2_cu_pred


###################################################################################################################################################################################################################################################
### Zn test 3

# choose subset of data to fit model to
sub.dat <- filter(dat, Experiment == 3 & Metal == "Zn")

# plot
ggplot(sub.dat, aes())

# allow for repeated measures by id
m1 <- glmer(cbind(live, dead) ~ dur.days:factor(Concentration) + 
              (dur.days|id), family = binomial(link = "log"), 
            data = sub.dat)
anova(m1)
summary(m1)

# make predictions
# first get combinations of duration and concentration to predict to
dur <- seq(0, 4, 0.1)
conc <- unique(sub.dat$Concentration)
newd <- expand.grid(dur, conc) 
names(newd) <- c("dur.days", "Concentration")
newd$id <- "new"

newd$pred <- predict(m1, newdata = newd, type = "response", allow.new.levels = T)
head(newd)

ggplot(sub.dat, aes(y = prop.live, x = dur.days)) +
  geom_line(aes(group = id), colour = "grey") +
  geom_point(aes(colour = factor(Concentration)), size = 3) +
  geom_line(data = newd, aes(y = pred, colour = factor(Concentration)), lwd = 1.1) +
  theme_bw() +
  theme(legend.position = "top")

######################################################################################
#-------------------------------------------------------------------------------------
# fit exponential model using JAGS

# data for JAGS model
n.live <- sub.dat$live
n.total <- sub.dat$live + sub.dat$dead
dur <- sub.dat$dur.days
id <- as.numeric(factor(sub.dat$id))

# numbers of observations / groups
N <- length(n.live)
N.id <- max(id)

# the concentration associated with each tube (id)
a <- table(sub.dat$id, sub.dat$Concentration)
a
conc <- apply(a, 1, function(x) which(x > 0))
N.conc <- max(conc)

# JAGS model
mod <- "model
  {
  # Likelihood
  for(i in 1:N) {
    n.live[i] ~ dbinom(mu[i], n.total[i])
    mu[i] <- int * exp(-lambda[id[i]] * dur[i])
  }
  
  for(i in 1:N.id) {
    lambda[i] ~ dnorm(mu.lambda[conc[i]], tau)
  }

  int ~ dunif(0, 1)
  for(i in 1:N.conc) {
    mu.lambda[i] ~ dnorm(0, 0.01) T(0, )
  }
  
  tau <- 1/(sigma * sigma)
  sigma ~ dunif(0, 10)


}"

# write model
write(mod, "survival.txt")

mod <- jags(model = "survival.txt",
            data = list(n.live = n.live, n.total = n.total, dur = dur, id = id,
                        conc = conc, N = N, N.id = N.id, N.conc = N.conc),
            param = c( "int", "mu.lambda", "sigma"),
            n.chains = 3,
            n.iter =11000,
            n.burnin = 1000,
            parallel = T)

all.sum <- mod$summary
all.sum 

#------------------------------------------------------------------------
# extract lambda coefficients and plot
lam <- ext(all.sum, "mu.lambda") %>%
  dplyr::rename(lambda = mean)

# add concentrations
lam$Concentration <- as.numeric(levels(factor(sub.dat$Concentration)))

# plot results
ggplot(lam, aes(y = lambda, x = Concentration)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0) +
  theme_bw()

#------------------------------------------------------------------------
# fit to data
dur <- seq(0, 4, 0.01)
conc <- lam$Concentration
newd <- expand.grid(dur, conc)
names(newd) <- c("dur.days", "Concentration")

# add the mean lambda value for each concentration
newd <- left_join(newd, dplyr::select(lam, lambda, Concentration)) 

# add the intercept term
newd$int <- ext(all.sum, "int")[1, 1]

# predict values and calculate hazard
newd <- newd %>%
  mutate(pred = int * exp(-lambda * dur.days),
         haz = lambda)

head(newd)

ggplot(sub.dat, aes(y = prop.live, x = dur.days)) +
  geom_line(aes(group = id), colour = "grey") +
  geom_point(aes(colour = factor(Concentration)), size = 3) +
  geom_line(data = newd, aes(y = pred, colour = factor(Concentration)), lwd = 1.1) +
  theme_bw() +
  theme(legend.position = "top")

#--------------------------------------------------------------------------
# plot the hazard functions

ggplot(newd, aes(y = haz, x = dur.days)) +
  geom_line(aes(colour = factor(Concentration)), lwd = 1.1) +
  theme_classic() +
  theme(legend.position = "top")

#---------------------------------------------------------------------------
# plot observed versus predicted
# add the int and lambda values for each concentraion to the data
pred.dat <- left_join(sub.dat, dplyr::select(lam, lambda, Concentration))
pred.dat$int <- ext(all.sum, "int")[1, 1]

pred.dat <- pred.dat %>%
  mutate(pred = int*(exp(-lambda*dur.days)))

glimpse(pred.dat)

ggplot(pred.dat, aes(y = prop.live, x = pred)) +
  geom_point(aes(colour = factor(Concentration)), size = 4) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic()


##############################################################################
#-----------------------------------------------------------------------------
# fit Weibull model using JAGS

# data for JAGS model
n.live <- sub.dat$live
n.total <- sub.dat$live + sub.dat$dead
dur <- sub.dat$dur.days
id <- as.numeric(factor(sub.dat$id))

# numbers of observations / groups
N <- length(n.live)
N.id <- max(id)

# the concentration associated with each tube (id)
a <- table(sub.dat$id, sub.dat$Concentration)
a
conc <- apply(a, 1, function(x) which(x > 0))

N.conc <- max(conc)

# JAGS model
mod <- "model
  {
  # Likelihood
  for(i in 1:N) {
    n.live[i] ~ dbinom(mu[i], n.total[i])
    mu[i] <- int * exp(-(lambda[id[i]] * dur[i])^k[id[i]])
  }
  
  for(i in 1:N.id) {
    lambda[i] ~ dnorm(mu.lambda[conc[i]], tau)
    k[i] ~ dnorm(mu.k[conc[i]], tau)
  }

  int ~ dunif(0, 1)
  for(i in 1:N.conc) {
    mu.lambda[i] ~ dnorm(0, 0.1) T(0, )
    mu.k[i] ~ dnorm(0, 0.1) T(0, )
  }
  
  tau <- 1/(sigma * sigma)
  sigma ~ dunif(0, 10)


}"

# write model
write(mod, "survival.txt")

mod <- jags(model = "survival.txt",
            data = list(n.live = n.live, n.total = n.total, dur = dur, id = id,
                        conc = conc, N = N, N.id = N.id, N.conc = N.conc),
            param = c( "int", "mu.k", "mu.lambda","sigma"),
            n.chains = 3,
            n.iter =11000,
            n.burnin = 1000,
            parallel = T)

all.sum <- mod$summary
all.sum 

#------------------------------------------------------------------------
# extract lambda coefficients and plot
lam <- ext(all.sum, "mu.lambda") %>%
  dplyr::rename(lambda = mean)

# add concentrations
lam$Concentration <- as.numeric(levels(factor(sub.dat$Concentration)))

x_breaks <- c(0, 0.02, 0.2, 2, 4, 6, 8)
x_labels <- c("0", "0.02", "0.2", "2", "4", "6", "8")


# plot results
lambda_zn3 <- ggplot(lam, aes(y = lambda, x = Concentration)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0) +scale_x_continuous(breaks = x_breaks, labels = x_labels) + theme_bw()+ xlab("Toxic unit") + ylab("Lambda") + ggeasy::easy_all_text_size(8) + ggeasy::easy_all_text_color("black")

lambda_zn3


view(lam)


#------------------------------------------------------------------------
# extract k coefficients and plot
k <- ext(all.sum, "mu.k") %>%
  dplyr::rename(k = mean)

# add concentrations
k$Concentration <- as.numeric(levels(factor(sub.dat$Concentration)))

# Define breaks and labels for the x-axis
x_breaks <- c(0.02, 0.2, 2, 4, 6, 8)
x_labels <- c("0.02", "0.2", "2", "4", "6", "8")

View(k)

# plot results
k_zn3 <- ggplot(k, aes(y = k, x = Concentration)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = c(0.2, 8)) +
  theme_bw() +
  xlab("Toxic unit") +
  ylab("K") +
  ggeasy::easy_all_text_size(8) +
  ggeasy::easy_all_text_color("black")

k_zn3


#------------------------------------------------------------------------
# fit to data
dur <- seq(0, 4, 0.01)
conc <- lam$Concentration
newd <- expand.grid(dur, conc)
names(newd) <- c("dur.days", "Concentration")

# add the mean lambda value for each concentration
newd <- left_join(newd, dplyr::select(lam, lambda, Concentration)) 

# add the k value for each concentration
newd <- left_join(newd, dplyr::select(k, k, Concentration)) 

# add the intercept term
newd$int <- ext(all.sum, "int")[1, 1]

# predict values and calculate the hazard
newd <- newd %>%
  mutate(pred = int * exp(-(lambda * dur.days)^k),
         haz = (lambda*k) * (lambda*dur.days)^k)

head(newd)

# Define breaks and labels for the y-axis

y_labels <- c("0", "0.25", "0.5", "0.75", "1")


exp_3_zn <- ggplot(sub.dat, aes(y = prop.live, x = dur.days)) +
  geom_line(aes(group = id), colour = "grey") +
  geom_point(aes(colour = factor(Concentration)), size = 2) +
  geom_line(data = newd, aes(y = pred, colour = factor(Concentration)), lwd = 1.1) +
  theme_classic() + scale_y_continuous(labels = y_labels, limits = c(0, 1))+
  theme(legend.position = "top")+
  theme_classic()+labs(x= "Exposure duration (day)",  y= "Survival probability") + ggeasy::easy_add_legend_title("Toxic unit") + ggeasy::easy_all_text_size(10) + ggeasy::easy_all_text_color("black") 

exp_3_zn

#------------------------------------------------------------------------------
# plot the hazard functions

ggplot(newd, aes(y = haz, x = dur.days)) +
  geom_line(aes(colour = factor(Concentration)), lwd = 1.1) +
  theme_classic() +
  theme(legend.position = "top")

#---------------------------------------------------------------------------
# plot observed versus predicted
# add the int, k and lambda values for each concentration to the data
pred.dat <- left_join(sub.dat, dplyr::select(lam, lambda, Concentration))
pred.dat <- left_join(pred.dat, dplyr::select(k, k, Concentration))
pred.dat$int <- ext(all.sum, "int")[1, 1]

pred.dat <- pred.dat %>%
  mutate(pred = int*(exp(-(lambda*dur.days)^k)))


x_breaks <- c(0, 0.25, 0.50, 0.75, 1.00)
y_breaks <- c(0, 0.25, 0.50, 0.75, 1.00)


exp3_zn_pred <- ggplot(pred.dat, aes(y = prop.live, x = pred)) +
  geom_point(aes(colour = factor(Concentration)), size = 2) +
  geom_abline(intercept = 0, slope = 1) + scale_x_continuous(breaks = x_breaks) + scale_y_continuous(breaks = y_breaks) +
  theme_classic()+labs(x= "Predicted probability",  y= "Observed probability") + ggeasy::easy_add_legend_title("Toxic unit") + ggeasy::easy_all_text_size(10) + ggeasy::easy_all_text_color("black") 

exp3_zn_pred


###################################################################################################################################################################################################################################################
### Zn test 2

# choose subset of data to fit model to
sub.dat <- filter(dat, Experiment == 2 & Metal == "Zn")

# plot
ggplot(sub.dat, aes())

# allow for repeated measures by id
m1 <- glmer(cbind(live, dead) ~ dur.days:factor(Concentration) + 
              (dur.days|id), family = binomial(link = "log"), 
            data = sub.dat)
anova(m1)
summary(m1)

# make predictions
# first get combinations of duration and concentration to predict to
dur <- seq(0, 4, 0.1)
conc <- unique(sub.dat$Concentration)
newd <- expand.grid(dur, conc) 
names(newd) <- c("dur.days", "Concentration")
newd$id <- "new"

newd$pred <- predict(m1, newdata = newd, type = "response", allow.new.levels = T)
head(newd)

ggplot(sub.dat, aes(y = prop.live, x = dur.days)) +
  geom_line(aes(group = id), colour = "grey") +
  geom_point(aes(colour = factor(Concentration)), size = 3) +
  geom_line(data = newd, aes(y = pred, colour = factor(Concentration)), lwd = 1.1) +
  theme_bw() +
  theme(legend.position = "top")

######################################################################################
#-------------------------------------------------------------------------------------
# fit exponential model using JAGS

# data for JAGS model
n.live <- sub.dat$live
n.total <- sub.dat$live + sub.dat$dead
dur <- sub.dat$dur.days
id <- as.numeric(factor(sub.dat$id))

# numbers of observations / groups
N <- length(n.live)
N.id <- max(id)

# the concentration associated with each tube (id)
a <- table(sub.dat$id, sub.dat$Concentration)
a
conc <- apply(a, 1, function(x) which(x > 0))
N.conc <- max(conc)

# JAGS model
mod <- "model
  {
  # Likelihood
  for(i in 1:N) {
    n.live[i] ~ dbinom(mu[i], n.total[i])
    mu[i] <- int * exp(-lambda[id[i]] * dur[i])
  }
  
  for(i in 1:N.id) {
    lambda[i] ~ dnorm(mu.lambda[conc[i]], tau)
  }

  int ~ dunif(0, 1)
  for(i in 1:N.conc) {
    mu.lambda[i] ~ dnorm(0, 0.01) T(0, )
  }
  
  tau <- 1/(sigma * sigma)
  sigma ~ dunif(0, 10)


}"

# write model
write(mod, "survival.txt")

mod <- jags(model = "survival.txt",
            data = list(n.live = n.live, n.total = n.total, dur = dur, id = id,
                        conc = conc, N = N, N.id = N.id, N.conc = N.conc),
            param = c( "int", "mu.lambda", "sigma"),
            n.chains = 3,
            n.iter =11000,
            n.burnin = 1000,
            parallel = T)

all.sum <- mod$summary
all.sum 

#------------------------------------------------------------------------
# extract lambda coefficients and plot
lam <- ext(all.sum, "mu.lambda") %>%
  dplyr::rename(lambda = mean)

# add concentrations
lam$Concentration <- as.numeric(levels(factor(sub.dat$Concentration)))

# plot results
ggplot(lam, aes(y = lambda, x = Concentration)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0) +
  theme_bw()

#------------------------------------------------------------------------
# fit to data
dur <- seq(0, 4, 0.01)
conc <- lam$Concentration
newd <- expand.grid(dur, conc)
names(newd) <- c("dur.days", "Concentration")

# add the mean lambda value for each concentration
newd <- left_join(newd, dplyr::select(lam, lambda, Concentration)) 

# add the intercept term
newd$int <- ext(all.sum, "int")[1, 1]

# predict values and calculate hazard
newd <- newd %>%
  mutate(pred = int * exp(-lambda * dur.days),
         haz = lambda)

head(newd)

ggplot(sub.dat, aes(y = prop.live, x = dur.days)) +
  geom_line(aes(group = id), colour = "grey") +
  geom_point(aes(colour = factor(Concentration)), size = 3) +
  geom_line(data = newd, aes(y = pred, colour = factor(Concentration)), lwd = 1.1) +
  theme_bw() +
  theme(legend.position = "top")

#--------------------------------------------------------------------------
# plot the hazard functions

ggplot(newd, aes(y = haz, x = dur.days)) +
  geom_line(aes(colour = factor(Concentration)), lwd = 1.1) +
  theme_classic() +
  theme(legend.position = "top")

#---------------------------------------------------------------------------
# plot observed versus predicted
# add the int and lambda values for each concentraion to the data
pred.dat <- left_join(sub.dat, dplyr::select(lam, lambda, Concentration))
pred.dat$int <- ext(all.sum, "int")[1, 1]

pred.dat <- pred.dat %>%
  mutate(pred = int*(exp(-lambda*dur.days)))

glimpse(pred.dat)

ggplot(pred.dat, aes(y = prop.live, x = pred)) +
  geom_point(aes(colour = factor(Concentration)), size = 4) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic()


##############################################################################
#-----------------------------------------------------------------------------
# fit Weibull model using JAGS

# data for JAGS model
n.live <- sub.dat$live
n.total <- sub.dat$live + sub.dat$dead
dur <- sub.dat$dur.days
id <- as.numeric(factor(sub.dat$id))

# numbers of observations / groups
N <- length(n.live)
N.id <- max(id)

# the concentration associated with each tube (id)
a <- table(sub.dat$id, sub.dat$Concentration)
a
conc <- apply(a, 1, function(x) which(x > 0))

N.conc <- max(conc)

# JAGS model
mod <- "model
  {
  # Likelihood
  for(i in 1:N) {
    n.live[i] ~ dbinom(mu[i], n.total[i])
    mu[i] <- int * exp(-(lambda[id[i]] * dur[i])^k[id[i]])
  }
  
  for(i in 1:N.id) {
    lambda[i] ~ dnorm(mu.lambda[conc[i]], tau)
    k[i] ~ dnorm(mu.k[conc[i]], tau)
  }

  int ~ dunif(0, 1)
  for(i in 1:N.conc) {
    mu.lambda[i] ~ dnorm(0, 0.1) T(0, )
    mu.k[i] ~ dnorm(0, 0.1) T(0, )
  }
  
  tau <- 1/(sigma * sigma)
  sigma ~ dunif(0, 10)


}"

# write model
write(mod, "survival.txt")

mod <- jags(model = "survival.txt",
            data = list(n.live = n.live, n.total = n.total, dur = dur, id = id,
                        conc = conc, N = N, N.id = N.id, N.conc = N.conc),
            param = c( "int", "mu.k", "mu.lambda","sigma"),
            n.chains = 3,
            n.iter =11000,
            n.burnin = 1000,
            parallel = T)

all.sum <- mod$summary
all.sum 

#------------------------------------------------------------------------
# extract lambda coefficients and plot
lam <- ext(all.sum, "mu.lambda") %>%
  dplyr::rename(lambda = mean)

# add concentrations
lam$Concentration <- as.numeric(levels(factor(sub.dat$Concentration)))

x_breaks <- c(0, 0.02, 0.2, 2, 4, 6, 8)
x_labels <- c("0", "0.02", "0.2", "2", "4", "6", "8")


# plot results
lambda_zn2 <- ggplot(lam, aes(y = lambda, x = Concentration)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0) +scale_x_continuous(breaks = x_breaks, labels = x_labels) + theme_bw()+ xlab("Toxic unit") + ylab("Lambda") + ggeasy::easy_all_text_size(8) + ggeasy::easy_all_text_color("black")

lambda_zn2

view(lam)



#------------------------------------------------------------------------
# extract k coefficients and plot
k <- ext(all.sum, "mu.k") %>%
  dplyr::rename(k = mean)

# add concentrations
k$Concentration <- as.numeric(levels(factor(sub.dat$Concentration)))

# Define breaks and labels for the x-axis
x_breaks <- c(0.02, 0.2, 2, 4, 6, 8)
x_labels <- c("0.02", "0.2", "2", "4", "6", "8")

View(k)

# plot results
k_zn2 <- ggplot(k, aes(y = k, x = Concentration)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = c(0.2, 8)) +
  theme_bw() +
  xlab("Toxic unit") +
  ylab("K") +
  ggeasy::easy_all_text_size(8) +
  ggeasy::easy_all_text_color("black")

k_zn2


#------------------------------------------------------------------------
# fit to data
dur <- seq(0, 4, 0.01)
conc <- lam$Concentration
newd <- expand.grid(dur, conc)
names(newd) <- c("dur.days", "Concentration")

# add the mean lambda value for each concentration
newd <- left_join(newd, dplyr::select(lam, lambda, Concentration)) 

# add the k value for each concentration
newd <- left_join(newd, dplyr::select(k, k, Concentration)) 

# add the intercept term
newd$int <- ext(all.sum, "int")[1, 1]

# predict values and calculate the hazard
newd <- newd %>%
  mutate(pred = int * exp(-(lambda * dur.days)^k),
         haz = (lambda*k) * (lambda*dur.days)^k)

head(newd)

# Define breaks and labels for the y-axis

y_labels <- c("0", "0.25", "0.5", "0.75", "1")


exp_2_zn <- ggplot(sub.dat, aes(y = prop.live, x = dur.days)) +
  geom_line(aes(group = id), colour = "grey") +
  geom_point(aes(colour = factor(Concentration)), size = 2) +
  geom_line(data = newd, aes(y = pred, colour = factor(Concentration)), lwd = 1.1) +
  theme_classic() + scale_y_continuous(labels = y_labels, limits = c(0, 1))+
  theme(legend.position = "top")+
  theme_classic()+labs(x= "Exposure duration (day)",  y= "Survival probability") + ggeasy::easy_add_legend_title("Toxic unit") + ggeasy::easy_all_text_size(10) + ggeasy::easy_all_text_color("black") 

exp_2_zn

#------------------------------------------------------------------------------
# plot the hazard functions

ggplot(newd, aes(y = haz, x = dur.days)) +
  geom_line(aes(colour = factor(Concentration)), lwd = 1.1) +
  theme_classic() +
  theme(legend.position = "top")

#---------------------------------------------------------------------------
# plot observed versus predicted
# add the int, k and lambda values for each concentration to the data
pred.dat <- left_join(sub.dat, dplyr::select(lam, lambda, Concentration))
pred.dat <- left_join(pred.dat, dplyr::select(k, k, Concentration))
pred.dat$int <- ext(all.sum, "int")[1, 1]

pred.dat <- pred.dat %>%
  mutate(pred = int*(exp(-(lambda*dur.days)^k)))


x_breaks <- c(0, 0.25, 0.50, 0.75, 1.00)
y_breaks <- c(0, 0.25, 0.50, 0.75, 1.00)


exp2_zn_pred <- ggplot(pred.dat, aes(y = prop.live, x = pred)) +
  geom_point(aes(colour = factor(Concentration)), size = 2) +
  geom_abline(intercept = 0, slope = 1) + scale_x_continuous(breaks = x_breaks) + scale_y_continuous(breaks = y_breaks) +
  theme_classic()+labs(x= "Predicted probability",  y= "Observed probability") + ggeasy::easy_add_legend_title("Toxic unit") + ggeasy::easy_all_text_size(10) + ggeasy::easy_all_text_color("black") 

exp2_zn_pred


###################################################################################################################################################################################################################################################
### Mixed test 3

# choose subset of data to fit model to
sub.dat <- filter(dat, Experiment == 3 & Metal == "Mixed")

# plot
ggplot(sub.dat, aes())

# allow for repeated measures by id
m1 <- glmer(cbind(live, dead) ~ dur.days:factor(Concentration) + 
              (dur.days|id), family = binomial(link = "log"), 
            data = sub.dat)
anova(m1)
summary(m1)

# make predictions
# first get combinations of duration and concentration to predict to
dur <- seq(0, 4, 0.1)
conc <- unique(sub.dat$Concentration)
newd <- expand.grid(dur, conc) 
names(newd) <- c("dur.days", "Concentration")
newd$id <- "new"

newd$pred <- predict(m1, newdata = newd, type = "response", allow.new.levels = T)
head(newd)

ggplot(sub.dat, aes(y = prop.live, x = dur.days)) +
  geom_line(aes(group = id), colour = "grey") +
  geom_point(aes(colour = factor(Concentration)), size = 3) +
  geom_line(data = newd, aes(y = pred, colour = factor(Concentration)), lwd = 1.1) +
  theme_bw() +
  theme(legend.position = "top")

######################################################################################
#-------------------------------------------------------------------------------------
# fit exponential model using JAGS

# data for JAGS model
n.live <- sub.dat$live
n.total <- sub.dat$live + sub.dat$dead
dur <- sub.dat$dur.days
id <- as.numeric(factor(sub.dat$id))

# numbers of observations / groups
N <- length(n.live)
N.id <- max(id)

# the concentration associated with each tube (id)
a <- table(sub.dat$id, sub.dat$Concentration)
a
conc <- apply(a, 1, function(x) which(x > 0))
N.conc <- max(conc)

# JAGS model
mod <- "model
  {
  # Likelihood
  for(i in 1:N) {
    n.live[i] ~ dbinom(mu[i], n.total[i])
    mu[i] <- int * exp(-lambda[id[i]] * dur[i])
  }
  
  for(i in 1:N.id) {
    lambda[i] ~ dnorm(mu.lambda[conc[i]], tau)
  }

  int ~ dunif(0, 1)
  for(i in 1:N.conc) {
    mu.lambda[i] ~ dnorm(0, 0.01) T(0, )
  }
  
  tau <- 1/(sigma * sigma)
  sigma ~ dunif(0, 10)


}"

# write model
write(mod, "survival.txt")

mod <- jags(model = "survival.txt",
            data = list(n.live = n.live, n.total = n.total, dur = dur, id = id,
                        conc = conc, N = N, N.id = N.id, N.conc = N.conc),
            param = c( "int", "mu.lambda", "sigma"),
            n.chains = 3,
            n.iter =11000,
            n.burnin = 1000,
            parallel = T)

all.sum <- mod$summary
all.sum 

#------------------------------------------------------------------------
# extract lambda coefficients and plot
lam <- ext(all.sum, "mu.lambda") %>%
  dplyr::rename(lambda = mean)

# add concentrations
lam$Concentration <- as.numeric(levels(factor(sub.dat$Concentration)))

# plot results
ggplot(lam, aes(y = lambda, x = Concentration)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0) +
  theme_bw()

#------------------------------------------------------------------------
# fit to data
dur <- seq(0, 4, 0.01)
conc <- lam$Concentration
newd <- expand.grid(dur, conc)
names(newd) <- c("dur.days", "Concentration")

# add the mean lambda value for each concentration
newd <- left_join(newd, dplyr::select(lam, lambda, Concentration)) 

# add the intercept term
newd$int <- ext(all.sum, "int")[1, 1]

# predict values and calculate hazard
newd <- newd %>%
  mutate(pred = int * exp(-lambda * dur.days),
         haz = lambda)

head(newd)

ggplot(sub.dat, aes(y = prop.live, x = dur.days)) +
  geom_line(aes(group = id), colour = "grey") +
  geom_point(aes(colour = factor(Concentration)), size = 3) +
  geom_line(data = newd, aes(y = pred, colour = factor(Concentration)), lwd = 1.1) +
  theme_bw() +
  theme(legend.position = "top")

#--------------------------------------------------------------------------
# plot the hazard functions

ggplot(newd, aes(y = haz, x = dur.days)) +
  geom_line(aes(colour = factor(Concentration)), lwd = 1.1) +
  theme_classic() +
  theme(legend.position = "top")

#---------------------------------------------------------------------------
# plot observed versus predicted
# add the int and lambda values for each concentraion to the data
pred.dat <- left_join(sub.dat, dplyr::select(lam, lambda, Concentration))
pred.dat$int <- ext(all.sum, "int")[1, 1]

pred.dat <- pred.dat %>%
  mutate(pred = int*(exp(-lambda*dur.days)))

glimpse(pred.dat)

ggplot(pred.dat, aes(y = prop.live, x = pred)) +
  geom_point(aes(colour = factor(Concentration)), size = 4) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic()


##############################################################################
#-----------------------------------------------------------------------------
# fit Weibull model using JAGS

# data for JAGS model
n.live <- sub.dat$live
n.total <- sub.dat$live + sub.dat$dead
dur <- sub.dat$dur.days
id <- as.numeric(factor(sub.dat$id))

# numbers of observations / groups
N <- length(n.live)
N.id <- max(id)

# the concentration associated with each tube (id)
a <- table(sub.dat$id, sub.dat$Concentration)
a
conc <- apply(a, 1, function(x) which(x > 0))

N.conc <- max(conc)

# JAGS model
mod <- "model
  {
  # Likelihood
  for(i in 1:N) {
    n.live[i] ~ dbinom(mu[i], n.total[i])
    mu[i] <- int * exp(-(lambda[id[i]] * dur[i])^k[id[i]])
  }
  
  for(i in 1:N.id) {
    lambda[i] ~ dnorm(mu.lambda[conc[i]], tau)
    k[i] ~ dnorm(mu.k[conc[i]], tau)
  }

  int ~ dunif(0, 1)
  for(i in 1:N.conc) {
    mu.lambda[i] ~ dnorm(0, 0.1) T(0, )
    mu.k[i] ~ dnorm(0, 0.1) T(0, )
  }
  
  tau <- 1/(sigma * sigma)
  sigma ~ dunif(0, 10)


}"

# write model
write(mod, "survival.txt")

mod <- jags(model = "survival.txt",
            data = list(n.live = n.live, n.total = n.total, dur = dur, id = id,
                        conc = conc, N = N, N.id = N.id, N.conc = N.conc),
            param = c( "int", "mu.k", "mu.lambda","sigma"),
            n.chains = 3,
            n.iter =11000,
            n.burnin = 1000,
            parallel = T)

all.sum <- mod$summary
all.sum 

#------------------------------------------------------------------------
# extract lambda coefficients and plot
lam <- ext(all.sum, "mu.lambda") %>%
  dplyr::rename(lambda = mean)

# add concentrations
lam$Concentration <- as.numeric(levels(factor(sub.dat$Concentration)))

x_breaks <- c(0, 0.22, 0.6, 3, 6, 10, 14)
x_labels <- c("0", "0.22", "0.6", "3", "6", "10", "14")


# plot results
lambda_mix3 <- ggplot(lam, aes(y = lambda, x = Concentration)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0) +scale_x_continuous(breaks = x_breaks, labels = x_labels) + theme_bw()+ xlab("Total toxic unit") + ylab("Lambda") + ggeasy::easy_all_text_size(8) + ggeasy::easy_all_text_color("black")

lambda_mix3

view(lam)



#------------------------------------------------------------------------
# extract k coefficients and plot
k <- ext(all.sum, "mu.k") %>%
  dplyr::rename(k = mean)

# add concentrations
k$Concentration <- as.numeric(levels(factor(sub.dat$Concentration)))

# Define breaks and labels for the x-axis
x_breaks <- c(0, 0.22, 0.6, 3, 6, 10, 14)
x_labels <- c("0", "0.22", "0.6", "3", "6", "10", "14")

View(k)

# plot results
k_mix3 <- ggplot(k, aes(y = k, x = Concentration)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = c(0.2, 14)) +
  theme_bw() +
  xlab("Total toxic unit") +
  ylab("K") +
  ggeasy::easy_all_text_size(8) +
  ggeasy::easy_all_text_color("black")

k_mix3


#------------------------------------------------------------------------
# fit to data
dur <- seq(0, 4, 0.01)
conc <- lam$Concentration
newd <- expand.grid(dur, conc)
names(newd) <- c("dur.days", "Concentration")

# add the mean lambda value for each concentration
newd <- left_join(newd, dplyr::select(lam, lambda, Concentration)) 

# add the k value for each concentration
newd <- left_join(newd, dplyr::select(k, k, Concentration)) 

# add the intercept term
newd$int <- ext(all.sum, "int")[1, 1]

# predict values and calculate the hazard
newd <- newd %>%
  mutate(pred = int * exp(-(lambda * dur.days)^k),
         haz = (lambda*k) * (lambda*dur.days)^k)

head(newd)

# Define breaks and labels for the y-axis

y_labels <- c("0", "0.25", "0.5", "0.75", "1")


exp_3_mix <- ggplot(sub.dat, aes(y = prop.live, x = dur.days)) +
  geom_line(aes(group = id), colour = "grey") +
  geom_point(aes(colour = factor(Concentration)), size = 2) +
  geom_line(data = newd, aes(y = pred, colour = factor(Concentration)), lwd = 1.1) +
  theme_classic() + scale_y_continuous(labels = y_labels, limits = c(0, 1))+
  theme(legend.position = "top")+
  theme_classic()+labs(x= "Exposure duration (day)",  y= "Survival probability") + ggeasy::easy_add_legend_title("Total toxic unit") + ggeasy::easy_all_text_size(10) + ggeasy::easy_all_text_color("black") 

exp_3_mix

#------------------------------------------------------------------------------
# plot the hazard functions

ggplot(newd, aes(y = haz, x = dur.days)) +
  geom_line(aes(colour = factor(Concentration)), lwd = 1.1) +
  theme_classic() +
  theme(legend.position = "top")

#---------------------------------------------------------------------------
# plot observed versus predicted
# add the int, k and lambda values for each concentration to the data
pred.dat <- left_join(sub.dat, dplyr::select(lam, lambda, Concentration))
pred.dat <- left_join(pred.dat, dplyr::select(k, k, Concentration))
pred.dat$int <- ext(all.sum, "int")[1, 1]

pred.dat <- pred.dat %>%
  mutate(pred = int*(exp(-(lambda*dur.days)^k)))


x_breaks <- c(0, 0.25, 0.50, 0.75, 1.00)
y_breaks <- c(0, 0.25, 0.50, 0.75, 1.00)


exp3_mix_pred <- ggplot(pred.dat, aes(y = prop.live, x = pred)) +
  geom_point(aes(colour = factor(Concentration)), size = 2) +
  geom_abline(intercept = 0, slope = 1) + scale_x_continuous(breaks = x_breaks) + scale_y_continuous(breaks = y_breaks) +
  theme_classic()+labs(x= "Predicted probability",  y= "Observed probability") + ggeasy::easy_add_legend_title("Total toxic unit") + ggeasy::easy_all_text_size(10) + ggeasy::easy_all_text_color("black") 

exp3_mix_pred

###################################################################################################################################################################################################################################################
### Mixed test 2

# choose subset of data to fit model to
sub.dat <- filter(dat, Experiment == 2 & Metal == "Mixed")

# plot
ggplot(sub.dat, aes())

# allow for repeated measures by id
m1 <- glmer(cbind(live, dead) ~ dur.days:factor(Concentration) + 
              (dur.days|id), family = binomial(link = "log"), 
            data = sub.dat)
anova(m1)
summary(m1)

# make predictions
# first get combinations of duration and concentration to predict to
dur <- seq(0, 4, 0.1)
conc <- unique(sub.dat$Concentration)
newd <- expand.grid(dur, conc) 
names(newd) <- c("dur.days", "Concentration")
newd$id <- "new"

newd$pred <- predict(m1, newdata = newd, type = "response", allow.new.levels = T)
head(newd)

ggplot(sub.dat, aes(y = prop.live, x = dur.days)) +
  geom_line(aes(group = id), colour = "grey") +
  geom_point(aes(colour = factor(Concentration)), size = 3) +
  geom_line(data = newd, aes(y = pred, colour = factor(Concentration)), lwd = 1.1) +
  theme_bw() +
  theme(legend.position = "top")

######################################################################################
#-------------------------------------------------------------------------------------
# fit exponential model using JAGS

# data for JAGS model
n.live <- sub.dat$live
n.total <- sub.dat$live + sub.dat$dead
dur <- sub.dat$dur.days
id <- as.numeric(factor(sub.dat$id))

# numbers of observations / groups
N <- length(n.live)
N.id <- max(id)

# the concentration associated with each tube (id)
a <- table(sub.dat$id, sub.dat$Concentration)
a
conc <- apply(a, 1, function(x) which(x > 0))
N.conc <- max(conc)

# JAGS model
mod <- "model
  {
  # Likelihood
  for(i in 1:N) {
    n.live[i] ~ dbinom(mu[i], n.total[i])
    mu[i] <- int * exp(-lambda[id[i]] * dur[i])
  }
  
  for(i in 1:N.id) {
    lambda[i] ~ dnorm(mu.lambda[conc[i]], tau)
  }

  int ~ dunif(0, 1)
  for(i in 1:N.conc) {
    mu.lambda[i] ~ dnorm(0, 0.01) T(0, )
  }
  
  tau <- 1/(sigma * sigma)
  sigma ~ dunif(0, 10)


}"

# write model
write(mod, "survival.txt")

mod <- jags(model = "survival.txt",
            data = list(n.live = n.live, n.total = n.total, dur = dur, id = id,
                        conc = conc, N = N, N.id = N.id, N.conc = N.conc),
            param = c( "int", "mu.lambda", "sigma"),
            n.chains = 3,
            n.iter =11000,
            n.burnin = 1000,
            parallel = T)

all.sum <- mod$summary
all.sum 

#------------------------------------------------------------------------
# extract lambda coefficients and plot
lam <- ext(all.sum, "mu.lambda") %>%
  dplyr::rename(lambda = mean)

# add concentrations
lam$Concentration <- as.numeric(levels(factor(sub.dat$Concentration)))

# plot results
ggplot(lam, aes(y = lambda, x = Concentration)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0) +
  theme_bw()

#------------------------------------------------------------------------
# fit to data
dur <- seq(0, 4, 0.01)
conc <- lam$Concentration
newd <- expand.grid(dur, conc)
names(newd) <- c("dur.days", "Concentration")

# add the mean lambda value for each concentration
newd <- left_join(newd, dplyr::select(lam, lambda, Concentration)) 

# add the intercept term
newd$int <- ext(all.sum, "int")[1, 1]

# predict values and calculate hazard
newd <- newd %>%
  mutate(pred = int * exp(-lambda * dur.days),
         haz = lambda)

head(newd)

ggplot(sub.dat, aes(y = prop.live, x = dur.days)) +
  geom_line(aes(group = id), colour = "grey") +
  geom_point(aes(colour = factor(Concentration)), size = 3) +
  geom_line(data = newd, aes(y = pred, colour = factor(Concentration)), lwd = 1.1) +
  theme_bw() +
  theme(legend.position = "top")

#--------------------------------------------------------------------------
# plot the hazard functions

ggplot(newd, aes(y = haz, x = dur.days)) +
  geom_line(aes(colour = factor(Concentration)), lwd = 1.1) +
  theme_classic() +
  theme(legend.position = "top")

#---------------------------------------------------------------------------
# plot observed versus predicted
# add the int and lambda values for each concentraion to the data
pred.dat <- left_join(sub.dat, dplyr::select(lam, lambda, Concentration))
pred.dat$int <- ext(all.sum, "int")[1, 1]

pred.dat <- pred.dat %>%
  mutate(pred = int*(exp(-lambda*dur.days)))

glimpse(pred.dat)

ggplot(pred.dat, aes(y = prop.live, x = pred)) +
  geom_point(aes(colour = factor(Concentration)), size = 4) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic()


##############################################################################
#-----------------------------------------------------------------------------
# fit Weibull model using JAGS

# data for JAGS model
n.live <- sub.dat$live
n.total <- sub.dat$live + sub.dat$dead
dur <- sub.dat$dur.days
id <- as.numeric(factor(sub.dat$id))

# numbers of observations / groups
N <- length(n.live)
N.id <- max(id)

# the concentration associated with each tube (id)
a <- table(sub.dat$id, sub.dat$Concentration)
a
conc <- apply(a, 1, function(x) which(x > 0))

N.conc <- max(conc)

# JAGS model
mod <- "model
  {
  # Likelihood
  for(i in 1:N) {
    n.live[i] ~ dbinom(mu[i], n.total[i])
    mu[i] <- int * exp(-(lambda[id[i]] * dur[i])^k[id[i]])
  }
  
  for(i in 1:N.id) {
    lambda[i] ~ dnorm(mu.lambda[conc[i]], tau)
    k[i] ~ dnorm(mu.k[conc[i]], tau)
  }

  int ~ dunif(0, 1)
  for(i in 1:N.conc) {
    mu.lambda[i] ~ dnorm(0, 0.1) T(0, )
    mu.k[i] ~ dnorm(0, 0.1) T(0, )
  }
  
  tau <- 1/(sigma * sigma)
  sigma ~ dunif(0, 10)


}"

# write model
write(mod, "survival.txt")

mod <- jags(model = "survival.txt",
            data = list(n.live = n.live, n.total = n.total, dur = dur, id = id,
                        conc = conc, N = N, N.id = N.id, N.conc = N.conc),
            param = c( "int", "mu.k", "mu.lambda","sigma"),
            n.chains = 3,
            n.iter =11000,
            n.burnin = 1000,
            parallel = T)

all.sum <- mod$summary
all.sum 

#------------------------------------------------------------------------
# extract lambda coefficients and plot
lam <- ext(all.sum, "mu.lambda") %>%
  dplyr::rename(lambda = mean)

# add concentrations
lam$Concentration <- as.numeric(levels(factor(sub.dat$Concentration)))

x_breaks <- c(0, 0.22, 0.6, 3, 6, 10, 14)
x_labels <- c("0", "0.22", "0.6", "3", "6", "10", "14")


# plot results
lambda_mix2 <- ggplot(lam, aes(y = lambda, x = Concentration)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0) +scale_x_continuous(breaks = x_breaks, labels = x_labels) + theme_bw()+ xlab("Total toxic unit") + ylab("Lambda") + ggeasy::easy_all_text_size(8) + ggeasy::easy_all_text_color("black")

lambda_mix2

view(lam)


#------------------------------------------------------------------------
# extract k coefficients and plot
k <- ext(all.sum, "mu.k") %>%
  dplyr::rename(k = mean)

# add concentrations
k$Concentration <- as.numeric(levels(factor(sub.dat$Concentration)))

# Define breaks and labels for the x-axis
x_breaks <- c(0, 0.22, 0.6, 3, 6, 10, 14)
x_labels <- c("0", "0.22", "0.6", "3", "6", "10", "14")

View(k)

# plot results
k_mix2 <- ggplot(k, aes(y = k, x = Concentration)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = c(0.2, 14)) +
  theme_bw() +
  xlab("Total toxic unit") +
  ylab("K") +
  ggeasy::easy_all_text_size(8) +
  ggeasy::easy_all_text_color("black")

k_mix2


#------------------------------------------------------------------------
# fit to data
dur <- seq(0, 4, 0.01)
conc <- lam$Concentration
newd <- expand.grid(dur, conc)
names(newd) <- c("dur.days", "Concentration")

# add the mean lambda value for each concentration
newd <- left_join(newd, dplyr::select(lam, lambda, Concentration)) 

# add the k value for each concentration
newd <- left_join(newd, dplyr::select(k, k, Concentration)) 

# add the intercept term
newd$int <- ext(all.sum, "int")[1, 1]

# predict values and calculate the hazard
newd <- newd %>%
  mutate(pred = int * exp(-(lambda * dur.days)^k),
         haz = (lambda*k) * (lambda*dur.days)^k)

head(newd)

# Define breaks and labels for the y-axis

y_labels <- c("0", "0.25", "0.5", "0.75", "1")


exp_2_mix <- ggplot(sub.dat, aes(y = prop.live, x = dur.days)) +
  geom_line(aes(group = id), colour = "grey") +
  geom_point(aes(colour = factor(Concentration)), size = 2) +
  geom_line(data = newd, aes(y = pred, colour = factor(Concentration)), lwd = 1.1) +
  theme_classic() + scale_y_continuous(labels = y_labels, limits = c(0, 1))+
  theme(legend.position = "top")+
  theme_classic()+labs(x= "Exposure duration (day)",  y= "Survival probability") + ggeasy::easy_add_legend_title("Total toxic unit") + ggeasy::easy_all_text_size(10) + ggeasy::easy_all_text_color("black") 

exp_2_mix

#------------------------------------------------------------------------------
# plot the hazard functions

ggplot(newd, aes(y = haz, x = dur.days)) +
  geom_line(aes(colour = factor(Concentration)), lwd = 1.1) +
  theme_classic() +
  theme(legend.position = "top")

#---------------------------------------------------------------------------
# plot observed versus predicted
# add the int, k and lambda values for each concentration to the data
pred.dat <- left_join(sub.dat, dplyr::select(lam, lambda, Concentration))
pred.dat <- left_join(pred.dat, dplyr::select(k, k, Concentration))
pred.dat$int <- ext(all.sum, "int")[1, 1]

pred.dat <- pred.dat %>%
  mutate(pred = int*(exp(-(lambda*dur.days)^k)))


x_breaks <- c(0, 0.25, 0.50, 0.75, 1.00)
y_breaks <- c(0, 0.25, 0.50, 0.75, 1.00)


exp2_mix_pred <- ggplot(pred.dat, aes(y = prop.live, x = pred)) +
  geom_point(aes(colour = factor(Concentration)), size = 2) +
  geom_abline(intercept = 0, slope = 1) + scale_x_continuous(breaks = x_breaks) + scale_y_continuous(breaks = y_breaks) +
  theme_classic()+labs(x= "Predicted probability",  y= "Observed probability") + ggeasy::easy_add_legend_title("Total toxic unit") + ggeasy::easy_all_text_size(10) + ggeasy::easy_all_text_color("black") 

exp2_mix_pred


############# Put all the plots together 

### Predicted vs. observed survival probability 

# Add labels to the plots
exp2_cu_pred <- exp2_cu_pred + labs(title = "Cu_test1")
exp3_cu_pred <- exp3_cu_pred + labs(title = "Cu_test2")
exp2_zn_pred <- exp2_zn_pred + labs(title = "Zn_test1")
exp3_zn_pred <- exp3_zn_pred + labs(title = "Zn_test2")
exp2_mix_pred <- exp2_mix_pred + labs(title = "Combined_test1")
exp3_mix_pred <- exp3_mix_pred + labs(title = "Combined_test2")


pred_obsrv <- grid.arrange(exp2_cu_pred, exp3_cu_pred, exp2_zn_pred, exp3_zn_pred, exp2_mix_pred, exp3_mix_pred, ncol = 2, nrow=3) 

ggsave("pred_obsrv.png", pred_obsrv, width = 8, height = 9, units = "in")

### Exponential survival probability with fitted models 

# Add labels to the plots
exp_2_cu <- exp_2_cu + labs(title = "Cu_test1")
exp_3_cu <- exp_3_cu + labs(title = "Cu_test2")
exp_2_zn <- exp_2_zn + labs(title = "Zn_test1")
exp_3_zn <- exp_3_zn + labs(title = "Zn_test2")
exp_2_mix <- exp_2_mix + labs(title = "Combined_test1")
exp_3_mix <- exp_3_mix + labs(title = "Combined_test2")

expo_survive <- grid.arrange(exp_2_cu, exp_3_cu, exp_2_zn, exp_3_zn, exp_2_mix, exp_3_mix, ncol = 2, nrow=3) 

ggsave("expo_survive.png", expo_survive, width = 8, height = 9, units = "in")



### Lambda and K values
## Test 2
# Add labels to the plots
lambda_cu2 <- lambda_cu2 + labs(title = "Cu_test1")
k_cu2 <- k_cu2 + labs(title = "Cu_test1")
lambda_zn2 <- lambda_zn2 + labs(title = "Zn_test1")
k_zn2 <- k_zn2 + labs(title = "Zn_test1")
lambda_mix2 <- lambda_mix2 + labs(title = "Combined_test1")
k_mix2 <- k_mix2 + labs(title = "Combined_test1")

lambda_k_2 <- grid.arrange(lambda_cu2, k_cu2, lambda_zn2, k_zn2, lambda_mix2, k_mix2, ncol = 2, nrow=3) 

ggsave("lambda_k_2.png", lambda_k_2, width = 9, height = 8, units = "in")


### Lambda and K values
## Test 3
# Add labels to the plots
lambda_cu3 <- lambda_cu3 + labs(title = "Cu_test2")
k_cu3 <- k_cu3 + labs(title = "Cu_test2")
lambda_zn3 <- lambda_zn3 + labs(title = "Zn_test2")
k_zn3 <- k_zn3 + labs(title = "Zn_test2")
lambda_mix3 <- lambda_mix3 + labs(title = "Combined_test2")
k_mix3 <- k_mix3 + labs(title = "Combined_test2")

lambda_k_3 <- grid.arrange(lambda_cu3, k_cu3, lambda_zn3, k_zn3, lambda_mix3, k_mix3, ncol = 2, nrow=3) 

ggsave("lambda_k_3.png", lambda_k_3, width = 9, height = 8, units = "in")








############################ Analyse the combined effects of Cu and Zn

### Cu test 3

view(dat)
head(dat)

lambda <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/Bd experiments data_new versions/lambda values for all tests.csv")

view(lambda)


# Load the necessary library
library(ggplot2)
library(dplyr)

# Create a new column for the sum of Cu and Zn lambda values
lambda <- lambda %>%
  mutate(Sum_Cu_Zn = ifelse(Metal %in% c("Cu", "Zn"), "Cu+Zn", Metal))

# Calculate the sum of Cu and Zn lambda values
lambda <- lambda %>%
  group_by(Experiment, Concentration, Sum_Cu_Zn) %>%
  summarise(lambda = sum(lambda))

view(lambda)

# Plot for Experiment 3
plot_exp_3 <- ggplot(lambda[lambda$Experiment == 3, ], aes(x = Concentration, y = lambda, color = Sum_Cu_Zn)) +
  geom_line(aes(group = Sum_Cu_Zn)) +
  geom_ribbon(data = lambda[lambda$Experiment == 3, ], aes(ymin = lcl, ymax = ucl, fill = Sum_Cu_Zn), alpha = 0.2) +
  labs(title = "Experiment 3", x = "Concentration", y = "Lambda") +
  scale_color_manual(values = c("Cu" = "red", "Zn" = "blue", "Mixed" = "green", "Cu+Zn" = "purple")) +
  scale_fill_manual(values = c("Cu" = "red", "Zn" = "blue", "Mixed" = "green", "Cu+Zn" = "purple")) +
  theme_minimal()

# Plot for Experiment 2
plot_exp_2 <- ggplot(lambda[lambda$Experiment == 2, ], aes(x = Concentration, y = lambda, color = Sum_Cu_Zn)) +
  geom_line(aes(group = Sum_Cu_Zn)) +
  geom_ribbon(data = lambda[lambda$Experiment == 2, ], aes(ymin = lcl, ymax = ucl, fill = Sum_Cu_Zn), alpha = 0.2) +
  labs(title = "Experiment 2", x = "Concentration", y = "Lambda") +
  scale_color_manual(values = c("Cu" = "red", "Zn" = "blue", "Mixed" = "green", "Cu+Zn" = "purple")) +
  scale_fill_manual(values = c("Cu" = "red", "Zn" = "blue", "Mixed" = "green", "Cu+Zn" = "purple")) +
  theme_minimal()

# Display the plots
plot_exp_3
plot_exp_2








