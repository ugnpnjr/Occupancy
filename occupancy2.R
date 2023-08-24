
# Single-season occupancy model (SSOM) (MacKenzie et al., 2002)
# April 2023
# Ugyen Penjor, Fauna & Flora

# Chimp data from Sapo NP, Liberia
# 2016 Dry season
# 56 sites, 4 replicates (temporal)
# Site covariates: Dist. village, dist. road, dist. river, elevation, treecover
# Observation covariates: human encounter frequency, transect length

# Load packages
#install.packages("unmarked", "AICcmodavg", "MuMIn")   # Do only once
library(unmarked)
library(AICcmodavg)
library(MuMIn)

### ------------------------------------------------------------------------ ###

# Load data 
dat <- read.csv("your_data.csv", header=T)
head(dat)

# Extract detection history
detn <- dat[, 2:5]
head(detn)

# Pull site covariates (occupancy) and standarise them
dvilS <- as.data.frame(scale(dat$dvillage))
droadS <- as.data.frame(scale(dat$droad))
drivS <- as.data.frame(scale(dat$driver))
treeS <- as.data.frame(scale(dat$treecover))
elevS <- as.data.frame(scale(dat$elevation))

tlength <- as.data.frame(scale(dat$tlength))
colnames(tlength) <- "tlength"

# Create a new data frame with standardised covariates
site.cov <- cbind(dvilS, droadS, drivS, treeS, elevS, tlength)
colnames(site.cov) <- c("dvillage", "droad", "driver", "treecover", "elevation", "tlength")
head(site.cov)

# Observation covariate (detection)
human <- dat[, 12:15]
humanS <- as.data.frame(scale(human))

### Calculate naive occupancy (sites with detection/total no. of sites)
naive.occ <- sum(ifelse(rowSums(detn, na.rm=T)>0, 1, 0))/nrow(detn)
naive.occ

# Prepare unmarked frame object 
sp.umf <- unmarkedFrameOccu(y=detn, siteCovs=site.cov, obsCovs=list(human=humanS))
str(sp.umf)

### -----------------------------------------------------------------------  ###

# First model detection probability
d1 <- occu(~1 ~dvillage+droad+driver+treecover+elevation, data=sp.umf)
d2 <- occu(~human ~dvillage+droad+driver+treecover+elevation, data=sp.umf)
d3 <- occu(~tlength ~dvillage+droad+driver+treecover+elevation, data=sp.umf)
d4 <- occu(~human+tlength ~dvillage+droad+driver+treecover+elevation, data=sp.umf)

det.mod <- list("p(.)"=d1, "p(HUM)"=d2, "p(LEN)"=d3, "p(HUM+LEN)"=d4)
aictab(cand.set=det.mod, second.ord=F) 

### -----------------------------------------------------------------------  ###

# Then model occupancy
globalmod <- occu(~human+tlength ~dvillage+droad+driver+treecover+elevation,
                  data=sp.umf)

### -----------------------------------------------------------------------  ###

# Goodness of fit
globmod.gof <- mb.gof.test(globalmod, nsim=1000, plot.hist=T)
globmod.gof

### -----------------------------------------------------------------------  ###

# Try different combination of covariates (but avoid dredging (!))
nullmod <- occu(~1 ~1, data=sp.umf)
mod1 <- occu(~human+tlength ~dvillage+treecover, data=sp.umf)
mod2 <- occu(~human+tlength ~dvillage+droad+treecover, data=sp.umf)
mod3 <- occu(~human+tlength ~treecover+elevation, data=sp.umf)
mod4 <- occu(~human+tlength ~dvillage+driver+treecover, data=sp.umf)
mod5 <- occu(~human+tlength ~dvillage+droad+driver+treecover, data=sp.umf)
mod6 <- occu(~human+tlength ~dvillage+droad+elevation, data=sp.umf)

# Prepare model list 
fmL <- fitList("psi(VIL+TREE)p(HUM+LEN)"              = mod1, 
               "psi(VIL+ROAD)p(HUM+LEN)"              = mod2, 
               "psi(TREE+ELE)p(HUM+LEN)"              = mod3,
               "psi(VIL+RIV+TREE)p(HUM+LEN)"          = mod4, 
               "psi(VIL+ROAD+RIV)p(HUM+LEN)"          = mod5, 
               "psi(VIL+ROAD+ELE)p(HUM+LEN)"          = mod6,
               "psi(.)p(.)"                           = nullmod, 
               "psi(VIL+ROAD+RIV+TREE+ELE)p(HUM+LEN)" = globalmod)

# Rank models by AIC
ms <- modSel(fmL)
ms

### -----------------------------------------------------------------------  ###

# Coefficients of the top model
coef(mod3)
confint(mod3, type="state")
confint(mod3, type="det")

# Export coefficients of all models 
coef(ms)
expCoef <- as(ms, "data.frame")

### -----------------------------------------------------------------------  ###

# unmarked can calculate Empirical Bayes estimate and confidence interval of the 
# number/proportion of study sites occupied. 
# This method uses the estimate of psi to generate a vector of 0/1 binary variables z 
# under the assumption that z~Bernoulli(psi); the sum of the mode, 
# lower, and upper quantiles and CI on the number/proportion of occupied sites.
# Technically, z is binary random effect in what can be called a Bernoulli-Bernoulli 
# random effects (or mixed) model, so to estimate it we can use the ranef() function.

# Prediction (proportion of area occupied)
top.mod <- occu(~human+tlength ~treecover+elevation, data=sp.umf)
re <- ranef(top.mod)
EBUP <- bup(re, stat="mean")
CI <- confint(re, level=0.95)
rbind(PAO=c(Estimate=sum(EBUP), colSums(CI))/nrow(detn))

# Convert log-odds of occupancy into probabilities
backTransform(linearComb(top.mod, coefficients=c(1, 0, 0), "state"))
backTransform(linearComb(top.mod, coefficients=c(1, 1, 0), "state"))
backTransform(linearComb(top.mod, coefficients=c(1, 1, 1), "state"))

# Do the same for detection probability
backTransform(linearComb(top.mod, coefficients=c(1, 1, 1), "det"))

# True number of occupied sites
sum(bup(ranef(top.mod), stat="mode"))

# Estimate occupancy probability from all models
psiAll <- predict(fmL, type="state")
cbind(mean=mean(psiAll$Predicted),
      se=mean(psiAll$SE),
      lower=mean(psiAll$lower),
      upper=mean(psiAll$upper))

# Detection probability
pAll <- predict(fmL, type="det")
cbind(mean=mean(pAll$Predicted),
      se=mean(pAll$SE),
      lower=mean(pAll$lower),
      upper=mean(pAll$upper))

### -----------------------------------------------------------------------  ###

# Model averaging
# Make model list (all models within delta AIC 2)
occu.model.list <- list(occ.1=mod3,
                        occ.2=mod6)

occ.avg <- model.avg(occu.model.list, fit=T)
coef(occ.avg)
confint(occ.avg, level=0.95)

### -----------------------------------------------------------------------  ###

# Predict the relationship between occupancy and a covariate

# Elevation 
# Set other covariates to 0. 
# Since we're using a standardised value, 0 is equivalent to setting 
# the raw cov at their mean.
newData <- data.frame(elevation=seq(-1.7, 1.2, by=0.1), treecover=0)
occ.prob <- predict(mod3, type="state", newdata=newData, appendData=T)
head(occ.prob)

plot(Predicted ~ elevation, occ.prob, type="l", lwd=2, ylim=c(0,1),
     xlab="Elevation (standardised)",
     ylab="Expected occupancy probability")                   # Mean
lines(lower ~ elevation, occ.prob, type="l", col=gray(0.5))   # Lower CI
lines(upper ~ elevation, occ.prob, type="l", col=gray(0.5))   # Upper CI

# Plot model-averaged prediction
newDataE1 <- data.frame(elevation=seq(-1.7, 1.2, by=0.1), 
                        treecover=0, dvillage=0, droad=0, driver=0)
occ.probE1 <- predict(fmL, type="state", newdata=newDataE1, appendData=T)
head(occ.probE1)

plot(Predicted ~ elevation, occ.probE1, type="l", lwd=2, ylim=c(0,1),
     xlab="Elevation (standardised)",
     ylab="Expected occupancy probability")
lines(lower ~ elevation, occ.probE1, type="l", col=gray(0.5))
lines(upper ~ elevation, occ.probE1, type="l", col=gray(0.5))

# Plot with x-axis on the original scale
plot(Predicted ~ elevation, occ.probE1, type="l", lwd=2, ylim=c(0,1), xaxt="n",
     xlab="Elevation (masl)",
     ylab="Expected occupancy probability")
xticks <- -1.5:1
xlabs <- xticks*sd(dat$elevation) + mean(dat$elevation)
axis(1, at=xticks, labels=round(xlabs, 1))
lines(lower ~ elevation, occ.probE1, type="l", col=gray(0.5))
lines(upper ~ elevation, occ.probE1, type="l", col=gray(0.5))

# Detection probability vs human encounter frequency
he <- seq(0, 5, length=51)
hes <- (he - mean(as.matrix(human)))/sd(as.matrix(human))
newdataD <- data.frame(human=hes, tlength=0)
ep <- predict(top.mod, type="det", newdata=newdataD, appendData=T)
ep$humanOrig <- seq(0, 5, length=51)
with(ep, {plot(humanOrig, Predicted, type="l", lwd=2, ylim=c(0, 1),
               xlab="Human encounter frequency", ylab="p", main="") 
  lines(humanOrig, lower, type="l", col="grey")
  lines(humanOrig, upper, type="l", col="grey")
}
)

### -----------------------------------------------------------------------  ###

# Occupancy map (detection-corrected species distribution)
library(raster)

# Load raster files
# All rasters must have same extent and projection
vil <- raster("DistVillage1.tif")
roa <- raster("DistRoad1.tif")
riv <- raster("DistRiver.tif")
frst <- raster("Treecover.tif")
ele <- raster("Elevation.tif")

# Standardise raster with the same covariate value used in modelling
vilS <- (vil - mean(dat$dvillage)) / sd(dat$dvillage)
roaS <- (roa - mean(dat$droad)) / sd(dat$droad)
#rivS <- (riv - mean(dat$driver)) / sd(dat$driver)
forS <- (frst - mean(dat$treecover)) / sd(dat$treecover)
eleS <- (ele - mean(dat$elevation)) / sd(dat$elevation)

# Extract coefficients
beta <- coef(occ.avg)

# Plug into the linear model with raster covariates
logitPsi <- beta[1] + (beta[2]*forS) + (beta[3]*eleS) + (beta[7]*vilS) + (beta[8]*roaS)

# Convert to the probability scale
psi <- exp(logitPsi)/(1+exp(logitPsi))   

# Plot the map
plot(psi, main="Chimpanzee occupancy 2016 (dry season)")

# Mean occupancy probability across all pixels (345268)
mean(values(psi), na.rm=T)  # check the estimate with PAO estimate
sd(values(psi), na.rm=T) 
quantile(values(psi), probs = c(0.25, 0.75), na.rm=T)

cellStats(psi, stat="mean", na.rm=T)  # using a function in raster package

### -----------------------------------------------------------------------  ###

# For those curious!

# Compare with conventional glm

# Collapse detection history matrix into single-column detection data
y.red <- apply(detn, 1, max, na.rm=T)

# Do the same for human encounter data
humanSum <- apply(human, 1, sum)
humanSum <- scale(humanSum)

# Bundle data and rename columns 
c2016 <- cbind(y.red, dvilS, droadS, drivS, treeS, elevS, humanSum, tlength)
colnames(c2016) <- c("detection", "dvillage", "droad", "driver", "treecover", 
                     "elevation", "humanEnc", "tlength")
head(c2016)

# Do conventional glm (global model)
glm.mod <- glm(detection ~ dvillage + droad + driver + treecover + elevation +
                 humanSum + tlength, data=c2016, family=binomial(link="logit"))
summary(glm.mod)

# To make it comparable with SSOM
glm.mod1 <- glm(detection ~ treecover + elevation + dvillage + droad, 
                data=c2016, family=binomial(link="logit"))
summary(glm.mod1)

# Produce species distribution map
beta.glm <- coef(glm.mod1)

# Linear model for prediction
logit.glm <- beta.glm[1] + (beta.glm[2]*forS) + (beta.glm[3]*eleS) + 
             (beta.glm[4]*vilS) + (beta.glm[5]*roaS)

# Convert to probability scale
glm.occ <- exp(logit.glm)/(1+exp(logit.glm))
plot(glm.occ, main="GLM")

mean(values(glm.occ), na.rm=T)  # compare with naive estimate and psi

# Compare ssom and glm estimates
round(cbind(SSOM=mean(values(psi), na.rm=T),
            GLM=mean(values(glm.occ), na.rm=T),
            NAIVE=naive.occ), 2)

# Compare the prediction maps from SSOM and GLM
par(mfrow=c(1, 2))
plot(psi, main="Occupancy model")
plot(glm.occ, main="GLM")
par(mfrow=c(1, 1))

# Bonus
# Plot the relationship between (apparent) occupancy and covariate

library(dplyr)

# I will demonstrate for elevation
# Prepare new data frame for prediction 
newdata <- data.frame(elevation=seq(-1.7, 1.2, by=0.1), treecover=0, dvillage=0, droad=0)

# Predict with lower and upper CIs
pred.glm <- predict(glm.mod1, type="link", newdata=newdata, se.fit=TRUE) %>%
  as.data.frame() %>%
  mutate(elevation=seq(-1.7, 1.2, by=0.1),
         lower=glm.mod1$family$linkinv(fit - 1.96*se.fit),
         mean=glm.mod1$family$linkinv(fit),
         upper=glm.mod1$family$linkinv(fit + 1.96*se.fit))

# Do the plot
plot(mean ~ elevation, pred.glm, type="l", lwd=2, ylim=c(0,1),
     xlab="Elevation (standardised)", ylab="Apparent occupancy probability")
lines(lower ~ elevation, pred.glm, type="l", col=gray(0.5))
lines(upper ~ elevation, pred.glm, type="l", col=gray(0.5))

# Using visreg package (produces beautiful ggplots)
visreg(fit=glm.mod1, xvar="elevation", scale="response", 
       ylab="Habitat use", xlab="Elevation (standardised)")


############################################################################
############################################################################

# AUC
# August 2023

# For all models within delta AIC 2
# Essentially your model-average

# Predict occupancy probability first
# Create data frame of all covariates in models within delta AIC 2
occ.pred <- predict(d2modList, 
                    newdata=data.frame(cbind("treecover" = SiteCov$treecover, 
                                             "elevation" = SiteCov$elevation,
                                             "dist_road" = siteCov$dist_road,
                                             "dist_vil"  = siteCov$dist_vil,
                                             "dist_riv"  = siteCov$dist_riv)), 
                    type="state")

# Extract required parameters to plot
psi <- occ.pred$Predicted

# A model list object - note here you have all the models within delta AIC 2 
# or any threshold of your choice
modList <- list(mod4, mod5, mod2, mod6, mod3)

# Create a function to incorporate ranef() and bup() functions to model list object
apply_bup <- function(model){
  random_effects <- ranef(model)
  bup(random_effects, stat="mean")
}

# Apply the custom function to your model list
bup_list <- purrr::map(modList, apply_bup)
str(bup_list)

# Calculate the mean from all the objects in the list
mean_b <- purrr::reduce(bup_list, `+`) / length(bup_list)
z <- round(mean_b)   # rounding off to 1s and 0s

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# For a single model - perhaps your top or best model
top.mod <- occu(~scale(effort) ~scale(elevation)+scale(treecover), data=sp.umf)
re <- ranef(top.mod)
str(re)
EBUP <- bup(re, stat="mean")
round(EBUP)
z <- round(EBUP)

# Predict occupancy probability
occ.pred <- predict(top.mod, 
                    newdata=data.frame(cbind("treecover"=SiteCov$treecover, 
                                             "elevation"=SiteCov$elevation)), 
                    type="state", inf.rm=T)

# Extract required parameters to plot
psi <- occ.pred$Predicted

################################################################################

# Make a plot
library(pROC)
proc.obj <- roc(z, psi, smoothed=T, 
                ci=T, ci.alpha=0.9, stratified=F, legacy.axes=FALSE,
                plot=T, auc.polygon=T, max.auc.polygon=T, grid=T, print.auc=T, show.thres=T)

##################################### END ##################################
