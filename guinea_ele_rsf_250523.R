
# Guinea elephant RSF
# May 2023

# Load packages
pcks <- c("move", "raster", "sf", "ggplot2", "plyr", "scales", "lubridate", "dplyr")
a <- sapply(pcks, require, character=TRUE)

# Load movement data, convert to simple feature
LJHRdata <- read.csv("20230504_DataP_LJHR.csv")
str(LJHRdata)

# Bind coordinates
datapointsLJHR <- SpatialPoints(cbind(LJHRdata$POINT_X,LJHRdata$POINT_Y))
crs(datapointsLJHR)

#3. Combining raster layers
elevation <- raster("elevation.tif")
treeheight <- raster("treecover2000_Hansen.tif")
rivers <- raster("riversorder2_distance.tif")
populated <- raster("OSMpopulated_distance.tif")
MODIS <- raster("MODIS.tif")
beehives <- raster("Beehives_kern_5000.tif")
dist_beehive <- raster("distance_to_beehive.tif")
# Project distance to beehive raster to match others
dist_beehive <- projectRaster(dist_beehive, crs=st_crs(treeheight)$proj4string)

# Make all rasters same size and extent
elevation.matched <- resample(elevation, treeheight) # , method="bilinear"
MODIS.matched <- resample(MODIS, treeheight)
rivers.matched <- resample(rivers, treeheight)
populated.matched <- resample(populated, treeheight)
beehives.matched <- resample(beehives, treeheight)
distbeehive.matched <- resample(dist_beehive, treeheight)

# Create raster brick
RSF.brick <- brick(elevation.matched, MODIS.matched, rivers.matched, populated.matched, 
                   beehives.matched, distbeehive.matched, treeheight)

# Plot to check
plot(RSF.brick[[7]])
plot(datapointsLJHR, col=alpha("blue", 0.25), pch=19, cex=0.5, add=TRUE)

# Compute MCP
ele.mcp <- adehabitatHR::mcp(datapointsLJHR, 100)
lines(ele.mcp, col="red", lwd=2)

# Random points within MCP
xy.obs <- LJHRdata[sample(1:nrow(LJHRdata), 200), c("POINT_X", "POINT_Y")]
xy.random <- spsample(ele.mcp, 1000, "random")@coords

plot(xy.random, asp=1, col=adjustcolor("darkblue", alpha=0.5), pch=19, cex=0.5)
points(xy.obs, pch=19, col="orange", cex=0.5)

# Stack used and available points and label them
data.rsf <- data.frame(X=c(xy.obs[, 1], xy.random[, 1]), Y=c(xy.obs[, 2], xy.random[, 2]),
                       Used=c(rep(TRUE, nrow(xy.obs)), rep(FALSE, nrow(xy.random))))

# Extract variables for used and available points
data.rsf$elevation <- raster::extract(RSF.brick[[1]], data.rsf[, 1:2])
data.rsf$Modis <- raster::extract(RSF.brick[[2]], data.rsf[, 1:2])
data.rsf$river <- raster::extract(RSF.brick[[3]], data.rsf[, 1:2])
data.rsf$settlement <- raster::extract(RSF.brick[[4]], data.rsf[, 1:2])
data.rsf$beehive_den <- raster::extract(RSF.brick[[5]], data.rsf[, 1:2])
data.rsf$beehive_dist <- raster::extract(RSF.brick[[6]], data.rsf[, 1:2])
data.rsf$treecover <- raster::extract(RSF.brick[[7]], data.rsf[, 1:2])

# Visualise used vs available points
boxplot(elevation ~ Used, data=data.rsf, main="Elevation")
boxplot(Modis ~ Used, data=data.rsf, main="Modis")
boxplot(river ~ Used, data=data.rsf, main="river")
boxplot(settlement ~ Used, data=data.rsf, main="settlement")
boxplot(beehive_den ~ Used, data=data.rsf, main="beehive_den")
boxplot(beehive_dist ~ Used, data=data.rsf, main="beehive_dist")
boxplot(treecover ~ Used, data=data.rsf, main="treecover")

# Fit an RSF model
rsf.fit <- glm(Used ~ scale(elevation) + scale(river) + scale(settlement) + scale(beehive_den) + scale(treecover),
               data=data.rsf, family="binomial")
summary(rsf.fit)

rsf.fit1 <- glm(Used ~ scale(elevation) + scale(river) + scale(settlement) + scale(sqrt(beehive_dist)) + scale(treecover),
               data=data.rsf, family="binomial")
summary(rsf.fit1)

require(sjPlot)
plot_model(rsf.fit1)

with(data.rsf, plot(elevation, settlement, col=factor(Used), cex=0.5, pch=21))
with(data.rsf, plot(beehive_dist, settlement, col=factor(Used), cex=0.5, pch=21))
with(data.rsf, plot(river, treecover, col=factor(Used), cex=0.5, pch=21))
with(data.rsf, plot(settlement, treecover, col=factor(Used), cex=0.5, pch=21))
with(data.rsf, plot(beehive_dist, treecover, col=factor(Used), cex=0.5, pch=21))

# Might want to check interaction terms if you think it is needed based on the plots above

rsf.fit2 <- glm(Used ~ scale(elevation) + scale(river) + scale(settlement) + scale(beehive_dist) + scale(treecover) +
                scale(settlement) * scale(beehive_dist),
                data=data.rsf, family="binomial")
sjPlot::plot_model(rsf.fit2)
summary(rsf.fit2)

# Positive interaction term indicates that further the distance away from beehive, 
# the greater the effect (+ve) of settlement on elephant resource use, meaning, elephants 
# tend to avoid areas closer to human settlement and beehives. 

(exp(0.77)-1)*100 # when dist to beehive = 0, dist. to sett is associated with x% decrease in resource use
(exp(0.45)-1)*100 # for every unit increase in dist to beehive, dist to sett effect increases by (interaction coef)

# Plot coefficients
require(broom); require(dplyr); require(ggplot2)
rsf.coef <- tidy(rsf.fit2, conf.int=T) |>
  mutate(type=ifelse(grepl("", term), "Continuous covariate"),
         term=gsub("", "", term),
         significance=cut(p.value, c(0, .01, .05, .1, 1))) |>
  mutate(term=factor(term, levels=term[order(estimate)]))
rsf.coef$term
  
coef_plot <- ggplot(rsf.coef, aes(estimate, term, xmin=conf.low, xmax=conf.high, col=significance)) +
  geom_errorbarh(height=.2) + 
  geom_point() + 
  geom_vline(xintercept=0, lty=3, lwd=1) + 
  xlab(expression(beta~coefficients)) + ylab("") +
  scale_y_discrete(limit=c("(Intercept)", "scale(elevation)", "scale(beehive_dist)", "scale(river)", 
                           "scale(treecover)", "scale(settlement)", "scale(settlement):scale(beehive_dist)"),
                   labels=c("Intercept", "Elevation", "Dist. beehive", "Dist. river", "Treecover", "Dist. settlement",
                            "Dist. settlement x Dist. beehive")) + # rename y-axis labels
  scale_colour_hue(l=40)

################################################################################

data.rsf2 <- data.rsf

# scale() changes the class of your columns
# Simply them to numeric vectors
# Use the dimension-stripping properties of c()

data.rsf2$elevation <- c(scale(data.rsf2$elevation))
data.rsf2$river <- c(scale(data.rsf2$river))
data.rsf2$settlement <- c(scale(data.rsf2$settlement))
data.rsf2$beehive_dist <- c(scale(data.rsf2$beehive_dist))
data.rsf2$treecover <- c(scale(data.rsf2$treecover))
str(data.rsf2)

rsf.fit3 <- glm(Used ~ elevation + river + settlement + beehive_dist + treecover +
                  settlement * beehive_dist, data=data.rsf2, family="binomial")
sjPlot::plot_model(rsf.fit3)
summary(rsf.fit3)

# Visualise RSF vs covariate 
library(visreg)

visreg(fit=rsf.fit3, xvar="elevation", 
       gg = TRUE, 
       scale="response") +
  labs(y = "Prob of use", 
       x = "Elevation (standardised)")

visreg(fit=rsf.fit3, xvar="river", scale="response", 
       xlab="Distance to river (standardised)", ylab="Prob of use")


visreg(fit=rsf.fit3, xvar="settlement", scale="response", 
       xlab="Distance to settlement (standardised)", ylab="Prob of use")


visreg(fit=rsf.fit3, xvar="beehive_dist", scale="response", 
       xlab="Distance to beehive (standardised)", ylab="Prob of use")

visreg(fit=rsf.fit3, xvar="treecover", scale="response", 
       xlab="Treecover (standardised)", ylab="Prob of use")

# Plot interaction
visreg(fit=rsf.fit3, xvar="beehive_dist", by="settlement", overlay=T, scale="response", 
       xlab="Distance to beehive (standardised)", ylab="Prob of use")

visreg(fit=rsf.fit3, xvar="settlement", by="beehive_dist", overlay=T, scale="response", 
       xlab="Distance to beehive (standardised)", ylab="Prob of use")

# 2D filled contour plots
visreg2d(rsf.fit3, "settlement", "beehive_dist")


################################################################################

# Model comparison
AIC(rsf.fit, rsf.fit1, rsf.fit2)

# RSF prediction map
p.rsf <- RSF.brick[[1]]

ele.vec <- RSF.brick[[1]]@data@values
riv.vec <- RSF.brick[[3]]@data@values
set.vec <- RSF.brick[[4]]@data@values
bee.vec <- RSF.brick[[6]]@data@values
tre.vec <- RSF.brick[[7]]@data@values

predict.df <- data.frame(elevation=ele.vec, river=riv.vec, settlement=set.vec,
                         beehive_dist=bee.vec, treecover=tre.vec)
str(predict.df)

p.rsf@data@values <- predict(rsf.fit2, newdata=predict.df)

mycol <- colorRampPalette(c("grey", "yellow", "orange", "red", "darkred"))

# Convert to probability scale
rsf_prob <- exp(p.rsf)/(1+exp(p.rsf))
plot(rsf_prob, main="RSF", col=mycol(100))
plot(datapointsLJHR, col=alpha("blue", 0.25), pch=19, cex=0.5, add=TRUE)

#################################################################################

# Using the dredge() function

library(MuMIn)
options(na.action=na.fail)

#scale(log(river)) + 
rsf.fit <- glm(Used ~ scale(elevation) + scale(river) + scale(settlement) + scale(beehive_den) + scale(treecover),
               data=data.rsf, family="binomial")
summary(rsf.fit)

(fm_dredge <- dredge(rsf.fit, rank="AICc", trace=TRUE))

# Sum of Aikake weights (Summed Model Weights (SMW))
sw(fm_dredge)

# Get top-most models, but fitted by REML (subset only models with SMW < 2)
(dd <- subset(fm_dredge, delta < 2))

# Generate confidence intervals
confset.95p <- get.models(fm_dredge, subset=delta < 2)
avgmod.95p <- model.avg(confset.95p)
summary(avgmod.95p)
confint(avgmod.95p)

#################################################################################