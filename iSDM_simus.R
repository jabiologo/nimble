set.seed(3)
library(terra)
library(sf)
library(dplyr)
library(dismo)
library(rgeos)


sarea <- raster(nrows = 100, ncols = 100, xmn = 0, xmx = 100, ymn = 0, ymx = 100)
# Distance to water point covariate
dwat <- scale(distanceFromPoints(sarea, c(30,30)))
# Tree cover covariate
tree <- raster(nrows = 10, ncols = 10, xmn = 0, xmx = 100, ymn = 0, ymx = 100)
tree[] <- runif(100, 1,10)
tree <- scale(disaggregate(tree,10, "bilinear"))

# Tree cover covariate
autoc <- raster(nrows = 10, ncols = 10, xmn = 0, xmx = 100, ymn = 0, ymx = 100)
autoc[] <- runif(100, 1,10)
autoc <- scale(disaggregate(autoc,10, "bilinear"))

# Lambda parameter for the Poisson distribution of the abundance will be 
# function from "distance to water point" and "tree cover" with the following
# coefficients
beta0 <- 2
beta1 <- -0.5
beta2 <- 0.3
beta3 <- 0.4
lambda <- exp(beta0 + beta1*(dwat) + beta2*(tree) + beta3*(autoc))

# Now we can fill each cell of our study area with a random number from a 
# Poisson distribution with a different lambda at each site/cell (IPPP)
for (i in 1:ncell(sarea)){
  sarea[i] <- rpois(1, lambda[i])
}

# Plot the different variables and the study area
par(mfrow = c(2,2))
plot(dwat, main = "Distance to the water point")
plot(tree, main = "Tree cover")
plot(autoc, main = "Spatial noise")
plot(lambda, main = "Lambda parameter of the IPPP")
plot(sarea, main = "Animal abundance per cell")

# Hunting ground simulation. Firstly, we create some municipalities using 
# random points and Voronoi polygons
l1 <- randomPoints(sarea, 150)
muni <- crop(voronoi(l1, ext = extent(sarea)),sarea)

# Here we transform our raster in polygons to work with dataframes
spoly <- rasterToPolygons(sarea)
names(spoly) <- "animals"
spoly$cellId <- 1:nrow(spoly)

# This ugly chunk is to assign each cell to a polygon
spoints <- rasterToPoints(sarea)
ex <- terra::extract(muni,spoints[,1:2])
ex <- ex[!duplicated(ex$id.y),]
spoly$muni <- ex$id

spoly_<-st_as_sf(spoly) %>% arrange(muni) %>% mutate(cellIdNew=1:nrow(spoly))

spoly_$dwat <- dwat[spoly_$cellId]
spoly_$tree <- tree[spoly_$cellId]
muni_ <- spoly_ %>% group_by(muni) %>% dplyr::summarise(animals=sum(animals),
                                                        minId=min(cellIdNew), 
                                                        maxId=max(cellIdNew),
                                                        dwat=mean(dwat),
                                                        tree=mean(tree))

plot(muni_["animals"])

muni_$press <- rnorm(nrow(muni_),0.35,0.1)
muni_$hy <- round(muni_$animals * muni_$press, 0)

plot(muni_["animals"])
plot(muni_["press"])
plot(muni_["hy"])

noise1 <- muni_$hy + rnorm(n= length(muni_$hy), mean = 0, sd = 500)
muni_$effCov1 <- scale(muni_$hy + rnorm(n= length(muni_$hy), mean = 0, sd = 500))
muni_$effCov2 <- scale(muni_$hy + rnorm(n= length(muni_$hy), mean = 0, sd = 200))
cor(muni_$hy, muni_$effCov1) # ~ 0.4
cor(muni_$hy, muni_$effCov2) # ~ 0.8

# For some municipalities we'll have a kind of good measure of animal abundance
set.seed(2)
nSamp <- 20
muni_sampled <- muni_ %>% sample_n(nSamp)
# Adding some noise
muni_sampled$abuSamp <- muni_sampled$animals + rnorm(nrow(muni_sampled),0,250)
muni_sampled$abuSamp[muni_sampled$abuSamp < 0] <- 0
cor(muni_sampled$animals, muni_sampled$abuSamp) # ~ 0.85

# Now we'll simulate for some cells latrine counts by using
# the equations in Cabezas and Virgós 2023
# y = a * (1-e^(-bx)) 
# log(abundance) = beta0L + beta1L*log(letrinas/km)
# Approximated values from Cabezas and Virgos
#letkm <- runif(20,1,100)
#abun <- exp(-7.5 + (2.2*log(letkm*2)))
#plot(abun, letkm, ylab = "let/2 km")
# From abundance to latrines
#exp((7.5+log(abun))/2.2)/2

# Selecting cells for latrine simulation
latSamp <- data.frame(ID = sample(1:ncell(sarea),100), latC = NA)
# Adding some noise
latSamp$latC <- (exp((7.5+log(extract(sarea,latSamp$ID)))/2.2)/2) +
  rnorm(nrow(latSamp),0,4)
latSamp$latC[latSamp$latC < 0] <- 0

# Plot the general issue
plot(sarea)
lines(muni_, lwd = 0.5)
lines(muni_sampled, lwd = 3, col = "darkblue")
points(xyFromCell(sarea,latSamp$ID), pch = 18)


















