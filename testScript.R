# Ugly R script for testing
set.seed(3)
library(terra)
library(sf)
library(dplyr)
library(dismo)
library(rgeos)

# Create a study area
sarea <- raster(nrows = 50, ncols = 50, xmn = 0, xmx = 50, ymn = 0, ymx = 50)
# Distance to water point covariate
dwat <- scale(distanceFromPoints(sarea, c(15,15)))
# Tree cover covariate
tree <- raster(nrows = 5, ncols = 5, xmn = 0, xmx = 50, ymn = 0, ymx = 50)
tree[] <- runif(25, 1,10)
tree <- scale(disaggregate(tree,10, "bilinear"))

# Lambda parameter for the Poisson distribution of the abundance will be 
# function from "distance to water point" and "tree cover" with the following
# coefficients
beta0 <- 2
beta1 <- -0.5
beta2 <- 0.3
lambda <- exp(beta0 + beta1*(dwat) + beta2*(tree))

# Now we can fill each cell of our study area with a random number from a 
# Poisson distribution with a different lambda at each site/cell (IPPP)
for (i in 1:ncell(sarea)){
  sarea[i] <- rpois(1, lambda[i])
}

# Plot the different variables and the study area
par(mfrow = c(2,2))
plot(dwat, main = "Distance to the water point")
plot(tree, main = "Tree cover")
plot(lambda, main = "Lambda parameter of the IPPP")
plot(sarea, main = "Animal abundance per cell")

# Firstly, we create some municipalities using random points and Voronoi 
# polygons
l1 <- randomPoints(sarea, 50)
muni <- crop(voronoi(l1, ext = extent(sarea)),sarea)

# Here we transform our raster in polygons to work with dataframes
spoly <- rasterToPolygons(sarea)
names(spoly) <- "animals"
spoly$cellId <- 1:nrow(spoly)

# This ugly chunk is to assign each cell to a polygon
spoints <- rasterToPoints(sarea)
ex <- terra::extract(muni,spoints[,1:2])
ex <- ex[!duplicated(ex$point.ID),]
spoly$muni <- ex$poly.ID

# Finally, another ugly code to name cell IDs in order inside each polygon
# To invoke Change of Support we need cells having an ascendant order inside
# each polygon, so cell IDs in polygon 1 will be from 1 to 9, cells in polygon
# 2 will be from 10 to 15, in polygon 3 will be from 16 to 22, etc.
# (Note that numbers are invented)
spoly_<-st_as_sf(spoly) %>% arrange(muni) %>% mutate(cellIdNew=1:nrow(spoly))
#plot(st_as_sf(spoly_)["cellId"])
#plot(st_as_sf(spoly_)["cellIdNew"])
spoly_$dwat <- dwat[spoly_$cellId]
spoly_$tree <- tree[spoly_$cellId]
# Here we count number of animals and we take the mean of predictor covariates
# for each municipality
muni_ <- spoly_ %>% group_by(muni) %>% dplyr::summarise(animals=sum(animals),
                                                        minId=min(cellIdNew), 
                                                        maxId=max(cellIdNew),
                                                        dwat=mean(dwat),
                                                        tree=mean(tree))
plot(sarea); lines(muni)
plot(muni_["animals"], pal = terrain.colors(13, rev = T))

######################################### HUNTING PRESSURE 1####################

# Including a categorical factor as a proxy of hunting pressure? (type of hunting ground)
# We will define here hunting pressure as the percentage of animals hunted in a season
# Following Vajas et al. 2021 (adapted)
# Hunting Pressure = catchability * effort


muni_$press <- runif(nrow(muni_),0,0.5)
muni_$hy <- round(muni_$animals * muni_$press, 0)

par(mfrow=c(1,3))
plot(muni_$animals, muni_$press, pch = 16, xlab = "ANIMALS", ylab = "HUNTING PRESSURE")
plot(muni_$animals, muni_$hy, pch = 16, xlab = "ANIMALS", ylab = "HUNTING YIELD")
plot(muni_$press, muni_$hy, pch = 16, xlab = "HUNTING PRESSURE", ylab = "HUNTING YIELD")
dev.off()
plot(muni_["animals"], pal = terrain.colors(13, rev = T), main = "Animals")
plot(muni_["press"],  pal = heat.colors(10, rev = T), main = "Hunting pressure")
plot(muni_["hy"], pal = terrain.colors(9, rev = T), main = "Hunting yield (animals * pressure)")


library(Hmisc)
library(fastDummies)
library(dplyr)
library(tidyverse)

library(nimble)
library(nimbleSCR)

ff <- cut2(muni_$press, g = 3)
levels(ff) <- c("1", "2", "3")

dd <- fastDummies::dummy_cols(ff) %>% as_tibble() #%>% 
  #dplyr::select(.data:.data_3) %>% 
names(dd) <- c("CatT", "Cat1", "Cat2", "Cat3")
muni_$catT <- as.numeric(as.character(dd$CatT))
muni_$cat1 <- dd$Cat1
muni_$cat2 <- dd$Cat2
muni_$cat3 <- dd$Cat3

plot(muni_["catT"], main = "Hunting ground Category")

constants <- list(ncell = nrow(spoly_),
                  nmuni = nrow(muni_),
                  low = muni_$minId,
                  high = muni_$maxId)

data <- list(animals = muni_$hy,
             dwat = spoly_$dwat,
             tree = spoly_$tree,
             cat1 = muni_$cat1,
             cat2 = muni_$cat2,
             cat3 = muni_$cat3)




simuCat <- nimble::nimbleCode( {
  # PRIORS
  
  b_intercept ~ dnorm(0, 2)
  b_dwat ~ dnorm(0, 2)
  b_tree ~ dnorm(0, 2)
  b_cat1 ~ dunif(0, 0.5)
  b_cat2 ~ dunif(0, 0.5)
  b_cat3 ~ dunif(0, 0.5)
  
  # LIKELIHOOD
  for(i in 1:ncell){
    log(lambda[i]) <- b_intercept + b_dwat*dwat[i] + b_tree*tree[i]
    n[i] ~ dpois(lambda[i])
  }
  
  # Sampling model. This is the part that changes respect the previous model
  # Here the counted animals per municipality is distributed following a
  # Poisson distribution with lambda = lambda_muni
  # lambda_muni is simply the sum of cell lambda in each municipality
  for(j in 1:nmuni){
    log(lambda_muni[j]) <- log(sum(lambda[low[j]:high[j]])) 
    animals[j] ~ dpois((b_cat1*cat1[j]*lambda_muni[j]) + (b_cat2*cat2[j]*lambda_muni[j]) + (b_cat3*cat3[j]*lambda_muni[j]))
  }
} )

# Once the model is defined, we should provide a function to get some random
# initial values for each of our parameters (sampled from an uniform 
# distribution, for example)

inits <- function() {
  base::list(n = rep(1, constants$ncell),
             b_intercept = runif(1, -1, 1),
             b_dwat = runif(1, -1, 1),
             b_tree = runif(1, -1, 1),
             b_cat1 = runif(1, 0, 1),
             b_cat2 = runif(1, 0, 1),
             b_cat3 = runif(1, 0, 1)
  )
}

# Set values we are interested in
keepers <- c("lambda", 'b_intercept', "b_dwat", "b_tree", "b_cat1", "b_cat2", "b_cat3")

# Finally we define the settings of our MCMC algorithm
nc <- 2 # number of chains
nb <- 1000 # number of initial MCMC iterations to discard
ni <- nb + 20000 # total number  of iterations

# Now he create the model
model <- nimble::nimbleModel(code = simuCat, 
                             data = data, 
                             constants = constants, 
                             inits = inits(),
                             calculate = FALSE)

# Check if everything is initialized (I understand this)
model$initializeInfo()

# Compile the model (I'm lost here. In general I understand, but I'm not able
# to modify any configuration right now)
c_model <- nimble::compileNimble(model)
model_conf <- nimble::configureMCMC(model,
                                    useConjugacy = FALSE)
model_conf$addMonitors(keepers)
model_mcmc <- nimble::buildMCMC(model_conf)
c_model_mcmc <- nimble::compileNimble(model_mcmc, project = model)

# Run the MCMC
samples <- nimble::runMCMC(c_model_mcmc, 
                           nburnin = nb, 
                           niter = ni, 
                           nchains = nc)

# We can use now the coda package to see MCMC results
samples_mcmc <- coda::as.mcmc.list(lapply(samples, coda::mcmc))

# Look at traceplots (3 chains) of the three parameters
par(mfrow=c(2,3))
coda::traceplot(samples_mcmc[, 1:6])
# Calculate Rhat convergence diagnostic for the three parameters
coda::gelman.diag(samples_mcmc[,1:6])

# extract mean for each parameter
samplesdf <- as.data.frame(rbind(samples_mcmc$chain1,samples_mcmc$chain2))
mValues <- colMeans(samplesdf)
# We can inspect the mean of posterior distributions for each parameter
# Remember that real values were: int=2; dwat=-0.5; tree=0.3
mValues[1:7]

# Now we can plot lambda predictions and SD for each cell
pred3 <- sarea
# Notice that we changed the cellID so we have to use old IDs
pred3[spoly_$cellId] <- mValues[7:length(mValues)]

par(mfrow = c(1,3))
plot(sarea, main = "Animal abundance per cell")
plot(pred3, main = "Predicted abundance per cell")
plot(pred3[], sarea[], pch = 16, cex = .8, xlab = "predicted", ylab = "observed", cex.lab = 1.5)
abline(a=1, b=1, col = "darkred", lwd = 2)
sum(sarea[])
sum(pred3[])


pressEst <- data.frame(cbind(c(samples_mcmc$chain2[,1],samples_mcmc$chain2[,2],samples_mcmc$chain2[,3]),
c(rep(1,length(samples_mcmc$chain2[,1])), rep(2,length(samples_mcmc$chain2[,1])), 
  rep(3,length(samples_mcmc$chain2[,1])))))
boxplot(X1 ~ X2, data = pressEst, xlab = "Hunting ground Category", ylab = "Estimated Pressure", cex.lab = 1.4,
        col = c("darkblue", "purple", "yellow"))

######################################### HUNTING PRESSURE 2###################
# Including a categorical factor as a proxy of hunting pressure? (type of hunting ground)
# We will define here hunting pressure as the percentage of animals hunted in a season
# Following Vajas et al. 2021 (adapted)
# Hunting Pressure = catchability * effort

muni_$press <- runif(nrow(muni_),0,0.5)
muni_$hy <- round(muni_$animals * muni_$press, 0)
muni_$area <- st_area(muni_)/1000000
muni_$hyDen <- muni_$hy/muni_$area

par(mfrow=c(1,3))
plot(muni_$animals, muni_$press, pch = 16, xlab = "ANIMALS", ylab = "HUNTING PRESSURE")
plot(muni_$animals, muni_$hy, pch = 16, xlab = "ANIMALS", ylab = "HUNTING YIELD")
plot(muni_$press, muni_$hy, pch = 16, xlab = "HUNTING PRESSURE", ylab = "HUNTING YIELD")
dev.off()
plot(muni_["animals"], pal = terrain.colors(13, rev = T), main = "Animals")
plot(muni_["press"],  pal = heat.colors(10, rev = T), main = "Hunting pressure")
plot(muni_["hy"], pal = terrain.colors(8, rev = T), main = "Hunting yield (animals * pressure)")


library(Hmisc)
library(fastDummies)
library(dplyr)
library(tidyverse)

library(nimble)
library(nimbleSCR)

# Let's fit hunting ground categories from hg rather than press
ff <- cut2(as.numeric(muni_$hyDen), g = 3)
levels(ff) <- c("1", "2", "3")

dd <- fastDummies::dummy_cols(ff) %>% as_tibble() #%>% 
#dplyr::select(.data:.data_3) %>% 
names(dd) <- c("CatT", "Cat1", "Cat2", "Cat3")
muni_$catT <- as.numeric(as.character(dd$CatT))
muni_$cat1 <- dd$Cat1
muni_$cat2 <- dd$Cat2
muni_$cat3 <- dd$Cat3

plot(muni_["catT"], main = "Hunting ground Category")

constants <- list(ncell = nrow(spoly_),
                  nmuni = nrow(muni_),
                  low = muni_$minId,
                  high = muni_$maxId)

data <- list(animals = muni_$hy,
             dwat = spoly_$dwat,
             tree = spoly_$tree,
             cat1 = muni_$cat1,
             cat2 = muni_$cat2,
             cat3 = muni_$cat3)

# 1) Vectorización
# log(lambda[1:ncell]) <- b_intercept + b_dwat*dwat[1:ncell] + b_tree*tree[1:ncell]
# dbinom_vector
# 2) Change parametrization following Olivier advice

dpois_vector <- function(){
  dbinom_vector(rep(5,10),size = rep(10,10), prob = rep(0.5,10))
}


simuCat <- nimble::nimbleCode( {
  # PRIORS
  
  b_intercept ~ dnorm(0, 2)
  b_dwat ~ dnorm(0, 2)
  b_tree ~ dnorm(0, 2)
  b_cat1 ~ dunif(0, 0.5)
  b_cat2 ~ dunif(0, 0.5)
  b_cat3 ~ dunif(0, 0.5)
  
  # LIKELIHOOD
  for(i in 1:ncell){
    log(lambda[i]) <- b_intercept + b_dwat*dwat[i] + b_tree*tree[i]
    n[i] ~ dpois(lambda[i])
  }
  
  # Sampling model. This is the part that changes respect the previous model
  # Here the counted animals per municipality is distributed following a
  # Poisson distribution with lambda = lambda_muni
  # lambda_muni is simply the sum of cell lambda in each municipality
  for(j in 1:nmuni){
    log(lambda_muni[j]) <- log(sum(lambda[low[j]:high[j]])) 
    animals[j] ~ dpois((b_cat1*cat1[j]*lambda_muni[j]) + (b_cat2*cat2[j]*lambda_muni[j]) + (b_cat3*cat3[j]*lambda_muni[j]))
  }
} )

# Once the model is defined, we should provide a function to get some random
# initial values for each of our parameters (sampled from an uniform 
# distribution, for example)

inits <- function() {
  base::list(n = rep(1, constants$ncell),
             b_intercept = runif(1, -1, 1),
             b_dwat = runif(1, -1, 1),
             b_tree = runif(1, -1, 1),
             b_cat1 = runif(1, 0, 1),
             b_cat2 = runif(1, 0, 1),
             b_cat3 = runif(1, 0, 1)
  )
}

# Set values we are interested in
keepers <- c("lambda", 'b_intercept', "b_dwat", "b_tree", "b_cat1", "b_cat2", "b_cat3")

# Finally we define the settings of our MCMC algorithm
nc <- 2 # number of chains
nb <- 1000 # number of initial MCMC iterations to discard
ni <- nb + 20000 # total number  of iterations

# Now he create the model
model <- nimble::nimbleModel(code = simuCat, 
                             data = data, 
                             constants = constants, 
                             inits = inits(),
                             calculate = FALSE)

# Check if everything is initialized (I understand this)
model$initializeInfo()

# Compile the model (I'm lost here. In general I understand, but I'm not able
# to modify any configuration right now)
c_model <- nimble::compileNimble(model)
model_conf <- nimble::configureMCMC(model,
                                    useConjugacy = FALSE)
model_conf$addMonitors(keepers)
model_mcmc <- nimble::buildMCMC(model_conf)
c_model_mcmc <- nimble::compileNimble(model_mcmc, project = model)

# Run the MCMC
samples <- nimble::runMCMC(c_model_mcmc, 
                           nburnin = nb, 
                           niter = ni, 
                           nchains = nc)

# We can use now the coda package to see MCMC results
samples_mcmc <- coda::as.mcmc.list(lapply(samples, coda::mcmc))

# Look at traceplots (3 chains) of the three parameters
par(mfrow=c(2,3))
coda::traceplot(samples_mcmc[, 1:6])
# Calculate Rhat convergence diagnostic for the three parameters
coda::gelman.diag(samples_mcmc[,1:6])

# extract mean for each parameter
samplesdf <- as.data.frame(rbind(samples_mcmc$chain1,samples_mcmc$chain2))
mValues <- colMeans(samplesdf)
# We can inspect the mean of posterior distributions for each parameter
# Remember that real values were: int=2; dwat=-0.5; tree=0.3
mValues[1:7]

# Now we can plot lambda predictions and SD for each cell
pred3 <- sarea
# Notice that we changed the cellID so we have to use old IDs
pred3[spoly_$cellId] <- mValues[7:length(mValues)]

par(mfrow = c(1,3))
plot(sarea, main = "Animal abundance per cell")
plot(pred3, main = "Predicted abundance per cell")
plot(pred3[], sarea[], pch = 16, cex = .8, xlab = "predicted", ylab = "observed", cex.lab = 1.5)
abline(a=1, b=1, col = "darkred", lwd = 2)
sum(sarea[])
sum(pred3[])


pressEst <- data.frame(cbind(c(samples_mcmc$chain2[,1],samples_mcmc$chain2[,2],samples_mcmc$chain2[,3]),
                             c(rep(1,length(samples_mcmc$chain2[,1])), rep(2,length(samples_mcmc$chain2[,1])), 
                               rep(3,length(samples_mcmc$chain2[,1])))))
boxplot(X1 ~ X2, data = pressEst, xlab = "Hunting ground Category", ylab = "Estimated Pressure", cex.lab = 1.4,
        col = c("darkblue", "purple", "yellow"))


###################### MVNormal Simulation #####################################
library(MASS)
library(matrixcalc)
library(lqmm)
# Create a study area
sarea <- raster(nrows = 4, ncols = 4, xmn = 0, xmx = 4, ymn = 0, ymx = 4)
# Distance to water point covariate
dwat <- scale(distanceFromPoints(sarea, c(0.5,0.5)))
coords <- xyFromCell(dwat,1:16)
dd <- dist(coords)

Cmat <- function(sigma2, distance, rho){
  C <- sigma2 * exp(-distance/rho)
  return(C)
}

C <- round(Cmat(0.1, dd, 1),3)

MVN <- mvrnorm(n = 1,mu = rep(0,16),Sigma = CCC)

ccc <- dwat
ccc[] <- MVN

lambda <- dwat

lambda[1:16] <- exp((0.3*dwat[1:16]) + MVN[1:16])
plot(lambda)

####################

expcov <- nimbleFunction(     
  run = function(dists = double(2), rho = double(0), sigma = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    sigma2 <- sigma*sigma   # calculate once
    for(i in 1:n)
      for(j in 1:n)
        result[i, j] <- sigma2*exp(-dists[i,j]/rho)
    return(result)
  })



expcov <- function(dists, rho, sigma){
  n <- dim(dists)[1]
  result <- matrix(nrow = n, ncol = n)
  sigma2 <- sigma*sigma
  for(i in 1:n){
    for(j in 1:n){
      result[i, j] <- sigma2*exp(-dists[i,j]/rho)
    }
  }
  return(result)
}


cov <- expcov(ddd, 1, 1)
N <- 4
mu <- c(1,1,1,1)
x[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])

dmnorm_chol(mu[1:N], dmnorm_chol = FALSE)

mvrnorm(n = 1,mu = rep(1,16),Sigma = cov)
