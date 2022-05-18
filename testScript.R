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

################################################################################

library(nimble)
library(nimbleSCR)

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
# To invoque Change of Support we need cells having an ascendent order inside
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

################################################################################
# Adding effort and hunting yields

muni_$eff <- runif(nrow(muni_),0,0.3)
muni_$hy <- round(muni_$animals * muni_$eff, 0)

plot(muni_$animals, muni_$eff)
cor(muni_$animals, muni_$eff)

plot(muni_$animals, muni_$hy)
cor(muni_$animals, muni_$hy)

plot(sarea)
lines(muni)
plot(muni_["animals"])
plot(muni_[c(2,9)])

constants <- list(ncell = nrow(spoly_),
                  nmuni = nrow(muni_),
                  low = muni_$minId,
                  high = muni_$maxId)

data <- list(animals = muni_$hy,
             dwat = spoly_$dwat,
             tree = spoly_$tree)

simuCoS <- nimble::nimbleCode( {
  # PRIORS
  
  b_intercept ~ dnorm(0, 2)
  b_dwat ~ dnorm(0, 2)
  b_tree ~ dnorm(0, 2)
  
  # LIKELIHOOD
  for(i in 1:ncell){
    log(lambda[i]) <- b_intercept + b_dwat*dwat[i] + b_tree*tree[i]
    n[i] ~ dpois(lambda[i])
  }
  
  # Sampling model. This is the part that changes respect the previous model
  # Here the counted animals per municipality is distributed following a
  # Poisson distribution with lambda = lamnda_muni
  # lambda_muni is simpy the summatory of cell lambda in each municipality
  for(j in 1:nmuni){
    log(lambda_muni[j]) <-log(sum(lambda[low[j]:high[j]])) 
    animals[j] ~ dpois(lambda_muni[j])
  }
} )

# Once the model is defined, we should provide a function to get some random
# initial values for each of our parameters (sampled from an uniform 
# distribution, for example)

inits <- function() {
  base::list(n = rep(1, constants$ncell),
             b_intercept = runif(1, -1, 1),
             b_dwat = runif(1, -1, 1),
             b_tree = runif(1, -1, 1)
  )
}

# Set values we are interested in
keepers <- c("lambda", 'b_intercept', "b_dwat", "b_tree")

# Finally we define the settings of our MCMC algorithm
nc <- 2 # number of chains
nb <- 1000 # number of initial MCMC iterations to discard
ni <- nb + 20000 # total number  of iterations

# Now he create the model
model <- nimble::nimbleModel(code = simuCoS, 
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
par(mfrow=c(1,3))
coda::traceplot(samples_mcmc[, 1:3])
# Calculate Rhat convergence diagnostic for the three parameters
coda::gelman.diag(samples_mcmc[,1:3])

# extract mean for each parameter
samplesdf <- as.data.frame(rbind(samples_mcmc$chain1,samples_mcmc$chain2))
mValues <- colMeans(samplesdf)
# We can inspect the mean of posterior distributions for each parameter
# Remember that real values were: int=2; dwat=-0.5; tree=0.3
mValues[1:4]

# Now we can plot lambda predictions and SD for each cell
pred1 <- sarea
# Notice that we changed the cellID so we have to use old IDs
pred1[spoly_$cellId] <- mValues[4:length(mValues)]

par(mfrow = c(1,3))
plot(sarea, main = "Animal abundance per cell")
plot(pred1, main = "Predicted abundance per cell")
plot(pred1[], sarea[], pch = 16, cex = .8)
abline(a=1, b=1, col = "darkred", lwd = 2)
sum(sarea[])
sum(pred1[])

# Including effort in modeling as a constant?

simuEff <- nimble::nimbleCode( {
  # PRIORS
  
  b_intercept ~ dnorm(0, 2)
  b_dwat ~ dnorm(0, 2)
  b_tree ~ dnorm(0, 2)
  effort ~ dunif(0, 0.5) # We think that effort could be somewhere between 0
                         # and 0.5, so we set the prior...
  
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
    log(lambda_muni[j]) <-log(sum(lambda[low[j]:high[j]])) 
    animals[j] ~ dpois(effort * lambda_muni[j]) # effort as a constant?
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
             effort = runif(1, 0, 0.5)
  )
}

# Set values we are interested in
keepers <- c("lambda", 'b_intercept', "b_dwat", "b_tree", "effort")

# Finally we define the settings of our MCMC algorithm
nc <- 2 # number of chains
nb <- 1000 # number of initial MCMC iterations to discard
ni <- nb + 20000 # total number  of iterations

# Now he create the model
model <- nimble::nimbleModel(code = simuEff, 
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
par(mfrow=c(1,3))
coda::traceplot(samples_mcmc[, 1:3])
# Calculate Rhat convergence diagnostic for the three parameters
coda::gelman.diag(samples_mcmc[,1:3])

# extract mean for each parameter
samplesdf <- as.data.frame(rbind(samples_mcmc$chain1,samples_mcmc$chain2))
mValues <- colMeans(samplesdf)
# We can inspect the mean of posterior distributions for each parameter
# Remember that real values were: int=2; dwat=-0.5; tree=0.3
mValues[1:4]

# Now we can plot lambda predictions and SD for each cell
pred2 <- sarea
# Notice that we changed the cellID so we have to use old IDs
pred2[spoly_$cellId] <- mValues[5:length(mValues)]

par(mfrow = c(1,3))
plot(sarea, main = "Animal abundance per cell")
plot(pred2, main = "Predicted abundance per cell")
plot(pred2[], sarea[], pch = 16, cex = .8)
abline(a=1, b=1, col = "darkred", lwd = 2)
sum(sarea[])
sum(pred2[])

par(mfrow = c(1,2))
plot(pred1[], sarea[], pch = 16, cex = .8)
abline(a=1, b=1, col = "darkred", lwd = 2)
plot(pred2[], sarea[], pch = 16, cex = .8)
abline(a=1, b=1, col = "darkred", lwd = 2)

# Including a categorical factor as a proxy of effort? (type of hunting ground)

muni_$eff <- runif(nrow(muni_),0,0.5)
muni_$hy <- round(muni_$animals * muni_$eff, 0)

plot(muni_$animals, muni_$eff, pch = 16, xlab = "ANIMALS", ylab = "HUNTING EFFORT")
cor(muni_$animals, muni_$eff)

plot(muni_$animals, muni_$hy, pch = 16, xlab = "ANIMALS", ylab = "HUNTING YIELD")
cor(muni_$animals, muni_$hy)

plot(muni_$eff, muni_$hy, pch = 16, xlab = "HUNTING EFFORT", ylab = "HUNTING YIELD")
cor(muni_$eff, muni_$hy)

plot(sarea)
lines(muni)
plot(muni_["animals"])
plot(muni_[c(2,9)])
plot(muni_["eff"],  pal = heat.colors(10, rev = T), main = "Hunting effort")
plot(muni_["hy"], pal = terrain.colors(10, rev = T), main = "Hunting yield (animals * effort)")

library(Hmisc)
library(fastDummies)
library(dplyr)
library(tidyverse)


ff <- cut2(muni_$eff, g = 3)
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
  #b_eff ~ dunif(0, 0.5)
  
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
             #eff = runif(1,0,0.5)
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

par(mfrow = c(1,2))
plot(pred1[], sarea[], pch = 16, cex = .8)
abline(a=1, b=1, col = "darkred", lwd = 2)
plot(pred3[], sarea[], pch = 16, cex = .8)
abline(a=1, b=1, col = "darkred", lwd = 2)

effEst <- data.frame(cbind(c(samples_mcmc$chain2[,1],samples_mcmc$chain2[,2],samples_mcmc$chain2[,3]),
c(rep(1,length(samples_mcmc$chain2[,1])), rep(2,length(samples_mcmc$chain2[,1])), 
  rep(3,length(samples_mcmc$chain2[,1])))))
boxplot(X1 ~ X2, data = effEst, xlab = "Hunting ground Category", ylab = "Estimated Effort", cex.lab = 1.4,
        col = c("darkblue", "purple", "yellow"))
