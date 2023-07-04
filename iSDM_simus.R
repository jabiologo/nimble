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
lambda <- exp(beta0 + beta1*(dwat) + beta2*(tree))# + beta3*(autoc))

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
muni_$pressCov <- runif(nrow(muni_),-2,2)
muni_$press <- plogis(-0.5 + 0.3*muni_$pressCov)    #rnorm(nrow(muni_),0.35,0.1)
muni_$hy <- round(muni_$animals * muni_$press, 0)

plot(muni_["animals"])
plot(muni_["press"])
plot(muni_["hy"])

muni_$effCov1 <- muni_$pressCov + rnorm(n= length(muni_$press), mean = 0, sd = 2)
muni_$effCov2 <- muni_$pressCov + rnorm(n= length(muni_$press), mean = 0, sd = 1)
cor(muni_$pressCov, muni_$effCov1) # ~ 0.5
cor(muni_$pressCov, muni_$effCov2) # ~ 0.7

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
dev.off()
plot(sarea)
lines(muni_, lwd = 0.5)
lines(muni_sampled, lwd = 3, col = "darkblue")
points(xyFromCell(sarea,latSamp$ID), pch = 18)

################################################################################
library(nimble)
library(coda)

# Declaring constants for all models
constantsOc <- list(ncell = ncell(sarea),
                  nhg = nrow(muni_),
                  nds = nrow(muni_sampled),
                  nlt = nrow(latSamp),
                  hgLow = muni_$minId,
                  hgHigh = muni_$maxId,
                  dsLow = muni_sampled$minId,
                  dsHigh = muni_sampled$maxId,
                  latID = latSamp$ID)

# Declaring data
dataOc <- list(harvest = muni_$hy,
             ds = muni_sampled$abuSamp,
             letri = latSamp$latC,
             dwat = dwat[],
             tree = tree[],
             eff = muni_$effCov2)#,
#constraint1 = rep(1,nrow(dcal)))

mOc <- nimbleCode( {
  # PRIORS
  b0 ~ dnorm(0, 10)
  b1 ~ dnorm(0, 10) 
  b2 ~ dnorm(0, 10)
  a0 ~ dnorm(0, 10)
  a1 ~ dnorm(0, 10)
  c0 ~ dnorm(0, 10)
  c1 ~ dnorm(0, 10)
  
  # LIKELIHOOD
  for(i in 1:ncell){
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- b0 + b1*dwat[i] + b2*tree[i]
  }
  # Obs process for hunting yield
  for(h in 1:nhg){
    harvest[h] ~ dpois(lambda_hg[h]*(S[h]))
    log(lambda_hg[h]) <- log(sum(lambda[hgLow[h]:hgHigh[h]]))
    logit(S[h]) <- a0 + a1*eff[h]
  }
  # # Obs process for distance sampling
  for(j in 1:nds){
     log(ds[j]) <- log(sum(lambda[dsLow[j]:dsHigh[j]]))
  }
  # Obs process for latrines
    for(l in 1:nlt){
      letri[l] <- (exp((c0+log(lambda[latID[l]]))/c1)/2)
    }
} )

keepers <- c("b0", "b1", "b2", "a0", "a1", "c0", "c1")

nc <- 1
nb <- 10000
ni <- nb + 50000

inits <- function() {
  base::list(N = rep(1, constantsOc$ncell),
             lambda = rep(1, constantsOc$ncell),
             b0 = runif(1, -1, 1),
             b1 = runif(1, -1, 1),
             b2 = runif(1, -1, 1),
             a0 = runif(1, -1, 1),
             a1 = runif(1, -1, 1)
  )
}

modelOc <- nimbleModel(code = mOc, data = dataOc, constants = constantsOc, calculate = FALSE)
c_modelOc <- compileNimble(modelOc)
model_confOc <- configureMCMC(modelOc, useConjugacy = FALSE)
model_confOc$addMonitors(keepers)
model_mcmcOc <- buildMCMC(model_confOc)
c_model_mcmcOc <- compileNimble(model_mcmcOc, project = modelOc)
samplesOc <- runMCMC(c_model_mcmcOc, nburnin = nb, niter = ni, nchains = nc)

# Trace plots
sampOc <- coda::mcmc(samplesOc)
#sampOc <- coda::as.mcmc.list(lapply(samplesOc, coda::mcmc))
par(mfrow=c(3,3))
traceplot(sampOc)
effectiveSize(sampOc)
gelman.diag(sampOc)

# Merge chains
sampSs <- rbind(sampSs$chain1, sampSs$chain2, sampSs$chain3)
#plotAuto(sampSs)
# Store coefficients for prediction (mean)
#a0Ss <- mean(sampSs[,1])
#a1Ss <- mean(sampSs[,2])
#a2Ss <- mean(sampSs[,3])
#a3Ss <- mean(sampSs[,4])
#a4Ss <- mean(sampSs[,5])
#b0Ss <- mean(sampSs[,6])
#b1Ss <- mean(sampSs[,7])
#b2Ss <- mean(sampSs[,8])
#print(c(a0Ss, a1Ss, a2Ss, a3Ss, a4Ss, b0Ss, b1Ss, b2Ss))
# Store coefficients for prediction (samples from posterior)
a0Ss <- sample(sampSs[,1], 10000)
a1Ss <- sample(sampSs[,2], 10000)
a2Ss <- sample(sampSs[,3], 10000)
a3Ss <- sample(sampSs[,4], 10000)
a4Ss <- sample(sampSs[,5], 10000)
b0Ss <- sample(sampSs[,6], 10000)
b1Ss <- sample(sampSs[,7], 10000)
b2Ss <- sample(sampSs[,8], 10000)

# Abundance at grid cell level (for coefficients = mean)
# grid_$SsN <- exp(b0Ss + b1Ss*grid_$cover + b2Ss*grid_$climate)
# Abundance at grid cell level (for coefficients = samples from posterior)
grid_SsN <- exp(b0Ss + b1Ss%o%grid_$cover + b2Ss%o%grid_$climate)
grid_$SsN_mean <- colMeans(grid_SsN)
ci_SsN <- apply(grid_SsN, MARGIN = 2, quantile, probs = c(0.05, 0.95))
grid_$SsN_low <- ci_SsN[1,]
grid_$SsN_high <- ci_SsN[2,]

# Abundance aggregation at hunting ground level
SsHgN_mean <- aggregate(x = grid_$SsN_mean, by = list(grid_$Hg), FUN = sum)
names(SsHgN_mean)[2] <- "SsHgN_mean"
SsHgN_low <- aggregate(x = grid_$SsN_low, by = list(grid_$Hg), FUN = sum)
names(SsHgN_low)[2] <- "SsHgN_low"
SsHgN_high <- aggregate(x = grid_$SsN_high, by = list(grid_$Hg), FUN = sum)
names(SsHgN_high)[2] <- "SsHgN_high"
hg_ <- left_join(hg_, SsHgN_mean, by = c("Hg" = "Group.1") )
hg_ <- left_join(hg_, SsHgN_low, by = c("Hg" = "Group.1") )
hg_ <- left_join(hg_, SsHgN_high, by = c("Hg" = "Group.1") )
# Hunting pressure / catch rate for each hunting ground (for coefficients = mean)
# hg_$SsP <- exp(a0Ss + a1Ss*hg_$Ax1Ss + a2Ss*hg_$Ax2Ss + a3Ss*hg_$CERRAMIENT + a4Ss*hg_$areaS) / (1 + exp(a0Ss + a1Ss*hg_$Ax1Ss + a2Ss*hg_$Ax2Ss + a3Ss*hg_$CERRAMIENT + a4Ss*hg_$areaS))
# Hunting pressure / catch rate for each hunting ground (for coefficients = samples from posterior)
hg_SsP <- exp(a0Ss + a1Ss%o%hg_$Ax1Ss + a2Ss%o%hg_$Ax2Ss + a3Ss%o%hg_$CERRAMIENT + a4Ss%o%hg_$areaS) / (1 + exp(a0Ss + a1Ss%o%hg_$Ax1Ss + a2Ss%o%hg_$Ax2Ss + a3Ss%o%hg_$CERRAMIENT + a4Ss%o%hg_$areaS))
hg_$SsP_mean <- colMeans(hg_SsP)
ci_SsP <- apply(hg_SsP, MARGIN = 2, quantile, probs = c(0.05, 0.95))
hg_$SsP_low <- ci_SsP[1,]
hg_$SsP_high <- ci_SsP[2,]

# Hunting yield / harvest prediction at hunting ground level
hg_$SsHarv_mean <- hg_$SsHgN_mean *hg_$SsP_mean
hg_$SsHarv_low <- hg_$SsHgN_low *hg_$SsP_low
hg_$SsHarv_high <- hg_$SsHgN_high *hg_$SsP_high
par(mfrow=c(1,1))
hist(hg_$SsP_mean, breaks = 50, col = "darkgrey", xlab = "p", ylab = "nº hunting grounds", main = "Wild boar")
hist(hg_$SsP_low, breaks = 50, add = TRUE, col = alpha("blue", 0.25))
hist(hg_$SsP_high, breaks = 50, add = TRUE, col = alpha("red", 0.25))
















