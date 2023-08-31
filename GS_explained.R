# Loading packages
library(plgp)
library(mvtnorm)
library(raster)
library(spBayes)
library(spOccupancy)

# References
# Gaussian process regression: https://bookdown.org/rbg/surrogates/chap5.html
# Jeff Dosser spOccupancy manual: https://www.jeffdoser.com/files/spoccupancy-web/articles/modelfitting#spPGOcc
#

# Setting seed
set.seed(1)

# Creating a 30x30 grid (study area)
nx <- 30
x <- seq(0, 3, length=nx)
X <- expand.grid(x, x)

# Setting sig2 and phi parameters. With this parameters we'll build the spatial 
# covariance matrix Sigma, following an exponential model (others are allowed)
sig2 <- 1 # This is the variance, value in the matrix diagonal
phi <- 5 # Spatial decay parameter, how fast covariance decay 
#help(mkSpCov)
Sigma <- mkSpCov(coords = as.matrix(X), K = as.matrix(sig2), 
                 Psi = as.matrix(0), theta = phi, cov.model = "exponential")

# To better understand the effects of sig2 and phi, you can change the values and
# plot the variance-covariance matrix Sigma
# Sig <- raster(nrows = 900, ncols = 900, xmn = 0, xmx = 900, ymn = 0, ymx = 900)
# Sig[] <- Sigma
# plot(Sig)

# Now we generate random numbers from a multivariate normal distribution with
# mean = 0 and covariance matrix = Sigma
Y <- rmvnorm(1, rep(0, 900), Sigma)

# This would be the "latent" or purely spatial component of the species distribution,
# not affected by any other covariate
r <- raster(nrows = 30, ncols = 30, xmn = 0, xmx = 3, ymn = 0, ymx = 3)
r[cellFromXY(r,X)] <- Y
plot(r)
# To better understand the effect of sig2 and phi you can change those values and
# plot also several r

# Now we will create an environmental covariate that affects the distribution
temper <- raster(nrows = 3, ncols = 3, xmn = 0, xmx = 3, ymn = 0, ymx = 3)
temper[] <- rnorm(9)
temper <- disaggregate(temper, 10, "bilinear")
names(temper) <- "temper"

# Finally, following the Jeff Dosser formulation we'll create the distribution
# probability with a logit link: 
# https://www.jeffdoser.com/files/spoccupancy-web/articles/modelfitting#spPGOcc
psi <- exp(r + (-1.2+2*temper))/(1+exp(r + (-1.2+2*temper)))

# Now we can draw the actual distribution thought a random binomial number using 
# this probability
distr <- psi
distr[] <- rbinom(900, 1, psi[])

# Plot
par(mfrow = c(2,2))
plot(r, main = "Pure spatial")
plot(temper, main = "temper")
plot(psi, main = "Psi")
plot(distr, main = "distribution")

################################################################################

# Now we'll simulate a sampling process: we'll visit 200 places with a fixed
# detection probability of 0.7
n <- 200 # Number of visited places
sampID <- sample(1:900, n) # Random selection of n places
coords <- xyFromCell(distr, sampID) # Getting coordinates of those places
occ.covs <- matrix(temper[sampID]) # Getting environmental covariates of those places
colnames(occ.covs) <- "temper"
y <- matrix(NA, n, 4) # Set of detections/non-detections 

for (j in 1:4){
  for (i in 1:length(sampID)){
    y[i,j] <- rbinom(1, extract(distr,sampID[i]), 0.7)
  }
}

# Finally we'll fit two models. One spatial + temperature and other non-spatial, 
# only with the temperature. We'll use spOccupancy.
# Setting data
data.list <- list(y = y, 
                  occ.covs = occ.covs, 
                  coords = coords)
# For the rest of the initial settings/priors, we'll follow the recomendations of 
# Jeff Doser at https://www.jeffdoser.com/files/spoccupancy-web/articles/modelfitting#spPGOcc
# He specifically stress that phi parameter is one of the most difficult to fit,
# as well as other related with spatial process (more for Matern cor model)

# Number of batches
n.batch <- 500
# Batch length
batch.length <- 25
n.iter <- n.batch * batch.length
# Pair-wise distances between all sites
dist.data.list <- dist(data.list$coords)
# Priors 
min.dist <- min(dist.data.list)
max.dist <- max(dist.data.list)
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                   alpha.normal = list(mean = 0, var = 2.72),
                   sigma.sq.ig = c(2, 2), 
                   phi.unif = c(3/max.dist, 3/min.dist)) #better c(3/max(dist.btw.loc),3/max(dist.btw.loc))

# Initial values
inits.list <- list(alpha = 0, beta = 0,
                   phi = 3 / mean(dist.data.list), # better 3/mean(dist.btw.loc) (Banerjee, Carlin, and Gelfand 2003)
                   sigma.sq = 2,
                   w = rep(0, n), # CHANGED
                   z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 1) 

mSpatial <- spPGOcc(occ.formula = ~ temper, 
               det.formula = ~ 1, 
               data = data.list, 
               inits = inits.list, 
               n.batch = n.batch, 
               batch.length = batch.length, 
               priors = prior.list,
               cov.model = "exponential", 
               tuning = tuning.list, 
               NNGP = TRUE, # Fitting a NNGP rather than GP
               n.neighbors = 5, # Fitting a NNGP rather than GP
               search.type = 'cb',# Fitting a NNGP rather than GP
               n.report = 100, 
               n.burn = 2000, 
               n.chains = 3)
mEnviron <- PGOcc(occ.formula = ~ temper, 
              det.formula = ~ 1, 
              data = data.list, 
              inits = inits.list, 
              n.samples = 5000,  
              priors = prior.list,
              n.report = 100, 
              n.burn = 2000, 
              n.chains = 3)
summary(mSpatial)
summary(mEnviron)
waicOcc(mSpatial)
waicOcc(mEnviron)

# Do prediction. 
library(stars)
pred.0 <- cbind(1, temper[]) 
coords.0 <- as.matrix(xyFromCell(distr, 1:900))
# Approx. run time: 6 min
out.sp.pred <- predict(mSpatial, pred.0, coords.0, n.omp.threads = 4, 
                            verbose = TRUE, n.report=100)
out.sp.pred2 <- predict(mEnviron, pred.0)
# Produce a species distribution map (posterior predictive means of occupancy)
plot.dat <- data.frame(x = coords.0[,1], 
                       y = coords.0[,2], 
                       mean.psi = apply(out.sp.pred$psi.0.samples, 2, mean), 
                       sd.psi = apply(out.sp.pred$psi.0.samples, 2, sd))
plot.dat2 <- data.frame(x = coords.0[,1], 
                        y = coords.0[,2], 
                        mean.psi = apply(out.sp.pred2$psi.0.samples, 2, mean), 
                        sd.psi = apply(out.sp.pred2$psi.0.samples, 2, sd))
mPred <- r
mPred[] <- plot.dat$mean.psi
mPred2 <- r
mPred2[] <- plot.dat2$mean.psi


par(mfrow = c(2,3))
plot(r, main = "Pure spatial effect")
plot(temper, main = "Temperature")
plot(psi, main = "True Psi")
plot(distr, main = "True distribution\n and sampling points")
points(xyFromCell(distr,sampID), pch = 16, cex = .6)
plot(mPred, main = "Prediction with\n spatial component")
plot(mPred2, main = "Prediction without\n spatial component")
