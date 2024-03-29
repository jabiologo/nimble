################
## 26/10/2023 ##
################
# Spatial random effects with spOccupancu and NIMBLE

# Script to use with spOccupancy and NIMBLE with spatial Gaussian Process for 
# spatial autocorrelation. We use adapted simulations from Guelat & Kery (2016): 
# Effects of spatial autocorrelation and imperfect detection on species distribution models
# Check also https://www.jeffdoser.com/files/spoccupancy-web/articles/modelfitting#spPGOcc
library(spOccupancy)
library(terra)
set.seed(2)

# Function to simulate random multivariate normal samples
rmvn <- function(n, mu = 0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V), p))))   stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}

# Parameters
nvisits <- 3 # Number of repeated visits
sampleSizes <- 500 # Number of sampled sites
detections <- 0.7 # Detectability
spCor <- 1 # Value of gamma from Guélat & Kery. Fixed to 1 since it's included in variance
niter <- 1 # Number of simulated datasets

# Set up a square lattice region
xdim <- 50
ydim <- 50
x.easting <- 1:xdim
x.northing <- 1:ydim
simgrid <- expand.grid(x.easting, x.northing)
n <- nrow(simgrid) # Number of pixels/cells in the study area

# Set up distance matrix between each cell in the study area
distance <- as.matrix(dist(simgrid))
# Covariate create a random covariate (eg. environmental fixed effects)
covx <- replicate(niter, rnorm(n), simplify = FALSE)
# Sampled sites
sites <- replicate(niter, sort(sample(1:n, size = sampleSizes[1])), simplify = FALSE)
# Spatial random effects covariate using multivariate normal (spatial Gaussian Process)
# Using the notation in https://www.jeffdoser.com/files/spoccupancy-web/articles/modelfitting#spPGOcc
# omega(s) ~ N(0, Sigma(s,s',theta)), where
# omega(s) = spatial random effect
# N= multivariate normal
# 0= mean
# Sigma= covariate matrix
# s,s'= distance between each pair of coordinates (distances)
# theta= parameters for the correlation function
# We use the exponential for the correlation function in Sigma. See help(mkSpCov) from spBayes package for
# an explanation of other correlation functions: Gaussian, Matern and Spherical.
# Exponential function has one parameter called decay parameter (phi). Here is set at 0.05
# Exponential function is multiply by variance (sig.sq). Here is set at 1.5
# Thus, Sigma(s,s',theta)= sig.sq * exp(-phi*distances) 
sig.sq <- 1.5
phi <- 0.05
rho <- replicate(niter, as.vector(rmvn(1, rep(0, n), sig.sq * exp(-(phi) * distance))), simplify = FALSE)
# Coefficients for intercept (alpha) and fixed effects/environmental covariate (beta)
alpha <- -0.8
beta <- 0.6

# Expected psi (occupancy). We use logit link
psi <- mapply(function(x, y) exp(alpha + beta*x + spCor*y) / (1 + exp(alpha + beta*x + spCor*y)), covx, rho, SIMPLIFY = FALSE)
# Realized presence/absence (random samples from a Benoulli distribution with probability = psi)
presence <- lapply(psi, function(x) rbinom(n = n, size = 1, p = x))
# Observed detections (repeated visits with detectability = detections)
Y <- lapply(presence, function(x) replicate(nvisits, rbinom(n = n, size = x, p = detections)))

# Produce rasters to easy plot
covxR <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
covxR[cellFromXY(covxR,simgrid-0.5)] <- covx[[1]]
rhoR <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
rhoR[cellFromXY(rhoR,simgrid-0.5)] <- rho[[1]] 
psiR <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
psiR[cellFromXY(psiR,simgrid-0.5)] <- psi[[1]] 
presenceR <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
presenceR[cellFromXY(presenceR,simgrid-0.5)] <- presence[[1]]  

# Plot (old school sorry, no ggplot)
par(mfrow = c(2,2))
plot(rhoR, main = "Pure spatial effect, rho")
plot(covxR, main = "Temperature")
plot(psiR, main = "True Psi")
plot(presenceR, main = "True presence")

###############
# spOccupancy #
###############
# Model fitting with spOccupancy 

# Prepara the data for spOccupancy
# Get coordinates of our sample sites
coords <- simgrid[sites[[1]][],]
# Get environmental covariate values of our sample sites
occ.covs <- matrix(covx[[1]][sites[[1]][]])
# We'll cal it temperature just to remember it
colnames(occ.covs) <- "temper"
# Get the detection/non-detection list for our sample sites
y <- Y[[1]][sites[[1]][],]

# Create a list
data.list <- list(y = y, occ.covs = occ.covs, coords = coords)

# A bunch of settings needed for spOccupancy
n.batch <- 500
batch.length <- 25
n.iter <- n.batch * batch.length
dist.data.list <- dist(data.list$coords)
tuning.list <- list(phi = 1) 

# Priors and initial values for parameters. 
# It's interesting to read https://www.jeffdoser.com/files/spoccupancy-web/articles/modelfitting#spPGOcc
# for a better explanation of priors settings, especially for phi and sigma.sq related to the spatial
# random effects... those parameters are usually the hardest to fit, so "good" priors and "good" initials 
# are important
min.dist <- min(dist.data.list)
max.dist <- max(dist.data.list)
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                   alpha.normal = list(mean = 0, var = 2.72),
                   sigma.sq.ig = c(2, 2), 
                   phi.unif = c(3/max.dist, 3/min.dist)) 
inits.list <- list(alpha = 0, beta = 0,
                   phi = 3 / mean(dist.data.list), 
                   sigma.sq = 2,
                   w = rep(0, 500), 
                   z = apply(y, 1, max, na.rm = TRUE))

mSpatial <- spPGOcc(occ.formula = ~ temper, # occupancy formula (here, only fixed effects covariate)
                    det.formula = ~ 1, # observational process formula (here detectability is fixed)
                    data = data.list, 
                    inits = inits.list, 
                    n.batch = n.batch, 
                    batch.length = batch.length, 
                    priors = prior.list,
                    cov.model = "exponential", # correlation function used for covariate model in spatial Gaussian Process
                    tuning = tuning.list, 
                    NNGP = TRUE, # Fitting a NNGP rather than GP (similar but faster)
                    n.neighbors = 5, # Fitting a NNGP rather than GP (similar but faster)
                    search.type = 'cb',# Fitting a NNGP rather than GP (similar but faster)
                    n.report = 100, 
                    n.burn = 2000, 
                    n.chains = 3)
# Summary of the model
summary(mSpatial)

# Prepare data for prediction
pred.0 <- cbind(1, covx[[1]][])
coords.0 <- as.matrix(simgrid-0.5)
# Predictions to the whole study area. Approx. run time: X min
out.sp.pred <- predict(mSpatial, pred.0, coords.0, n.omp.threads = 4, 
                       verbose = TRUE, n.report=100)
# Produce a species distribution map (posterior predictive means of occupancy)
plot.dat <- data.frame(x = coords.0[,1], 
                       y = coords.0[,2], 
                       mean.psi = apply(out.sp.pred$psi.0.samples, 2, mean), 
                       sd.psi = apply(out.sp.pred$psi.0.samples, 2, sd))

# Produce raster
predR <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
predR[cellFromXY(predR,simgrid-0.5)] <- plot.dat$mean.psi

# Some plots of Predicted vs Observed Occupancy
par(mfrow = c(2,3))
plot(rhoR, main = "Spatial random effects")
plot(covxR, main = "Temperature")
plot(psiR, main = "True Psi")
plot(presenceR, main = "True presence")
plot(predR, main = "Predicted presence")
plot(predR[], psiR[])
abline(a=0, b=1, col = "darkred")

##########
# NIMBLE #
##########
# Trying to reproduce spatial Gaussian Process model in NIMBLE
library(nimble)
library(MCMCvis)

# First, we'll create a compiled function in NIMBLE to compute the exponential covariance model
# Function is taken from here: https://r-nimble.org/nimbleExamples/gaussian_process.html
expcov <- nimbleFunction(
  run = function(dists = double(2), phi = double(0), sigma2 = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n)
      for(j in 1:n)
        result[i, j] <- sigma2*exp(-dists[i,j]*phi)
    return(result)
  })
cExpcov <- compileNimble(expcov)

# Prepare the data for NIMBLE
dists <- as.matrix(dist(coords))
# Normalize to max distance of 1. We will not do this now, but it would be necessary 
# in real life examples. I have some problems for backtransformation after the model 
# fitting, so by now we'll use raw distances.
#dists <- dists / max(dists)  # normalize to max distance of 1

# Compute this values for the phi prior (check in the model)
min.dists <- 3/unique(dists[order(dists)])[2] # min(dist) after 0
max.dists <- 3/max(dists)

# Constants for NIMBLE
constants <- list(N = nrow(y),
                  nsurveys = 3,
                  dists = dists, 
                  ones = rep(1, nrow(y)))

# Data for NIMBLE
data <- list(y = y,
             z = as.numeric(occ.covs))

code <- nimbleCode({
  # Priors
  alpha ~ dnorm(mean = 0, sd = 2.72) # The same as spOccupancy
  beta ~ dnorm(mean = 0, sd = 2.72) # The same as spOccupancy
  sigma2 ~ dinvgamma(shape = 2, scale = 2) # The same as spOccupancy
  phi ~ dunif(3, 203) # Using here the values computed before dunif(max.dists, min.dists)
  p ~ dbeta(1,1) # Non informative
  
  # Latent spatial process
  mu[1:N] <- 0*ones[1:N] # mean is fixed to 0 in our model
  cov[1:N, 1:N] <- expcov(dists[1:N, 1:N], phi, sigma2) # Exponentail function for covariance (Sigma)
  s[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N]) # Multivariate normal distribution with mean and Sigma
  # Likelihood for a site-occupancy model with spatial random effects in Occupancy process
  for (i in 1:N) {
    Y[i] ~ dbern(psi[i]) # True occupancy status
    # Occupancy (psi) is fixed effects plus spatial random effects "s[i]"
    logit(psi[i]) <- alpha + beta * z[i] + s[i] 
    for (j in 1:nsurveys) {
      y[i,j] ~ dbern(Y[i] * p) # Observed data
    }
  }
  
})

# Initial values. Remember that they are important for convergence
inits <- list(alpha = 0, beta = 1, sigma2 = 1, phi = 3/mean(dists), 
              p = 0.5, Y = apply(y, 1, max, na.rm = TRUE))
## setup initial spatial random effects initial values
inits$cov <- expcov(dists, inits$phi, inits$sigma)
inits$s <-  t(chol(inits$cov)) %*% rnorm(nrow(y))
inits$s <- inits$s[ , 1]  # so can give to NIMBLE a vector rather than one-column matrix

model <- nimbleModel(code, constants = constants, data = data, inits = inits)
cModel <- compileNimble(model)
conf <- configureMCMC(model)
conf$addMonitors('s')

# You can "ignore" next chunk of code.
# I just tried to improve MCMC performance but I had problems... I'll check it in the future
# Changing a tunable parameter in the adaptation of RW_block makes a big difference.
# adaptFactorExponent. Exponent controling the rate of decay of the scale adaptation factor. See Shaby and Wells, 2011, for details. (default = 0.8)
# IT'S TOO SLOW. ALSO PROBLEMS STORING s()
#conf$removeSamplers('s[1:500]')
#conf$addSampler('s[1:500]', 'RW_block', control = list(adaptFactorExponent = 0.25))
# PROBLEMS STORING s[] CHECK KEEPERS
# IMPORTANT!!!! WE SHOULD BACK-TRNSFORM (phi and sigma are for standarized distance)
# [Note] Assigning an RW_block sampler to nodes with very different scales can 
# result in low MCMC efficiency.  If all nodes assigned to RW_block are not on 
# a similar scale, we recommend providing an informed value for the "propCov" 
# control list argument, or using the AFSS sampler instead.

MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC, project = cModel)
samples <- runMCMC(cMCMC, niter = 100000, nburnin = 10000, thin = 10, nchains = 2)

# Saving/loading MCMC runs
# save(samples, file = "samples_guelat_kery_nimble.RData")
# load("samples_guelat_kery_nimble.RData")

# MCMC results
MCMCsummary(samples, round = 2,
            params = c("alpha", "beta", "sigma2", "phi", "p"))
MCMCtrace(samples,
          pdf = FALSE,
          ind = TRUE,
          Rhat = TRUE,
          n.eff = TRUE,
          params = c("alpha", "beta", "sigma2", "phi", "p")) 

###############################
# Predicting in new locations #
###############################

# It's necessary to make some additional steps to predict in new locations, 
# mainly because you have to compute spatial random effects in new places

# Merging chains
sampJoint <- rbind(samples$chain1, samples$chain2)

# Preparing data for predictions at new locations
newlocs <- simgrid # New loc coordinates (in our case, the whole study area)
dist11 <- fields::rdist(coords) # Matrix of distances among our sampled sites
dist21 <- fields::rdist(newlocs, coords) # Distances among new locations and our sampled sites
dist22 <- fields::rdist(newlocs) # Distances among the new locations

# Now we'll define two functions that are needed to compute preditions in new locations
# This is explained a bit in the attached PDF from the NIMBLE course. This is out of my
# understanding because matrix operations... 

sample_xstar <- nimbleFunction(
  # need types for inputs
  run = function(s = double(1), mu = double(1), muNew = double(1), sigma = double(0), phi = double(0), 
                 dist11 = double(2), dist22 = double(2), dist21 = double(2)) {
    returnType(double(1))
    n <- length(muNew)
    sigma2 <- sigma*sigma
    C22 <- sigma2 * exp(-dist22 * phi) # CHANGED
    C11 <- sigma2 * exp(-dist11 * phi) # CHANGED
    C21 <- sigma2 * exp(-dist21 * phi) # CHANGED
    # Note that this could be made a bit more efficient by using the Cholesky
    # decomposition rather than solve().
    xstar <- muNew + (C21 %*% solve(C11, s - mu))[,1]
    # THIS IS NOT WORKING IF PIVOT = FALSE
    # Formula in the attached PDF
    xstar <- xstar + (t(chol(C22 - C21 %*% solve(C11, t(C21)), pivot = TRUE)) %*% rnorm(n))[,1] # CHANGED
    return(xstar)
  })

get_samples <- nimbleFunction(
  run = function(samples = double(2), z = double(1), zNew = double(1), 
                 dist11 = double(2), dist22 = double(2), dist21 = double(2)) {
    returnType(double(2))
    m <- dim(samples)[1]
    nstar <- dim(dist21)[1]
    output <- matrix(0, nrow = m, ncol = nstar)
    for(i in 1:m) {
      # Extract parameter values from the input matrix based on numeric column indexes.
      # Inelegant because hard-coded based on looking at column names of samples matrix.
      
      # IMPORTANT: Here you should match the columns from the MCMC samples
      alpha <- samples[i, 1] # CHANGED
      beta <- samples[i, 2] # CHANGED
      #mu0 <- samples[i, 3] # CHANGED
      phi <- samples[i, 4] # CHANGED
      sigma <- samples[i, 505] # CHANGED
      s <- samples[i, 5:504] # CHANGED
      mu <- alpha + beta*z # CHANGED
      muNew <- alpha + beta*zNew # CHANGED
      # Get other parameters
      output[i, ] <- sample_xstar(s, mu, muNew, sigma, phi, dist11, dist22, dist21)
      print(i)
    }
    return(output)
  }
)

# Obtain spatial random effects predictions
xstar_samples <- get_samples(samples = sampJoint[seq(1,8000,100),], z = data$z, zNew = covx[[1]], 
                             dist11 = dist11, dist22 = dist22, dist21 = dist21)
# You could compile this function in NIMBLE, but is not much faster
#cget_samples <- compileNimble(get_samples)
#xstar_samples2 <- cget_samples(samples = samples2[seq(1,90000,100),], z = data$z, zNew = temper[], 
#                               dist11 = dist11, dist22 = dist22, dist21 = dist21)

# Obtaining parameter means for predictions
alphaEst <- mean(sampJoint[,1])
betaEst <- mean(sampJoint[,2])
# Storing in a raster spatial random effects
spatR <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
spatR[cellFromXY(spatR,simgrid)] <- colMeans(xstar_samples)

# Model prediction (backtransform from logit)
predRNim <- exp(alpha + beta*covxR + spatR) / (1 + exp(alpha + beta*covxR + spatR))

# Some summary plots comparing predictions from spOccupancy and NIMBLE
par(mfrow = c(2,4))
plot(covxR, main = "Temperature")
plot(rhoR, main = "Pure spatial effect, rho")
plot(psiR, main = "True Psi")
plot(presenceR, main = "True presence")
plot(predR, main = "Prediction with spOccupancy")
plot(predR[], psiR[])
plot(predRNim, main = "Prediction with NIMBLE")
plot(predRNim[], psiR[])
    






