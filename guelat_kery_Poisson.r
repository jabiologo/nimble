### Guelat & Kery: Effects of spatial autocorrelation and imperfect detection on species distribution models
### Adapted script to use Poisson and Gaussian model for spatial autocorrelation in NIMBLE
library(terra)
set.seed(2)

# Random multivariate normal
rmvn <- function(n, mu = 0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V), p))))   stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}

# Parameters
nvisits <- 1 # single visit
sampleSizes <- 500
detections <- 1 # perfect detectability
spCor <- 1

# Number of simulated datasets
niter <- 1

# Set up a square lattice region
xdim <- 50
ydim <- 50
x.easting <- 1:xdim
x.northing <- 1:ydim
simgrid <- expand.grid(x.easting, x.northing)
n <- nrow(simgrid)

# Set up distance matrix
distance <- as.matrix(dist(simgrid))
# Covariate
#covx <- replicate(niter, rnorm(n), simplify = FALSE)
covx <- replicate(niter, as.vector(rmvn(1, rep(0, n), 1.5 * exp(-2.1 * distance))), simplify = FALSE)
# Sampled sites
sites <- replicate(niter, sort(sample(1:n, size = sampleSizes[1])), simplify = FALSE)
# Patchy random fields
rho <- replicate(niter, as.vector(rmvn(1, rep(0, n), 1.5 * exp(-0.05 * distance))), simplify = FALSE)
# Coefficients 
alpha <- -0.6
beta <- 0.8

# Expected lambda
lambda <- mapply(function(x, y) exp(alpha + beta*x + spCor*y), covx, rho, SIMPLIFY = FALSE)
# Realized abun
abun <- lapply(lambda, function(x) rpois(n = n, lambda = x))
# Observed counts
Y <- lapply(abun, function(x) replicate(nvisits, rbinom(n = n, size = x, p = detections)))

covxR <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
covxR[cellFromXY(covxR,simgrid-0.5)] <- covx[[1]]
rhoR <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
rhoR[cellFromXY(rhoR,simgrid-0.5)] <- rho[[1]] 
lambdaR <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
lambdaR[cellFromXY(lambdaR,simgrid-0.5)] <- lambda[[1]] 
abunR <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
abunR[cellFromXY(abunR,simgrid-0.5)] <- abun[[1]]  

par(mfrow = c(2,2))
plot(rhoR, main = "Pure spatial effect, rho")
plot(covxR, main = "Temperature")
plot(lambdaR, main = "True Lambda")
plot(abunR, main = "True abundance")

################################################################################
# NIMBLE approach
library(nimble)
library(MCMCvis)

coords <- simgrid[sites[[1]][],]
occ.covs <- matrix(covx[[1]][sites[[1]][]])
colnames(occ.covs) <- "temper"
y <- Y[[1]][sites[[1]][],]

dists <- as.matrix(dist(coords))
dists <- dists / max(dists)  # normalize to max distance of 1

min.dists <- 3/unique(dists[order(dists)])[2] # min(dist) after 0
max.dists <- 3/max(dists)

constants <- list(N = length(y),  # nrow(y),
                  dists = dists, 
                  ones = rep(1, length(y)))  # nrow(y)

data <- list(y = y,
             z = as.numeric(occ.covs))

expcov <- nimbleFunction(     
  run = function(dists = double(2), phi = double(0), sigma2 = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    #sigma2 <- sigma*sigma
    for(i in 1:n)
      for(j in 1:n)
        result[i, j] <- sigma2*exp(-dists[i,j]*phi)
    return(result)
  })
cExpcov <- compileNimble(expcov)

code <- nimbleCode({
  # Priors
  alpha ~ dnorm(mean = 0, sd = 2.72)
  beta ~ dnorm(mean = 0, sd = 2.72)
  #mu0 ~ dnorm(mean = 0, sd = 0.5) # We'll use a zero-mean spatial Gaussian Process (as spOccupancy)
  sigma2 ~ dinvgamma(shape = 2, scale = 1.5) # shape = 2, scale = is the mean
  phi ~ dunif(3, 203) #max.dists, min.dists
  # p ~ dbeta(1,1)
  
  # Latent spatial process
  mu[1:N] <- 0*ones[1:N]
  cov[1:N, 1:N] <- expcov(dists[1:N, 1:N], phi, sigma2)
  s[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dpois(lam[i]) 
    log(lam[i]) <- alpha + beta * z[i] + s[i] 
  }
})

inits <- list(alpha = 0, beta = 0, sigma2 = 1, phi = 3/mean(dists))#, 
              #p = 0.5, Y = apply(y, 1, max, na.rm = TRUE))

## setup initial spatially-correlated latent process values
inits$cov <- expcov(dists, inits$phi, inits$sigma)
inits$s <-  t(chol(inits$cov)) %*% rnorm(length(y))
inits$s <- inits$s[ , 1]  # so can give nimble a vector rather than one-column matrix

model <- nimbleModel(code, constants = constants, data = data, inits = inits)
cModel <- compileNimble(model)
conf <- configureMCMC(model)
conf$addMonitors('s')
#conf$removeSamplers('s[1:500]')
# Changing a tunable parameter in the adaptation of RW_block makes a big difference.
# adaptFactorExponent. Exponent controling the rate of decay of the scale 
# adaptation factor. See Shaby and Wells, 2011, for details. (default = 0.8)
# IT'S TOO SLOW. ALSO PROBLEMS STORING s()
# PROBLEMS STORING s[] CHECK KEEPERS
# IMPORTANT!!!! WE SHOULD BACK-TRNSFORM (phi and sigma are for standarized distance)
#conf$addSampler('s[1:500]', 'RW_block', control = list(adaptFactorExponent = 0.25))
## Changing a tunable parameter in the adaptation of RW_block makes a big difference.
# [Note] Assigning an RW_block sampler to nodes with very different scales can 
# result in low MCMC efficiency.  If all nodes assigned to RW_block are not on 
# a similar scale, we recommend providing an informed value for the "propCov" 
# control list argument, or using the AFSS sampler instead.


MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC, project = cModel)

samples <- runMCMC(cMCMC, niter = 10000, nburnin = 5000, thin = 10, nchains = 2)
#save(samples, file = "samples_guelat_kery_nimble_poisson_standarized.RData")
#load("samples_guelat_kery_nimble_poisson.RData")

MCMCsummary(samples, round = 2,
            params = c("alpha", "beta", "sigma2", "phi"))
MCMCtrace(samples,
          pdf = FALSE,
          ind = TRUE,
          Rhat = TRUE,
          n.eff = TRUE,
          params = c("alpha", "beta", "sigma2", "phi")) 

###############################
# Predicting in new locations #
###############################

sampJoint <- rbind(samples$chain1, samples$chain2)

newlocs <- simgrid
#dist11 <- fields::rdist(coords)
#dist21 <- fields::rdist(newlocs, coords)
#dist22 <- fields::rdist(newlocs)

# For standardized distances
dist11 <- fields::rdist(coords) / max(dists)
dist21 <- fields::rdist(newlocs, coords) / max(dists)
dist22 <- fields::rdist(newlocs) / max(dists)

sample_xstar <- nimbleFunction(
  # need types for inputs
  run = function(s = double(1), mu = double(1), muNew = double(1), sigma = double(0), phi = double(0), 
                 dist11 = double(2), dist22 = double(2), dist21 = double(2)) {
    returnType(double(1))
    n <- length(muNew)
    sigma2 <- sigma#*sigma # Check parameterization
    C22 <- sigma2 * exp(-dist22 * phi) # CHANGED
    C11 <- sigma2 * exp(-dist11 * phi) # CHANGED
    C21 <- sigma2 * exp(-dist21 * phi) # CHANGED
    # Note that this could be made a bit more efficient by using the Cholesky
    # decomposition rather than solve().
    xstar <- muNew + (C21 %*% solve(C11, s - mu))[,1]
    # THIS IS NOT WORKING IF PIVOT = FALSE 
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
      
      alpha <- samples[i, 1] # CHANGED
      beta <- samples[i, 2] # CHANGED
      #mu0 <- samples[i, 3] # CHANGED
      phi <- samples[i, 3] # CHANGED
      sigma <- samples[i, 504] # CHANGED
      s <- samples[i, 4:503] # CHANGED
      mu <- alpha + beta*z # CHANGED
      muNew <- alpha + beta*zNew # CHANGED
      # Get other parameters
      output[i, ] <- sample_xstar(s, mu, muNew, sigma, phi, dist11, dist22, dist21)
      print(i)
    }
    return(output)
  }
)

# Prediction of s[] at new locations
# This function could be compiled in NIMBLE (see NIMBLE manual) 
xstar_samples <- get_samples(samples = sampJoint[seq(1,1000,10),], z = data$z, zNew = covx[[1]], 
                             dist11 = dist11, dist22 = dist22, dist21 = dist21)

# Obtaining estimates
alphaEst <- mean(sampJoint[,1])
betaEst <- mean(sampJoint[,2])
spatR <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
spatR[cellFromXY(spatR,simgrid-0.5)] <- colMeans(xstar_samples)
# Abundance predictions
predRNim <- exp(alphaEst + betaEst*covxR + spatR)

# Some plots
par(mfrow = c(2,4))
plot(rhoR, main = "True spatial effect")
plot(covxR, main = "Temperature")
plot(lambdaR, main = "True Lambda")
plot(abunR, main = "True abundance")
plot(spatR, main = "Predicted spatial effect")
plot(spatR[], rhoR[], main = "Spatial effects")
abline(a=0, b=1)
plot(predRNim, main = "Prediction with NIMBLE")#,  range = c(0,30))
plot(predRNim[], abunR[], main = "Abundance")
abline(a=0, b=1)

plot(predRNim[], abunR[])
abline(a=0, b=1)
cor(na.omit(as.numeric(predRNim[])), na.omit(as.numeric(abunR[])))

