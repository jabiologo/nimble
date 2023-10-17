### Guelat & Kery: Effects of spatial autocorrelation and imperfect detection on species distribution models
### Adapted script to use with spOccupancy
library(spOccupancy)
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
nvisits <- 3
sampleSizes <- 500
detections <- 0.7
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
covx <- replicate(niter, rnorm(n), simplify = FALSE)
#covx <- replicate(niter, as.vector(rmvn(1, rep(0, n), 1.5 * exp(-0.05 * distance))), simplify = FALSE)
# Sampled sites
sites <- replicate(niter, sort(sample(1:n, size = sampleSizes[1])), simplify = FALSE)
# Patchy random fields
rho <- replicate(niter, as.vector(rmvn(1, rep(0, n), 1.5 * exp(-0.05 * distance))), simplify = FALSE)
# Coefficients 
alpha <- -0.8
beta <- 0.6

# Expected psi
psi <- mapply(function(x, y) exp(alpha + beta*x + spCor*y) / (1 + exp(alpha + beta*x + spCor*y)), covx, rho, SIMPLIFY = FALSE)
# Realized presence
presence <- lapply(psi, function(x) rbinom(n = n, size = 1, p = x))
# Observed detections
Y <- lapply(presence, function(x) replicate(nvisits, rbinom(n = n, size = x, p = detections)))

covxR <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
covxR[cellFromXY(covxR,simgrid)] <- covx[[1]]
rhoR <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
rhoR[cellFromXY(rhoR,simgrid)] <- rho[[1]] 
psiR <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
psiR[cellFromXY(psiR,simgrid)] <- psi[[1]] 
presenceR <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
presenceR[cellFromXY(presenceR,simgrid)] <- presence[[1]]  

par(mfrow = c(2,2))
plot(rhoR, main = "Pure spatial effect, rho")
plot(covxR, main = "Temperature")
plot(psiR, main = "True Psi")
plot(presenceR, main = "True presence")

###############
# spOccupancy #
###############

coords <- simgrid[sites[[1]][],]
occ.covs <- matrix(covx[[1]][sites[[1]][]])
colnames(occ.covs) <- "temper"
y <- Y[[1]][sites[[1]][],]

data.list <- list(y = y, occ.covs = occ.covs, coords = coords)

n.batch <- 500
batch.length <- 25
n.iter <- n.batch * batch.length
dist.data.list <- dist(data.list$coords)

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


summary(mSpatial)

pred.0 <- cbind(1, covx[[1]][])
coords.0 <- as.matrix(simgrid)
# Approx. run time: X min
out.sp.pred <- predict(mSpatial, pred.0, coords.0, n.omp.threads = 4, 
                       verbose = TRUE, n.report=100)
# Produce a species distribution map (posterior predictive means of occupancy)
plot.dat <- data.frame(x = coords.0[,1], 
                       y = coords.0[,2], 
                       mean.psi = apply(out.sp.pred$psi.0.samples, 2, mean), 
                       sd.psi = apply(out.sp.pred$psi.0.samples, 2, sd))



predR <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
predR[cellFromXY(predR,simgrid)] <- plot.dat$mean.psi

par(mfrow = c(2,3))
plot(rhoR, main = "Pure spatial effect, rho")
plot(covxR, main = "Temperature")
plot(psiR, main = "True Psi")
#points(xyFromCell(Y,sampID), pch = 16, cex = .6)
plot(presenceR, main = "True presence")
plot(predR, main = "Prediction with\n spatial component")
plot(predR[], psiR[])

################################################################################
# NIMBLE approach
library(nimble)
library(MCMCvis)

dists <- as.matrix(dist(coords))
dists <- dists / max(dists)  # normalize to max distance of 1

min.dists <- 3/unique(dists[order(dists)])[2] # min(dist) after 0
max.dists <- 3/max(dists)

constants <- list(N = nrow(y),
                  nsurveys = 3,
                  dists = dists, 
                  ones = rep(1, nrow(y)))

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
  sigma2 ~ dinvgamma(shape = 2, scale = 2) #################################
  phi ~ dunif(3, 203) # (max.dists, min.dists)
  p ~ dbeta(1,1)
  
  # Latent spatial process
  mu[1:N] <- 0*ones[1:N]
  cov[1:N, 1:N] <- expcov(dists[1:N, 1:N], phi, sigma2)
  s[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
  # Likelihood
  for (i in 1:N) {
    Y[i] ~ dbern(psi[i]) # True occupancy status
    logit(psi[i]) <- alpha + beta * z[i] + s[i] 
    for (j in 1:nsurveys) {
      y[i,j] ~ dbern(Y[i] * p) # Observed data
    }
  }
  
  # for (m in 1:M) {
  #   logit(psi[m]) <- beta[1] + beta[2] * elev[m] + beta[3] * forest[m]
  #   z[m] ~ dbern(psi[m])
  #   for (j in 1:J) {
  #     logit(p[m, j]) <- alpha[1] + alpha[2] * wind[m, j]
  #     y[m, j] ~ dbern(z[m] * p[m, j])
  #   }
  # }
  
  
  
})

inits <- list(alpha = 0, beta = 0, sigma2 = 1, phi = 3/mean(dists), 
              p = 0.5, Y = apply(y, 1, max, na.rm = TRUE))

## setup initial spatially-correlated latent process values
inits$cov <- expcov(dists, inits$phi, inits$sigma)
inits$s <-  t(chol(inits$cov)) %*% rnorm(nrow(y))
inits$s <- inits$s[ , 1]  # so can give nimble a vector rather than one-column matrix

model <- nimbleModel(code, constants = constants, data = data, inits = inits)
cModel <- compileNimble(model)
conf <- configureMCMC(model)
#conf$addMonitors('s')
conf$removeSamplers('s[1:500]')
# Changing a tunable parameter in the adaptation of RW_block makes a big difference.
# adaptFactorExponent. Exponent controling the rate of decay of the scale 
# adaptation factor. See Shaby and Wells, 2011, for details. (default = 0.8)
# IT'S TOO SLOW. ALSO PROBLEMS STORING s()

# PROBLEMS STORING s[] CHECK KEEPERS
# IMPORTANT!!!! WE SHOULD BACK-TRNSFORM (phi and sigma are for standarized distance)
conf$addSampler('s[1:500]', 'RW_block', control = list(adaptFactorExponent = 0.25))
## Changing a tunable parameter in the adaptation of RW_block makes a big difference.
# [Note] Assigning an RW_block sampler to nodes with very different scales can 
# result in low MCMC efficiency.  If all nodes assigned to RW_block are not on 
# a similar scale, we recommend providing an informed value for the "propCov" 
# control list argument, or using the AFSS sampler instead.


MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC, project = cModel)

samples <- runMCMC(cMCMC, niter = 50000, nburnin = 10000, thin = 10, nchains = 2)
#save(samples, file = "samples_guelat_kery_nimble.RData")
load("samples_guelat_kery_nimble.RData")

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

sampJoint <- rbind(samples$chain1, samples$chain2)

newlocs <- simgrid
dist11 <- fields::rdist(coords)
dist21 <- fields::rdist(newlocs, coords)
dist22 <- fields::rdist(newlocs)

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

#zNew <- runif(ncell(abu))  # fake new z values

## Now compare uncompiled and compiled speed of carrying out the sampling.
xstar_samples <- get_samples(samples = sampJoint[seq(1,8000,10),], z = data$z, zNew = covx[[1]], 
                             dist11 = dist11, dist22 = dist22, dist21 = dist21)
cget_samples <- compileNimble(get_samples)
set.seed(1)
xstar_samples2 <- cget_samples(samples = samples2[seq(1,90000,100),], z = data$z, zNew = temper[], 
                               dist11 = dist11, dist22 = dist22, dist21 = dist21)

alphaEst <- mean(sampJoint[,1])
betaEst <- mean(sampJoint[,2])
spatR <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
spatR[cellFromXY(spatR,simgrid)] <- colMeans(xstar_samples)
predRNim <- exp(alpha + beta*covxR + spatR) / (1 + exp(alpha + beta*covxR + spatR))

par(mfrow = c(2,4))
plot(covxR, main = "Temperature")
plot(rhoR, main = "Pure spatial effect, rho")
plot(psiR, main = "True Psi")
#points(xyFromCell(Y,sampID), pch = 16, cex = .6)
plot(presenceR, main = "True presence")
plot(predR, main = "Prediction with spOccupancy")
plot(predR[], psiR[])
plot(predRNim, main = "Prediction with NIMBLE")
plot(predRNim[], psiR[])







