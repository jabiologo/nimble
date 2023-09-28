# Loading packages
library(plgp)
library(mvtnorm)
library(raster)
library(spBayes)
library(nimble)

# Trying to fit a spatial Poisson model (abundance) with Gaussian Process in NIMBLE

# References
# Gaussian process regression: https://bookdown.org/rbg/surrogates/chap5.html
# Jeff Dosser spOccupancy manual: https://www.jeffdoser.com/files/spoccupancy-web/articles/modelfitting#spPGOcc
# https://r-nimble.org/nimbleExamples/gaussian_process.html

# Setting seed
set.seed(1)

# Creating a 30x30 grid (study area)
nx <- 50
x <- seq(0, 50, length=nx)
X <- expand.grid(x, x)

# Setting sig2 and phi parameters. With this parameters we'll build the spatial 
# covariance matrix Sigma, following an exponential model (others are allowed)
sig2 <- 1.5 # This is the variance, value in the matrix diagonal
phi <- 0.05 # Spatial decay parameter, how fast covariance decay 
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
Y <- rmvnorm(n=1, mean=rep(0, 2500), sigma=Sigma)

# This would be the "latent" or purely spatial component of the species distribution,
# not affected by any other covariate
r <- raster(nrows = 50, ncols = 50, xmn = 0, xmx = 50, ymn = 0, ymx = 50)
r[cellFromXY(r,X)] <- Y
#plot(r)
# To better understand the effect of sig2 and phi you can change those values and
# plot also several r

# Now we will create an environmental covariate that affects the distribution
#temper <- raster(nrows = 3, ncols = 3, xmn = 0, xmx = 3, ymn = 0, ymx = 3)
#temper[] <- rnorm(9)
#temper <- disaggregate(temper, 10, "bilinear")
temper <- raster(nrows = 50, ncols = 50, xmn = 0, xmx = 50, ymn = 0, ymx = 50)
temper[] <- rnorm(2500)
names(temper) <- "temper"

# Finally, following the Jeff Dosser formulation we'll create the distribution
# probability with a logit link: 
# https://www.jeffdoser.com/files/spoccupancy-web/articles/modelfitting#spPGOcc
lambda <- exp(1.5*r + (-1.5 +0.8*temper))

# Now we can draw the actual distribution thought a random binomial number using 
# this probability
abu <- lambda
abu[] <- rpois(2500, lambda[])

# Plot
par(mfrow = c(2,2), mai = c(0.5, 0.5, 0.5, 0.5), oma = c(1, 1, 1, 1))
plot(r, main = "Pure spatial")
plot(temper, main = "temper")
plot(lambda, main = "Lambda")
plot(abu, main = "abundance")

################################################################################

# Now we'll simulate a sampling process: we'll visit 200 places with a fixed
# detection probability of 0.7
n <- 500 # Number of visited places
sampID <- sample(1:2500, n) # Random selection of n places
locs <- xyFromCell(abu, sampID) # Getting coordinates of those places
points(locs, pch = 16, cex = 0.5)
z <- temper[sampID] # Getting environmental covariates of those places
# We'll fit a simple Poisson model with perfect detectability
y <- matrix(NA, n, 1) # Set of detections/non-detections 

for (j in 1:1){
  for (i in 1:length(sampID)){
    y[i,j] <- rbinom(1, extract(abu,sampID[i]), 1)
  }
}

# NIMBLE

dists <- as.matrix(dist(locs))
dists <- dists / max(dists)  # normalize to max distance of 1

constants <- list(N = n, 
                  dists = dists, 
                  ones = rep(1, n))

data <- list(y = as.numeric(y),
             z = z)

expcov <- nimbleFunction(     
  run = function(dists = double(2), phi = double(0), sigma = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    sigma2 <- sigma*sigma
    for(i in 1:n)
      for(j in 1:n)
        result[i, j] <- sigma2*exp(-dists[i,j]*phi)
    return(result)
  })
cExpcov <- compileNimble(expcov)


code <- nimbleCode({
  # Priors
  beta0 ~ dnorm(0, 5)
  beta1 ~ dnorm(0, 5)
  g1 ~ dnorm(0, 5)
  mu0 ~ dnorm(0, 5) # We'll use a zero-mean spatial Gaussian Process (as spOccupancy)
  sigma ~ dunif(0, 10)
  phi ~ dunif(0, 10)
  
  # Latent spatial process
  mu[1:N] <- mu0*ones[1:N]
  cov[1:N, 1:N] <- expcov(dists[1:N, 1:N], phi, sigma)
  s[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
  # Likelihood
  # vectorized
  log(lambda[1:N]) <- beta0 + beta1 * z[1:N] + g1 * s[1:N]
  #y[1:N] ~ dpois_vec(lambda[1:N])
   for(i in 1:N) {
     # It could be also vectorized but we should write the function
     y[i] ~ dpois(lambda[i])
   }
})

inits <- list(mu0 = 0, beta1 = 0, sigma = 1, phi = 1, g1 = 0)

set.seed(1)

## setup initial spatially-correlated latent process values
inits$cov <- expcov(dists, inits$phi, inits$sigma)
inits$s <-  t(chol(inits$cov)) %*% rnorm(n)
inits$s <- inits$s[ , 1]  # so can give nimble a vector rather than one-column matrix

model <- nimbleModel(code, constants = constants, data = data, inits = inits)
cModel <- compileNimble(model)
conf <- configureMCMC(model)
conf$addMonitors('s')
#conf$removeSamplers('s[1:n]')
# Changing a tunable parameter in the adaptation of RW_block makes a big difference.
# adaptFactorExponent. Exponent controling the rate of decay of the scale 
# adaptation factor. See Shaby and Wells, 2011, for details. (default = 0.8)
# IT'S TOO SLOW. ALSO PROBLEMS STORING s()
# conf$addSampler('s[1:n]', 'RW_block', control = list(adaptFactorExponent = 0.25))

MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC, project = cModel)

samples <- runMCMC(cMCMC, niter = 50000, nburnin = 10000, thin = 10, nchains = 1)
#save(samples, file = "/home/javifl/iSDM/simus/samplesSp.RData")
library(MCMCvis)
MCMCsummary(samples, round = 2,
            params = c("beta1", "sigma", "phi", "mu0", "g1"))
MCMCtrace(samples,
          pdf = FALSE,
          ind = TRUE,
          Rhat = TRUE,
          n.eff = TRUE,
          params = c("beta1", "sigma", "phi", "mu0", "g1")) 

samples2 <- rbind(samples$chain1,samples$chain2)

colnames(samples2)
b1 <- mean(samples2[,1])
b1 <- mean(samples2[,2])
s <- colMeans(samples2[,5:204])

pred <- exp(b0 + b1*data$z + s)
plot(pred,data$y)
cor(pred,data$y)

# Predicting in new locations
newlocs <- xyFromCell(abu, 1:ncell(abu))
dist11 <- fields::rdist(locs)
dist21 <- fields::rdist(newlocs, locs)
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
      
      b0 <- samples[i, 1] # CHANGED
      b1 <- samples[i, 2] # CHANGED
      mu0 <- samples[i, 3] # CHANGED
      phi <- samples[i, 4] # CHANGED
      sigma <- samples[i, 205] # CHANGED
      s <- samples[i, 5:204] # CHANGED
      mu <- mu0 + b0 + b1*z # CHANGED
      muNew <- mu0 + b0 + b1*zNew # CHANGED
      # Get other parameters
      output[i, ] <- sample_xstar(s, mu, muNew, sigma, phi, dist11, dist22, dist21)
      print(i)
    }
    return(output)
  }
)

#zNew <- runif(ncell(abu))  # fake new z values

## Now compare uncompiled and compiled speed of carrying out the sampling.
set.seed(1)
xstar_samples <- get_samples(samples = samples2[seq(1,90000,50),], z = data$z, zNew = temper[], 
                             dist11 = dist11, dist22 = dist22, dist21 = dist21)
cget_samples <- compileNimble(get_samples)
set.seed(1)
xstar_samples2 <- cget_samples(samples = samples2[seq(1,90000,100),], z = data$z, zNew = temper[], 
                               dist11 = dist11, dist22 = dist22, dist21 = dist21)

save(xstar_samples, file = "/home/javifl/iSDM/simus/xstar_samples.RData")

xstar <- abu
xstar[] <- colMeans(xstar_samples)
predAbu <- exp(b0 + b1*temper + xstar)

dev.off()
par(mfrow = c(2,3))
plot(temper)
plot(r)
plot(lambda)
plot(abu)
plot(xstar)
plot(predAbu)


