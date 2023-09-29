library(plgp)
library(mvtnorm)
library(raster)
library(spBayes)
library(spOccupancy)

# Random multivariate normal from Guelat & Kery 2018
rmvn <- function(n, mu = 0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V), p))))   stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}

distance <- as.matrix(dist(X))
rmnor <- rmvn(1, rep(0, n), 1.5 * exp(-0.05 * distance))
r[] <- rmnor
plot(r)
Y <- rmvnorm(n=1, mean=rep(0, 2500), sigma=Sigma)
r[] <- Y
plot(r)

df <- data.frame(simSig2 = rep(1.5, 50),
                 simPhi = rep(0.5, 50),
                 simB0 = rep(-1.2, 50),
                 simB1 = rep(2, 50),
                 simP = rep(0.7, 50),
                 estSig2 = rep(NA, 50),
                 estPhi = rep(NA, 50),
                 estB0 = rep(NA, 50),
                 estB1 = rep(NA, 50),
                 estP = rep(NA, 50))



for (q in 2:50){
  
  print(q)
  
  nx <- 50
  x <- seq(0, 50, length=nx)
  X <- expand.grid(x, x)
  
  sig2 <- 1.5 
  phi <- 0.05 
  
  Sigma <- mkSpCov(coords = as.matrix(X), K = as.matrix(sig2), 
                   Psi = as.matrix(0), theta = phi, cov.model = "exponential")
  
  Y <- rmvnorm(n=1, mean=rep(0, 2500), sigma=Sigma)
  
  r <- raster(nrows = 50, ncols = 50, xmn = 0, xmx = 50, ymn = 0, ymx = 50)
  r[cellFromXY(r,X)] <- Y
  
  temper <- raster(nrows = 50, ncols = 50, xmn = 0, xmx = 50, ymn = 0, ymx = 50)
  temper[] <- rnorm(2500)
  names(temper) <- "temper"
  
  psi <- exp(r + (-1.2+2*temper))/(1+exp(r + (-1.2+2*temper)))
  
  # Now we can draw the actual distribution thought a random binomial number using 
  # this probability
  distr <- psi
  distr[] <- rbinom(2500, 1, psi[])
  
  # Plot
  # par(mfrow = c(2,2))
  # plot(r, main = "Pure spatial")
  # plot(temper, main = "temper")
  # plot(psi, main = "Psi")
  # plot(distr, main = "distribution")
  
  n <- 500 
  sampID <- sample(1:2500, n)
  coords <- xyFromCell(distr, sampID) 
  occ.covs <- matrix(temper[sampID])
  colnames(occ.covs) <- "temper"
  y <- matrix(NA, n, 4)
  
  for (j in 1:4){
    for (i in 1:length(sampID)){
      y[i,j] <- rbinom(1, extract(distr,sampID[i]), 0.7)
    }
  }
  
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
                     w = rep(0, n), 
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
  df$estB0[q] <- mean(mSpatial$beta.samples[,1])
  df$estB1[q] <- mean(mSpatial$beta.samples[,2])
  df$estP[q] <- mean(mSpatial$alpha.samples)
  df$estSig2[q] <- mean(mSpatial$theta.samples[,1])
  df$estPhi[q] <- mean(mSpatial$theta.samples[,2])
  
}
save(df, file = "/home/javifl/iSDM/simus/spOccupancy_simus.RData")




























