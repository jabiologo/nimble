# Preliminary not-spatial simple simulation for data integration
library(nimble)
library(coda)

################################################################################
# Latent process (abundance). Two environmental predictor covariates (x1 and x2)
# b0 = 2; b1 = -0.5; b2 = 0.3
x1 <- runif(1000, -1,1)
x2 <- runif(1000, -1,1)
lam <- exp(2 -0.5*x1 + 0.3*x2)
N <- rpois(1000, lam)

################################################################################
# Dataset 1: Distance-sampling: abundances at cell level
samp1 <- round(runif(100,1,1000),0)
y1 <- N[samp1] + rnorm(length(samp1), 0, 2)

# Dataset 2: latrine counts at level cell
samp2 <- round(runif(500,1,1000),0)
# Following Cabezas-Díaz & Virgós (2023) https://doi.org/10.1016/j.ecolind.2022.109684
# c0 = 7.5; c1 = 2
y2 <- exp((7.5+log(N[samp2]))/2.2)/2 + rnorm(length(samp2), 0, 5)

# Dataset 3: Harvest data at hunting state level following an effort covariate
# Harvest = sum(abundance_cell) * thinning_parameter (thinned Poisson point model)
# logit(thinning_parameter) = a0 + a1*effort_covariate
# a0 = -1.4; a1 = 0.5

hg <- data.frame(N=N, hgID = rep(1:50, each = 20))
samp3 <- aggregate(hg$N, by = list(hg$hgID), FUN = sum)
names(samp3) <- c("ID", "N")
samp3$eff <- rnorm(nrow(samp3))
samp3$pres <- plogis(-1.4 + 0.5*samp3$eff)
samp3$y3 <- round(samp3$N * samp3$pres,0)
samp3$ini <- seq(from = 1, to = 1000-19, by = 20)
samp3$fin <- seq(from = 20, to = 1000, by = 20)
head(samp3)

################################################################################
# Full integrated model
# Some errors when trying to estimate c1 and c2 with non-informative priors
# Obtaining "warning: logProb of data blabla"
# Better estimates when used slightly informative priors (we know the slope c1
# should be > 0...). But with a weird behavior since abundance parameters (b0 and 
# b1) are well estimated, but c0 and c1 are not... Interesting related to the 
# question about how different data sources contributes to inform latent variable 
# estimates...
# Estimates improved (obviously?) when using very informative priors such as
# c0 ~ dnorm(7.5, 1)
# c1 ~ dnorm(2.2, 1)


constants <- list( ncell = length(N),
                   nsamp1 = length(samp1),
                   nsamp2 = length(samp2),
                   nsamp3 = nrow(samp3),
                   sampID1 = samp1,
                   sampID2 = samp2,
                   samp3ini = samp3$ini,
                   samp3fin = samp3$fin)

data <- list(y1 = y1,
             y2 = y2,
             y3 = samp3$y3,
             x1 = x1,
             x2 = x2,
             eff = samp3$eff)

m5 <- nimbleCode( {
  # PRIORS
  b0 ~ dnorm(0, 2)
  b1 ~ dnorm(0, 2) 
  b2 ~ dnorm(0, 2)
  c0 ~ dnorm(0, 10) 
  c1 ~ dunif(0, 10) # light informative priors (we know the slope should be >0)
  a0 ~ dnorm(0, 2)
  a1 ~ dnorm(0, 2)
  
  # LIKELIHOOD
  for(i in 1:ncell){
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- b0 + b1*x1[i] + b2*x2[i]
  }
  for(j in 1:nsamp1){
    y1[j] ~ dnorm(mean=lambda[sampID1[j]], sd=2)
  }
  for(k in 1:nsamp2){
    y2[k] ~ dnorm(mean=y2_1km[k]/2, sd=5) # I didn't know how to remove exp with this /2 so I divided it in two lines
    log(y2_1km[k]) <- c0+log(lambda[sampID2[k]])/c1
  }
  for(l in 1:nsamp3){
    y3[l] ~ dpois(lambdaHg[l]*S[l])
    lambdaHg[l] <- sum(lambda[samp3ini[l]:samp3fin[l]])
    logit(S[l]) <- a0 + a1*eff[l]
  }
} )

keepers <- c("b0", "b1", "b2", "c0", "c1", "a0", "a1")

nc <- 1 # tested with 3
nb <- 10000 # tested with 100000
ni <- nb + 100000 # tested with 1000000

model <- nimbleModel(code = m5, data = data, constants = constants, calculate = FALSE)
c_model <- compileNimble(model)
model_conf <- configureMCMC(model, useConjugacy = FALSE)
model_conf$addMonitors(keepers)
model_mcmc <- buildMCMC(model_conf)
c_model_mcmc <- compileNimble(model_mcmc, project = model)
samples <- runMCMC(c_model_mcmc, nburnin = nb, niter = ni, nchains = nc)

# Trace plots
samp <- coda::mcmc(samples)
#samp <- coda::as.mcmc.list(lapply(samples, coda::mcmc))
par(mfrow=c(3,3))
traceplot(samp)










