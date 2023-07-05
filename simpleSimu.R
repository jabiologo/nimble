library(nimble)
library(coda)

x1 <- runif(1000, -1,1)
x2 <- runif(1000, -1,1)
lam <- exp(2 -0.5*x1 + 0.3*x2)
N <- rpois(1000, lam)

################################################################################
# Dataset 1

samp1 <- round(runif(100,1,1000),0)
y1 <- N[samp1] + rnorm(length(samp1), 0, 2)

constants <- list( ncell = length(N),
                    nsamp1 = length(samp1),
                    sampID1 = samp1)
data <- list(y1 = y1,
             x1 = x1,
             x2 = x2)

m1 <- nimbleCode( {
  # PRIORS
  b0 ~ dnorm(0, 10)
  b1 ~ dnorm(0, 10) 
  b2 ~ dnorm(0, 10)
  
  # LIKELIHOOD
  for(i in 1:ncell){
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- b0 + b1*x1[i] + b2*x2[i]
  }
  for(j in 1:nsamp1){
    y1[j] ~ dnorm(mean=lambda[sampID1[j]], sd=2)
  }
} )

keepers <- c("b0", "b1", "b2")

nc <- 1
nb <- 10000
ni <- nb + 50000

model <- nimbleModel(code = m1, data = data, constants = constants, calculate = FALSE)
c_model <- compileNimble(model)
model_conf <- configureMCMC(model, useConjugacy = FALSE)
model_conf$addMonitors(keepers)
model_mcmc <- buildMCMC(model_conf)
c_model_mcmc <- compileNimble(model_mcmc, project = model)
samples <- runMCMC(c_model_mcmc, nburnin = nb, niter = ni, nchains = nc)

# Trace plots
samp <- coda::mcmc(samples)
#sampOc <- coda::as.mcmc.list(lapply(samplesOc, coda::mcmc))
par(mfrow=c(3,3))
traceplot(samp)


################################################################################
# Dataset 2

# Include sd maybe doesn't make sense... it's just an error that we can obtain from
# the data like in Distance sampling?
# Initial values https://oliviergimenez.github.io/nimble-workshop/#114
# Does't converge when including c0 and c1 and sd is fixed
# It kind of improves if we include informative priors (not really)
# If we fix c1 or c0, does converge but with poor fit
# It makes sense the model fail... Actually we are asking many things from a poor
# informative dataset. We are trying to obtain information about 2 different process:
# 1) How does the letrines relates to abundance?
# 2) Hoe does the abundance relates to predictor covariates?
# Only data about latrine have not enough info to describe those relations using this model
# Maybe we could relate latrines directly with predictor covariates...?


samp2 <- round(runif(500,1,1000),0)
y2 <- exp((7.5+log(N[samp2]))/2.2)/2 + rnorm(length(samp2), 0, 5)

constants <- list( ncell = length(N),
                   nsamp2 = length(samp2),
                   sampID2 = samp2)
data <- list(y2 = y2,
             x1 = x1,
             x2 = x2)

m2 <- nimbleCode( {
  # PRIORS
  b0 ~ dnorm(0, 10)
  b1 ~ dnorm(0, 10) 
  b2 ~ dnorm(0, 10)
  #c0 ~ dnorm(0, 10)
  #c1 ~ dnorm(0, 10)
  c0 ~ dnorm(7.5, 1)
  c1 ~ dnorm(2.2, 1)
  
  # LIKELIHOOD
  for(i in 1:ncell){
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- b0 + b1*x1[i] + b2*x2[i]
  }
  for(j in 1:nsamp2){
    y2[j] ~ dnorm(mean=exp((c0+log(lambda[sampID2[j]]))/c1)/2, sd=5)
    #y2[j] ~ dpois(log(exp((c0+log(lambda[sampID2[j]]))/c1)/2))
  }
} )

keepers <- c("b0", "b1", "b2", "c0", "c1")

nc <- 1
nb <- 100000
ni <- nb + 500000

model <- nimbleModel(code = m2, data = data, constants = constants, calculate = FALSE)
c_model <- compileNimble(model)
model_conf <- configureMCMC(model, useConjugacy = FALSE)
model_conf$addMonitors(keepers)
model_mcmc <- buildMCMC(model_conf)
c_model_mcmc <- compileNimble(model_mcmc, project = model)
samples <- runMCMC(c_model_mcmc, nburnin = nb, niter = ni, nchains = nc)

# Trace plots
samp <- coda::mcmc(samples)
#sampOc <- coda::as.mcmc.list(lapply(samplesOc, coda::mcmc))
par(mfrow=c(3,3))
traceplot(samp)


################################################################################
# Integrated model
# warning: logProb of data strongly depends on initial values
# Initial values https://oliviergimenez.github.io/nimble-workshop/#114
# It works very good if we fix c0 or c1, but doesn't converge if we estimate both

constants <- list( ncell = length(N),
                   nsamp1 = length(samp1),
                   nsamp2 = length(samp2),
                   sampID1 = samp1,
                   sampID2 = samp2)
data <- list(y1 = y1,
             y2 = y2,
             x1 = x1,
             x2 = x2)

m3 <- nimbleCode( {
  # PRIORS
  b0 ~ dnorm(0, 10)
  b1 ~ dnorm(0, 10) 
  b2 ~ dnorm(0, 10)
  c0 ~ dnorm(0, 10)
  c1 ~ dnorm(0, 10)
  
  # LIKELIHOOD
  for(i in 1:ncell){
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- b0 + b1*x1[i] + b2*x2[i]
  }
  for(j in 1:nsamp1){
    y1[j] ~ dnorm(mean=lambda[sampID1[j]], sd=2)
  }
  for(k in 1:nsamp2){
    y2[k] ~ dnorm(mean=exp((c0+log(lambda[sampID2[k]]))/c1)/2, sd=5)
  }
} )

keepers <- c("b0", "b1", "b2", "c0", "c1")

nc <- 1
nb <- 50000
ni <- nb + 100000

model <- nimbleModel(code = m3, data = data, constants = constants, calculate = FALSE)
c_model <- compileNimble(model)
model_conf <- configureMCMC(model, useConjugacy = FALSE)
model_conf$addMonitors(keepers)
model_mcmc <- buildMCMC(model_conf)
c_model_mcmc <- compileNimble(model_mcmc, project = model)
samples <- runMCMC(c_model_mcmc, nburnin = nb, niter = ni, nchains = nc)

# Trace plots
samp <- coda::mcmc(samples)
#sampOc <- coda::as.mcmc.list(lapply(samplesOc, coda::mcmc))
par(mfrow=c(3,3))
traceplot(samp)

################################################################################
# Dataset 3

hg <- data.frame(N=N, hgID = rep(1:50, each = 20))
samp3 <- aggregate(hg$N, by = list(hg$hgID), FUN = sum)
names(samp3) <- c("ID", "N")
samp3$eff <- rnorm(nrow(samp3))
samp3$pres <- plogis(-1.4 + 0.5*samp3$eff)
samp3$y3 <- round(samp3$N * samp3$pres,0)
samp3$ini <- seq(from = 1, to = 1000-19, by = 20)
samp3$fin <- seq(from = 20, to = 1000, by = 20)
head(samp3)

constants <- list( ncell = length(N),
                   nsamp3 = nrow(samp3),
                   samp3ini = samp3$ini,
                   samp3fin = samp3$fin)

data <- list(y3 = samp3$y3,
             x1 = x1,
             x2 = x2,
             eff = samp3$eff)

m4 <- nimbleCode( {
  # PRIORS
  b0 ~ dnorm(0, 10)
  b1 ~ dnorm(0, 10) 
  b2 ~ dnorm(0, 10)
  a0 ~ dnorm(0, 10)
  a1 ~ dnorm(0, 10)
  
  # LIKELIHOOD
  for(i in 1:ncell){
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- b0 + b1*x1[i] + b2*x2[i]
  }
  for(l in 1:nsamp3){
    y3[l] ~ dpois(lambdaHg[l]*S[l])
    lambdaHg[l] <- sum(lambda[samp3ini[l]:samp3fin[l]])
    logit(S[l]) <- a0 + a1*eff[l]
  }
} )

keepers <- c("b0", "b1", "b2", "a0", "a1")

nc <- 1
nb <- 10000
ni <- nb + 50000

model <- nimbleModel(code = m4, data = data, constants = constants, calculate = FALSE)
c_model <- compileNimble(model)
model_conf <- configureMCMC(model, useConjugacy = FALSE)
model_conf$addMonitors(keepers)
model_mcmc <- buildMCMC(model_conf)
c_model_mcmc <- compileNimble(model_mcmc, project = model)
samples <- runMCMC(c_model_mcmc, nburnin = nb, niter = ni, nchains = nc)

# Trace plots
samp <- coda::mcmc(samples)
#sampOc <- coda::as.mcmc.list(lapply(samplesOc, coda::mcmc))
par(mfrow=c(3,3))
traceplot(samp)

################################################################################
# Full integrated model
# Some errors when trying to estimate c1 and c2
# Estimates improved when using very informative priors c0 ~ dnorm(7.5, 1)
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
  c0 ~ dnorm(7.5, 1)
  c1 ~ dnorm(2.2, 1)
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
    y2[k] ~ dnorm(mean=exp((c0+log(lambda[sampID2[k]]))/c1)/2, sd=5)
  }
  for(l in 1:nsamp3){
    y3[l] ~ dpois(lambdaHg[l]*S[l])
    lambdaHg[l] <- sum(lambda[samp3ini[l]:samp3fin[l]])
    logit(S[l]) <- a0 + a1*eff[l]
  }
} )

keepers <- c("b0", "b1", "b2", "c0", "c1", "a0", "a1")

nc <- 1
nb <- 100000
ni <- nb + 1000000

model <- nimbleModel(code = m5, data = data, constants = constants, calculate = FALSE)
c_model <- compileNimble(model)
model_conf <- configureMCMC(model, useConjugacy = FALSE)
model_conf$addMonitors(keepers)
model_mcmc <- buildMCMC(model_conf)
c_model_mcmc <- compileNimble(model_mcmc, project = model)
samples <- runMCMC(c_model_mcmc, nburnin = nb, niter = ni, nchains = nc)

# Trace plots
samp <- coda::mcmc(samples)
#sampOc <- coda::as.mcmc.list(lapply(samplesOc, coda::mcmc))
par(mfrow=c(3,3))
traceplot(samp)

################################################################################

# Playing with latrine abundance relationship

N <- seq(1:50)
let1 <- exp((7.5+log(N))/2.2)/2# + rnorm(length(samp2), 0, 5)
plot(N, let1)
let2 <- N/0.106
lines(N,let2)


