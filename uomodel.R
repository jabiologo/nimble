library(nimble)
library(nimbleSCR)

# Data simulation for harvest model with harvest over-reporting

# x1 is a random predictor of abundance
x1 <- scale(rnorm(100,20,5))[1:100]
# lambda simulation
lambda <- exp(5 + 0.7*(x1))
# Number of rabbits at each site
N <- rpois(100,lambda)
# Number of harvest rabbits (we capture 33% from each place)
H <- rpois(100, N*0.33)
# Land use: 1 = crops; 0 = others
luse <- rbinom(100,1,0.5)
# Under/Over (in percentage) reporting depends on land use
# Sites with crops reports 50% more harvest rabbits
# (They want to kill them all so they report more than they actually hunt)
a <- 0.5
UO = 1 + a*luse
# "Trustability" is a normal random number with mean in UO
Tr <- rnorm(100,UO,0.01)
# Reported is based in harvest times "trustability"
R <- rpois(100, H*Tr)

simdata <- data.frame(cbind(x1, lambda, N,H,luse,UO,Tr,H*Tr,R))
#head(simdata)
ourdata <- simdata[,c(1,5,9)]
head(ourdata)

# Beta moments
betaM <- function(m,var){
  alpha <- ((1-m)/var) - (1/m)
  beta <- alpha * ((1/m) - 1)
  return(c(alpha, beta))
}

betaM(0.33,0.01)
################################################################################
# NIMBLE model

constants <- list(n = 100,
                  x1 = x1,
                  luse = luse)
data <- list(R = R)

report <- nimble::nimbleCode( {
  # PRIORS
  b0 ~ dnorm(0, 10)
  b1 ~ dnorm(0, 10)
  h ~ dbeta(1, 1)
  #h ~ dbeta(5, 5) 
  #h ~ dbeta(63.9697, 129.8779) # very informative prior LOL
  a1 ~ dnorm(0, 10) 
  
  # LIKELIHOOD
  for(i in 1:n){
    log(lambda[i]) <- b0 + b1*x1[i]
    N[i] ~ dpois(lambda[i])
  }
  for(j in 1:n){
    H[j] ~ dpois(h*N[j])
  }
  for(k in 1:n){
    R[k] ~ dpois((1+a1*luse[k])*H[k])
    }
} )

inits <- function() {
  base::list(N =  rep(1, constants$n),
             h = runif(1, 0, 1),
             a0 = runif(1, -1, 1),
             a1 = runif(1, -1, 1)
             
  )
}

keepers <- c("a1", "h", "N", 'H', "b0", "b1")

nc <- 3 
nb <- 10000 
ni <- nb + 100000

# Now he create the model
model <- nimble::nimbleModel(code = report, 
                             data = data, 
                             constants = constants, 
                             inits = inits(),
                             calculate = FALSE)
model$initializeInfo()
c_model <- nimble::compileNimble(model)
model_conf <- nimble::configureMCMC(model,
                                    useConjugacy = FALSE)
model_conf$addMonitors(keepers)
model_mcmc <- nimble::buildMCMC(model_conf)
c_model_mcmc <- nimble::compileNimble(model_mcmc, project = model)
samples_ou <- nimble::runMCMC(c_model_mcmc, 
                           nburnin = nb, 
                           niter = ni, 
                           nchains = nc)



samples_mcmc_ou <- coda::as.mcmc.list(lapply(samples_ou, coda::mcmc))
par(mfrow=c(2,2))
coda::traceplot(samples_mcmc_ou[,201:204])
# Calculate Rhat convergence diagnostic for the three parameters
coda::gelman.diag(samples_mcmc_ou[,201:204])


################################################################################
