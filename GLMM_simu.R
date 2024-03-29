# GLMM simulation and fitting (random intercept)
# 29-03-2024

library(fastDummies)

df <- data.frame(x1 = rnorm(1300), 
                 gr = as.factor(sample(LETTERS, 1300, replace = TRUE)))
df$grN <- as.numeric(df$gr)

r1 <- rnorm(26, mean = 0, sd = 1.9)
r1D <- dummy_columns(df$gr)[,-1]
beta1 <- 0.6

mu <- (r1[1]*r1D[,1] + 
         r1[2]*r1D[,2] +
         r1[3]*r1D[,3] +
         r1[4]*r1D[,4] +
         r1[5]*r1D[,5] +
         r1[6]*r1D[,6] +
         r1[7]*r1D[,7] +
         r1[8]*r1D[,8] +
         r1[9]*r1D[,9] +
         r1[10]*r1D[,10] +
         r1[11]*r1D[,11] +
         r1[12]*r1D[,12] +
         r1[13]*r1D[,13] +
         r1[14]*r1D[,14] +
         r1[15]*r1D[,15] +
         r1[16]*r1D[,16] +
         r1[17]*r1D[,17] +
         r1[18]*r1D[,18] +
         r1[19]*r1D[,19] +
         r1[20]*r1D[,20] +
         r1[21]*r1D[,21] +
         r1[22]*r1D[,22] +
         r1[23]*r1D[,23] +
         r1[24]*r1D[,24] + 
         r1[25]*r1D[,25] + 
         r1[26]*r1D[,26] + beta1 * df$x1)

df$Y <- rnorm(1300, mean = mu, sd = 1)

m1 <- glm(Y ~ gr + x1 - 1, data = df)
summary(m1)

library(lme4)
m2 <- lmer(Y ~ x1 + (1|gr) -1 , data = df)
summary(m2)
plot(r1, m1$coefficients[1:26], xlab = "simulados", ylab = "predichos",
     xlim = c(-2.1,2.1), ylim = c(-2.1,2.1))
points(r1, m2@u, xlab = "simulados", ylab = "predichos", col = "red", cex = 0.5, pch = 19)
abline(a = 0, b = 1)

# NIMBLE

library(nimble)

constants <- list(nsamp = 1300,
                  ngr = 26,
                  group = df$grN)
data <- list(Y = df$Y,
             x1 = df$x1)

glmm1 <- nimbleCode({
  # PRIORS
  sd.y ~ dunif(0,100)
  b ~ dnorm(0, 100)
  # Prior para el factor aleatorio
  for(j in 1:ngr){
    a[j] ~ dnorm(mean = mu.a, sd = sd.a)
  }
  mu.a ~ dnorm(0, 100)
  sd.a ~ dunif(0,100)
  
  # LIKELIHOOD
  for(i in 1:nsamp){
    Y[i] ~ dnorm(mean = mu[i], sd.y)
    mu[i] <- a[group[i]] + b*x1[i]
  }
  })

inits <- function() {
  base::list(b = runif(1, -1, 1),
             sd.y = runif(1, 0, 100),
             mu.a = runif(1, -1, 1),
             sd.a = runif(1, 0, 100)
  )
}

keepers <- c("b", 'sd.y', "mu.a", "sd.a")

# Finally we define the settings of our MCMC algorithm
nc <- 1 # number of chains
nb <- 5000 # number of initial MCMC iterations to discard
ni <- nb + 50000 # total number  of iterations

model <- nimble::nimbleModel(code = glmm1, 
                             data = data, 
                             constants = constants, 
                             inits = inits(),
                             calculate = FALSE)

c_model <- nimble::compileNimble(model)
model_conf <- nimble::configureMCMC(model,
                                    useConjugacy = FALSE)
model_conf$addMonitors(keepers)
model_mcmc <- nimble::buildMCMC(model_conf)
c_model_mcmc <- nimble::compileNimble(model_mcmc, project = model)

# Run the MCMC
samples <- nimble::runMCMC(c_model_mcmc, 
                           nburnin = nb, 
                           niter = ni, 
                           nchains = nc)

MCMCvis::MCMCplot(samples)
MCMCvis::MCMCsummary(samples)
