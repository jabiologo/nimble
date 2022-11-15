library(sf)
library(nimble)
library(dplyr)
library(coda)
library(ade4)
library(Hmisc)
library(RCurl)

x <- getURL("https://raw.github.com/aronlindberg/latent_growth_classes/master/LGC_data.csv")
y <- read.csv(text = x)

hg_ <- read_sf("/home/javifl/margaritaSalas/ciencia/clmHunting/gis/hgNew.shp")
hg_$CERRAMIENT <- as.numeric(hg_$CERRAMIENT)
grid_ <- read_sf("/home/javifl/margaritaSalas/ciencia/clmHunting/gis/grNew.shp")
hg_$areaS <- scale(hg_$area)[,1]
#hg_$area2 <- scale(hg_$area*hg_$area)[,1]

# Data used for model calibration
dcal <- slice_sample(hg_, prop = 0.7)

constants <- list(ncell = nrow(grid_),
                  nhg = nrow(dcal),
                  low = dcal$minId,
                  high = dcal$maxId)

data <- list(harvest = dcal$JABALI,
             X1 = grid_$cover,
             X2 = grid_$climate,
             W1 = dcal$Axis1,
             W2 = dcal$CERRAMIENT,
             W3 = dcal$areaS)#,
#constraint1 = rep(1,nrow(dcal)))

m1 <- nimbleCode( {
  
  # PRIORS
  b0 ~ dnorm(0, 100)
  b1 ~ dnorm(0, 100) 
  b2 ~ dnorm(0, 100)
  a0 ~ dnorm(0, 100)
  a1 ~ dnorm(0, 100)
  a2 ~ dnorm(0, 100)
  a3 ~ dnorm(0, 100)
  
  # LIKELIHOOD
  for(i in 1:ncell){
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- b0 + b1*X1[i] + b2*X2[i]
  }
  for(h in 1:nhg){
    harvest[h] ~ dpois(lambda_hg[h]*(S[h]))
    log(lambda_hg[h]) <- log(sum(lambda[low[h]:high[h]]))
    logit(S[h]) <- a0 + a1*W1[h] + a2*W2[h] + a3*W3[h]
    #constraint1[h] ~ dconstraint(logit(S[h]) < 0.6)
  }
} )

keepers <- c("b0", "b1", "b2", "a0", "a1", "a2", "a3")

nc <- 3
nb <- 10000 
ni <- nb + 100000

model <- nimbleModel(code = m1, 
                     data = data, 
                     constants = constants,
                     calculate = FALSE)
c_model <- compileNimble(model)
model_conf <- configureMCMC(model, useConjugacy = FALSE)
model_conf$addMonitors(keepers)
model_mcmc <- buildMCMC(model_conf)
c_model_mcmc <- compileNimble(model_mcmc, project = model)
samplesCoS <- runMCMC(c_model_mcmc, nburnin = nb, niter = ni, nchains = nc)

# Trace plots
#sampCoS <- coda::mcmc(samplesCoS)
sampCoS <- coda::as.mcmc.list(lapply(samplesCoS, coda::mcmc))
par(mfrow=c(3,3))
traceplot(sampCoS)
# Merge chains
sampCoS <- rbind(sampCoS$chain1, sampCoS$chain2, sampCoS$chain3)

# Store coefficients for prediction (mean)
a0 <- mean(sampCoS[,1])
a1 <- mean(sampCoS[,2])
a2 <- mean(sampCoS[,3])
a3 <- mean(sampCoS[,4])
b0 <- mean(sampCoS[,5])
b1 <- mean(sampCoS[,6])
b2 <- mean(sampCoS[,7])
print(c(a0, a1, a2, a3, b0, b1, b2))

# Abundance at grid cell level
grid_$m1 <- exp(b0 + b1*grid_$cover + b2*grid_$climate)

# Abundance aggregation at hunting ground level
ag <- aggregate(x = grid_$m1, by = list(grid_$Hg), FUN = sum)
names(ag)[2] <- "aggHg"

# Hunting pressure / catch rate for each hunting ground
hg_$p1 <- exp(a0 + a1*hg_$Axis1 + a2*hg_$CERRAMIENT + a3*hg_$areaS) / (1 + exp(a0 + a1*hg_$Axis1 + a2*hg_$CERRAMIENT + a3*hg_$areaS))
hg_ <- left_join(hg_, ag, by = c("Hg" = "Group.1") )
# Hunting yield / harvest prediction at hunting ground level 
hg_$m1Pred <- hg_$aggHg*hg_$p1
#par(mfrow=c(1,1))
#hist(hg_$p1, breaks = 50)

# Validation dataset
dval <- hg_[!(hg_$Hg %in% dcal$Hg),]

# Validations with harvest data
par(mfrow=c(1,1))
plot(dval$m1Pred, dval$JABALI, xlim = c(0,350), ylim = c(0,350), pch = 16, cex = .5)
abline(a=0,b=1, col = "red")
cor(dval$m1Pred, dval$JABALI)

s_class<-cut2(dval$m1Pred, g=9) # defining bins (percentiles) on the predicted
s_mean<-as.matrix(tapply(dval$m1Pred, s_class,mean)) # calculating mean values for predicted
y_mean<-as.matrix(tapply(dval$JABALI,s_class,mean)) # calculating mean values for observed

plot(s_mean,y_mean,pch=19, # the calibration plot
     ylab="Observed hunting bag",
     xlab="Predicted hunting bag",
     xlim=c(0,max(s_mean, y_mean)),ylim=c(0,max(y_mean, s_mean)),cex=1.2,
     cex.lab=1.3,cex.axis=1.2)
abline(a=0,b=1,col="black",cex=2) # identity

# Prediction to shapefiles
#write_sf(grid_, "/home/javifl/margaritaSalas/ciencia/clmHunting/gis/grid_res.shp")
#write_sf(hg_, "/home/javifl/margaritaSalas/ciencia/clmHunting/gis/hg_res.shp")

