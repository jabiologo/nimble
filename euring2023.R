library(sf)
library(nimble)
library(dplyr)
library(coda)
library(ade4)
library(Hmisc)
library(RCurl)

set.seed(1)
# Load the data
hg <- getURL("https://raw.githubusercontent.com/jabiologo/nimble/main/data/hgNew.csv")
gr <- getURL("https://raw.githubusercontent.com/jabiologo/nimble/main/data/grNew.csv")
hg_ <- read.csv(text = hg)
grid_ <- read.csv(text = gr)

# Data used for model calibration
dcal <- slice_sample(hg_, prop = 0.7)

constants <- list(ncell = nrow(grid_),
                  nhg = nrow(dcal),
                  low = dcal$minId,
                  high = dcal$maxId)

data <- list(harvest = dcal$JABALI, # hunting ground level: wild boar harvest
             X1 = grid_$cover,      # grid cell level: land cover
             X2 = grid_$climate,    # grid cell level: climate
             W1 = dcal$Axis1,       # hunting ground level: continuous type
             W2 = dcal$CERRAMIENT,  # hunting ground level: fenced/not-fenced
             W3 = dcal$areaS)#,     # hunting ground level: area
#constraint1 = rep(1,nrow(dcal)))  # using dconstrain()

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
    #constraint1[h] ~ dconstraint(logit(S[h]) < 0.6)        # using dconstrain()
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
sampCoS <- coda::as.mcmc.list(lapply(samplesCoS, coda::mcmc))
par(mfrow=c(3,3))
traceplot(sampCoS)
# Calculate Rhat convergence diagnostic
gelman.diag(sampCoS)
effectiveSize(sampCoS)
# Merge chains
sampCoS <- mcmc(rbind(sampCoS$chain1, sampCoS$chain2, sampCoS$chain3))

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
#abline(v=0.16, col = "red")
#abline(v=0.69, col = "red")
#abline(v=0.33, col = "blue")

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

################################################################################
# Uncertainty
# Some diagnostic plots (posterior probs)
library(plotMCMC)
library(bayestestR)
library(plotrix)
plotAuto(sampCoS)
plotDens(sampCoS)
plotQuant(sampCoS)

samp2 <- sample_n(data.frame(sampCoS), 10000)

a0 <- samp2[,1]
a1 <- samp2[,2]
a2 <- samp2[,3]
a3 <- samp2[,4]
b0 <- samp2[,5]
b1 <- samp2[,6]
b2 <- samp2[,7]

# Abundance at grid cell level
predDens <- data.frame(exp(b0 + b1%o%grid_$cover + b2%o%grid_$climate))
predDensCI <- sapply(predDens, ci)
predDensSE <- sapply(predDens, std.error)
grid_$m1low <- as.numeric(predDensCI[2,])
grid_$m1high <- as.numeric(predDensCI[3,])
grid_$m1se <- predDensSE

# Abundance aggregation at hunting ground level
agLow <- aggregate(x = grid_$m1low, by = list(grid_$Hg), FUN = sum)
names(agLow)[2] <- "agLow"
agHigh <- aggregate(x = grid_$m1high, by = list(grid_$Hg), FUN = sum)
names(agHigh)[2] <- "agHigh"
hg_ <- left_join(hg_, agLow, by = c("Hg" = "Group.1"))
hg_ <- left_join(hg_, agHigh, by = c("Hg" = "Group.1"))

# Hunting pressure / catch rate for each hunting ground
pDens <- data.frame(exp(a0 + a1%o%hg_$Axis1 + a2%o%hg_$CERRAMIENT + a3%o%hg_$areaS) / (1 + exp(a0 + a1%o%hg_$Axis1 + a2%o%hg_$CERRAMIENT + a3%o%hg_$areaS)))
pDensCI <- sapply(pDens, ci)
hg_$p1low <- as.numeric(pDensCI[2,])
hg_$p1high <- as.numeric(pDensCI[3,])

# Hunting yield / harvest prediction at hunting ground level 
hg_$m1PredLow <- hg_$agLow*hg_$p1low
hg_$m1PredHigh <- hg_$agHigh*hg_$p1high
#par(mfrow=c(1,1))
#hist(hg_$p1, breaks = 50)

# Validation dataset
dval <- hg_[!(hg_$Hg %in% dcal$Hg),]

s_class<-cut2(dval$m1Pred, g=9) # defining bins (percentiles) on the predicted
s_mean<-as.matrix(tapply(dval$m1Pred, s_class,mean)) # calculating mean values for predicted
x_mean<-as.matrix(tapply(dval$JABALI,s_class,mean)) # calculating mean values for observed

s_meanL<-as.matrix(tapply(dval$m1PredLow, s_class,mean)) # calculating mean values for predicted
s_meanH<-as.matrix(tapply(dval$m1PredHigh, s_class,mean)) # calculating mean values for predicted

dfVal <- data.frame(x_mean, s_mean, s_meanL, s_meanH)

plot(x_mean,s_mean,pch=19, # the calibration plot
     xlab="Observed hunting bag",
     ylab="Predicted hunting bag",
     xlim=c(0,max(s_meanH, x_mean)),ylim=c(0,max(x_mean, s_meanH)),cex=1.2,
     cex.lab=1.3,cex.axis=1.2)
abline(a=0,b=1,col="black",cex=2) # identity
points(x_mean,s_meanL)
points(x_mean,s_meanH)

# Basic scatter plot
ggplot(dfVal, aes(x=x_mean, y=s_mean)) + 
  geom_point(size = 3) +
  xlim(0,max(s_meanH)) + 
  ylim(0,max(s_meanH)) +
  geom_pointrange(aes(ymin = s_meanL, ymax = s_meanH)) + 
  geom_abline(intercept = 0, slope = 1) +
  labs(title="", x="Observed harvest", y = "Predicted harvest") +
  theme(axis.text=element_text(size=14),
      axis.title=element_text(size=18, margin = margin(t = 0, r = 15, b = 0, l = 0)),
      axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
      axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
      legend.position="none", 
      plot.title=element_text(hjust=0.5,size=27,face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white",
                                      colour = "black",
                                      size = 1, linetype = "solid"))

sum(grid_$m1low)
sum(grid_$m1)
sum(grid_$m1high)
