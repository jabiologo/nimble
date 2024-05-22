library(terra)
library(spAbundance)
# Download the Mona Lisa paint from https://github.com/jabiologo/nimble/blob/main/gio.tif 
# Load the raster
r <- rast("/mypath/gio.tif")
# Set number of samples
n <- 1000
# Sample the RGB values
samp <- sample(1:ncell(r), n)
# Plot the paint and the sampled pixels
plot(r)
points(xyFromCell(r,samp), cex = 0.3, col = "red")

#############################
# Multivariate Spatial GLMM #
#############################
# Prepare "count" data
counts <- t(as.matrix(cbind(r$gio_1[samp],r$gio_2[samp],r$gio_3[samp])))
# Explore species "counts" for the 10 first sampled pixels
counts[1:3,1:10]
# Obtain coordinates from sampled pixels
coords <- xyFromCell(r,samp)
# No fixed covariates, only intercept
covs <- list(int = rep(1,n))
# Wrap all input data
data.list <- list(y = counts, covs = covs, coords = coords)

# Setting-up priors 
#(check https://www.jeffdoser.com/files/spabundance-web/articles/glmm#spatial-factor-multivariate-glmms)
prior.list <- list(beta.comm.normal = list(mean = 0, var = 100),
                   phi.unif = list(a = 0.01, b = 3 / 0.5),
                   tau.sq.beta.ig = list(a = .1, b = .1))
# Setting-up inits
inits.list <- list(beta.comm = 0, beta = 0,
                   tau.sq.beta = 1, phi = 3 / 1)

# Running the model, it will take for a while... be patient!
out <- sfMsAbund(formula = ~ 1, data = data.list, n.batch = 10,
                 inits = inits.list, priors = prior.list,
                 NNGP = TRUE, cov.model = 'exponential',
                 n.neighbors = 15, n.factors = 3, batch.length = 1000,
                 n.omp.threads = 3, n.report = 1, n.burn = 1000,
                 n.thin = 100, n.chains = 1)

#summary(out)

# Setting-up data to predict 
coords.0 <- xyFromCell(r, 1:ncell(r))
X.0 <- matrix(rep(1,ncell(r)), ncell(r),1)
colnames(X.0) <- "int"
# Doing predictions, it will take for a while...  be patient!
out.pred <- predict(out, X.0, coords.0, verbose = TRUE)

# Sorting predictions
plot.dat <- data.frame(x = coords.0[,1], 
                       y = coords.0[,2], 
                       mean.abu1 = apply(out.pred$y.0.samples[,1,,], 2, mean),
                       mean.abu2 = apply(out.pred$y.0.samples[,2,,], 2, mean),
                       mean.abu3 = apply(out.pred$y.0.samples[,3,,], 2, mean))
# Create raster with three layers to store and plot predicted "abundances"
r.pred <- r
r.pred$gio_1[] <- plot.dat$mean.abu1
r.pred$gio_2[] <- plot.dat$mean.abu2
r.pred$gio_3[] <- plot.dat$mean.abu3
# Creating RGB colors from predicted abundances 
RGB(r.pred) <- c(1,2,3)

# Plot observed vs predicted RGB colors
par(mfrow = c(1,2))
plot(r)
text(200,30, "Sampled pixels (input)")
points(xyFromCell(r,samp), cex = 0.3, col = "red")
plot(r.pred)
text(200,30, "Predicted abundance (merged RGB)")