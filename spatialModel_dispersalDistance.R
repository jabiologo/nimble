### Guelat & Kery: Effects of spatial autocorrelation and imperfect detection on species distribution models
### Adapted script to use with spOccupancy
library(terra)
library(spatstat)
library(dagR)
set.seed(3)

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
xdim <- 100
ydim <- 100
x.easting <- 1:xdim
x.northing <- 1:ydim
simgrid <- expand.grid(x.easting, x.northing)
n <- nrow(simgrid)
simgrid <- simgrid-0.5

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

# Expected lambda
lam <- mapply(function(x, y) exp(alpha + beta*x + spCor*y), covx, rho, SIMPLIFY = FALSE)
# Realized abundance
abun <- lapply(lam, function(x) rpois(n = n, lambda = x))

covxR <- rast(nrows=100, ncols=100, xmin=0, xmax=100, ymin = 0, ymax = 100)
covxR[cellFromXY(covxR,simgrid)] <- covx[[1]]
rhoR <- rast(nrows=100, ncols=100, xmin=0, xmax=100, ymin = 0, ymax = 100)
rhoR[cellFromXY(rhoR,simgrid)] <- rho[[1]] 
lamR <- rast(nrows=100, ncols=100, xmin=0, xmax=100, ymin = 0, ymax = 100)
lamR[cellFromXY(lamR,simgrid)] <- lam[[1]] 
abunR <- rast(nrows=100, ncols=100, xmin=0, xmax=100, ymin = 0, ymax = 100)
abunR[cellFromXY(abunR,simgrid)] <- abun[[1]]  

par(mfrow = c(2,2))
plot(rhoR, main = "Pure spatial effect, rho")
plot(covxR, main = "Temperature")
plot(lamR, main = "True Lambda")
plot(abunR, main = "True abundance")

lamI <- as.im(matrix(flip(lamR)[],100, 100, byrow = T),
              W = owin(xrange = c(0, 100), yrange = c(0, 100)))
pp <- rpoispp(lamI) 
plot(lamR)
points(pp, cex = 0.1)

move <- function(x, rate){
  df <- data.frame(x = NA, y = NA)
  for(i in 1:x$n){
    df[i,] <- anglePoint(c(x$x[i],x$y[i]), runif(1,0,360), rexp(1, rate = rate))
  }
  return(df)
}

ppM <- move(pp, 0.5)

plot(abunR)
points(pp, cex = 0.1)
points(ppM, cex = 0.1, col = "darkred")

count <- table(cellFromXY(abunR,ppM))

abunRM <- abunR
abunRM[] <- 0
abunRM[as.numeric(names(count))] <- as.numeric(count)
plot(abunR)
plot(abunRM)

sampID <- sample(1:10000, 1000)
coord <- xyFromCell(abunR, sampID)

df <- data.frame(cbind(abunR[sampID], abunRM[sampID], covxR[sampID], rhoR[sampID]))
colnames(df) <- c("abOr", "abMo", "cov", "rho")
m1 <- glm(abOr ~ cov + rho, data = df, family = poisson(link = "log"))
m2 <- glm(abMo ~ cov + rho, data = df, family = poisson(link = "log"))
summary(m1)
summary(m2)

# Extract the standardized residuals
# http://rstudio-pubs-static.s3.amazonaws.com/9688_a49c681fab974bbca889e3eae9fbb837.html
# https://gwenantell.com/exploring-spatial-autocorrelation-in-r/
library(sp)

rs <- rstandard(m2)
# Create a spatial object (SpatialPointsDataFrame)
spdata <- data.frame(resid = rs, x = coord[,1], y = coord[,2])
coordinates(spdata) <- c("x", "y")

bubble(spdata, "resid", col = c("blue", "orange"), main = "Residuals", xlab = "X-coordinates", 
       ylab = "Y-coordinates")

library(ncf)
fit1 <- correlog(x = spdata$x, y = spdata$y, z = spdata$resid, increment = 2, 
                 resamp = 50, quiet = TRUE)
plot(fit1)


variocloud <- variogram(resid ~ 1, spdata, cloud = TRUE)
plot(variocloud)
vario <- variogram(resid ~ 1, spdata, cutoff = 30)
plot(vario, pch = 16, cex = 1.5)
