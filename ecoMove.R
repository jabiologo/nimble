# ecoMove simulation
# Simulation of habitat suitability. Deployment of several activity centers based
# on habitat suitability. Finally, animal movement simulation under an OU model
# by using the activity center as the optimum.

# In a second stage, we should deploy camera trap stations in a regular grid or
# in different designs to estimate animal abundance by using REM or other methods

library(terra)
set.seed(1)

# Random multivariate normal
rmvn <- function(n, mu = 0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V), p))))   stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}


# Set up a square lattice region
xdim <- 50
ydim <- 50
x.easting <- 1:xdim
x.northing <- 1:ydim
simgrid <- expand.grid(x.easting, x.northing)
n <- nrow(simgrid)

# Set up distance matrix
distance <- as.matrix(dist(simgrid))

# Covariates
clim <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
phi_cl <- 0.05
var_cl <- 1.5
clim_values <- as.vector(rmvn(1, rep(0, n), var_cl * exp(-phi_cl * distance)))
clim[cellFromXY(clim,simgrid-0.5)] <- clim_values

cover <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
phi_co <- 2
var_co <- 1.5
cover_values <- as.vector(rmvn(1, rep(0, n), var_co * exp(-phi_co * distance)))
cover[cellFromXY(cover,simgrid-0.5)] <- cover_values

# Psi (habitat suitability = occupancy)
b_cl <- 3.7
b_co <- -1.5
psi <- exp(-3 + b_cl*clim + b_co*cover)/(1 + exp(b_cl*clim + b_co*cover))

# Realized abundance
suit <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
suit[] <- psi
plot(suit)

acCent <- rbinom(2500,1,suit[])

centers <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
centers[] <- acCent
plot(centers)
sum(acCent)

centXY <- xyFromCell(centers, cells(centers, 1)$lyr.1)
plot(suit)

for(j in 1:nrow(centXY)){
  print(j)
  pp <- data.frame(x = centXY[j,1], y = centXY[j,2])
  for(i in 1:200){
    pp[i+1,1] <- pp[i,1] + rnorm(1) + (0.2 * (centXY[j,1] - pp[i,1]))
    pp[i+1,2] <- pp[i,2] + rnorm(1) + (0.2 * (centXY[j,2] - pp[i,2]))
  }
  lines(pp)
}




