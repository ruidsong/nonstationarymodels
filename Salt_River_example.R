#### Data preparation ####
# install.packages('SSN')
# install.packages('geoR')
library(SSN)
library(geoR)

## Import the data from the .ssn directory and create a SpatialStreamNetwork object
SaltWQ <- importSSN("/mnt/nfs/clasnetappvm/fs2/ruidsong/Github/TroutDensity_BlockKrige.ssn")

## Create distance matrix among network observation sites
createDistMat(SaltWQ, o.write = TRUE)

## Remove corresponding rows and columns of four outliers
half_dist_mat <- as.matrix(getStreamDistMat(SaltWQ, Name = 'obs')[[1]])
half_dist_mat_4less <- half_dist_mat[-c(59,64,87,88), -c(59,64,87,88)]

file_handle = file(file.path(SaltWQ@path, "distance",
                              "obs", "dist.net1.RData"), open = "wb")
serialize(half_dist_mat_4less, file_handle, ascii = FALSE)
close(file_handle)

# A quick check
dim(as.matrix(getStreamDistMat(SaltWQ, Name = 'obs')[[1]])) == c(104, 104)

## Remove corresponding rows of the four outliers from all data sources 
SaltWQ@obspoints@SSNPoints[[1]]@network.point.coords <- 
  SaltWQ@obspoints@SSNPoints[[1]]@network.point.coords[-c(59,64,87,88), ]
SaltWQ@obspoints@SSNPoints[[1]]@point.coords <- 
  SaltWQ@obspoints@SSNPoints[[1]]@point.coords[-c(59,64,87,88), ]
SaltWQ@obspoints@SSNPoints[[1]]@point.data <- 
  SaltWQ@obspoints@SSNPoints[[1]]@point.data[-c(59,64,87,88), ]


#### Plots in explanatory analysis ####
## Plot Salt River network and locations of 104 trout density observations
plot(SaltWQ, lwdLineCol = "afvArea", lwdLineEx = 5, lineCol = "blue",
     pch = 19, xlab = "x-coordinate (m)",
     ylab = "y-coordinate (m)", asp = 1)

## Plot values of 104 trout density observations
brks <- plot(SaltWQ, "trout_100m", lwdLineCol = "afvArea",
             lwdLineEx = 5, lineCol = "black", xlab = "x-coordinate" ,
             ylab = "y-coordinate", asp=1)

## Scatterplot of resid vs. upDist from OLS 
formula <- trout_100m ~ SLOPE + CANOPY + S1_93_11

OLS_fit <- glmssn(formula, SaltWQ, CorModels = NULL, EstMeth = "ML")
resid_ssn <- residuals(OLS_fit)
resid <- getSSNdata.frame(resid_ssn)[,'_resid_']

upstream_dist <- getSSNdata.frame(resid_ssn)[, 'upDist']

par(mfrow = c(2, 2))
plot(upstream_dist, resid, xlab = 'Upstream distance', ylab = 'Residual')

## Nearest-neighbor scatterplot for 1st, 2nd, 3rd upDist portions
L3rd_ind <- which(upstream_dist <= quantile(upstream_dist, 1/3))
M3rd_ind <- which((upstream_dist > quantile(upstream_dist, 1/3)) & 
                    (upstream_dist <= quantile(upstream_dist, 2/3)))
H3rd_ind <- which(upstream_dist > quantile(upstream_dist, 2/3))

# A quick check
sort(c(L3rd_ind, M3rd_ind, H3rd_ind)) == 1:104

mat_dist <- half_dist_mat_4less + t(half_dist_mat_4less)
isSymmetric(mat_dist)

NN_ind <- rep(NA, 104)
for (i in 1:104) {
  NN_ind[i] <- which(mat_dist[i,] == sort(mat_dist[i,])[2])
}

NN_resid <- resid[NN_ind] 

plot(resid[L3rd_ind], NN_resid[L3rd_ind], xlim = c(-45, 60), ylim = c(-45, 60), 
     xlab = 'Residual at originating site', ylab = 'Residual at nearest neighbor')
plot(resid[M3rd_ind], NN_resid[M3rd_ind], xlim = c(-45, 60), ylim = c(-45, 60), 
     xlab = 'Residual at originating site', ylab = 'Residual at nearest neighbor')
plot(resid[H3rd_ind], NN_resid[H3rd_ind], xlim = c(-45, 60), ylim = c(-45, 60), 
     xlab = 'Residual at originating site', ylab = 'Residual at nearest neighbor')

par(mfrow = c(1, 1))


#### Nonstationary model fitting ####
## Load functions
source('/mnt/nfs/clasnetappvm/fs2/ruidsong/Github/nonstationary_function.R')

## Define multiple starting values for nonstationary parameters 
## range1 = psi in elastic models and range1=alpha1 in SVMA models
gamma1_seq <- seq(-2, 0, length.out = 5)
range1_seq <- seq(0, 1, length.out = 5)
phi1_seq <- seq(-2, 0, length.out = 5)
kappa_seq <- seq(0.5, 1, length.out = 5)

## Define attributes
formula <- trout_100m ~ S1_93_11
ssn.object <- SaltWQ
dichvar <- "upDist"
family <- "Gaussian"
use.nugget <- TRUE
use.anisotropy <- FALSE
addfunccol <- "afvArea"
trialscol <- NULL
useTailDownWeight <- FALSE
trans.power <- NULL
trans.shift <- 0
updist_scaler <- min(SaltWQ@obspoints@SSNPoints[[1]]@point.data$upDist)

## IMPORTANT: copy and paste the downloaded SSN folder and name it "TroutDensity_BlockKrige.assist.ssn"
## Load the assistant Salt River data
ssn.assist <- importSSN("/mnt/nfs/clasnetappvm/fs2/ruidsong/Github/TroutDensity_BlockKrige.assist.ssn")

## We use elastic tail-down exponential model as an example, readers may replace 
## CorModels = 'EL.Exponential.taildown' with other models, possible models are
## EL.Exponential.taildown, EL.Exponential.tailup
## EL.LinearSill.taildown, EL.LinearSill.tailup
## EL.Spherical.taildown, EL.Spherical.tailup
## SVMA.Exponential.taildown, SVMA.Exponential.tailup
## SVMA.LinearSill.taildown, SVMA.LinearSill.tailup
## SVMA.Spherical.taildown, SVMA.Spherical.tailup
## Exponential.isotropic, Matern.isotropic

EL_EXP_td <- c()

## Step 1: automated process
for (i in 1:length(gamma1_seq)) {
  for (j in 1:length(range1_seq)) {
    for (k in 1:length(phi1_seq)) {
      fit <- try(glmssn_nonstationary(formula, ssn.object, dichvar, CorModels = 'EL.Exponential.taildown', 
                                      EstMeth = "ML", init = c(gamma1_seq[i], range1_seq[j], phi1_seq[k])))
      if (class(fit) == "glmssn") {
        out <- InfoCritCompare(list(fit))
        out <- cbind(out, parsil1_seq[i], range1_seq[j], nugget1_seq[k], t(as.vector(fit$estimates$theta)), 
                     t(fit$estimates$betahat), t(c(fit$estimates$covb[1,1], fit$estimates$covb[2,2])))
        EL_EXP_td <- rbind(EL_EXP_td, out)
      }
      print(c(i,j,k))
    }
  }
}
EL_EXP_td

## Step 2: find the smallest AIC in EL_EXP_td and record the corresponding parameter estimates,
## evaluate the log-likelihood on a grid overlaid on this vector estimate
grid_search(formula, ssn.object, dichvar, CorModels = 'EL.Exponential.taildown', EstMeth = "ML",
            lower.bound = c(6.68, 0.0042, 6.59, -0.093, 0.015, -0.425), 
            upper.bound = c(6.70, 0.0044, 6.61, -0.091, 0.017, -0.423), grid_depth = 10)


#### Isotropic models ####
## for the shape parameter kappa, we need to constrain it in (0,1) for power exponential isotropic model,
## and in (0,1/2) for matern isotropic model
glmssn_nonstationary(formula, ssn.object, dichvar, CorModels = 'Matern.isotropic', EstMeth = "ML",
                     init = c(-1,0.5,-1, 0.5), upper.bound = c(Inf, Inf, Inf, Inf, Inf, Inf, 0.5),
                     lower.bound = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 0))
