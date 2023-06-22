
library(mev)
library(xts)
library(lubridate)
data(frwinds, package = "mev")
lyon <- with(frwind,
             xts(x = S2, order.by = date))
# Create series of yearly maximum
ymax <- apply.yearly(lyon, max)


################################################
###             Block maximum          #########
################################################

# Fit a GEV distribution via maximum likelihood
opt_gev <- mev::fit.gev(xdat = ymax, show = TRUE)
# Extract MLE
mle <- coef(opt_gev)
# Check convergence manually
isTRUE(all.equal(rep(0,3),
                 mev::gev.score(par = mle, dat = ymax),
                 check.attributes = FALSE,
                 tolerance = 1e-5))
# Score should be zero at MLE

# Compute observed information matrix
jmat <- mev::gev.infomat(par = mle, dat = ymax)
# Compute standard errors
se_mle <- sqrt(diag(solve(jmat)))
# Compare with
isTRUE(all.equal(se_mle, opt_gev$std.err))

# To see other methods, query
methods(class = "mev_gev")
# For example, to fit a Gumbel model and compare the fit
opt_gumb <- mev::fit.gev(xdat = ymax,
                         fpar = list(shape = 0))
# Likelihood ratio test for nested models
anova(opt_gev, opt_gumb)

# PP and QQ plots
par(mfrow = c(1,2))
plot(opt_gev)
dev.off()
# Invariance: compute a risk summary
# MLE of expectation of maximum of 50 blocks
gev.mle(xdat = ymax, args = "Nmean", N = 50)
# Compute profile log-likelihood
prof <- mev::gev.pll(param = "Nmean", dat = ymax, N = 50)
# Extract confidence intervals
(confint(prof))

################################################
###        Threshold exceedances       #########
################################################


with(frwind, plot(
  x = yday(date),
  y = S2,
  xlab = "day of year",
  ylab = "wind at Lyon St-Exupéry (in km/h)"))
# Data are not stationary, so we pick a period
# where this seems more or less the case
windlyon <- with(frwind,
                 S2[month(date) <= 4 | month(date) >= 9])
# Keep only 100 exceedances - more about threshold
# selection later

# Should always have at least 20 observations
#  for fit of GP to be reliable
#  ML also exhibits strong small sample bias
qulev <- 1-100/nrow(windlyon)
u <- quantile(windlyon, 1-100/length(windlyon))

# Fit generalized Pareto distribution
opt_gp <- mev::fit.gpd(
  xdat = windlyon,
  threshold = u,
  show = TRUE)

# Threshold stability plots
useq <- quantile(windlyon, seq(0.95, 0.99, by = 0.005))
tstab.gpd(windlyon,
          method = "profile",
          thresh = useq)
# Go backward and try to see if higher values
# fall within confidence intervals
W.diag(xdat = windlyon, u = useq)
# Includes threshold stability plots
# (but simultaneous intervals)
# White noise sequence from increments:
# consider normal centered at zero as H0

# Northrop-Coleman diagnostics
# want p-value to go above 5% or so
NC.diag(xdat = windlyon, u = useq)


## So far, contradictory results...
# Pick a threshold
u <- useq[7] #98%

# How many periods (year)?
nyear <- diff(range(year(frwind$date)))
# Fit the inhomogeneous Poisson point process
opt_pp <- fit.pp(
  xdat = windlyon,
  threshold = u,
  np = nyear)
# Can specify either number of periods of data (np)
# or average number of exceedances per period

# Point estimator and variance (delta-method based)
Nyr_ipp <- gev.Nyr(
  par = coef(opt_pp),
  N = 50,
  type = "mean",
  nobs = nyear)
# Compare point estimates between IPP and GEV
c(GEV = prof$mle['Nmean'], IPP = Nyr_ipp$est)


# Computing loss functions
# The EVA 2023 Data challenge featured
# a loss function for return levels, of the form


# We can compute the analog of a profile
# For each value of q, find the best parameter combo
loss_fn <- function(psi, data, Nyr = 50){
  data <- as.numeric(data)
  # Define loss function
  loss <- function(qhat, q){
    mean(ifelse(0.99*q > qhat,
                0.99*(0.99*q-qhat),
                ifelse(1.01*q < qhat,
                       0.1*(qhat-1.01*q),
                       0)))
  }
  retlev_fn <- function(par, Nyr){
    loc <- par[1]
    scale <- par[2]
    shape <- par[3]
    qgev(p = 0.368,
         loc = loc + scale*(Nyr^shape-1)/shape,
         scale = scale*Nyr^shape,
         shape = shape)
  }
  # Create optimization wrapper
  opt_fun <- function(par, qi, data = data, Nyr = Nyr){
  log(loss(qhat = retlev_fn(par, Nyr = Nyr),
       q = qi)) -
    mean(dgev(x = data,
         loc = par[1],
         scale = par[2],
         shape = par[3],
         log = TRUE))
  }
  mle <- coef(fit.gev(xdat = data))
  loss_pt <- vector(mode = "numeric", length = length(psi))
  param_mat <- matrix(nrow = length(psi), ncol = 3L)
  init <- mle
  for(i in seq_along(psi)){
    opt <- Rsolnp::solnp(
       fun = opt_fun,
       data = data,
       Nyr = Nyr,
      # Use previous optimum as starting value
       pars = init,
       qi = psi[i],
       ineqfun = function(par, qi, data, Nyr){
         c(par[2],
           par[2] + par[3]*(range(data)-par[1]),
           par[3])
       },
       # Put sensible values for upper bounds
       ineqLB = c(rep(0, 3), -1),
       ineqUB = c(rep(1e4, 3), 2),
       control = list(trace = 0))
    param_mat[i,] <- opt$pars
    init <- opt$pars
    loss_pt[i] <- tail(opt$values, 1)
  }
  return(list(psi = psi, loss = loss_pt, pars = param_mat))
}

# Compute loss function pointwise for posterior
loss <- lossfun(qhat = (qu <- seq(180, 250, by = 0.01)),
                q = retlev_gp)
qu[which.min(loss)]
