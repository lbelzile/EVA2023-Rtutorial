# Load packages and data
library(mev)
library(xts)
library(lubridate)
data(frwinds, package = "mev")
# Extract time series of mean wind speed for Lyon
lyon <- with(frwind, 
             xts(x = S2, order.by = date))
# Create series of yearly maximum
ymax <- apply.yearly(lyon, max)

# Fit a GEV distribution to annual max
# Using maximum likelihood
# 'show' = TRUE prints the result to the console
opt_gev <- mev::fit.gev(xdat = ymax, show = TRUE)
# Extract maximum likelihood estimates
mle <- coef(opt_gev)

# check which S3 methods are available
methods(class = "mev_gev")
# Quantile-quantile and probability-probability plots
par(mfrow = c(1,2))
plot(opt_gev)

# Parameter invariance - compute MLE, 
# MLE of expectation of maximum of 50 blocks
gev.mle(xdat = ymax, args = "Nmean", N = 50)


# Check that gradient of log likelihood is zero at MLE
mev::gev.score(par = mle, dat = ymax) 
# Compute observed information matrix
jmat <- mev::gev.infomat(par = mle, dat = ymax)
# Compute standard errors
sqrt(diag(solve(jmat)))
# Compare with opt$std.err
# Not valid if xi < -0.5

# Confidence intervals for functional of interest
# via profile likelihood
prof <- mev::gev.pll(
   param = "Nmean", # functional of interest
   dat = ymax, 
   N = 50) # number of 'years'
# Obtain confidence intervals based on chi-square (1) 
# asymptotic distribution
(confint(prof))

prof2 <- mev::gev.pll(
   param = "Nquant",
   dat = ymax, 
   q = 0.5, # probability level (here median)
   N = 50) # number of 'years'
confint(prof2, print = TRUE)
# Peaks-over-threshold approach
# First plot the time series
plot(lyon)
# Plot as a function of the day of the year
plot(x = yday(frwind$date),
     y = frwind$S2,
     xlab = "day of year",
     ylab = "mean wind speed at Lyon (km/h)",
     pch = 20)

# Extract data for some months to deal with seasonality
windlyon <- with(frwind, S2[month(date) <= 4 | month(date) >= 9])
# Pick a quantile level - there won't be 100 obs due to ties
qulev <- 1-100/nrow(windlyon)
# Extract (interpolated) empirical quantile for threshold
u <- quantile(windlyon, 1-100/length(windlyon))

# Fit a generalized Pareto distribution to exceedances above u
opt_gp <- mev::fit.gpd(
  xdat = windlyon, 
  threshold = u, 
  show = TRUE)

# In practice, we need a threshold such that we have
# approximate threshold stability (needed for extrapolation)
useq <- quantile(windlyon, seq(0.9, 0.99, by = 0.01))
tstab <- tstab.gpd(windlyon,
          method = "profile",
          thresh = useq)
# Read from right: find lowest threshold such that later estimates
# are included inside confidence intervals
par(mfrow = c(1,1))
plot(tstab, which = 2, sub = "", main = "")
# Wadsworth white noise challenge
W.diag(xdat = windlyon, u = useq)
# Northrop-Coleman score test for penultimate approx
NC.diag(xdat = windlyon, u = useq)


# No clear choice for threshold
# Pick a high quantile level
u <- quantile(windlyon, 0.99)
# Fit inhomogeneous Poisson likelihood
opt.pp <- fit.pp(
  xdat = windlyon, 
  threshold = u, 
  show = TRUE,
  # 'np' = number of periods of data
  np = diff(range(lubridate::year(frwind$date))))

# Compute 50-year return level
gev.Nyr(par = coef(opt.pp),
        nobs = opt.pp$nat, # number above threshold 
        type = "retlev", 
        p = 1/50)
# Standard errors are based on delta method