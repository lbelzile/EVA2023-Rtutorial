# Load packages and data
library(evgam)
library(lubridate)
library(ggplot2)
data(frwind, package = "mev")
# Create time series of wind speed at Lyon
lyon <- with(frwind,
             xts::xts(x = S2, order.by = date))
# Extract annual maximum
ymax <- xts::apply.yearly(lyon, max)
ymax <- ymax[-length(ymax)] # Remove 2023 (incomplete year)

# Fit a GEV model using 'evgam'
# cubic spline for location based on year
# constant scale and shape
opt_gev <- evgam::evgam(
  # data frame with response and covariates
  data = data.frame(
    year = year(ymax),
    ymax = ymax),
  # formula is list of formulae for loc, scale, shape
  formula = list(ymax ~ s(year, bs = "cr"),
                 ~ 1, 
                 ~ 1),
  family = "gev") # extreme value distribution
# see ?evgam::evgam

## Summary with coefficient estimates
summary(opt_gev)
# plot of fitted spline for location parameter 
# of generalized extreme value distribution as 
# a function of day of year (scaled to unit interval)."
plot(opt_gev)

# Careful! - default is type='link' (log scale)
# Response here is mu, sigma, xi
# Based on normal approximation to posterior
sim_post <- simulate(object = opt_gev,
                     nsim = 1000L,
                     type = "response")
# Simulated parameter for each combo
# averaged across observed values of
# covariates in original data frame.
retlev <- with(sim_post,
  evgam::qev(p = 1-1/50,
             loc = location,
             scale = scale,
             shape = shape,
            family = "gev"))
# Compute point estimate
mean(retlev)

# More straightforwardly, simulate directly quantiles
mean(simulate(opt_gev, 
              nsim = 1000L, 
              type = "quantile", 
              probs = 1-1/50))
# Another model, this time restricting the number of knots
opt_gev_spl <- evgam::evgam(
  data = data.frame(
    syear = scale(year(ymax)),
    ymax = ymax),
  # k controls the knots
  formula = list(ymax ~ s(syear, k = 5, bs = "cr"),
                 ~ 1, 
                 ~ 1),
  family = "gev")
## Summary with coefficients
post_sim <- simulate(opt_gev_spl, nsim = 1000L, seed = 2023)
# Plot posterior of location function 
# note the varying uncertainty as we move 
# towards higher mean value
ggplot(
  data = data.frame(
    location = c(post_sim$location),
    year = factor(rep(1:length(years), 
                      length.out = prod(dim(post_sim$location))))),
  mapping = aes(x = location,
                color = year,
                group = year)) +
  geom_density() +
  theme_minimal() +
  viridis::scale_color_viridis(discrete = TRUE) +
  theme(legend.position = "none")



## Fit a nonstationary threshold using 
# asymmetric Laplace distribution to get
# a working likelihood with a modified check function
qreg <- evgam::evgam(
  # Cubic cyclic splines
  formula = list(wind ~ s(tday, bs = "cc"), 
                 ~ s(tday, bs = "cc")), 
  data = data.frame(wind = as.numeric(lyon),
                    tday = yday(lyon)),
  family = "ald", 
  ald.args = list(tau = 0.95)) # quantile level

plot(qreg)

# Now extract the fitted quantiles (correspond to loc)
# to use as threshold
u <- fitted(qreg)[,'location']
# Create a data frame containing
# exceedance, thresholds and covariates
df <- data.frame(exc = as.numeric(lyon - u),
                 yday = yday(lyon),
                 u = u)
# Extract exceedances from the latter
df <- df[lyon >u, ]

# Fit a generalized Pareto distribution 
opt_gpd <- evgam::evgam(
  formula = list(exc ~ 1, ~1), # scale shape formulaes - constant
  data = df,
  family = "gpd")
# Simulate parameters of generalized Pareto
post_sim <- simulate(opt_gpd, 
                     type = "response")

# Compute average number of exceedances per year
npy <- diff(range(date(lyon))) / nrow(df)
# Compute return levels - must be high enough
# for empirical probability component to be zero
retlev_gpd <- qev(
  p = 1-1/50, # quantile level for return level
  m = as.numeric(npy), # number of observations per year (average)
  loc = df$u, # threshold
  scale = post_sim$scale,
  shape = post_sim$shape,
  family = "gpd")
# Compute average return level
mean(retlev_gpd)
