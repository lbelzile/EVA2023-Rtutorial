library(evd)
library(exdex)
library(extRemes)
library(extremogram)
library(lite)
library(POT)
data(frwind, package = "mev")
? frwind

# Time-series extremes ####

# We here consider the wind speed observations from the first station. We show the time series data for the first year of observation:
obs = frwind$S1
plot(obs[1:365],
     type = "h",
     xlab = "Day",
     ylab = "Observation")


# Next, we explore the extremogram (also called auto tail dependence function, as implemented in the package extRemes. We set a probability level for the threshold corresponding to correspond to approximately one exceedance per month.
prob = 1 - 1 / 30
maxlag = 7
atdf(obs, u = prob, lag.max = maxlag)
# Clearly, extremal dependance decays quite fast, and the chi-bar coefficient (here called "rhobar") points towards asymptotic independence for all positive time lags.
# Nevertheless, using the standard methods assuming asymptotic dependence can still be useful but may overestimate the strength of extremal correlation at very high quantile levels.

# Chi-bar analysis using tsdep.plot in POT package ####
# If necessary, we can reset graphics device with dev.off().
dev.off()
tsdep.plot(obs, u, lag.max = 7) # POT package
# The confidence bounds are relatively large but they tend to be bounded away from 1 (implying asymptotic independence).

# A similar analysis can be done with the extremogram package.
# We here estimate the upper-tail extremogram.
# Again, we reset graphics device if necessary:
# dev.off()
# We first produce an extremogram plot. Confidence bounds can also be obtained.
ext = extremogram1(obs,
                   quant = prob,
                   type = 1,
                   maxlag = maxlag)
bootconf1(
  obs,
  R = 100,
  l = 30,
  maxlag = maxlag,
  quant = prob,
  type = 1,
  par = 1,
  alpha = 0.05
)
points(1:(maxlag - 1), ext[-1], pch = 19, col = "blue")
# We note that there is still room for improvement of the output of these functions.

# It is also possible to calculate cross-extremogram between two measurement stations, and corresponding confidence bounds.
obs1 = frwind$S2
obs2 = frwind$S3
ext_cross = extremogram2(
  cbind(obs1, obs2),
  quant1 = prob,
  quant2 = prob,
  type = 1,
  maxlag = maxlag
)
bootconf2(
  cbind(obs1, obs2),
  R = 100,
  l = 30,
  maxlag = maxlag,
  quant1 = prob,
  quant2 = prob,
  type = 1,
  par = 1,
  alpha = 0.05
)
points(1:(maxlag - 1), ext_cross[-1], pch = 19, col = "blue")

# We find that the extremal dependence across the two series is very weak here.

# Estimation of the extremal index ####

# A large number of estimation approaches are available in the exdex package (type "?exdex" for the help package); see also the package vignette for more detailed information:
# https://cran.r-project.org/web/packages/exdex/vignettes/exdex-vignette.html

# We first explore estimation based on maxima for sliding blocks.
# Here, we use block sizes of approximately one month. We further obtain and visualize symmetric confidence intervals.
theta = spm(obs, 30)
summary(theta)
conf = confint(theta)
conf = confint(theta, interval_type = "lik")
plot(conf)

# The estimation of the extremal index can be sensitive to the block size. We check this for longer and shorter blocks.
theta2 = spm(obs, 14)
theta3 = spm(obs, 60)
cbind(theta2$theta_sl, theta$theta_sl, theta3$theta_sl)
# There is some moderate variability in the estimates.

# Moreover, there is a tool to help choosing an optimal block size. The execution of the following command can take more than one minute.
# Here, we explore block sizes between one and ten weeks. Note that certain small b-values may be too small for the sampling variance of the sliding blocks
# estimator to be estimated.
b_vals = 7 * 1:10
res = choose_b(obs, b_vals)
plot(res, ylim = c(0, 1))
# Here, even for a relatively small block size (b=7), the estimate is not significantly different from results for larger blocks, and it makes sense to use the estimate obtained for b=7.

# Another possibility is to estimate the extremal index based on threshold exceedances. We here explore two possible estimator, the K-gaps and D-gaps estimators.
prob = 1 - 1 / 30
u = quantile(obs, prob)
theta = kgaps(obs, u, k = 1)
summary(theta)
theta = dgaps(obs, u, D = 3)
summary(theta)

# Again, for the sake of choosing optimal tuning parameters, we have visual tool for the K-gaps estimator to explore different thresholds and K-values.
res = choose_uk(obs, u = quantile(obs, 1 - 1 / c(7, 14, 30, 60, 90, 180)), k = 1:7)
plot(res, y = "theta")
# The estimates look more stable only at relatively high thresholds.
# In practice, this lack of asymptotic stability makes threshold choice difficult.
# However, sometimes, there may be a specific application-relevant threshold that provides a natural choice.
# For example, the monthly return level with K=3 leads to estimates similar to the blocks estimates.

# Another implementation of an extremal-index estimator comes from the extRemes package.
extremalindex(obs, u, method = "runs", run.length = 2)
extremalindex(obs, u, method = "intervals")

# In general, it is recommended to check sensitivity of estimates to  different methods/thresholds.

# Peaks-over-threshold with declustering ####

# Here, we exemplify the approach using the implementation in the extRemes package. We also provide a comparison to estimation without any prior declustering.

obs_decluster = decluster(obs, u, method = "runs", run.length = 2)
fit = fevd(
  obs_decluster,
  threshold = u,
  type = "GP",
  time.units = "365/year"
)
fit
plot(fit)
fit2 = fevd(obs,
            threshold = u,
            type = "GP",
            time.units = "365/year")
fit2
plot(fit2)
# As suggested by theory, the estimated values are relatively close. The standard error estimates without declustering are lower but may be biased towards too low values.

# Likelihood-based estimation of extremal index and GPD without explicit declustering ####
# The lite package uses all exceedances for likelihood-based estimation with an adjustment to avoid biased uncertainty estimates. This avoids wasting information when we keep only cluster maxima.
# The K-gaps estimator is used for the extremal index. Remember that there is the function choose_uk(...) from above to help with the choice of K and the threshold. Finally, we can produce unbiased confidence intervals for parameters.
fit = flite(obs, u = u, k = 3)
summary(fit)
confint(fit)

# Joint modeling of margins and dependence of time-series extremes ####

# We show how a first-order Markov chain model can be fitted using a censored likelihood thanks to the POT package.
# Here, we use the bivariate logistic distribution to define the Markov chain transitions densities. Some graphical output is available to diagnostic the fitted model. We first illustrate the estimation for simulated data where the true parameters are known.

set.seed(2)
obs_simulated = simmc(
  n = 5000,
  alpha = .5,
  model = "log",
  margin = "gumbel"
)
fit_mc_simulated = fitmcgpd(obs_simulated,
                            threshold = quantile(obs_simulated, 0.95),
                            model = "log")
fit_mc_simulated
plot(fit_mc_simulated)
# Estimates are quite close to simulated parameters (scale = 1, shape = 0, alpha = 0.5)

# Next, we apply the method to our windspeed data. We also compare its estimation of marginal parameters to a marginal fit that does not take into account the dependence.
fit_mc = fitmcgpd(obs, u, model = "log")
fit_mc
plot(fitmc)
fit_gpd = fitgpd(obs, u, est = "mle")
fit_gpd
# Good news: The marginal parameters estimates are very similar.

# Exercice
# Perform a similar analysis for one or several of the other stations (windspeeds) or for temperature data, and interpret results and compare them across stations.
