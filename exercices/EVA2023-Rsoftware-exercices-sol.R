
library(mev)
nrep <- 1000L
ests <- matrix(0, nrow = nrep, ncol = 2)
for(i in 1:nrep){
  set.seed(i)
  samp <- rgp(n = 15, scale = 1, shape = -0.1)
  ests[i,] <- coef(fit.gpd(samp))
}
# Strong small-sample bias, negative for shape
summary(ests)
# In small samples, empirical coefficient of variation
# may be less than 1
mean(ests[,2] == -1)

data(frwind, package = 'mev')
# Create time series
montelimar_ts <- with(frwind,
  xts::xts(S4, order.by = date))
# Extract 3-largest observations per year
tlarg <- xts::apply.yearly(
  montelimar_ts, 
  function(x){
    # need to cast to numeric, as 'sort' does not work as intended
    head(sort(as.numeric(x), decreasing = TRUE), n = 3)
  })
# Fit models with GEV and r-largest
mle_gev <- fit.gev(xdat = tlarg[,1])
mle_rlarg <- fit.rlarg(xdat = tlarg)
# Maximum likelihood estimates
(mle_gev)
(mle_rlarg)
# Variance reduction factor here is 
(det(vcov(mle_rlarg))/det(vcov(mle_gev)))^(1/3)

null_gev <- fit.gev(xdat = tlarg[,1], fpar = list(shape = 0))

null_score <- gev.score(par = null_gev$param, dat = as.numeric(tlarg[,1]))
null_info <- gev.infomat(par = null_gev$param, dat = as.numeric(tlarg[,1]), method = "exp")
# Compute score statistic
score_test <- t(null_score) %*% solve(null_info) %*% null_score
pval <- pchisq(score_test, df = 1, lower.tail = FALSE)
as.numeric(pval)
# Compare with likelihood ratio test
anova(object2 = null_gev, mle_gev)
# Profile and confidence intervals
prof <- gev.pll(psi = seq(27, 34, length.out = 51),
                param = "quant", 
                p = 1-1/50, 
                dat = tlarg[,1])
confint(prof, level = 0.5)
prof <- gev.pll(psi = seq(27, 34, length.out = 51),
                param = "Nquant", 
                q = 0.368,
                N = 50,
                dat = tlarg[,1])
confint(prof, level = 0.5)


## Bayesian and nonstationary models

library(evgam)
library(lubridate)
library(dplyr)
data <- frwind|>
  filter(month(date) <= 4 | month(date) >= 9) |>
  select(S1:S4)
thexc <- apply(data, 2, function(x){
  ui <- quantile(x, 0.98, na.rm = TRUE)
  x[x>ui] - ui
})
nexc <- sapply(thexc, length)
exc_data <- data.frame(
  exc = unlist(thexc),
  group = factor(
    rep(1:4, times = sapply(thexc, length))
  ))
gpd_mult <- evgam(formula = list(exc ~ group, ~1),
             family = "gpd",
             data = exc_data)
gpd_mult$logLik
# Fit the model separately to each series 
# (different shape parameters
gpd_sep <- sapply(thexc, function(x){
  logLik(fit.gpd(x))
})
# Likelihood ratio test
pchisq(q = 2*(sum(gpd_sep) - gpd_mult$logLik),
       df = 3, 
       lower.tail = FALSE)  
# In parametrization, the first site gets 
# the intercept for log-scale
# the others are different relative to that
# Since all other sites have lower scale
# the first site is the one with the largest return levels
coef(gpd_mult)
# Fit the Bayesian model
library(revdbayes)
# Extract data for montelimar
S4_post_samp <- 
  rpost_rcpp(
    n = 1e4L, 
    thresh = quantile(data$S4, 0.98),
    model = "bingp",
    data = data$S4,
    prior = set_prior(prior = "mdi", model = "gp"),
    bin_prior = set_bin_prior(prior = "beta"))
plot(S4_post_samp$sim_vals)

S4_evgam <- evgam(
  formula = list(exc ~ 1, ~1),
  data = exc_data |> filter(group == 4),
  family = "gpd")
S4_evgam_post <- simulate(
  object = S4_evgam,
  nsim = 1000,
  type = "response")
points(x = S4_evgam_post$scale, 
       y = S4_evgam_post$shape, 
       col = scales::alpha(colour = 2, alpha = 0.1))
# Doesn't quite capture the shape and the range of values of xi...

gpd_retlev <- function(zetau, sigma, xi, N = 50, npy, thresh = 0){
  thresh + sigma/xi*((npy*N*zetau)^xi-1)
}
# Number of exceedances per year on average
npy <- nexc[4]/(diff(range(year(frwind$date)))-1)
th <- quantile(data$S4, 0.98)
retlev_rou <- with(S4_post_samp,
     gpd_retlev(zetau = bin_sim_vals,
                sigma = sim_vals[,1],
                xi = sim_vals[,2],
                N = 50,
                npy = npy,
                thresh = th))
quantile(retlev_rou, probs = c(0.05, 0.95))
# Do the same, but fix the binomial probability
retlev_evgam <- with(S4_evgam_post,
     gpd_retlev(zetau = 0.02,
                sigma = scale,
                xi = shape,
                N = 50,
                npy = npy,
                thresh = th))
quantile(retlev_evgam, probs = c(0.05, 0.95))
