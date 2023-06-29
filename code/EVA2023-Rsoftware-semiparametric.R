# Load data
data(frwind, package = "mev")
# Extract Lyon series and sort in decreasing order
lyon <- sort(frwind$S2, decreasing = TRUE)
# random block maximum of Wager 
remotes::install_github("lbelzile/rbm", quiet = TRUE)
# Fit Hill's estimator
hill_est <- rbm::hill(data_array = lyon,
                      idx = 1:400)
# Hill plot
plot(hill_est)
# Shifter staircase patterns are a result of rounding

# Compute a smoothed Hill plot (log scale, etc.)
evmix::hillplot(data = lyon[1:300],
                hill.type = "SmooHill", 
                r = 3L, 
                x.theta = TRUE)

# Use the random block max estimator
# Restrict to largest 300 exceedances since we wouldn't
# choose a lower value - also ensures that the computational
# burden is limited
rbm_est <- rbm::rbm.point_estimate(lyon[1:300])
# Plot the smooth curve with the optimal order statistic
# selected by the algorithm
plot <- rbm::rbm.plot(lyon[1:300])

# The 'tea' package includes many different methods
# for semiparametric threshold selection estimation
# Minimization of asymptotic mean squared error of Hill
# assuming second-order regular variation 
est_AMSE <- tea::dAMSE(lyon[1:300])

# Compute return level by extrapolating
# using Weissman estimator
qweissman <- function(n, p, k, thresh, shape){
  thresh * ((k + 1) / (1 - p)/ (n + 1))^shape
}
# Compute quantiles at different probability levels
quants <- qweissman(
  n = length(lyon), 
  p = seq(0.99, 0.999, length.out = 101), 
  thresh = est_AMSE$threshold, # threshold selected
  k = est_AMSE$k0, # number of points above threshold
  shape = 1/est_AMSE$tail.index) # shape xi

