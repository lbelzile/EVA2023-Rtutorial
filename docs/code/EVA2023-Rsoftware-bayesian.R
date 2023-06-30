# Load packages and data
library(revdbayes)
library(ggplot2)
data("frwind", package = "mev")
# Extract time series for Lyon
lyon <- with(frwind,
             xts::xts(x = S2, order.by = date))
# Create series of yearly maximum
ymax <- as.numeric(xts::apply.yearly(lyon, max))
# Set priors: maximal data information, Martins-Stedinger
# and a trivariate normal prior for mu, log(sigma), xi
prior1 <- set_prior(prior = "mdi", model = "gev")
prior2 <- set_prior(prior = "beta", model = "gev")
prior3 <- set_prior(prior = "norm", 
                model = "gev", 
                mean = c(mean(ymax), log(sd(ymax)), 0), 
                cov = diag(c(1000, 1000, 1)))
# Use `revdbayes` to obtain samples from the posterior 
# Using ratio-of-uniform method
# GEV distribution fit to annual max
post_1 <- revdbayes::rpost_rcpp(
  n = 1e4L, # number of simulations
  model = "gev", # extreme value family
  data = ymax,
  prior = prior1, #prior function
  nrep = 100) # number of replications from the posterior predictive
# The latter are used for validation purposes later

# Extract posterior samples of theta
post_samp1 <- post_1$sim_vals
# Same, but with a different prior
post_samp2 <- revdbayes::rpost_rcpp(
  n = 1e4L, 
  model = "gev",
  data = ymax,
  prior = prior2)$sim_vals
post_samp3 <- revdbayes::rpost_rcpp(
  n = 1e4L, 
  model = "gev", # extreme value distribution
  data = ymax, # vector of yearly maximum
  prior = prior3)$sim_vals
# Draw marginal posterior for shape
# This amounts to only selecting the column for the shape
ggplot(data = data.frame(
  shape = c(post_samp1[,'xi'],
            post_samp2[,'xi'],
            post_samp3[,'xi']),
  prior = rep(c("mdi", "beta", "normal"), 
              each = 1e4L)),
  mapping = aes(x = shape,
                col = prior, 
                group = prior)) +
  geom_density() +
  scale_color_viridis_d() + 
  theme_minimal() +
  theme(legend.position = "bottom")

# Create a function to compute the N-year maximum median
# Apply is equivalent to a for-loop, so for each posterior draw
# Compute function with those arguments
post_gev_Nmed <- apply(post_1$sim_vals, 1, function(x){
# *gev from 'revdbayes' has an argument 'm' that changes
# the arguments so that the parameters correspond to m-max
  revdbayes::qgev(p = 0.5, loc = x[1], scale = x[2], 
                  shape = x[3], m = 50)
})
# Posterior mean of median function
mean(post_gev_Nmed)
# To get a 95% credible interval, simply compute quantiles
quantile(post_gev_Nmed, c(0.025, 0.975))


# Posterior predictive samples by simulating
# new GEV with the parameters drawn
postpred <- revdbayes::rgev(
  n = nrow(post_1$sim_vals), 
  loc = post_1$sim_vals[,1],
  scale = post_1$sim_vals[,2],
  shape = post_1$sim_vals[,3],
  m = 50)


# Density of posterior predictive (black) and posterior median (grey) of 50 year maximum.
ggplot(data = data.frame(
  postpred = postpred,
  postmed = post_gev_Nmed)) +
  geom_density(mapping = aes(x = postmed), color = "grey") +
  geom_density(mapping = aes(x = postpred), color = "black", linewidth = 2) +
  labs(x = "average wind speed (in km/h)") + 
  theme_classic() +
  theme(legend.position = "none") 

# Posterior predictive checks
# Can use any summary statistic
# Need 'nrep' argument larger than 1
# the more, the better the inference.
pp_check(post, stat = median)

# Different package that allows for regression  model
library(texmex)
post_reg <- texmex::evm(
  y = ymax, 
  data = data.frame(
    syear = scale(1976:2023), 
    ymax = ymax),
  family = texmex::gev,
  method = "simulate",
  burn = 1000L,
  chains = 4L,
  iter  = 1.1e4,
  proposal.dist = "gaussian",
  mu = ~ syear,
  verbose = FALSE)


# Summary of MCMC
summary(post_reg)

# Default 'plot' method for texmex objects:
# density, traceplots and correlograms plots
library(gridExtra)
ggplot(post_reg)


# There are methods for the output, but the 'coda' package
# gives more choice
chains <- coda::as.mcmc.list(
  lapply(post_reg$chains, coda::as.mcmc))
# effective sample size is sufficient here
coda::effectiveSize(chains)


## Default 'plot' method for texmex objects:
# density, traceplots and correlograms plots
plot(chains, density = FALSE)


# Plot correlogram for Markov chains of different parameters
coda::acfplot(chains)


# Obtain summary statistics (mean, quartiles, range) for chains
summary(chains[[1]])


# Loss function for wind speed data, based on fitting a 
# generalized extreme value distribution to annual maxima.
gev_retlev <- function(par, N, p = 0.368){
  # Map parameters via GEV max-stability
  mu <- par[1] + par[2]*(N^par[3]-1)/par[3]
  sigma <- par[2]*N^par[3]; 
  xi <- par[3]
  # quantile of N-block maximum
  mev::qgev(p = p, loc = mu, scale = sigma, shape = xi)
}

# Loss function
loss <- function(qhat, q){
    mean(ifelse(0.99*q > qhat,
           0.99*(0.99*q-qhat),
           ifelse(1.01*q < qhat,
                  0.1*(qhat-1.01*q),
                  0)))
}
# Compute the posterior of the return levels
retlev_post <- apply(post_samp1, 1, gev_retlev, N = 50)

# Create a grid of values over which to estimate the risk
retlev_psi <- seq(
  from = quantile(retlev_post, 0.2),
  to = quantile(retlev_post, 0.99), 
  length.out = 101)
  
# Create a container to store results
risk <- numeric(length = length(retlev_psi))
for(i in seq_along(risk)){
  # Compute integral (Monte Carlo average over draws)
 risk[i] <- loss(q = retlev_post, qhat = retlev_psi[i])
}
# Plot loss function
ggplot(data = data.frame(
  loss = risk, 
  retlev = retlev_psi), 
  mapping = aes(x = retlev, y = loss)) +
  geom_line() +
  geom_vline(xintercept = mean(retlev_post)) +
  labs(x = "return level") +
  theme_minimal()

