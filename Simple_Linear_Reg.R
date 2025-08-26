####Basic lm in JAGS

set.seed(42) # Set a random seed for reproducibility of the simulation

samplesize <- 30 # Number of data points
b_length <- sort(rnorm(samplesize)) # Body length (explanatory variable)

int_true <- 30 # True intercept
slope_true <- 10 # True slope
mu <- int_true + slope_true * b_length # True means of normal distributions
sigma <- 5 # True standard deviation of normal distributions

b_mass <- rnorm(samplesize, mean = mu, sd = sigma) # Body mass (response variable)

snakes1 <- data.frame(b_length = b_length, b_mass = b_mass)
head(snakes1) ## these will have some negatives despite not being realistic 

library(R2jags)



#convert to different format (vector for JAGS)
jagsdata_s1 <- with(snakes1, list(b_mass = b_mass, b_length = b_length, N = length(b_mass)))


# similar to base r but jags is dif (~ is the signifier for random variable distribution)
# this first part specifies the likelihood of data given the lm
# values randomly drawn from norm
#
lm1_jags <- function(){
  for (i in 1:N) {
    b_mass[i]~dnorm(mu[i],tau) #tau being precision
    mu[i] <- alpha+beta*b_length[i]
  }
  #define the priors 
  # Priors:
  alpha ~ dnorm(0, 0.01) # intercept
  beta ~ dnorm(0, 0.01) # slope
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS  
}


# MCMC sampler 
init_values <-  function(){
  list(alpha = rnorm(1),beta=rnorm(1),sigma=runif(1))
}

params <- c("alpha", "beta", "sigma")

#run the MCMC chain with 12000 iterations and has burn in of 2000, every 10th is saved or "thinned"
fit_lm1 <- jags(data = jagsdata_s1, inits = init_values, parameters.to.save = params, model.file = lm1_jags,
                n.chains = 3, n.iter = 12000, n.burnin = 2000, n.thin = 10, DIC = F)


## For each parameter, n.eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
fit_lm1


traceplot(fit_lm1, mfrow = c(2, 2), ask = F) 

# visualization of posterior
lm1_mcmc <- as.mcmc(fit_lm1)
plot(lm1_mcmc)



# plot data together with model prediction 
nvalues <- 100 #100 values

#new seq of predictor values
b_length_new <- seq(min(snakes1$b_length), max(snakes1$b_length), length.out=nvalues)


# 3 chains into one to use the MCMC samples for prediction
lm1_mcmc_combi <- as.mcmc(rbind(lm1_mcmc[[1]],lm1_mcmc[[2]],lm1_mcmc[[3]]))

# calc the expected value of body mass for each new body length
pred_mean_mean <- mean(lm1_mcmc_combi[, "alpha"]) + b_length_new * mean(lm1_mcmc_combi[, "beta"])


# uncertainty about the true parameter values
pred_mean_dist <- matrix(NA, nrow = nrow(lm1_mcmc_combi), ncol = nvalues)

for (i in 1:nrow(pred_mean_dist)){
  pred_mean_dist[i,] <- lm1_mcmc_combi[i,"alpha"] + b_length_new * lm1_mcmc_combi[i,"beta"]
}
credible_lower <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
credible_upper <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)


# incorporates both true parameter uncertainty, and stochasticity of relationship
lm1_mcmc_combi_rep <- do.call(rbind, rep(list(lm1_mcmc_combi), 50)) # replication

# Draw random values for all parameter combinations (rows) and body length values (columns):
pred_data_dist <- matrix(NA, nrow = nrow(lm1_mcmc_combi_rep), ncol = nvalues)
for (i in 1:nrow(pred_data_dist)){
  pred_data_dist[i,] <- lm1_mcmc_combi_rep[i,"alpha"] + b_length_new * lm1_mcmc_combi_rep[i,"beta"] +
    rnorm(nvalues, mean = 0, sd = lm1_mcmc_combi_rep[i, "sigma"])
}

# Calculate quantiles:
uncertain_lower <- apply(pred_data_dist, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper <- apply(pred_data_dist, MARGIN = 2, quantile, prob = 0.975)



plot(b_mass ~ b_length, data = snakes1)
lines(b_length_new, pred_mean_mean)
lines(b_length_new, credible_lower, lty = 2)
lines(b_length_new, credible_upper, lty = 2)
lines(b_length_new, uncertain_lower, lty = 2, col = "red")
lines(b_length_new, uncertain_upper, lty = 2, col = "red")





#### With categorical predictor 
set.seed(42)

samplesize <- 50 # Larger sample size because we're fitting a more complex model
b_length <- sort(rnorm(samplesize)) # Body length
sex <- sample(c(0, 1), size = samplesize, replace = T) # Sex (0: female, 1: male), all cat jags data must be numeric

int_true_f <- 30 # Intercept of females
int_true_m_diff <- 5 # Difference between intercepts of males and females
slope_true_f <- 10 # Slope of females
slope_true_m_diff <- -3 # Difference between slopes of males and females

mu <- int_true_f + sex * int_true_m_diff + (slope_true_f + sex * slope_true_m_diff) * b_length # True means
sigma <- 5 # True standard deviation of normal distributions

b_mass <- rnorm(samplesize, mean = mu, sd = sigma) # Body mass (response variable)

# Combine into a data frame:
snakes2 <- data.frame(b_length = b_length, b_mass = b_mass, sex = sex)
head(snakes2)

# reformat data
jagsdata_s2 <- with(snakes2, list(b_mass = b_mass, b_length = b_length, sex = sex, N = length(b_mass)))


#rewrite func to account for sex and body length relationship 
lm2_jags <- function(){
  # Likelihood:
  for (i in 1:N){
    b_mass[i] ~ dnorm(mu[i], tau) # tau is precision (1 / variance)
    mu[i] <- alpha[1] + sex[i] * alpha[2] + (beta[1] + beta[2] * sex[i]) * b_length[i]
  }
  # Priors:
  for (i in 1:2){
    alpha[i] ~ dnorm(0, 0.01)
    beta[i] ~ dnorm(0, 0.01)
  }
  sigma ~ dunif(0, 100)
  tau <- 1 / (sigma * sigma)
}




init_values <- function(){
  list(alpha = rnorm(2), beta = rnorm(2), sigma = runif(1))
}

params <- c("alpha", "beta", "sigma")

#tell what data to use, where to start gen data, file spec, chains, iterations, burnins, thinning, no deviance info criterion 
fit_lm2 <- jags(data = jagsdata_s2, inits = init_values, parameters.to.save = params, model.file = lm2_jags,
                n.chains = 3, n.iter = 12000, n.burnin = 2000, n.thin = 10, DIC = F)
fit_lm2


traceplot(fit_lm2, mfrow = c(2, 2), ask = F) 









#### Linear mixed model 
set.seed(42)

samplesize <- 200
nsites  <-  10 # number of sites in this model, assign total number places sampled 
b_length <-  sort(rnorm(samplesize)) 
sites <-  sample(1:10, samplesize,replace = T) #site or grouping variable
table(sites) #see how many samples get assigned to each site 


int_true_mean <- 45 # True mean intercept
int_true_sigma <- 10 # True standard deviation of intercepts
int_true_sites <- rnorm(n = nsites, mean = int_true_mean, sd = int_true_sigma) # True intercept of each site


# Intercept of each snake individual (depending on the site where it was captured):
sitemat <-  matrix(0,nrow = samplesize, ncol=nsites)

for (i in 1:nrow(sitemat)) sitemat [i,sites[i]] <-1
int_true <- sitemat %*% int_true_site

slope_true <-  10
mu <-  int_true + slope_true * b_length
sigma <-  5 #true stdev
b_mass <-  rnorm(samplesize,mean=mu, sd= sigma)
snakes3 <-  data.frame(b_length=b_length, b_mass=b_mass, site=sites)
head(snakes3)

plot(b_mass ~ b_length, col = site, data = snakes3)


#  include number of sites in list
Nsites <- length(levels(as.factor(snakes3$site)))

jagsdata_s3 <- with(snakes3, list(b_mass = b_mass, b_length = b_length, site = site,
                                  N = length(b_mass), Nsites = Nsites))


lm3_jags <- function(){
  # Likelihood:
  for (i in 1:N){
    b_mass[i] ~ dnorm(mu[i], tau) # tau is precision (1 / variance)
    mu[i] <- alpha + a[site[i]] + beta * b_length[i] # Random intercept for site
  }
  # Priors:
  alpha ~ dnorm(0, 0.01) # intercept
  sigma_a ~ dunif(0, 100) # standard deviation of random effect (variance between sites)
  tau_a <- 1 / (sigma_a * sigma_a) # convert to precision
  for (j in 1:Nsites){
    a[j] ~ dnorm(0, tau_a) # random intercept for each site
  }
  beta ~ dnorm(0, 0.01) # slope
  sigma ~ dunif(0, 100) # standard deviation of fixed effect (variance within sites)
  tau <- 1 / (sigma * sigma) # convert to precision
}

# each site gets its own init values
init_values <- function(){
  list(alpha = rnorm(1), sigma_a = runif(1), beta = rnorm(1), sigma = runif(1))
}

params <- c("alpha","beta","sigma","sigma_a")
fit_lm3 <-  jags( data=jagsdata_s3,inits=init_values,parameters.to.save = params, model.file = lm3_jags,
                  n.chains = 3, n.iter = 20000, n.burnin = 5000, n.thin = 10, DIC = F)
fit_lm3
