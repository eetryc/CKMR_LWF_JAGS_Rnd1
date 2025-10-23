# Use data from LWF_data_simulation test to get sample information (cohort years,sampling etc)


# Blanket survival rate
surv ~ dbeta(1, 1)

#Male proportion (inverse = female proportion)

propM ~ dbeta(1, 1)

# Adult numbers (can be split by sex) (uninformative prior)



#### JAGS model ####
# Set up JAGs model defining what the equations are and what variables need to be provided 
cat("model{
  surv ~ dunif(0.5, 0.95)

  mu ~ dunif(1,100000)
  sd ~ dunif(1,100000)

  # Adult numbers
  
	for(i in 1:POP_years) {
	
	  Nadult[i] ~ dnorm(mu,1/(sd^2)) T(0,)

		POPTrue[i] ~ dbinom((2 / Nadult[i])*(surv^AgeDif[i]), POP_viable_count[i])

	}


}", file = "POPsingleY.jags")






#### Data to define for JAGS model ####


# Viable POP pair (0/1)
kinships$POP_viable = as.numeric((kinships$Cohort_1 + (5)) <= kinships$Cohort_2) # this is essentially POP_candidate, but real iterations should be able to be put into this pipeline without set up

# Extract years to estimate for (using cohort of the 2nd individual)
POP_Cohort_years <- kinships %>% 
  group_by(Cohort_2) %>% 
  summarise(POP_viable_count = sum(POP_viable > 0, na.rm = TRUE), .groups = "drop") %>% 
  filter(POP_viable_count > 0) %>% 
  arrange(Cohort_2)


# Count number of true pairs (here, they are assigned based on observation probability)
kinships <-  kinships %>% 
  mutate(TruePairs = if_else(
    POP_candidate == 1, # only assign values to possible pairs. Should be changed to POP_viable in the future
    rbinom(n(), size = 1, prob = 0.001),0)) # generate the true pairs since we don't have this, at a rate of 1%

# Number of cohort years to estimate for 
POP_years <- length(unique(POPonly$Cohort_2))  # number of cohorts




#### Run for single year ----
data = list( POP_years = POP_years,  # number of cohorts,
             AgeDif=kinships$AgeDif,
             POPTrue = kinships$TruePairs,
             POP_viable_count = POP_Cohort_years$POP_viable_count)

# Initial values
inits = function() {
  list(
    mu = runif(1, 1, 100000),
    sd = runif(1, 1, 100000),
    Nadult = rep(10000, POP_years)   # vector of length POP_years
  )
}

# Parameters to follow
params = c("Nadult")

nburn <- 3000
nchains <- 3
niter <- 5000
n.cores = 3


Out = jagsUI::jags(data, inits, params, "POPsingleY.jags", n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = T, verbose = T)
