# Use data from LWF_data_simulation test to get sample information (cohort years,sampling etc)



#### JAGS model ####
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
kinships$POP_viable = as.numeric((kinships$Cohort_1 + (5)) <= kinships$Cohort_2) # essentially POP_candidate, you can do either but this needs to be adjusted 


set.seed(777)
kinships <- kinships %>%
  mutate(TruePairs = if_else(
    POP_candidate == 1, #replace with pop viable in future
    rbinom(n(), size = 1, prob = 0.001),
    0
  ))

POP_Cohort_years <- kinships %>% 
  group_by(Cohort_2) %>% 
  summarise(
    POP_viable_count = sum(POP_candidate > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  filter(POP_viable_count > 0) %>%     # keep only cohorts with >0
  arrange(Cohort_2)

POP_years <- nrow(POP_Cohort_years)


#### Run for single year ----
data = list( POP_years = POP_years,  # number of cohorts,
             POPTrue = kinships$TruePairs,
             AgeDif=kinships$AgeDif,
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
