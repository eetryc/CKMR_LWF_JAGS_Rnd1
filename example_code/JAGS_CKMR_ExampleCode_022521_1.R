#### Calculate mortality and whether potential pairs are viable ----

#Get number of years that prev-sampled adult had to survive to spawn during juvenile cohort year
Fish.dt$ParMortYears = Fish.dt$OffBYear - Fish.dt$ParYear
#If negative (sampled adult after juvenile cohort year) set to zero
Fish.dt[which(Fish.dt$ParMortYears < 0),]$ParMortYears = 0

#Calculate whether the pair has potential to be parent-offspring based on birth years
Fish.dt$OldEno = as.numeric((Fish.dt$ParBYear + (age.maturity)) <= Fish.dt$OffBYear)
#Multiply by number of chances (two if unknown sex, one if known)
Fish.dt$OldEno = Fish.dt$OldEno * 2

#### Single year JAGs model ----
cat("model{

	# Adult numbers (uninformative prior)

	for(i in 1:nad) {
		Nadult[i] ~ dnorm(0,1.0E-6)
		POPs[i] ~ dbinom(2 / Nadult[i], OldEno[i])

	}


}", file = "C:\\Users\\Ben\\Google Drive\\Research\\HAL CKMR\\Estimation Model\\JAGS_Test\\CKMRsingleY.jags")



#### Multi-year JAGs model ----
cat("model{

	# Adult numbers (uninformative prior)
	for(i in 1:nad) {
		Nadult[i] ~ dnorm(0,1.0E-6)
	}

	M ~ dbeta(1,1)

	for(j in 1:ncom) {
	  POPs[j] ~ dbern((OldEno[j] * ((M)^ParYear[j]))/(Nadult[BYear[j]])) # Repro chances reduced by chance that adult died prior to spawning
	}

}", file = "C:\\Users\\Ben\\Google Drive\\Research\\HAL CKMR\\Estimation Model\\JAGS_Test\\CKMRmultiY.jags")

#### Set up input data ----
BYear = Fish.dt$OffBYear #Birth years
BYear = rank(unique(BYear))[match(BYear,unique(BYear))] #Reduce to index

POPs = Fish.dt$POP #Observed POPs. Vector of 0/1 for multi-year, per-cohort sum of 0/1 for single year (binomial faster than repeated bernoullis when feasible)
OldEno = Fish.dt$OldEno #Indicator of viable combination. Vector of 0/2 for multi-year, per-cohort vector length for single year

ParYear = Fish.dt$ParMortYears #Number of years between adult sampling and subsequent juvenile sampling for multi-year model

nad = length(unique(Fish.dt$OffBYear)) #Number of cohorts
ncom = length(POPs) #Number of pairs to test

#### Run for single year ----
data = list(nad = nad,
            POPs = POPs,
            ROs = ROs)

# Initial values
inits = function() {
  list(Nadult = rep(1000,nad))
}

# Parameters to follow
params = c("Nadult")

Out = jags(data, inits, params, "C:\\Users\\Ben\\Google Drive\\Research\\HAL CKMR\\Estimation Model\\JAGS_Test\\CKMRsingleY.jags", 
           n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = parallel, n.cores = n.cores, verbose = verbose)

#### Or run for multi-year ----

# Data for JAGS
data = list(nad = nad,
            ncom = ncom,
            BYear = BYear,
            POPs = POPs,
            OldEno = ROs,
            ParYear = ParYear)

# Initial values
inits = function() {
  list(Nadult = rep(1000,nad),
       M = 0.5)
}

# Parameters to follow
params = c("Nadult", "M")

Out = jags(data, inits, params, "C:\\Users\\Ben\\Google Drive\\Research\\HAL CKMR\\Estimation Model\\JAGS_Test\\CKMRmultiY.jags", 
           n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = parallel, n.cores = n.cores, verbose = verbose)

