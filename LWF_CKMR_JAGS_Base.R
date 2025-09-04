#### First attempt 

## use the data frame created by CKMR_LWF_Data_Sim




#### Calculate mortality and whether potential pairs are viable ----
kinships 
POPonly <- kinships %>% 
  filter(POP_candidate == 1)

#Get number of years that prev-sampled adult had to survive to spawn during juvenile cohort year
POPonly$ParMortYears = POPonly$Cohort_2 - POPonly$SampleYear_1
#If negative (sampled adult after juvenile cohort year) set to zero
POPonly[which(POPonly$ParMortYears < 0),]$ParMortYears = 0

#Calculate whether the pair has potential to be parent-offspring based on birth years
POPonly$OldEno = as.numeric((POPonly$Cohort_1 + (5)) <= POPonly$Cohort_2)
#Multiply by number of chances (two if unknown sex, one if known)
POPonly$OldEno = POPonly$OldEno * 2


cohort_counts <- POPonly %>% 
  group_by(Cohort_2) %>% 
  summarise(OldEno_count = sum(OldEno > 0, na.rm = TRUE), .groups = "drop") %>% 
  arrange(Cohort_2)


#### Single year JAGs model ----
cat("model{

	# Adult numbers (uninformative prior)

	for(i in 1:nad) {
		Nadult[i] ~ dnorm(100000, 1.0E-6) T(0,)
		POPs[i] ~ dbinom(2 / Nadult[i], OldEno_count[i])

	}


}", file = "CKMRsingleY.jags")
# move nad out of loop (1 nad being est), no i, take i off po data (every year nad is the same)







#### Set up input data ----
BYear = POPonly$Cohort_2 #Birth years
BYear <- rank(unique(POPonly$Cohort_2))[match(POPonly$Cohort_2, unique(POPonly$Cohort_2))] #Reduce to index

OldEno = POPonly$OldEno #Indicator of viable combination. Vector of 0/2 for multi-year, per-cohort vector length for single year

set.seed(777)
ncom <- nrow(POPonly) 
POPonly <-  POPonly %>% 
  mutate(TruePairs = rbinom(n = ncom, size = 1, prob = 0.001)) # generate the true pairs since we don't have this, at a rate of 1%

ParYear = POPonly$ParMortYears #Number of years between adult sampling and subsequent juvenile sampling for multi-year model

nad <- length(unique(POPonly$Cohort_2))  # number of cohorts

ncom = length(ROs) #Number of pairs to test


#### Run for single year ----
data = list( nad = length(unique(POPonly$Cohort_2)),  # number of cohorts,
            ncom = nrow(POPonly),
            POPs = POPonly$TruePairs,
            BYear = BYear,
            OldEno_count = cohort_counts$OldEno_count)

# Initial values
inits = function() {
  list(Nadult = rep(10000,nad))
}

# Parameters to follow
params = c("Nadult")

nburn <- 3000
nchains <- 3
niter <- 5000
n.cores = 3


Out = jagsUI::jags(data, inits, params, "CKMRsingleY.jags", n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = T, verbose = T)

























Outdata = list(
  ncombs = dim(POPonly)[1],
  PotPairs = POPonly$POP_candidate,
  TruePairs = POPonly$TruePairs,
  RObase = POPonly$RObasePOP,
  PairMortYears = POPonly$ParMortYears
)






CKMR_Nsimple = function()
{
  mu ~ dunif(1,1000000)
  sd ~ dunif(1,1000000)
  Nadult_all ~ dnorm(mu, 1/(sd^2))
  surv ~ dunif(0.05, 0.95)
  
  for(j in 1:ncombs)
  {
    TruePairs[j] ~ dbinom((RObase[j] * surv^PairMortYears[j]) / 
                            Nadult_all,PotPairs[j])
  }
}



jags_inits = function(nc) {
  inits = list()
  for(c in 1:nc){
    inits[[c]] = list(
      surv = 0.5,
      Nadult_all = runif(1,1000000)
    )
  }
  
  return(inits)
}


post.Nsimple = jagsUI::jags(data = data,  
                            Nsimple.file, 
                            inits = jags_inits(jags_dims["nc"]),parameters.to.save = c("Nadult_all", "surv"), 
                            n.burnin = 50000, n.chains = 3, 
                            n.iter = 100000, parallel = T,  verbose = T)

