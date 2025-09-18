## Setup ----
library(rjags)
library(runjags)
library(postpack)

## Convert input data to JAGS format ----

#Total everything up
true.kin = rbind(
  true.PO %>%
    #mutate(SpawnSpecific = 0) %>%
    group_by(Cohort_2,nyears) %>%
    summarize(PotPairs = n(), TruePairs = sum(isPO)) %>%
    rename(Year = Cohort_2, PairMortYears = nyears) %>%
    mutate(RObase = 2, Group = "PO"),
  true.HS %>%
    group_by(Cohort_1,nyears) %>%
    summarize(PotPairs = n(), TruePairs = sum(isHS)) %>%
    rename(Year = Cohort_1, PairMortYears = nyears) %>%
    mutate(RObase = 2, Group = "HS") %>%
    arrange(desc(Year), PairMortYears))

true.kin = true.kin %>% filter(PairMortYears <= 5)

# Notes to Chloe:
#Cohort_1 = cohort of older half of pair, Cohort_2 = younger
#isPO and isHS is a logical 0/1 indicating whether genotyping confirmed kinship
#PairMortYears is different depending on relationship type
#For POPs it is years between pot. parent being sampled and pot. offspring being conceived with all negative values set to zero
#For HSPs it is the difference in years between the older and younger cohorts (always at least one because same-cohort comparisons are removed)
#RObase is there for flexibility. Allows single sex models (1 chance of POP per indiv., RObase=1), accounting for fecundity/capture prob interactions, etc.

## Bring the data into JAGs ----

## Add to JAGS list
jags_data = list(
  ncombs = dim(true.kin)[1],
  PotPairs = true.kin$PotPairs,
  TruePairs = true.kin$TruePairs,
  RObase = true.kin$RObase,
  PairMortYears = true.kin$PairMortYears,
)


## Specify the model and write to file
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


Nsimple.file = "D:\\Google Drive\\Research\\Postdoc\\Projects\\LT_Genetics\\CKMR_Nsimple.txt"
write_model(CKMR_Nsimple, Nsimple.file)

## Run JAGS and collect output
post.Nsimple = jagsUI::jags(data = jags_data,  
                    Nsimple.file, 
                    inits = jags_inits(jags_dims["nc"]),parameters.to.save = c("Nadult_all", "surv"), 
                    n.burnin = 50000, n.chains = 3, 
                    n.iter = 100000, parallel = T,  verbose = T)


