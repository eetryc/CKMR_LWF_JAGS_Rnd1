#### Set up data frame for analysis with both HSPs and POPs ####

library(tidyverse)

# load in data table with full information 
original=read.csv("LWF_Info.csv",h=T)

original_clean <- original %>% 
  drop_na(FinalAge)


# select columns and add cohort year column

LWF_final <- original_clean %>%
  mutate(Cohort = SampleYear - FinalAge) %>% 
  select(SampleID,Cohort,SampleYear,SampleDate,Sex) %>% 
  mutate(SampleDate = mdy(SampleDate),
         Month = month(SampleDate))

### POPs filter 
POpairs <- LWF_final %>%
  rename_with(~ paste0(.x, "_1"), everything()) %>% # add the columns fro each indiv
  cross_join(LWF_final %>% rename_with(~ paste0(.x, "_2"), everything())) %>% # remove self-pairings
  mutate(conception = ymd(paste0(Cohort_2-1,"-10-01"))) %>% ###### Change as needed to fiddle with conception date
  # mutate to add the POP possibility (y/n) with the year and reproductive age filters
  mutate(HSPPOP_candidate = ifelse(Cohort_1 < Cohort_2 & 
                                     (Cohort_2 - Cohort_1 >= 5) & 
                                     (SampleYear_1 >= Cohort_2) & 
                                     (SampleDate_1 >= conception),"1","0")) %>% 
  #get rid of the ones that are the same comp, eg 1 and 2 vs 2 and 1 
  filter(SampleID_1 != SampleID_2) %>% 
  mutate(PairID = paste(SampleID_1, SampleID_2, sep = "_")) %>%
  mutate(older_first=ifelse(Cohort_1 < Cohort_2 |(Cohort_1 == Cohort_2 & SampleID_1 < SampleID_2), T,F) ) %>% 
  filter(older_first) %>%
  select(-older_first) %>% 
  mutate(AgeDif = Cohort_2 - Cohort_1) %>% 
  dplyr::rename(Individual_1 = SampleID_1,Individual_2 = SampleID_2) 



# count numnber of yes or no to get a better idea of how they are changing
POpairs %>% 
  dplyr::count(HSPPOP_candidate)


POPonly <- POpairs %>% 
  filter(HSPPOP_candidate == 1) %>% 
  mutate(RObase=2)


### HSPs and combine ###

HSpairs <- LWF_final %>%
  rename_with(~ paste0(.x, "_1"), everything()) %>% # sort into col for each indiv
  cross_join(LWF_final %>% rename_with(~ paste0(.x, "_2"), everything())) %>% # remove self-pairings
  mutate(conception = ymd(paste0(Cohort_2-1,"-10-01"))) %>% ###### Change as needed to fiddle with conception date
  mutate(HSPPOP_candidate = ifelse(Cohort_1 < Cohort_2 & (Cohort_1 != Cohort_2),"1","0")) %>% 
  #get rid of the ones that are the same comp, eg 1 and 2 vs 2 and 1 
  filter(SampleID_1 != SampleID_2) %>% 
  mutate(PairID = paste(SampleID_1, SampleID_2, sep = "_")) %>%
  mutate(older_first=ifelse(Cohort_1 < Cohort_2 |(Cohort_1 == Cohort_2 & SampleID_1 < SampleID_2), T,F) ) %>% 
  filter(older_first) %>%
  select(-older_first) %>% 
  mutate(AgeDif = Cohort_2 - Cohort_1) %>% 
  dplyr::rename(Individual_1 = SampleID_1,Individual_2 = SampleID_2) 

# count numnber of yes or no to get a better idea of how they are changing
HSpairs %>% 
  dplyr::count(HSPPOP_candidate)

HSPonly <- HSpairs %>% 
  filter(HSPPOP_candidate == 1) %>% 
  mutate(RObase=4)

# join the two tables
library(plyr)
kinships <- rbind(data.frame=POPonly,
                  data.frame=HSPonly) 

















#### JAGS model ####
cat("model{
	
	surv ~ dunif(0.5, 0.95)

  mu ~ dunif(1,100000)
  sd ~ dunif(1,100000)

  # Adult numbers
  
	for(i in 1:years) {

	  Nadult[i] ~ dnorm(mu, 1/(sd^2)) T(0, 1e+09)

		TruePairs[i] ~ dbinom((RObase[i]*(surv^AgeDif[i]))/(Nadult[i]), Pair_viable_count[i])

	}


}", file = "HSPPOPsingleY.jags")




#### Data to define for JAGS model ####

# Real POP/HSP - here assign randomly based on probability 

set.seed(777)
kinships <- kinships %>%
  mutate(TruePairs = case_when(
    RObase == 2 ~ rbinom(nrow(.), size = 1, prob = 0.001),
    RObase == 4 ~ rbinom(nrow(.), size = 1, prob = 0.0001),
    TRUE ~ 0
  ))

# Alternative to assign based on age difference likelihoods 
kinships <- kinships %>%
  mutate(prob = case_when(
    RObase == 2 ~ 0.01 * exp(-0.06 * AgeDif),
    RObase == 4 ~ 0.0001 * exp(-0.06 * AgeDif),
    TRUE ~ 0
  ),
  TruePairs = rbinom(nrow(.), size = 1, prob = prob))



true_per_year <- kinships %>% 
  group_by(Cohort_2,AgeDif) %>% 
  dplyr::summarise(
    Pair_true_count = as.integer(sum(TruePairs > 0, na.rm = TRUE))
  ) %>% 
  filter(Pair_true_count > 0) %>%     # keep only cohorts with >0
  arrange(Cohort_2)




Cohort_years <- kinships %>% 
  group_by(Cohort_2) %>% 
  dplyr::summarise(
    Pair_viable_count = as.integer(sum(HSPPOP_candidate > 0, na.rm = TRUE))
  ) %>% 
  filter(Pair_viable_count > 0) %>%     # keep only cohorts with >0
  arrange(Cohort_2)

years <- nrow(Cohort_years)

# baseline probability (4 for unsexed HSPs, 2 for unsexed POPs)

#### Run for single year ----
data = list( years = years,  # number of cohorts,
             TruePairs = kinships$TruePairs,
             AgeDif=kinships$AgeDif,
             Pair_viable_count = Cohort_years$Pair_viable_count,
             RObase = kinships$RObase)

# Initial values
inits = function() {
  list(
    mu = runif(1, 1, 100000),
    sd = runif(1, 1, 100000),
    Nadult = rep(10000, years),   # vector of length POP_years
    surv=rep(.6)
  )
}

# Parameters to follow
params = c("Nadult","surv")

nburn <- 3000
nchains <- 3
niter <- 5000
n.cores = 3


Out = jagsUI::jags(data, inits, params, "HSPPOPsingleY.jags", n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = T, verbose = T)

# intial, unedited run with this worked great, only bout 5 minutes to run and very good Rhat (1.000-1.002 for the 28 years)
