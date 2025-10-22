#### Set up data frame for analysis with both HSPs and POPs ####

library(tidyverse)

setwd("~/GitHub/CKMR_LWF_JAGS_Rnd1")

# load in data table with full information 
original=read.csv("LWF_Info.csv",h=T)

original_clean <- original %>% 
  drop_na(FinalAge)


# select columns and add cohort year column

LWF_final <- original_clean %>%
  mutate(Cohort = SampleYear - FinalAge) %>% 
  select(SampleID,Cohort,SampleYear,SampleDate,Sex,Length_mm) %>% 
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
	

  mu ~ dunif(1,100000)
  sd ~ dunif(1,100000)

  # Adult numbers
  
	for(i in 1:years) {

	  surv[i] ~ dunif(0.5, 0.95)

	  Nadult[i] ~ dnorm(mu, 1/(sd^2)) T(0, 1e+09)

		TruePairs[i] ~ dbinom((RObase[i]*(surv[i]^AgeDif[i]))/(Nadult[i]), Pair_viable_count[i])

	}


}", file = "HSPPOPsingleYsurv.jags")




#### Data to define for JAGS model ####

# Real POP/HSP - here assign randomly based on probability 

set.seed(777)
kinships <- kinships %>%
  mutate(TruePairs = case_when(
    RObase == 2 ~ rbinom(nrow(.), size = 1, prob = 0.001),
    RObase == 4 ~ rbinom(nrow(.), size = 1, prob = 0.0001),
    TRUE ~ 0
  ))

Cohort_years <- kinships %>% 
  group_by(Cohort_2) %>% 
  dplyr::summarise(
    Pair_viable_count = sum(HSPPOP_candidate > 0, na.rm = TRUE)
  ) %>% 
  filter(Pair_viable_count > 0) %>%     # keep only cohorts with >0
  arrange(Cohort_2)

years <- nrow(Cohort_years)


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
    surv=rep(.6,years)
  )
}

# Parameters to follow
params = c("Nadult","surv")

nburn <- 3000
nchains <- 3
niter <- 5000
n.cores = 3


Out = jagsUI::jags(data, inits, params, "HSPPOPsingleYsurv.jags", n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = T, verbose = T)

# inital attempt worked well, took about 5 minutes and produced great rhat values for both

























#### Grouped or binned years

cat("model{
	
  mu ~ dunif(1,100000)
  sd ~ dunif(1,100000)
  
  # Grouped survival estimate
  for(g in 1:ngroups) {
    surv[g] ~ dunif(0.5,0.95)
  }


  # Adult numbers
	for(i in 1:years) {
	  Nadult[i] ~ dnorm(mu, 1/(sd^2)) T(1, 1e+09)
	}
	  
	for (j in 1:nobs) {
		TruePairs[j] ~ dbinom((RObase[j]*(surv[group_index[j]]^AgeDif[j]))/(Nadult[year_index[j]]), Pair_viable_count[j])
	}


}", file = "HSPPOPsingleYsurvBin.jags")




#### Data to define for JAGS model ####


## Cohort years for the i loop
Cohort_years <- kinships %>% 
  group_by(Cohort_2) %>% 
  dplyr::summarise(
    Pair_viable_count = sum(HSPPOP_candidate > 0, na.rm = TRUE)
  ) %>% 
  filter(Pair_viable_count > 0) %>%     # keep only cohorts with >0
  arrange(Cohort_2)

years <- nrow(Cohort_years)


kinships <- kinships %>% 
  mutate(
    sex_num = case_when(
      Sex_1 == "U" ~ NA_real_, # unknown to na
      Sex_1 == "M" ~ 1, # male is 1
      Sex_1 == "F" ~ 0 # female is 0
    )
  )

# add the counts to each observation
Cohort_years_cohort_age <- kinships %>% group_by(Cohort_2, AgeDif) %>% dplyr::summarise( Pair_viable_count_per_agedif = sum(HSPPOP_candidate > 0, na.rm = TRUE))

kinships_obs <- kinships %>%
  left_join(Cohort_years_cohort_age, by = c("Cohort_2", "AgeDif"))

# year index to link kinship to years (just the number to match it, not the actual year)
year_index <- factor(kinships_obs$Cohort_2, levels = Cohort_years$Cohort_2)
kinships_obs <- kinships_obs %>%
  mutate(
    year_index = as.integer(year_index)
  )




# Count the ones with known sex (U = 0)
# Count females (F = 1)
kinships_obs <- kinships_obs %>%
  mutate(
    knownSex = as.numeric(!is.na(sex_num)),
    IsFemale = as.numeric(!is.na(sex_num) & (sex_num == 0))  # if sex_num encoded 0 = female
  )

# Real POP/HSP - here assign randomly based on probability 

set.seed(777)
kinships_obs <- kinships_obs %>% 
  mutate(TruePairs = case_when(
    RObase == 1 ~ rbinom(nrow(.), size = 1, prob = 0.001),
    RObase == 2 ~ rbinom(nrow(.), size = 1, prob = 0.001),
    RObase == 4 ~ rbinom(nrow(.), size = 1, prob = 0.0001),
    TRUE ~ 0
  ))




#### Now collapse based on cohort 2, age difference, RObase, if sex is known, and female or male 
group_vars <- c("Cohort_2", "AgeDif", "knownSex", "IsFemale", "RObase")

collapsed <- kinships_obs %>%
  group_by(year_index,across(all_of(group_vars))) %>%
  dplyr::summarise(
    Pair_viable_count_per_agedif = n(),                   
    TruePairs = sum(TruePairs, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  arrange(Cohort_2, AgeDif)


## Grouped years to measure survival across time 
# use a group index
collapsed <- collapsed %>%
  mutate(YearGroup = case_when(
    Cohort_2 <= 2000 ~ 1,
    Cohort_2 <= 2010 ~ 2,
    TRUE ~ 3
  ))


group_levels <- sort(unique(collapsed$YearGroup))
ngroups <- length(group_levels)

group_index <- as.integer(factor(collapsed$YearGroup, levels = group_levels))


#### Run for single year ----
data = list( years = years,  # number of cohorts,
             TruePairs = collapsed$TruePairs,
             AgeDif=collapsed$AgeDif,
             Pair_viable_count = collapsed$Pair_viable_count_per_agedif,
             RObase = collapsed$RObase,
             ngroups = ngroups,
             nobs = nrow(collapsed),
             group_index = group_index,
             year_index = collapsed$year_index
)

# Initial values
inits = function() {
  list(
    mu = runif(1, 1, 100000),
    sd = runif(1, 1, 100000),
    Nadult = rep(10000, years),   # vector of length POP_years
    surv=rep(.6,ngroups)
  )
}

# Parameters to follow
params = c("Nadult","surv")

nburn <- 3000
nchains <- 3
niter <- 5000
n.cores = 3


Out = jagsUI::jags(data, inits, params, "HSPPOPsingleYsurvBin.jags", n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = T, verbose = T)
