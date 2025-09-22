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
  dplyr::rename(Individual_1 = SampleID_1,Individual_2 = SampleID_2) %>% 
  mutate(nyears = Cohort_2 - Cohort_1)





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
  dplyr::rename(Individual_1 = SampleID_1,Individual_2 = SampleID_2) %>% 
  mutate(nyears = Cohort_2 - Cohort_1)

HSpairs <- HSpairs %>%
  mutate(nyears = Cohort_2 - Cohort_1)

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

kinships <- kinships %>% 
  mutate(
    sex_num = case_when(
      Sex_1 == "U" ~ NA_real_, # unknown to na
      Sex_1 == "M" ~ 1, # male is 1
      Sex_1 == "F" ~ 0 # female is 0
    )
  )













#### JAGS model ####
cat("model{
  
  # survival uniform 50-95%
  
  surv ~ dunif(0.5, 0.95)

  # proportion male in population (inverse of female, 1-propM)
  propM ~ dbeta(1,1)

  # Set up mean and stdev

  mu ~ dunif(1,100000)
  sd ~ dunif(1,100000)
  
  # Adult numbers
  
	for(i in 1:years) {
	
	  Nadult[i] ~ dnorm(mu, 1/(sd^2)) T(0, 1e+09)
  
	}
	
  for(j in 1:nobs) {
  
    TruePairs[j] ~ dbinom(
      (RObase[j] * (surv ^ AgeDif[j])) /
      ((Nadult[j] * knownSex[j] * abs(propM - IsFemale[j])) +
      (Nadult[j] * (1 - knownSex[j]))),
      Pair_viable_count[j]
    )
  }


}", file = "HSPPOPsingleYsex.jags")




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

knownSex <- kinships %>% 
  as.numeric(!is.na(kinships$sex_num))


# group together true by years
ntrue <- kinships %>%
  group_by(AgeDif,Cohort_2) %>%
  dplyr::summarise(
    ntrue = sum(TruePairs, na.rm = TRUE))

  
nobs = 

KnownSex = as.numeric(!is.na(kinships$sex_num))

IsFemale <- as.numeric(!is.na(kinships$sex_num) & kinships$sex_num == "0")


# baseline probability (4 for unsexed HSPs, 2 for unsexed POPs)

#### Run for single year ----
data = list( years = years,  # number of cohorts,
             PairTrue = kinships$TruePairs,
             AgeDif = kinships$AgeDif,
             Pair_viable_count = Cohort_years$Pair_viable_count,
             RObase = kinships$RObase,
             knownSex = KnownSex,
             nobs = nobs,
             IsFemale = IsFemale
             )

# Initial values
inits = function() {
  list(
    mu = runif(1, 1, 100000),
    sd = runif(1, 1, 100000),
    propM = runif(1,0.01, 0.99),
    Nadult = rep(10000, years)   # vector of length POP_years
  )
}

# Parameters to follow
params = c("Nadult","surv","propM")

nburn <- 3000
nchains <- 3
niter <- 5000
n.cores = 3


Out = jagsUI::jags(data, inits, params, "HSPPOPsingleYsex.jags", n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = T, verbose = T)
