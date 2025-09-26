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


# select only potential pairs, add baseline prob (2 if unknown, 1 if sex known)
POPonly <- POpairs %>% 
  filter(HSPPOP_candidate == 1) %>% 
  mutate(RObase=ifelse(Sex_1 == "U",2,1)) # change known sex to 1 (not compensating for there being two possible parents because we know it's mom or dad)


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
kinshipsSex <- rbind(data.frame=POPonly,
                  data.frame=HSPonly) 

kinshipsSex <- kinshipsSex %>% 
  mutate(
    sex_num = case_when(
      Sex_1 == "U" ~ NA_real_, # unknown to na
      Sex_1 == "M" ~ 1, # male is 1
      Sex_1 == "F" ~ 0 # female is 0
    )
  )













#### JAGS model ####
cat("model{
  
  # survival uniform 50-95% (can be put in the j loop to estimate yearly survival)
  
  # proportion male in population (inverse of female, 1-propM)
  propM ~ dbeta(1,1)

  # Set up mean and stdev

  mu ~ dunif(1,100000)
  sd ~ dunif(1,100000)
  
  # Adult numbers
  
	for(i in 1:years) {
	
	  Nadult[i] ~ dnorm(mu, 1/(sd^2)) T(0, 1e+09)
  
    surv[i] ~ dunif(0.5, 0.95)

	}
	
  for(j in 1:nobs) {
    
    TruePairs[j] ~ dbinom(
      (RObase[j] * (surv[year[j]]^ AgeDif[j])) /
      ((Nadult[year[j]] * knownSex[j] * abs(propM - IsFemale[j])) +
      (Nadult[year[j]] * (1 - knownSex[j]))),
      Pair_viable_count_per_agedif[j]
    )
  }


}", file = "HSPPOPsingleYsex.jags")




#### Data to define for JAGS model ####



# years being estimated (years that actually have potential)
Cohort_years <- kinshipsSex %>% 
  group_by(Cohort_2) %>% 
  dplyr::summarise(
    Pair_viable_count = sum(HSPPOP_candidate > 0, na.rm = TRUE)
  ) %>% 
  filter(Pair_viable_count > 0) %>%     # keep only cohorts with >0
  arrange(Cohort_2)

# extract years (count to give to be i in JAGS)
years <- nrow(Cohort_years)



# group together by cohort year and age dif (differing probabilities)
Cohort_years_cohort_age <- kinshipsSex %>%
  group_by(Cohort_2, AgeDif) %>%
  dplyr::summarise(
    Pair_viable_count_per_agedif = sum(HSPPOP_candidate > 0, na.rm = TRUE))





# add the counts to each observation
kinshipsSex_obs <- kinshipsSex %>%
  left_join(Cohort_years_cohort_age, by = c("Cohort_2", "AgeDif"))

# year index to link kinship to years (just the number to match it, not the actual year)
year_index <- factor(kinshipsSex_obs$Cohort_2, levels = Cohort_years$Cohort_2)
kinshipsSex_obs <- kinshipsSex_obs %>%
  mutate(
    year_index = as.integer(year_index)
  )




# Count the ones with known sex (U = 0)
# Count females (F = 1)
kinshipsSex_obs <- kinshipsSex_obs %>%
  mutate(
    knownSex = as.numeric(!is.na(sex_num)),
    IsFemale = as.numeric(!is.na(sex_num) & (sex_num == 0))  # if sex_num encoded 0 = female
  )

# Real POP/HSP - here assign randomly based on probability 

set.seed(777)
kinshipsSex_obs <- kinshipsSex_obs %>% 
  mutate(TruePairs = case_when(
    RObase == 1 ~ rbinom(nrow(.), size = 1, prob = 0.001),
    RObase == 2 ~ rbinom(nrow(.), size = 1, prob = 0.001),
    RObase == 4 ~ rbinom(nrow(.), size = 1, prob = 0.0001),
    TRUE ~ 0
  ))




#### Now collapse based on cohort 2, age difference, RObase, if sex is known, and female or male 
group_vars <- c("Cohort_2", "AgeDif", "knownSex", "IsFemale", "RObase")

collapsed <- kinshipsSex_obs %>%
  group_by(year_index,across(all_of(group_vars))) %>%
  dplyr::summarise(
    Pair_viable_count_per_agedif = n(),                   
    TruePairs = sum(TruePairs, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  arrange(Cohort_2, AgeDif)



#### Run for single year ----
# use as.numeric if any are giving issues
# make sure all (except years for i loop) are equal to nobs! Check with length first...
data = list( years = years,  # number of cohorts,
             TruePairs = collapsed$TruePairs,
             AgeDif = collapsed$AgeDif,
             Pair_viable_count = collapsed$Pair_viable_count_per_agedif,
             RObase = collapsed$RObase,
             knownSex = as.numeric(collapsed$knownSex),
             nobs = nrow(collapsed),
             IsFemale = collapsed$IsFemale,
             Pair_viable_count_per_agedif = collapsed$Pair_viable_count_per_agedif,
             year=collapsed$year_index
             )

# Initial values
inits = function() {
  list(
    mu = runif(1, 1, 100000),
    sd = runif(1, 1, 100000),
    propM = runif(1,0.01, 0.99),
    Nadult = rep(10000, years),   # vector of length POP_years
    surv = (rep(.7,years))
  )
}

# Parameters to follow
params = c("Nadult","surv","propM")

nburn <- 3000
nchains <- 3
niter <- 5000
n.cores = 3


Out = jagsUI::jags(data, inits, params, "HSPPOPsingleYsex.jags", n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = T, verbose = T)
# initial did not have fantastic rhat values, specifically for the years with lowered potential and actual pairs (which makes sense, especially given the stochastic way I assigned the true pairs)
# but it did run in under 2 minutes
# Same deal with addition of survival to every year. these rhat values were pretty solid. however, again maybe not super accurate to do a per year estimate for this 

