#### Set up data frame for analysis with both HSPs and POPs ####
setwd("~/GitHub/CKMR_LWF_JAGS_Rnd1")

library(tidyverse)

# load in data table with full information 
original=read.csv("LWF_Info.csv",h=T)

original_clean <- original %>% 
  drop_na(FinalAge)


# select columns and add cohort year column
set.seed(102)
library(dplyr)
library(lubridate)

set.seed(123)

LWF_final_GMZ_mixed <- original_clean %>%
  mutate(
    Cohort = SampleYear - FinalAge,
    SampleDate = mdy(SampleDate),
    Month = month(SampleDate)
  ) %>%
  select(SampleID, Cohort, SampleYear, SampleDate, Sex, Month,Length_mm) %>%
  mutate(
    Stock = case_when(
      Month %in% 9:12 ~ sample(
        c("GMZ1N", "GMZ1S", "GMZ2", "GMZ3", "GMZ4", "GMZ5"),
        size = nrow(.), replace = TRUE,
        prob = c(0.02, 0.02, 0.02, 0.9, 0.02, 0.02)
      ),
      Month %in% 1:8 ~ sample(
        c("GMZ1N", "GMZ1S", "GMZ2", "GMZ3", "GMZ4", "GMZ5"),
        size = nrow(.), replace = TRUE,
        prob = c(0.05, 0.00, 0.35, 0.3, 0.15, 0.15)
      ),
      TRUE ~ NA_character_
    )
  )


# test for a little dispersal
table(unlist(LWF_final_GMZ_mixed$Stock))

# visualize stock assignment distribution 
ggplot(LWF_final_GMZ_mixed, aes(x = Stock)) +
  geom_bar() +
  theme_minimal() +
  xlab("Source") +
  ylab("Percent of Individuals Caught in WFM-06") 


# well mixed
set.seed(123)

LWF_final_well_mixed <- original_clean %>%
  mutate(
    Cohort = SampleYear - FinalAge,
    SampleDate = mdy(SampleDate),
    Month = month(SampleDate)
  ) %>%
  select(SampleID, Cohort, SampleYear, SampleDate, Sex, Month) %>%
  mutate(
    Stock = case_when(
      Month %in% 9:12 ~ sample(
        c("GMZ1N", "GMZ1S", "GMZ2", "GMZ3", "GMZ4", "GMZ5"),
        size = nrow(.), replace = TRUE,
        prob = c(0.167, 0.167, 0.167, 0.167, 0.167, 0.167)
      ),
      Month %in% 1:8 ~ sample(
        c("GMZ1N", "GMZ1S", "GMZ2", "GMZ3", "GMZ4", "GMZ5"),
        size = nrow(.), replace = TRUE,
        prob = c(0.167, 0.167, 0.167, 0.167, 0.167, 0.167)
      ),
      TRUE ~ NA_character_
    )
  )


# test for a little dispersal
table(unlist(LWF_final_well_mixed$Stock))




# visualize stock assignment distribution 
ggplot(LWF_final_well_mixed, aes(x = Stock),colo) +
  geom_bar() +
  theme_minimal() +
  xlab("Source") +
  ylab("Percent of Individuals Caught in WFM-06") 



### POPs filter 
POpairs <- LWF_final_GMZ_mixed %>%
  rename_with(~ paste0(.x, "_1"), everything()) %>% # add the columns fro each indiv
  cross_join(LWF_final_GMZ_mixed %>% rename_with(~ paste0(.x, "_2"), everything())) %>% # remove self-pairings
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
  mutate(StockPair = paste(Stock_1, Stock_2, sep = "_")) %>% 
  mutate(nyears = Cohort_2 - Cohort_1)




# count numnber of yes or no to get a better idea of how they are changing
POpairs %>% 
  dplyr::count(HSPPOP_candidate)


POPonly <- POpairs %>% 
  filter(HSPPOP_candidate == 1) %>% 
  mutate(RObase=ifelse(Sex_1 == "U",2,1)) %>% 
  mutate(Stock_1 = unlist(Stock_1),
         Stock_2 = unlist(Stock_2),
         StockWeight = ifelse(Stock_1 == Stock_2, 1, 0.1)) # likllihoods for if they are from the same place vs not


### HSPs and combine ###

HSpairs <- LWF_final_GMZ_mixed %>%
  rename_with(~ paste0(.x, "_1"), everything()) %>% # sort into col for each indiv
  cross_join(LWF_final_GMZ_mixed %>% rename_with(~ paste0(.x, "_2"), everything())) %>% # remove self-pairings
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
  mutate(StockPair = paste(Stock_1, Stock_2, sep = "_"))%>% 
  mutate(nyears = Cohort_2 - Cohort_1)


# count numnber of yes or no to get a better idea of how they are changing
HSpairs %>% 
  dplyr::count(HSPPOP_candidate)

HSPonly <- HSpairs %>% 
  filter(HSPPOP_candidate == 1) %>% 
  mutate(RObase=4) %>% 
  mutate(Stock_1 = unlist(Stock_1),
         Stock_2 = unlist(Stock_2),
         StockWeight = ifelse(Stock_1 == Stock_2, 1, 0.1)) # liklihoods for if they are from the same place vs not

# join the two tables
library(plyr)
kinshipsSexStock <- rbind(data.frame=POPonly,
                       data.frame=HSPonly)  %>% 
  mutate(
    sex_num = case_when(
      Sex_1 == "U" ~ NA_real_, # unknown to na
      Sex_1 == "M" ~ 1, # male is 1
      Sex_1 == "F" ~ 0 # female is 0
    )
  )





























#### JAGS model ####
cat("model{

  # proportion male in population (inverse of female, 1-propM)
  propM ~ dbeta(1,1)
  
  
  mu ~ dunif(1,100000)
  sd ~ dunif(1,100000)
  
  # Adult numbers
	for(i in 1:years) {
	
	  Nadult[i] ~ dnorm(mu, 1/(sd^2)) T(1, 1e+09)

    surv[i] ~ dunif(0.5, 0.95)

	}
	
	
  for(j in 1:nobs) {
    
    TruePairs[j] ~ dbinom(
      (RObase[j] * (surv[year[j]]^ AgeDif[j])*StockWeight[j]) /
      ((Nadult[year[j]] * knownSex[j] * abs(propM - IsFemale[j])) +
      (Nadult[year[j]] * (1 - knownSex[j]))),
      Pair_viable_count_per_agedif[j]
    )
  }	

}", file = "HSPPOP.sex.stocks.singleY.jags")




#### Data to define for JAGS model ####


# years being estimated (years that actually have potential)
Cohort_years <- kinshipsSexStock %>% 
  group_by(Cohort_2) %>% 
  dplyr::summarise(
    Pair_viable_count = sum(HSPPOP_candidate > 0, na.rm = TRUE)
  ) %>% 
  filter(Pair_viable_count > 0) %>%     # keep only cohorts with >0
  arrange(Cohort_2)

# extract years (count to give to be i in JAGS)
years <- nrow(Cohort_years)


# group together by cohort year and age dif (differing probabilities)
Cohort_years_cohort_age <- kinshipsSexStock %>%
  group_by(Cohort_2, AgeDif) %>%
  dplyr::summarise(
    Pair_viable_count_per_agedif = sum(HSPPOP_candidate > 0, na.rm = TRUE))


# add the counts to each observation
kinshipsSex_obs <- kinshipsSexStock %>%
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

# parameters for decay rate of pair likelihood with age
HS_true <- 0.0001 # edit this number to change the total percentage of pairs that are true
HS_beta <- 0.3 # decay rate of half siblings based on survival of shared parent 
PO_true <- 0.001 # edit this number to change the total percentage of pairs that are true



# assign true pairs systematically 
kinshipsSex_obs <- kinshipsSex_obs %>% 
  mutate(trueProb = case_when(
    RObase == 1 ~ PO_true,
    RObase == 2 ~ PO_true,
    RObase == 4 ~ HS_true*exp(-HS_beta*abs(AgeDif)),
    TRUE ~ 0),
    TruePairs = rbinom(nrow(.), size = 1, prob = trueProb))


#### Now collapse based on cohort 2, age difference, RObase, if sex is known, and female or male 
group_vars <- c("Cohort_2", "AgeDif", "knownSex", "IsFemale", "RObase","StockWeight")

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
             year=collapsed$year_index,
             StockWeight = collapsed$StockWeight
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


Out = jagsUI::jags(data, inits, params, "HSPPOP.sex.stocks.singleY.jags", n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = T, verbose = T)
# worked initially will a few funky Rhat values for some Nhat years where Neff was low. 

# if issues occur, run un parallelized
Out <- jagsUI::jags(
  data = data,
  inits = inits,
  parameters.to.save = params,
  model.file = "HSPPOP.sex.stocks.singleY.jags",
  n.burnin = nburn,
  n.chains = nchains,
  n.iter = niter,
  parallel = FALSE,
  verbose = TRUE
)













# Extract posterior summary as a data frame
nhat <- as.data.frame(Out$summary)

# Keep only the Nadult rows
nhat <- nhat[grep("^Nadult", rownames(nhat)), ]

# Add a Year index (1, 2, 3, ...)
nhat$Year.index <- 1:nrow(nhat)
nhat$Year <-1993:2020
# Select the useful columns
nhat_df <- nhat %>%
  dplyr::select(Year, mean, `2.5%`, `50%`, `97.5%`)


library(ggplot2)

ggplot(nhat_df, aes(x = Year, y = mean)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue") +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, fill = "lightblue") +
  labs(
    x = "Year (Cohort index)",
    y = expression(hat(N)[adult]),
    title = "Posterior estimates of adult abundance (N-hat) 
  in a mixed fishery"
  ) +
  theme_minimal(base_size = 14)







# Extract posterior summary as a data frame
Rhat <- as.data.frame(Out$summary)

# Keep only the Rhat rows
Rhat <- Rhat[grep("^Nadult", rownames(Rhat)), ]

# Add a Year index (1, 2, 3, ...)
Rhat$Year <- 1:nrow(Rhat)

# Select the useful columns
Rhat_df <- Rhat %>%
  dplyr::select(Year, mean, `2.5%`, `50%`, `97.5%`,Rhat)
ggplot(Rhat_df, aes(x = Year, y = Rhat)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue") +
  labs(
    x = "Year (Cohort index)",
    y = "Convergence (R-hat)",
    title = "Posterior estimates of Convergence (R-hat) 
  in a well mixed fishery"
  ) +
  theme_minimal(base_size = 14)


















#################### Test mixing effects on total, GMZ3 estimates of abundance


#### Now test with only the fish in GMZ3



kinshipsGMZ3 <- kinshipsSexStock



# Real POP/HSP - here assign randomly based on probability 

set.seed(777)

# parameters for decay rate of pair likelihood with age
HS_trueGMZ3 <- 0.0001 # edit this number to change the total percentage of pairs that are true
HS_betaGMZ3 <- 0.3 # decay rate of half siblings based on survival of shared parent 
PO_trueGMZ3 <- 0.001 # edit this number to change the total percentage of pairs that are true



# assign true pairs systematically 
kinshipsSex_obsGMZ3 <- kinshipsGMZ3 %>% 
  mutate(trueProb = case_when(
    RObase == 1 ~ PO_true,
    RObase == 2 ~ PO_true,
    RObase == 4 ~ HS_true*exp(-HS_beta*abs(AgeDif)),
    TRUE ~ 0),
    TruePairs = rbinom(nrow(.), size = 1, prob = trueProb))


kinshipsSex_obsGMZ3 <- kinshipsSex_obsGMZ3 %>% 
  filter(StockPair == "GMZ3_GMZ3")



# extract years (count to give to be i in JAGS)
Cohort_yearsGMZ3 <- kinshipsSex_obsGMZ3 %>% 
  group_by(Cohort_2) %>% 
  dplyr::summarise(
    Pair_viable_count = sum(HSPPOP_candidate > 0, na.rm = TRUE)
  ) %>% 
  filter(Pair_viable_count > 0) %>%     # keep only cohorts with >0
  arrange(Cohort_2)
yearsGMZ3 <- nrow(Cohort_yearsGMZ3)


# group together by cohort year and age dif (differing probabilities)
Cohort_years_cohort_ageGMZ3 <- kinshipsSex_obsGMZ3 %>%
  group_by(Cohort_2, AgeDif) %>%
  dplyr::summarise(
    Pair_viable_count_per_agedif = sum(HSPPOP_candidate > 0, na.rm = TRUE))


# add the counts to each observation
kinshipsSex_obsGMZ3 <- kinshipsSex_obsGMZ3 %>%
  left_join(Cohort_years_cohort_ageGMZ3, by = c("Cohort_2", "AgeDif"))





# Count the ones with known sex (U = 0)
# Count females (F = 1)
kinshipsSex_obsGMZ3 <- kinshipsSex_obsGMZ3 %>%
  mutate(
    knownSex = as.numeric(!is.na(sex_num)),
    IsFemale = as.numeric(!is.na(sex_num) & (sex_num == 0))  # if sex_num encoded 0 = female
  )






# year index to link kinship to years (just the number to match it, not the actual year)
year_indexGMZ3 <- factor(kinshipsSex_obsGMZ3$Cohort_2, levels = Cohort_years$Cohort_2)
kinshipsSex_obsGMZ3 <- kinshipsSex_obsGMZ3 %>%
  mutate(
    year_index = as.integer(year_indexGMZ3)
  )

#### Now collapse based on cohort 2, age difference, RObase, if sex is known, and female or male 
group_varsGMZ3 <- c("Cohort_2", "AgeDif", "knownSex", "IsFemale", "RObase","StockWeight")

collapsedGMZ3 <- kinshipsSex_obsGMZ3 %>%
  group_by(year_indexGMZ3,across(all_of(group_varsGMZ3))) %>%
  dplyr::summarise(
    Pair_viable_count_per_agedif = n(),                   
    TruePairs = sum(TruePairs, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  arrange(Cohort_2, AgeDif)



#### Run for single year ----
# use as.numeric if any are giving issues
# make sure all (except years for i loop) are equal to nobs! Check with length first...
GMZ3data = list( years = yearsGMZ3,  # number of cohorts,
             TruePairs = collapsedGMZ3$TruePairs,
             AgeDif = collapsedGMZ3$AgeDif,
             RObase = collapsedGMZ3$RObase,
             knownSex = as.numeric(collapsedGMZ3$knownSex),
             nobs = nrow(collapsedGMZ3),
             IsFemale = collapsedGMZ3$IsFemale,
             Pair_viable_count_per_agedif = collapsedGMZ3$Pair_viable_count_per_agedif,
             year=as.numeric(collapsedGMZ3$year_indexGMZ3),
             StockWeight = collapsedGMZ3$StockWeight
)

# Initial values
initsGMZ3 = function() {
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


OutGMZ3 = jagsUI::jags(GMZ3data, initsGMZ3, params, "HSPPOP.sex.stocks.singleY.jags", n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = T, verbose = T)
# worked initially will a few funky Rhat values for some Nhat years where Neff was low. 

# if issues occur, run un parallelized
OutGMZ3 <- jagsUI::jags(
  data = GMZ3data,
  inits = initsGMZ3,
  parameters.to.save = params,
  model.file = "HSPPOP.sex.stocks.singleY.jags",
  n.burnin = nburn,
  n.chains = nchains,
  n.iter = niter,
  parallel = FALSE,
  verbose = TRUE
)













# Extract posterior summary as a data frame
nhatGMZ3 <- as.data.frame(OutGMZ3$summary)
nhatGMZ3 <- nhatGMZ3[grep("^Nadult", rownames(nhatGMZ3)), ]
# Add a Year index (1, 2, 3, ...)
nhat_dfGMZ3$Year.index <- 1:nrow(nhat)
nhat_dfGMZ3$Year <-1993:2020

# Select the useful columns
nhat_dfGMZ3 <- nhat_dfGMZ3 %>%
  dplyr::select(Year, mean, `2.5%`, `50%`, `97.5%`)


library(ggplot2)

ggplot(nhat_dfGMZ3, aes(x = Year, y = mean)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue") +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, fill = "lightblue") +
  labs(
    x = "Year (Cohort index)",
    y = expression(hat(N)[adult]),
    title = "Posterior estimates of adult abundance (N-hat) in GMZ3
  in a mixed fishery"
  ) +
  theme_minimal(base_size = 14)



















#### Compare annual results to a well mixed fishery 
# well mixed
set.seed(123)

LWF_final_well_mixed <- original_clean %>%
  mutate(
    Cohort = SampleYear - FinalAge,
    SampleDate = mdy(SampleDate),
    Month = month(SampleDate)
  ) %>%
  select(SampleID, Cohort, SampleYear, SampleDate, Sex, Month) %>%
  mutate(
    Stock = case_when(
      Month %in% 9:12 ~ sample(
        c("GMZ1N", "GMZ1S", "GMZ2", "GMZ3", "GMZ4", "GMZ5"),
        size = nrow(.), replace = TRUE,
        prob = c(0.167, 0.167, 0.167, 0.167, 0.167, 0.167)
      ),
      Month %in% 1:8 ~ sample(
        c("GMZ1N", "GMZ1S", "GMZ2", "GMZ3", "GMZ4", "GMZ5"),
        size = nrow(.), replace = TRUE,
        prob = c(0.167, 0.167, 0.167, 0.167, 0.167, 0.167)
      ),
      TRUE ~ NA_character_
    )
  )


# test for a little dispersal
table(unlist(LWF_final_well_mixed$Stock))




# visualize stock assignment distribution 
ggplot(LWF_final_well_mixed, aes(x = Stock),colo) +
  geom_bar() +
  theme_minimal() +
  xlab("Source") +
  ylab("Percent of Individuals Caught in WFM-06") 



### POPs filter 
POpairs_well_mixed <- LWF_final_well_mixed %>%
  rename_with(~ paste0(.x, "_1"), everything()) %>% # add the columns fro each indiv
  cross_join(LWF_final_well_mixed %>% rename_with(~ paste0(.x, "_2"), everything())) %>% # remove self-pairings
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
  mutate(StockPair = paste(Stock_1, Stock_2, sep = "_")) %>% 
  mutate(nyears = Cohort_2 - Cohort_1)




# count numnber of yes or no to get a better idea of how they are changing
POpairs_well_mixed %>% 
  dplyr::count(HSPPOP_candidate)


POPonly_well_mixed <- POpairs_well_mixed %>% 
  filter(HSPPOP_candidate == 1) %>% 
  mutate(RObase=ifelse(Sex_1 == "U",2,1)) %>% 
  mutate(Stock_1 = unlist(Stock_1),
         Stock_2 = unlist(Stock_2),
         StockWeight = ifelse(Stock_1 == Stock_2, 1, 0.1)) # likllihoods for if they are from the same place vs not


### HSPs and combine ###

HSpairs_well_mixed <- LWF_final_well_mixed %>%
  rename_with(~ paste0(.x, "_1"), everything()) %>% # sort into col for each indiv
  cross_join(LWF_final_well_mixed %>% rename_with(~ paste0(.x, "_2"), everything())) %>% # remove self-pairings
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
  mutate(StockPair = paste(Stock_1, Stock_2, sep = "_"))%>% 
  mutate(nyears = Cohort_2 - Cohort_1)


# count numnber of yes or no to get a better idea of how they are changing
HSpairs_well_mixed %>% 
  dplyr::count(HSPPOP_candidate)

HSPonly_well_mixed <- HSpairs_well_mixed %>% 
  filter(HSPPOP_candidate == 1) %>% 
  mutate(RObase=4) %>% 
  mutate(Stock_1 = unlist(Stock_1),
         Stock_2 = unlist(Stock_2),
         StockWeight = ifelse(Stock_1 == Stock_2, 1, 0.1)) # liklihoods for if they are from the same place vs not

# join the two tables
library(plyr)
kinshipsSexStock_well_mixed <- rbind(data.frame=POPonly_well_mixed,
                          data.frame=HSPonly_well_mixed)  %>% 
  mutate(
    sex_num = case_when(
      Sex_1 == "U" ~ NA_real_, # unknown to na
      Sex_1 == "M" ~ 1, # male is 1
      Sex_1 == "F" ~ 0 # female is 0
    )
  )




# Real POP/HSP - here assign randomly based on probability 

set.seed(777)

# parameters for decay rate of pair likelihood with age
HS_true_well_mixed <- 0.0001 # edit this number to change the total percentage of pairs that are true
HS_beta_well_mixed <- 0.3 # decay rate of half siblings based on survival of shared parent 
PO_true_well_mixed <- 0.001 # edit this number to change the total percentage of pairs that are true



# assign true pairs systematically 
kinshipsSex_obs_well_mixed <- kinshipsSexStock_well_mixed %>% 
  mutate(trueProb = case_when(
    RObase == 1 ~ PO_true_well_mixed,
    RObase == 2 ~ PO_true_well_mixed,
    RObase == 4 ~ HS_true_well_mixed*exp(-HS_beta_well_mixed*abs(AgeDif)),
    TRUE ~ 0),
    TruePairs = rbinom(nrow(.), size = 1, prob = trueProb))





# extract years (count to give to be i in JAGS)
Cohort_years_well_mixed <- kinshipsSex_obs_well_mixed %>% 
  group_by(Cohort_2) %>% 
  dplyr::summarise(
    Pair_viable_count = sum(HSPPOP_candidate > 0, na.rm = TRUE)
  ) %>% 
  filter(Pair_viable_count > 0) %>%     # keep only cohorts with >0
  arrange(Cohort_2)
years_well_mixed <- nrow(Cohort_years_well_mixed)


# group together by cohort year and age dif (differing probabilities)
Cohort_years_cohort_age_well_mixed <- kinshipsSex_obs_well_mixed %>%
  group_by(Cohort_2, AgeDif) %>%
  dplyr::summarise(
    Pair_viable_count_per_agedif = sum(HSPPOP_candidate > 0, na.rm = TRUE))


# add the counts to each observation
kinshipsSex_obs_well_mixed <- kinshipsSex_obs_well_mixed %>%
  left_join(Cohort_years_cohort_age_well_mixed, by = c("Cohort_2", "AgeDif"))





# Count the ones with known sex (U = 0)
# Count females (F = 1)
kinshipsSex_obs_well_mixed <- kinshipsSex_obs_well_mixed %>%
  mutate(
    knownSex = as.numeric(!is.na(sex_num)),
    IsFemale = as.numeric(!is.na(sex_num) & (sex_num == 0))  # if sex_num encoded 0 = female
  )






# year index to link kinship to years (just the number to match it, not the actual year)
year_index_well_mixed <- factor(kinshipsSex_obs_well_mixed$Cohort_2, levels = Cohort_years$Cohort_2)
kinshipsSex_obs_well_mixed <- kinshipsSex_obs_well_mixed %>%
  mutate(
    year_index = as.integer(year_index_well_mixed)
  )

#### Now collapse based on cohort 2, age difference, RObase, if sex is known, and female or male 
group_vars_well_mixed <- c("Cohort_2", "AgeDif", "knownSex", "IsFemale", "RObase","StockWeight")

collapsed_well_mixed <- kinshipsSex_obs_well_mixed %>%
  group_by(year_index_well_mixed,across(all_of(group_vars_well_mixed))) %>%
  dplyr::summarise(
    Pair_viable_count_per_agedif = n(),                   
    TruePairs = sum(TruePairs, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  arrange(Cohort_2, AgeDif)



#### Run for single year ----
# use as.numeric if any are giving issues
# make sure all (except years for i loop) are equal to nobs! Check with length first...
well_mixeddata = list( years = years_well_mixed,  # number of cohorts,
                 TruePairs = collapsed_well_mixed$TruePairs,
                 AgeDif = collapsed_well_mixed$AgeDif,
                 RObase = collapsed_well_mixed$RObase,
                 knownSex = as.numeric(collapsed_well_mixed$knownSex),
                 nobs = nrow(collapsed_well_mixed),
                 IsFemale = collapsed_well_mixed$IsFemale,
                 Pair_viable_count_per_agedif = collapsed_well_mixed$Pair_viable_count_per_agedif,
                 year=as.numeric(collapsed_well_mixed$year_index_well_mixed),
                 StockWeight = collapsed_well_mixed$StockWeight
)

# Initial values
inits_well_mixed = function() {
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


Out_well_mixed = jagsUI::jags(well_mixeddata, inits_well_mixed, params, "HSPPOP.sex.stocks.singleY.jags", n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = T, verbose = T)
# worked initially will a few funky Rhat values for some Nhat years where Neff was low. 

# if issues occur, run un parallelized
Out_well_mixed <- jagsUI::jags(
  data = well_mixeddata,
  inits = inits_well_mixed,
  parameters.to.save = params,
  model.file = "HSPPOP.sex.stocks.singleY.jags",
  n.burnin = nburn,
  n.chains = nchains,
  n.iter = niter,
  parallel = FALSE,
  verbose = TRUE
)













# Extract posterior summary as a data frame
nhat_well_mixed <- as.data.frame(Out_well_mixed$summary)
nhat_well_mixed <- nhat_well_mixed[grep("^Nadult", rownames(nhat_well_mixed)), ]
# Add a Year index (1, 2, 3, ...)
nhat_well_mixed$Year.index <- 1:nrow(nhat)
nhat_well_mixed$Year <-1993:2020

# Select the useful columns
nhat_well_mixed <- nhat_well_mixed %>%
  dplyr::select(Year, mean, `2.5%`, `50%`, `97.5%`)


library(ggplot2)

ggplot(nhat_well_mixed, aes(x = Year, y = mean)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue") +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, fill = "lightblue") +
  labs(
    x = "Year (Cohort index)",
    y = expression(hat(N)[adult]),
    title = "Posterior estimates of adult abundance (N-hat) in
  in a well mixed fishery"
  ) +
  theme_minimal(base_size = 14)

















#### Now test with only the fish in GMZ3



kinshipsGMZ3_well_mixed <- kinshipsSexStock_well_mixed



# Real POP/HSP - here assign randomly based on probability 

set.seed(777)

# parameters for decay rate of pair likelihood with age
HS_trueGMZ3_well_mixed <- 0.0001 # edit this number to change the total percentage of pairs that are true
HS_betaGMZ3_well_mixed <- 0.3 # decay rate of half siblings based on survival of shared parent 
PO_trueGMZ3_well_mixed <- 0.001 # edit this number to change the total percentage of pairs that are true



# assign true pairs systematically 
kinshipsSex_obsGMZ3_well_mixed <- kinshipsGMZ3_well_mixed %>% 
  mutate(trueProb = case_when(
    RObase == 1 ~ PO_true,
    RObase == 2 ~ PO_true,
    RObase == 4 ~ HS_true*exp(-HS_beta*abs(AgeDif)),
    TRUE ~ 0),
    TruePairs = rbinom(nrow(.), size = 1, prob = trueProb))


kinshipsSex_obsGMZ3_well_mixed <- kinshipsSex_obsGMZ3_well_mixed %>% 
  filter(StockPair == "GMZ3_GMZ3")



# extract years (count to give to be i in JAGS)
Cohort_yearsGMZ3_well_mixed <- kinshipsSex_obsGMZ3_well_mixed %>% 
  group_by(Cohort_2) %>% 
  dplyr::summarise(
    Pair_viable_count = sum(HSPPOP_candidate > 0, na.rm = TRUE)
  ) %>% 
  filter(Pair_viable_count > 0) %>%     # keep only cohorts with >0
  arrange(Cohort_2)
yearsGMZ3_well_mixed <- nrow(Cohort_yearsGMZ3_well_mixed)


# group together by cohort year and age dif (differing probabilities)
Cohort_years_cohort_ageGMZ3_well_mixed <- kinshipsSex_obsGMZ3_well_mixed %>%
  group_by(Cohort_2, AgeDif) %>%
  dplyr::summarise(
    Pair_viable_count_per_agedif = sum(HSPPOP_candidate > 0, na.rm = TRUE))


# add the counts to each observation
kinshipsSex_obsGMZ3_well_mixed <- kinshipsSex_obsGMZ3_well_mixed %>%
  left_join(Cohort_years_cohort_ageGMZ3_well_mixed, by = c("Cohort_2", "AgeDif"))





# Count the ones with known sex (U = 0)
# Count females (F = 1)
kinshipsSex_obsGMZ3_well_mixed <- kinshipsSex_obsGMZ3_well_mixed %>%
  mutate(
    knownSex = as.numeric(!is.na(sex_num)),
    IsFemale = as.numeric(!is.na(sex_num) & (sex_num == 0))  # if sex_num encoded 0 = female
  )






# year index to link kinship to years (just the number to match it, not the actual year)
year_indexGMZ3_well_mixed <- factor(kinshipsSex_obsGMZ3_well_mixed$Cohort_2, levels = Cohort_yearsGMZ3_well_mixed$Cohort_2)
kinshipsSex_obsGMZ3_well_mixed <- kinshipsSex_obsGMZ3_well_mixed %>%
  mutate(
    year_index = as.integer(year_indexGMZ3_well_mixed)
  )

#### Now collapse based on cohort 2, age difference, RObase, if sex is known, and female or male 
group_varsGMZ3_well_mixed <- c("Cohort_2", "AgeDif", "knownSex", "IsFemale", "RObase","StockWeight")

collapsedGMZ3_well_mixed <- kinshipsSex_obsGMZ3_well_mixed %>%
  group_by(year_indexGMZ3_well_mixed,across(all_of(group_varsGMZ3_well_mixed))) %>%
  dplyr::summarise(
    Pair_viable_count_per_agedif = n(),                   
    TruePairs = sum(TruePairs, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  arrange(Cohort_2, AgeDif)



#### Run for single year ----
# use as.numeric if any are giving issues
# make sure all (except years for i loop) are equal to nobs! Check with length first...
GMZ3_well_mixeddata = list( years = yearsGMZ3_well_mixed,  # number of cohorts,
                 TruePairs = collapsedGMZ3_well_mixed$TruePairs,
                 AgeDif = collapsedGMZ3_well_mixed$AgeDif,
                 RObase = collapsedGMZ3_well_mixed$RObase,
                 knownSex = as.numeric(collapsedGMZ3_well_mixed$knownSex),
                 nobs = nrow(collapsedGMZ3_well_mixed),
                 IsFemale = collapsedGMZ3_well_mixed$IsFemale,
                 Pair_viable_count_per_agedif = collapsedGMZ3_well_mixed$Pair_viable_count_per_agedif,
                 year=as.numeric(collapsedGMZ3_well_mixed$year_indexGMZ3_well_mixed),
                 StockWeight = collapsedGMZ3_well_mixed$StockWeight
)

# Initial values
initsGMZ3_well_mixed = function() {
  list(
    mu = runif(1, 1, 100000),
    sd = runif(1, 1, 100000),
    propM = runif(1,0.01, 0.99),
    Nadult = rep(10000, yearsGMZ3_well_mixed),   # vector of length POP_years
    surv = (rep(.7,yearsGMZ3_well_mixed))
  )
}

# Parameters to follow
params = c("Nadult","surv","propM")

nburn <- 3000
nchains <- 3
niter <- 5000
n.cores = 3


OutGMZ3_well_mixed = jagsUI::jags(GMZ3_well_mixeddata, initsGMZ3_well_mixed, params, "HSPPOP.sex.stocks.singleY.jags", n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = T, verbose = T)
# worked initially will a few funky Rhat values for some Nhat years where Neff was low. 

# if issues occur, run un parallelized
OutGMZ3_well_mixed <- jagsUI::jags(
  data = GMZ3_well_mixeddata,
  inits = initsGMZ3_well_mixed,
  parameters.to.save = params,
  model.file = "HSPPOP.sex.stocks.singleY.jags",
  n.burnin = nburn,
  n.chains = nchains,
  n.iter = niter,
  parallel = FALSE,
  verbose = TRUE
)













# Extract posterior summary as a data frame
nhatGMZ3_well_mixed <- as.data.frame(OutGMZ3_well_mixed$summary)
nhatGMZ3_well_mixed <- nhatGMZ3_well_mixed[grep("^Nadult", rownames(nhatGMZ3_well_mixed)), ]
# Add a Year index (1, 2, 3, ...)
nhatGMZ3_well_mixed$Year.index <- 1:nrow(nhatGMZ3_well_mixed)
nhatGMZ3_well_mixed$Year <-1994:2019

# Select the useful columns
nhatGMZ3_well_mixed <- nhatGMZ3_well_mixed %>%
  dplyr::select(Year, mean, `2.5%`, `50%`, `97.5%`)


library(ggplot2)

ggplot(nhatGMZ3_well_mixed, aes(x = Year, y = mean)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue") +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, fill = "lightblue") +
  labs(
    x = "Year (Cohort index)",
    y = expression(hat(N)[adult]),
    title = "Posterior estimates of adult abundance (N-hat) in GMZ3
  in a well mixed fishery"
  ) +
  theme_minimal(base_size = 14)

