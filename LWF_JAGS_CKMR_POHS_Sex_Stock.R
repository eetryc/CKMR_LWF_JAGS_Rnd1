#### Set up data frame for analysis with both HSPs and POPs ####

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
  select(SampleID, Cohort, SampleYear, SampleDate, Sex, Month) %>%
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
        prob = c(0.05, 0.00, 0.4, 0.25, 0.15, 0.15)
      ),
      TRUE ~ NA_character_
    )
  )


# test for a little dispersal
table(unlist(LWF_final_GMZ_mixed$Stock))

# visualize stock assignment distribution 
ggplot(LWF_final_GMZ_mixed, aes(x = Stock),colo) +
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
                       data.frame=HSPonly) 
































#### JAGS model ####
cat("model{
  # survival
  surv ~ dunif(0.5, 0.95) 

  mu ~ dunif(1,100000)
  sd ~ dunif(1,100000)
  
  # Adult numbers
	for(i in 1:years) {
	
	  Nadult[i] ~ dnorm(mu, 1/(sd^2)) T(0, 1e+09)

		TruePairs[i] ~ dbinom((RObase[i]*(surv^AgeDif[i])*StockWeight[i])/(Nadult[i]), Pair_viable_count[i])

	}


}", file = "HSPPOP.stocks.singleY.jags")




#### Data to define for JAGS model ####

# Real POP/HSP - here assign randomly based on probability 

set.seed(777)
kinshipsStock <- kinshipsStock %>%
  mutate(TruePairs = case_when(
    RObase == 2 ~ rbinom(nrow(.), size = 1, prob = 0.001),
    RObase == 4 ~ rbinom(nrow(.), size = 1, prob = 0.0001),
    TRUE ~ 0
  ))

Cohort_years <- kinshipsStock %>% 
  group_by(Cohort_2) %>% 
  dplyr::summarise(
    Pair_viable_count = sum(HSPPOP_candidate > 0, na.rm = TRUE)
  ) %>% 
  filter(Pair_viable_count > 0) %>%     # keep only cohorts with >0
  arrange(Cohort_2)

years <- nrow(Cohort_years)

# baseline probability (4 for unsexed HSPs, 2 for unsexed POPs)

#### Run for single year ----
data = list( years = years,  # number of cohorts,
             TruePairs = kinshipsStock$TruePairs,
             AgeDif=kinshipsStock$AgeDif,
             Pair_viable_count = Cohort_years$Pair_viable_count,
             RObase = kinshipsStock$RObase,
             StockWeight = kinshipsStock$StockWeight)

# convert to numeric if not already (try this first if node issues occur)
data <- list(
  years = as.numeric(years),
  TruePairs = as.numeric(kinshipsStock$TruePairs),
  AgeDif = as.numeric(kinshipsStock$AgeDif),
  Pair_viable_count = as.numeric(Cohort_years$Pair_viable_count),
  RObase = as.numeric(kinshipsStock$RObase),
  StockWeight = as.numeric(kinshipsStock$StockWeight)
)




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


Out = jagsUI::jags(data, inits, params, "HSPPOP.stocks.singleY.jags", n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = T, verbose = T)


# if issues occur, run un parallelized
Out <- jagsUI::jags(
  data = data,
  inits = inits,
  parameters.to.save = params,
  model.file = "HSPPOP.stocks.singleY.jags",
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
nhat$Year <- 1:nrow(nhat)

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
  in a well mixed fishery"
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
