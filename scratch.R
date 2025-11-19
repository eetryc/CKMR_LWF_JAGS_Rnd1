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
      (RObase[j] * (surv[year[j]]^ AgeDif[j])) /
      ((Nadult[year[j]] * knownSex[j] * abs(propM - IsFemale[j])) +
      (Nadult[year[j]] * (1 - knownSex[j]))),
      Pair_viable_count_per_agedif[j]
    )
  }	

}", file = "HSPPOP.sex.stocks.singleY.jags")




data = list( years = years,  # number of cohorts,
             TruePairs = collapsed$TruePairs,
             AgeDif = collapsed$AgeDif,
             RObase = collapsed$RObase,
             knownSex = as.numeric(collapsed$knownSex),
             nobs = nrow(collapsed),
             IsFemale = collapsed$IsFemale,
             Pair_viable_count_per_agedif = collapsed$Pair_viable_count_per_agedif,
             year=collapsed$year_index)

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
    title = "Posterior estimates of adult abundance (N-hat)"  ) +
  theme_minimal(base_size = 14)





















#### Different adult dist
cat("model{

  # proportion male in population (inverse of female, 1-propM)
  propM ~ dbeta(1,1)
  
  
  mu ~ dunif(1,100000)
  sd ~ dunif(1,100000)
  
  # Adult numbers
  for (i in 1:years) {
    r_N[i] <- (mu^2) / (sd^2 - mu)
    p_N[i] <- r_N[i] / (r_N[i] + mu)
   Nadult[i] ~ dnbinom(p_N[i], r_N[i]) T(1, 1e+09)
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

}", file = "HSPPOP.sex.stocks.singleY.dnbinom.jags")




data = list( years = years,  # number of cohorts,
             TruePairs = collapsed$TruePairs,
             AgeDif = collapsed$AgeDif,
             RObase = collapsed$RObase,
             knownSex = as.numeric(collapsed$knownSex),
             nobs = nrow(collapsed),
             IsFemale = collapsed$IsFemale,
             Pair_viable_count_per_agedif = collapsed$Pair_viable_count_per_agedif,
             year=collapsed$year_index)

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


Out.dnbinom = jagsUI::jags(data, inits, params, "HSPPOP.sex.stocks.singleY.dnbinom.jags", n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = T, verbose = T)



# Extract posterior summary as a data frame
nhat.dnbinom <- as.data.frame(Out.dnbinom$summary)

# Keep only the Nadult rows
nhat.dnbinom <- nhat.dnbinom[grep("^Nadult", rownames(nhat)), ]

# Add a Year index (1, 2, 3, ...)
nhat.dnbinom$Year.index <- 1:nrow(nhat.dnbinom)
nhat.dnbinom$Year <-1993:2020
# Select the useful columns
nhat.dnbinom_df <- nhat.dnbinom %>%
  dplyr::select(Year, mean, `2.5%`, `50%`, `97.5%`)


library(ggplot2)

ggplot(nhat.dnbinom_df, aes(x = Year, y = mean)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue") +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, fill = "lightblue") +
  labs(
    x = "Year (Cohort index)",
    y = expression(hat(N)[adult]),
    title = "Posterior estimates of adult abundance (N-hat)"  ) +
  theme_minimal(base_size = 14)



# Extract posterior summary as a data frame
Rhatdnbinom <- as.data.frame(Out.dnbinom$summary)

# Keep only the Rhat rows
Rhatdnbinom <- Rhatdnbinom[grep("^Nadult", rownames(Rhatdnbinom)), ]

# Add a Year index (1, 2, 3, ...)
Rhatdnbinom$Year <- 1:nrow(Rhatdnbinom)

# Select the useful columns
Rhat_dfdnbinom <- Rhatdnbinom %>%
  dplyr::select(Year, mean, `2.5%`, `50%`, `97.5%`,Rhat)
ggplot(Rhat_dfdnbinom, aes(x = Year, y = Rhat)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue") +
  labs(
    x = "Year (Cohort index)",
    y = "Convergence (R-hat)",
    title = "Posterior estimates of Convergence (R-hat) 
  in a well mixed fishery"
  ) +
  theme_minimal(base_size = 14) +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "darkred")









#### poisson distribution

cat("model{

  # proportion male in population (inverse of female, 1-propM)
  propM ~ dbeta(1,1)
  
  
  mu ~ dunif(1,100000)
  sd ~ dunif(1,100000)
  
  # Adult numbers
  for (i in 1:years) {

   Nadult[i] ~ dpois(mu)
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

}", file = "HSPPOP.sex.stocks.singleY.dpois.jags")




data = list( years = years,  # number of cohorts,
             TruePairs = collapsed$TruePairs,
             AgeDif = collapsed$AgeDif,
             RObase = collapsed$RObase,
             knownSex = as.numeric(collapsed$knownSex),
             nobs = nrow(collapsed),
             IsFemale = collapsed$IsFemale,
             Pair_viable_count_per_agedif = collapsed$Pair_viable_count_per_agedif,
             year=collapsed$year_index)

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


Out.dpois = jagsUI::jags(data, inits, params, "HSPPOP.sex.stocks.singleY.dpois.jags", n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = T, verbose = T)



# Extract posterior summary as a data frame
nhat.dpois <- as.data.frame(Out.dpois$summary)

# Keep only the Nadult rows
nhat.dpois <- nhat.dpois[grep("^Nadult", rownames(nhat.dpois)), ]

# Add a Year index (1, 2, 3, ...)
nhat.dpois$Year.index <- 1:nrow(nhat.dpois)
nhat.dpois$Year <-1993:2020
# Select the useful columns
nhat_df.dpois <- nhat.dpois %>%
  dplyr::select(Year, mean, `2.5%`, `50%`, `97.5%`)


library(ggplot2)

ggplot(nhat_df.dpois, aes(x = Year, y = mean)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue") +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, fill = "lightblue") +
  labs(
    x = "Year (Cohort index)",
    y = expression(hat(N)[adult]),
    title = "Posterior estimates of adult abundance (N-hat)"  ) +
  theme_minimal(base_size = 14)


# Extract posterior summary as a data frame
Rhatdpois <- as.data.frame(Out.dpois$summary)

# Keep only the Rhat rows
Rhatdpois <- Rhatdpois[grep("^Nadult", rownames(Rhatdnbinom)), ]

# Add a Year index (1, 2, 3, ...)
Rhatdpois$Year <- 1:nrow(Rhatdpois)

# Select the useful columns
Rhat_dfpois <- Rhatdpois %>%
  dplyr::select(Year, mean, `2.5%`, `50%`, `97.5%`,Rhat)

ggplot(Rhat_dfpois, aes(x = Year, y = Rhat)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue") +
  labs(
    x = "Year (Cohort index)",
    y = "Convergence (R-hat)",
    title = "Posterior estimates of Convergence (R-hat) 
  in a well mixed fishery"
  ) +
  theme_minimal(base_size = 14) +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "darkred")









kinshipsSex_obspair <- kinshipsSex_obs %>% 
  mutate(matched = ifelse(Stock_1 == Stock_2, 1, 0.1))



ggplot(kinshipsSex_obspair,aes(x=StockPair)) +
  geom_bar(aes(fill=matched)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


























GMZ3data = list( years = yearsGMZ3,  # number of cohorts,
                 TruePairs = collapsedGMZ3$TruePairs,
                 AgeDif = collapsedGMZ3$AgeDif,
                 RObase = collapsedGMZ3$RObase,
                 knownSex = as.numeric(collapsedGMZ3$knownSex),
                 nobs = nrow(collapsedGMZ3),
                 IsFemale = collapsedGMZ3$IsFemale,
                 Pair_viable_count_per_agedif = collapsedGMZ3$Pair_viable_count_per_agedif,
                 year=as.numeric(collapsedGMZ3$year_indexGMZ3))

# Initial values
initsGMZ3 = function() {
  list(
    mu = runif(1, 1, 25000),
    sd = runif(1, 1, 25000),
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


#### Graph

# Extract posterior summary as a data frame
nhatGMZ3 <- as.data.frame(OutGMZ3$summary)
nhatGMZ3 <- nhatGMZ3[grep("^Nadult", rownames(nhatGMZ3)), ]
# Add a Year index (1, 2, 3, ...)
nhatGMZ3$Year.index <- 1:nrow(nhat)
nhatGMZ3$Year <-1993:2020

# Select the useful columns
nhatGMZ3 <- nhatGMZ3 %>%
  dplyr::select(Year, mean, `2.5%`, `50%`, `97.5%`)



ggplot(nhatGMZ3, aes(x = Year, y = mean)) +
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




















#### Lumping years with limited observations ####

group_vars <- c("Cohort_2", "AgeDif", "knownSex", "IsFemale", "RObase","StockWeight")

collapsed <- kinshipsSex_obs %>%
  group_by(year_index,across(all_of(group_vars))) %>%
  dplyr::summarise(
    Pair_viable_count_per_agedif = n(),                   
    TruePairs = sum(TruePairs, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  arrange(Cohort_2, AgeDif)


pair_counts_collapsed <- kinshipsSex_obs %>%
  group_by(year_index) %>%
  dplyr::summarise(
    Pair_viable_count_per_agedif = n(),                   
    TruePairs = sum(TruePairs, na.rm=TRUE),
    .groups = "drop"
  ) 
pair_counts_collapsed %>% 
  ggplot(aes(x = year_index, y = All_pairs)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1000, linetype = "dashed", color = "darkred")


pair_counts_collapsed %>% 
  ggplot(aes(x = year_index, y = TruePairs)) +
  geom_line() +
  geom_point() 
  # geom_hline(yintercept = 1000, linetype = "dashed", color = "darkred")


# group collapsed values

pair_counts_collapsed <- pair_counts_collapsed %>% 
  mutate(group_index = case_when(
    Pair_viable_count_per_agedif < 5000 ~ 1,
    TRUE ~ year_index
  )) %>% 
  group_by(group_index,year_index,TruePairs) %>% 
  dplyr::summarise(All_pairs = sum(Pair_viable_count_per_agedif))


collapsed_grouped <- left_join(collapsed,pair_counts_collapsed, by = "year_index")

collapsed_grouped <- collapsed %>%
  mutate(group_index = if_else(
    Pair_viable_count_per_agedif < 5000, 
    1L,
    as.integer(year_index)
  ))

years_grouped = nrow(pair_counts_collapsed)

data_grouped = list( years = years_grouped,  # number of cohorts,
             TruePairs = collapsed_grouped$TruePairs,
             AgeDif = collapsed_grouped$AgeDif,
             RObase = collapsed_grouped$RObase,
             knownSex = as.numeric(collapsed_grouped$knownSex),
             nobs = nrow(collapsed_grouped),
             IsFemale = collapsed_grouped$IsFemale,
             Pair_viable_count_per_agedif = collapsed_grouped$Pair_viable_count_per_agedif,
             year=collapsed_grouped$group_index,
             StockWeight = collapsed_grouped$StockWeight)


inits_grouped = function() {
  list(
    mu = runif(1, 1, 100000),
    sd = runif(1, 1, 100000),
    propM = runif(1,0.01, 0.99),
    Nadult = rep(10000, years_grouped),   # vector of length POP_years
    surv = (rep(.7,years_grouped))
  )
}
Out_grouped = jagsUI::jags(data_grouped, inits_grouped, params, "HSPPOP.sex.stocks.singleY.jags", n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = T, verbose = T)





