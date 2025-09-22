#Load script-specific libraries
library(CKMRsim)
library(rjags)
library(runjags)
library(jagsUI)
library(postpack)

#Load the combined data file ----
combdata = read.csv("path to data file")

## Run CKMRsim ----
unique(combdata$SampGroup)
genos.pop = combdata %>% filter(SampGroup %in% c("LC-Cap", "KIN-Off", "KIN-Par", "LC-Spawn", "LC-Matt", "LC-LIP") & !is.na(Allele1))

#Reformat for ckmrsim
genos.freq = rbind(data.frame(Indiv = genos.pop$NovoID, Locus = genos.pop$Locus,
                              gene_copy = 1, Allele = genos.pop$Allele1, stringsAsFactors = F),
                   data.frame(Indiv = genos.pop$NovoID, Locus = genos.pop$Locus,
                              gene_copy = 2, Allele = genos.pop$Allele2, stringsAsFactors = F))

genos.freq.pop = genos.freq

#Calculate allele frequencies
loci = unique(genos.freq.pop$Locus)

alle_freqs <- genos.freq.pop %>%
  filter(!is.na(Allele)) %>%
  group_by(Locus, Allele) %>%
  summarize(n = n()) %>%
  group_by(Locus) %>%
  mutate(Freq = n / sum(n),
         Chrom = "Unk",
         Pos = as.integer(factor(Locus, levels = loci))) %>%
  ungroup() %>%
  select(Chrom, Pos, Locus, Allele, Freq) %>%
  arrange(Pos, desc(Freq)) %>%
  mutate(AlleIdx = NA,
         LocIdx = NA) %>%
  filter(!is.na(Allele))

badalleles = alle_freqs %>% filter(Freq < 0.01)

genos.pop = genos.pop %>% filter(!(paste0(Locus, Allele1) %in% paste0(badalleles$Locus, badalleles$Allele) |
                                     paste0(Locus, Allele2) %in% paste0(badalleles$Locus, badalleles$Allele)))

#Reformat for ckmrsim
genos.freq = rbind(data.frame(Indiv = genos.pop$NovoID, Locus = genos.pop$Locus,
                              gene_copy = 1, Allele = genos.pop$Allele1, stringsAsFactors = F),
                   data.frame(Indiv = genos.pop$NovoID, Locus = genos.pop$Locus,
                              gene_copy = 2, Allele = genos.pop$Allele2, stringsAsFactors = F))

genos.freq.pop = genos.freq

#Calculate allele frequencies
loci = unique(genos.freq.pop$Locus)

alle_freqs <- genos.freq.pop %>%
  filter(!is.na(Allele)) %>%
  group_by(Locus, Allele) %>%
  summarize(n = n()) %>%
  group_by(Locus) %>%
  mutate(Freq = n / sum(n),
         Chrom = "Unk",
         Pos = as.integer(factor(Locus, levels = loci))) %>%
  ungroup() %>%
  select(Chrom, Pos, Locus, Allele, Freq) %>%
  arrange(Pos, desc(Freq)) %>%
  mutate(AlleIdx = NA,
         LocIdx = NA) %>%
  filter(!is.na(Allele))


master.error = 0.2 #Can play with this, but 0.2 is default

sim_ckmr = create_ckmr(reindex_markers(as_tibble(alle_freqs)),
                       ge_mod_assumed = ge_model_TGIE,
                       ge_mod_true = ge_model_TGIE,
                       ge_mod_assumed_pars_list = list(epsilon = master.error),
                       ge_mod_true_pars_list = list(epsilon = master.error)) #Format markers and create object

#Quick check for duplicates...
dups = (genos.freq.pop %>%
          dplyr::group_by(Indiv, Locus, gene_copy) %>%
          dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
          dplyr::filter(n > 1L))

unique(dups$Indiv)
genos.freq.pop = genos.freq.pop %>% filter(!Indiv %in% dups$Indiv)

# Calculate kinship log-likelihoods
ratios.PO = pairwise_kin_logl_ratios(
  genos.freq.pop,
  genos.freq.pop,
  sim_ckmr,
  numer = "PO",
  denom = "U",
  keep_top = NULL,
  num_cores = 30) %>%
  arrange(desc(logl_ratio))

ratios.HS = pairwise_kin_logl_ratios(
  genos.freq.pop,
  genos.freq.pop,
  sim_ckmr,
  numer = "HS",
  denom = "U",
  keep_top = NULL,
  num_cores = 30) %>%
  arrange(desc(logl_ratio))

ratios.FS = pairwise_kin_logl_ratios(
  genos.freq.pop,
  genos.freq.pop,
  sim_ckmr,
  numer = "FS",
  denom = "U",
  keep_top = NULL,
  num_cores = 30) %>%
  arrange(desc(logl_ratio))

#Pull out the known-kinship sample
ratios.PO.kin = ratios.PO %>% filter(D1_indiv %in% (combdata %>% filter(grepl("KIN", SampGroup)))$NovoID & 
                                       D2_indiv %in% (combdata %>% filter(grepl("KIN", SampGroup)))$NovoID)
ratios.HS.kin = ratios.HS %>% filter(D1_indiv %in% (combdata %>% filter(grepl("KIN", SampGroup)))$NovoID & 
                                       D2_indiv %in% (combdata %>% filter(grepl("KIN", SampGroup)))$NovoID)
ratios.FS.kin = ratios.FS %>% filter(D1_indiv %in% (combdata %>% filter(grepl("KIN", SampGroup)))$NovoID & 
                                       D2_indiv %in% (combdata %>% filter(grepl("KIN", SampGroup)))$NovoID)
ratios.kin = rbind(cbind(data.frame(Comp = "PO"),ratios.PO.kin),
                   cbind(data.frame(Comp = "HS"),ratios.HS.kin),
                   cbind(data.frame(Comp = "FS"),ratios.FS.kin))

#Assign mothers and fathers to each combo
ratios.kin$Mom1 = paste0("KK-F", substring(ratios.kin$D1_indiv,4,4))
ratios.kin$Dad1 = paste0("KK-M", substring(ratios.kin$D1_indiv,5,5))
ratios.kin$Mom2 = paste0("KK-F", substring(ratios.kin$D2_indiv,4,4))
ratios.kin$Dad2 = paste0("KK-M", substring(ratios.kin$D2_indiv,5,5))

ratios.kin[which(ratios.kin$D1_indiv %in% c(ratios.kin$Mom2, ratios.kin$Dad2)),c("Mom1", "Dad1")] = NA
ratios.kin[which(ratios.kin$D2_indiv %in% c(ratios.kin$Mom1, ratios.kin$Dad1)),c("Mom2", "Dad2")] = NA

ratios.kin$TrueRelat = "U"

ratios.kin[which(ratios.kin$D1_indiv == ratios.kin$Mom2 | ratios.kin$D1_indiv == ratios.kin$Dad2 |
                   ratios.kin$D2_indiv == ratios.kin$Mom1 | ratios.kin$D2_indiv == ratios.kin$Dad1),]$TrueRelat = "PO"
ratios.kin[which(ratios.kin$Mom1 == ratios.kin$Mom2 | ratios.kin$Dad1 == ratios.kin$Dad2),]$TrueRelat = "HS"
ratios.kin[which(ratios.kin$Mom1 == ratios.kin$Mom2 & ratios.kin$Dad1 == ratios.kin$Dad2),]$TrueRelat = "FS"

plot.KK.PO = ggplot(ratios.kin %>% filter(Comp == "PO"))+
  geom_density(aes(x = logl_ratio, fill = TrueRelat), alpha = 0.2)+
  scale_x_continuous(breaks = seq(-100,100, by = 10))+
  ggtitle("Known Kin - PO ratios")

plot.KK.HS = ggplot(ratios.kin %>% filter(Comp == "HS"))+
  geom_density(aes(x = logl_ratio, fill = TrueRelat), alpha = 0.2)+
  scale_x_continuous(breaks = seq(-100,100, by = 10))+
  ggtitle("Known Kin - HS ratios")

plot.KK.FS = ggplot(ratios.kin %>% filter(Comp == "FS"))+
  geom_density(aes(x = logl_ratio, fill = TrueRelat), alpha = 0.2)+
  scale_x_continuous(breaks = seq(-100,100, by = 10))+
  ggtitle("Known Kin - FS ratios")

plot.KK.PO + plot.KK.HS + plot.KK.FS + plot_layout(nrow = 3)

#Use plots to decide cutoff thresholds
PO.thresh = "Parent-offspring threshold"
HS.thresh = "Half-sibling threshold"
FS.thresh = "Full-sibling threshold"

(plot.KK.PO + geom_vline(aes(xintercept = PO.thresh), lty = 3)) + 
  (plot.KK.HS + geom_vline(aes(xintercept = HS.thresh), lty = 3)) + 
  (plot.KK.FS + geom_vline(aes(xintercept = FS.thresh), lty = 3)) + 
     plot_layout(nrow = 3)


#Get false-negative rates for later adjustment
PO.FNrate = dim((ratios.kin %>% filter(Comp == "PO", TrueRelat == "PO", logl_ratio < PO.thresh)))[1] /
  dim((ratios.kin %>% filter(Comp == "PO", TrueRelat == "PO")))[1]

HS.FNrate = dim((ratios.kin %>% filter(Comp == "HS", TrueRelat == "HS",  logl_ratio < HS.thresh)))[1] /
  dim((ratios.kin %>% filter(Comp == "HS", TrueRelat == "HS")))[1]

cat(paste0("\nThresholds result in false-negative rates of ", round(PO.FNrate,3), " for POPs and ",
           round(HS.FNrate,3), " for HSPs\n"))

#Infer relationships

ratios.comb = rbind(cbind(data.frame(Comp = "PO"),ratios.PO),
                    cbind(data.frame(Comp = "HS"),ratios.HS),
                    cbind(data.frame(Comp = "FS"),ratios.HS))

ratios.comb = ratios.comb %>% left_join(ratios.kin)

#Bind HS and PO ratios together
ratios.pred = ratios.PO %>%
  rename(PO_logl = logl_ratio, PO_nloc = num_loc) %>%
  full_join(ratios.HS %>%
              rename(HS_logl = logl_ratio, HS_nloc = num_loc)) %>%
  full_join(ratios.FS %>%
              rename(FS_logl = logl_ratio, FS_nloc = num_loc))

#Remove the known-kin dataset
ratios.pred = ratios.pred %>% filter(!grepl("KK", D1_indiv) & !grepl("KK", D2_indiv))

# If the D2_indiv is older (or if no age, larger) than the D1_indiv, switch them (order matters for potential parent-offspring)
ratios.pred$ParIndiv = "D1"
ratios.pred[which(combdata$TL_mm[match(ratios.pred$D1_indiv, combdata$NovoID)] < 
                    combdata$TL_mm[match(ratios.pred$D2_indiv, combdata$NovoID)] & 
                    (is.na(combdata$Cohort[match(ratios.pred$D1_indiv, combdata$NovoID)]) &
                       is.na(combdata$Cohort[match(ratios.pred$D2_indiv, combdata$NovoID)]))),]$ParIndiv = "D2"
ratios.pred[which(combdata$Cohort[match(ratios.pred$D1_indiv, combdata$NovoID)] > 
                    combdata$Cohort[match(ratios.pred$D2_indiv, combdata$NovoID)]),]$ParIndiv = "D2"

temp.indiv1 = ratios.pred[which(ratios.pred$ParIndiv == "D2"),]$D2_indiv
temp.indiv2 = ratios.pred[which(ratios.pred$ParIndiv == "D2"),]$D1_indiv
ratios.pred[which(ratios.pred$ParIndiv == "D2"),c("D1_indiv", "D2_indiv")] = cbind(temp.indiv1, temp.indiv2)

ratios.pred = ratios.pred %>% select(!ParIndiv)

#Append additional info
ratios.pred$Year_2 = combdata$Year[match(ratios.pred$D2_indiv, combdata$NovoID)]
ratios.pred$TL_2 = combdata$TL_mm[match(ratios.pred$D2_indiv, combdata$NovoID)]
ratios.pred$Origin_2 = combdata$StockOrigin[match(ratios.pred$D2_indiv, combdata$NovoID)]
ratios.pred$Region_2 = combdata$Region[match(ratios.pred$D2_indiv, combdata$NovoID)]
ratios.pred$PopAssign_2 = combdata$AssignPop[match(ratios.pred$D2_indiv, combdata$NovoID)]
ratios.pred$Cohort_2 = combdata$Cohort[match(ratios.pred$D2_indiv, combdata$NovoID)]
ratios.pred$Sex_2 = combdata$Sex[match(ratios.pred$D2_indiv, combdata$NovoID)]

ratios.pred$Year_1 = combdata$Year[match(ratios.pred$D1_indiv, combdata$NovoID)]
ratios.pred$TL_1 = combdata$TL_mm[match(ratios.pred$D1_indiv, combdata$NovoID)]
ratios.pred$Origin_1 = combdata$StockOrigin[match(ratios.pred$D1_indiv, combdata$NovoID)]
ratios.pred$Region_1 = combdata$Region[match(ratios.pred$D1_indiv, combdata$NovoID)]
ratios.pred$PopAssign_1 = combdata$AssignPop[match(ratios.pred$D1_indiv, combdata$NovoID)]
ratios.pred$Cohort_1 = combdata$Cohort[match(ratios.pred$D1_indiv, combdata$NovoID)]
ratios.pred$Sex_1 = combdata$Sex[match(ratios.pred$D1_indiv, combdata$NovoID)]

#Put together all the possible potential combinations
allcombs = expand.grid(ID1 = unique(combdata$NovoID), ID2 = unique(combdata$NovoID))
dim(allcombs)
#Remove self combos
allcombs = allcombs %>% filter(ID1 != ID2)

# Associate age and length data for each fish
allcombs[,c("CapYear1","Cohort1","StockOrigin1", "Length1", "LethalSamp1", "SpawnSurvey1", "AssignPop1")] = 
  combdata[match(allcombs$ID1, combdata$NovoID),c("Year","Cohort","StockOrigin", "TL_mm", "LethalSamp", "SpawnSurvey", "AssignPop")]
allcombs[,c("CapYear2","Cohort2","StockOrigin2", "Length2", "LethalSamp2", "SpawnSurvey2", "AssignPop2")] = 
  combdata[match(allcombs$ID2, combdata$NovoID),c("Year","Cohort","StockOrigin", "TL_mm", "LethalSamp", "SpawnSurvey", "AssignPop")]

allcombs$CapYear1 = as.numeric(allcombs$CapYear1)
allcombs$CapYear2 = as.numeric(allcombs$CapYear2)

# Filter a bit more (some data needs are non-negotiable)
allcombs = allcombs %>% filter(!is.na(CapYear1) & !is.na(CapYear2) & !is.na(StockOrigin1) & !is.na(StockOrigin2))

# Keep only the half of each double comparison where ind. 1 is older/larger 
allcombs$ID1older = (allcombs$Cohort1 < allcombs$Cohort2)
allcombs$ID1larger = allcombs$Length1 > allcombs$Length2

# Deal with same cohort comparisons (cut down to one comparison for each combo)
allcombs$SameCohort = !is.na(allcombs$Cohort1) & (allcombs$Cohort1 == allcombs$Cohort2)
allcombs$SameCohortKeep = NA
allcombs[which(allcombs$SameCohort & allcombs$Length1 > allcombs$Length2),]$SameCohortKeep = TRUE
allcombs[which(allcombs$SameCohort & allcombs$Length1 == allcombs$Length2 & 
                 allcombs$CapYear1 > allcombs$CapYear2),]$SameCohortKeep = TRUE

temp.data = allcombs[which(allcombs$SameCohort & 
                             allcombs$Length1 == allcombs$Length2 & 
                             allcombs$CapYear1 == allcombs$CapYear2),]
temp.checked = NULL
temp.combos = paste(temp.data$ID1,temp.data$ID2)
for(i in seq(length(temp.combos)))
{
  cat(paste0(i,"    \r"))
  if(temp.combos[i] %in% temp.checked)
  {
    temp.data[i,]$SameCohortKeep = FALSE
    next
  }
  temp.data[i,]$SameCohortKeep = TRUE
  temp.checked = c(temp.checked, paste(temp.data$ID2[i],temp.data$ID1[i]))
}

allcombs[which(allcombs$SameCohort & 
                 allcombs$Length1 == allcombs$Length2 & 
                 allcombs$CapYear1 == allcombs$CapYear2),] = temp.data

allcombs = allcombs[which(allcombs$ID1older | (is.na(allcombs$ID1older) & allcombs$ID1larger) | allcombs$SameCohortKeep),]

## Calculate whether each is a valid PO and/or HS pair
ageofmat = 5
# Valid half-sib combos are both wild fish and don't share a cohort year (for independence reasons, could relax)
# Also cut out fish born way before we noted wild recruitment as they're likely just missed clips
allcombs$ValidHS = (allcombs$StockOrigin1 == "wild" & allcombs$StockOrigin2 == "wild" & !is.na(allcombs$Cohort1) & !is.na(allcombs$Cohort2))

# Valid parent-offspring combos are ones where the youngest indiv is wild is the age gap between the two is equal to or great than maturity (assuming spawning the year prior to the "cohort" year)
# Also make sure the potential parent wasn't sampled lethally prior to the potential offspring's birth
allcombs$ValidPO = 
  (allcombs$StockOrigin2 == "wild" & ((allcombs$Cohort2 - 1) - allcombs$Cohort1) >= ageofmat & (allcombs$CapYear1 > allcombs$Cohort2 | !allcombs$LethalSamp1)) |
  ((is.na(allcombs$Cohort1) | is.na(allcombs$Cohort2)) & allcombs$StockOrigin2 == "wild" & abs(allcombs$CapYear2 - allcombs$CapYear1) <= 2 & ((allcombs$Length1 - allcombs$Length2) > 500) & (allcombs$CapYear1 > allcombs$Cohort2 | !allcombs$LethalSamp1))

# If a combo is potentially either relationship type, avoid double-counting and prioritize parent-offspring
allcombs[which(allcombs$ValidPO),]$ValidHS = FALSE

# Filter out the invalid comparisons
kinship.PO = ratios.pred %>% filter((paste(D1_indiv,D2_indiv) %in% paste(allcombs[which(allcombs$ValidPO),]$ID1,allcombs[which(allcombs$ValidPO),]$ID2)) | 
                                      (paste(D1_indiv,D2_indiv) %in% paste(allcombs[which(allcombs$ValidPO),]$ID2,allcombs[which(allcombs$ValidPO),]$ID1)))

kinship.HS = ratios.pred %>% filter((paste(D1_indiv,D2_indiv) %in% paste(allcombs[which(allcombs$ValidHS),]$ID1,allcombs[which(allcombs$ValidHS),]$ID2)) | 
                                      (paste(D1_indiv,D2_indiv) %in% paste(allcombs[which(allcombs$ValidHS),]$ID2,allcombs[which(allcombs$ValidHS),]$ID1)))

#Get true pairs (based on thresholds you set)
true.PO = kinship.PO %>%
  filter(!is.na(Cohort_2)) %>%
  #filter(!D1_indiv %in% (combdata %>% filter(SpawnSurvey))$NovoID)
  mutate(nyears = Cohort_2 - Year_1, isPO = PO_logl >= PO.thresh) %>%
  mutate(nyears = ifelse(nyears < 0, 0, nyears))

true.PO$SpawnSpecific = as.numeric(allcombs$SpawnSurvey1[match(true.PO$D1_indiv, allcombs$ID1)])

true.PO.sum = true.PO %>%
  mutate(ParOrigin = ifelse(Origin_1 == 'wild', 'wild', PopAssign_1)) %>%
  filter(!is.na(ParOrigin)) %>%
  group_by(nyears, Cohort_2, ParOrigin, SpawnSpecific, Sex_1) %>% 
  summarize(ncombs = n(), ntrue = sum(isPO)) %>%
  replace_na(list(ntrue = 0)) %>% arrange(SpawnSpecific, Cohort_2, ParOrigin)

true.HS = kinship.HS %>%
  filter(!is.na(Cohort_1) & !is.na(Cohort_2)) %>%
  filter(Cohort_1 != Cohort_2) %>%
  mutate(nyears = Cohort_2 - Cohort_1, isHS = (HS_logl >= HS.thresh))


#Remove lipid fish (chance they're mistakenly matching with other fish?)
true.HS = true.HS %>% filter(!grepl("LIP", D2_indiv) & !grepl("LIP", D1_indiv))

true.HS.sum = true.HS %>%
  group_by(nyears, Cohort_1) %>%
  summarize(ncombs = n(), ntrue = sum(isHS)) %>%
  replace_na(list(ntrue = 0))


#Adjust for false-negative rates
true.PO.sum$Adj_ntrue = round(true.PO.sum$ntrue / (1 - PO.FNrate),1)
true.HS.sum$Adj_ntrue = round(true.HS.sum$ntrue / (1 - HS.FNrate),1)

as.data.frame(true.PO.sum)
as.data.frame(true.HS.sum)

cat(paste0("\nDetected a total of ",sum(true.PO.sum$ntrue), " POPs (adj: ",sum(true.PO.sum$Adj_ntrue),") and ",
           sum(true.HS.sum$ntrue), " HSPs (adj: ",sum(true.HS.sum$Adj_ntrue),") across all years\n"))

## Run through JAGs model

true.PO.sum %>% 
  group_by(Cohort_2, nyears) %>% 
  summarize(ncombs = sum(ncombs), 
            ntrue = sum(ntrue), 
            Adj_ntrue = sum(Adj_ntrue)) %>%
  mutate(Est = (2 * ncombs) / Adj_ntrue)

sum(true.PO.sum$ntrue)


true.HS.sum %>% 
  group_by(Cohort_1, nyears) %>% 
  summarize(ncombs = sum(ncombs), 
            ntrue = sum(ntrue), 
            Adj_ntrue = sum(Adj_ntrue)) %>%
  mutate(Est = (4 * ncombs)  / (Adj_ntrue * surv^nyears)) %>%
  filter(!is.infinite(Est))

### Find PBT matches for VT feral/domestic experiment

ratios.pred %>% 
  filter((Cohort_2 == 2021 &
            Origin_2 == "stocked" &
            (D1_indiv %in% 
               c(combdata %>% filter(Year == "2020", Site == "Gordon Landing"))$NovoID | 
               grepl("EW-", D1_indiv))) |
           (Cohort_1 == 2021  &
              Origin_1 == "stocked"  &
              (D2_indiv %in% 
                 c(combdata %>% filter(Year == "2020", Site == "Gordon Landing"))$NovoID | 
                 grepl("EW-", D2_indiv)))) %>%
  arrange(desc(PO_logl)) %>%
  mutate(YoungIndiv = ifelse((is.na(Cohort_2) | Cohort_2 != 2021), D1_indiv, D2_indiv)) %>%
  mutate(OldIndiv = ifelse((is.na(Cohort_2) | Cohort_2 != 2021), D2_indiv, D1_indiv)) %>%
  mutate(isPar = ifelse(PO_logl >= 4, TRUE, FALSE)) %>% #Use a relaxed threshold compared to CKMR because influence of false-positives vs. false-negatives much more even
  group_by(YoungIndiv) %>%
  summarize(isPar = max(as.numeric(isPar)), maxPO = max(PO_logl), highestMatch = first(OldIndiv)) %>%
  arrange(desc(maxPO)) %>%
  group_by(isPar) %>%
  summarize(n = n())

ratios.pred %>% 
  filter((Cohort_2 == 2021 &
            Origin_2 == "stocked" &
            (D1_indiv %in% 
               c(combdata %>% filter(Year == "2020", Site == "Gordon Landing"))$NovoID | 
               grepl("EW-", D1_indiv))) |
           (Cohort_1 == 2021  &
              Origin_1 == "stocked"  &
              (D2_indiv %in% 
                 c(combdata %>% filter(Year == "2020", Site == "Gordon Landing"))$NovoID | 
                 grepl("EW-", D2_indiv)))) %>%
  arrange(desc(PO_logl)) %>%
  mutate(YoungIndiv = ifelse((is.na(Cohort_2) | Cohort_2 != 2021), D1_indiv, D2_indiv)) %>%
  mutate(OldIndiv = ifelse((is.na(Cohort_2) | Cohort_2 != 2021), D2_indiv, D1_indiv)) %>%
  mutate(isPar = ifelse(PO_logl >= PO.thresh, TRUE, FALSE)) %>%
  group_by(YoungIndiv) %>%
  summarize(isPar = max(as.numeric(isPar)), maxPO = max(PO_logl), highestMatch = first(OldIndiv)) %>%
  arrange(desc(maxPO)) %>%
  mutate(isPar = PO_logl) %>%
  group_by(isPar) %>%
  summarize(n = n())
  
  ### Run CKMR model using previously identified relationships ----

  
  JAGS.input = true.HS.sum %>% 
    filter(Cohort_1 >= 2012) %>%
    mutate(Cohort_1 = "All") %>%
    group_by(Cohort_1, nyears) %>% 
    summarize(ncombs = sum(ncombs), 
              ntrue = sum(ntrue), 
              Adj_ntrue = sum(Adj_ntrue)) %>%
    rename(EstYear = Cohort_1, 
           PairMortYears = nyears) %>%
    mutate(RObase = 4)
  
  
  JAGS.input = true.PO.sum %>% 
    mutate(Cohort_2 = "All") %>%
    group_by(Cohort_2, nyears, Sex_1) %>% 
    summarize(ncombs = sum(ncombs), 
              ntrue = sum(ntrue), 
              Adj_ntrue = sum(Adj_ntrue)) %>%
    rename(EstYear = Cohort_2, 
           PairMortYears = nyears) %>%
    mutate(RObase = 2)
  
  
  JAGS.input = rbind((true.HS.sum %>% 
                        filter(Cohort_1 >= 2014) %>%
                        #filter(nyears <=3) %>%
                        mutate(Cohort_1 = "All") %>%
                        group_by(Cohort_1, nyears) %>% 
                        summarize(ncombs = sum(ncombs), 
                                  ntrue = sum(ntrue), 
                                  Adj_ntrue = sum(Adj_ntrue)) %>%
                        rename(EstYear = Cohort_1, 
                               PairMortYears = nyears) %>%
                        mutate(RObase = 4, Sex_1 = NA, SpawnSamp = 1, KinError = (1 - HS.FNrate))),
                     (true.PO.sum %>% 
                        filter(Cohort_2 >= 2012) %>%
                        mutate(Cohort_2 = "All") %>%
                        group_by(Cohort_2, nyears, Sex_1,SpawnSpecific) %>% 
                        summarize(ncombs = sum(ncombs), 
                                  ntrue = sum(ntrue), 
                                  Adj_ntrue = sum(Adj_ntrue)) %>%
                        rename(EstYear = Cohort_2, 
                               PairMortYears = nyears) %>%
                        mutate(RObase = ifelse(is.na(Sex_1),2,1), SpawnSamp = 0, 
                               KinError = (1 - PO.FNrate)))) %>%
    arrange(EstYear)
  
  #Remove any years with true zeroes
  JAGS.input.Ns = JAGS.input %>% 
    group_by(EstYear) %>% 
    summarize(ntrue = sum(ntrue))
  
  JAGS.input = JAGS.input %>% 
    filter(EstYear %in% (JAGS.input.Ns %>% filter(ntrue >= 3))$EstYear)
  
  as.data.frame(JAGS.input)
  
  jags_data = list(
    nestyears = length(unique(JAGS.input$EstYear)),
    nobs = length(JAGS.input$ntrue),
    estyear = match(JAGS.input$EstYear,unique(JAGS.input$EstYear)),
    PairMortYears = JAGS.input$PairMortYears,
    KnownSex = as.numeric(!is.na(JAGS.input$Sex_1)),
    IsFemale = as.numeric(!is.na(JAGS.input$Sex_1) & JAGS.input$Sex_1 == "F"),
    SpawnSamp = JAGS.input$SpawnSamp,
    RObase = JAGS.input$RObase,
    PotPairs = JAGS.input$ncombs,
    KinError = JAGS.input$KinError,
    TruePairs = round(JAGS.input$Adj_ntrue,0))
  
  estyears.out = unique(JAGS.input$EstYear)
  
  
  CKMR_est = function(){
    
    #Blanket survival rate
    surv ~ dbeta(1, 1)
    
    #Male proportion (inverse = female proportion)
    propM ~ dbeta(1, 1)
    
    #Proportion successfully spawning
    propSpawn ~ dbeta(1, 1)
    
    #One adult estimate (all origins combined) for each year/year combo included in input
    for(i in 1:nestyears)
    {
      Nadult_all[i] ~ dunif(100,1000000)
    }
    
    #Match to lines in observation
    for(j in 1:nobs)
    {
      TruePairs[j] ~ dbinom((RObase[j] * surv^PairMortYears[j] * KinError[j]) / 
                              (((Nadult_all[estyear[j]] * (1 - KnownSex[j])) +
                                  (Nadult_all[estyear[j]] * KnownSex[j] * abs(propM - IsFemale[j]))) * propSpawn^SpawnSamp[j]),
                            PotPairs[j])
    }
    
  }
  
  
  jags_inits = function(nc) {
    inits = list()
    for(c in 1:nc){
      inits[[c]] = list(
        surv = runif(1,0.01, 0.99),
        propM = runif(1,0.01, 0.99),
        propSpawn = runif(1,0.01, 0.99),
        Nadult_all = runif(length(unique(JAGS.input$EstYear)),100, 1000000)
        
      )
    }
    
    return(inits)
  }
  
  
  
  # Write model
  jags_file = "Model location\\CKMR_est.txt"
  write_model(CKMR_est, jags_file)
  
  jags_params = c("Nadult_all", "surv", "propM", "propSpawn")
  n_params = length(jags_params) #used to autofill dataframe later
  
  
  #------------- STEP 5: SET MCMC DIMENSIONS ---------------#
  jags_dims = c(
    ni = 5000000,  # number of post-burn-in samples per chain
    nb = 1500000,  # number of burn-in samples
    nt = 5,     # thinning rate
    nc = 3      # number of chains
  )
  
  
  MCMC.settings <- paste0("thin", jags_dims[names(jags_dims) == "nt"], "_draw", jags_dims[names(jags_dims) == "ni"], "_burn", jags_dims[names(jags_dims) == "nb"])
  
  #---------------- STEP 6: RUN JAGS ---------------#
  post = jagsUI::jags(data = jags_data, inits = jags_inits(jags_dims["nc"]), 
                      parameters.to.save = jags_params, jags_file, 
                      n.burnin = jags_dims["nb"], n.chains = jags_dims["nc"], 
                      n.iter = sum(jags_dims[c("ni", "nb")]), parallel = T,  verbose = T)
  
  
  post
  
  post.sum = as.data.frame(post$summary)
  row.names(post.sum) = c(paste0("Nadult-", estyears.out),row.names(post.sum)[which(!grepl("Nadult",row.names(post.sum)))])
  post.sum$year = NA
  post.sum[which(grepl("Nadult", row.names(post.sum))),]$year = estyears.out
  
  ggplot(post.sum %>% filter(!is.na(year)))+
    geom_pointrange(aes(x = year, y = `50%`, ymin = `2.5%`, ymax = `97.5%`))+
    scale_y_continuous("Estimated adult abundance", 
                       breaks = seq(0,100000, by = 10000),
                       labels = format(seq(0,100000, by = 10000), big.mark = ","))+
    #scale_x_continuous("",breaks = seq(2012,2023))+
    coord_cartesian(ylim = c(0,100000))
  
 
