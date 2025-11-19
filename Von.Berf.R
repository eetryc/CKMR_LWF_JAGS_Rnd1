library(FSA)     # for vbFuns(), vbStarts(), confint.bootCase()
library(nlstools)
library(plotrix)

female <- filter(original_clean,Sex=="F")
f.agesum <- female %>%
  summarize(minage=min(FinalAge),maxage=max(FinalAge))
f.agesum


male <- filter(original_clean,Sex=="M") 
m.agesum <- male %>%
  summarize(minage=min(FinalAge),maxage=max(FinalAge))
m.agesum



vbmod <- Length_mm ~ Linf * (1 - exp(-K * (FinalAge - t0)))


starts <- vbStarts(formula = Length_mm~FinalAge,data=female)


# fit model function using least squares
mymod <- nls(vbmod, data = female, start = starts)

summary(mymod)

pred <- predict(mymod)

vbO <- vbFuns("typical")

vb_fit <- nls(Length_mm~vbO(FinalAge,Linf,K, t0), data=female, start=starts)

boot_fit <- nlsBoot(vb_fit)

boot_preds <- data.frame(
  predict(boot_fit, vbO, t = sort(unique(female$FinalAge))))


names(boot_preds) <- c("FinalAge", "fit", "lwr", "upr")

female_preds <- merge(female, boot_preds, by = "FinalAge")




# plot
ggplot(female_preds, aes(x = FinalAge, y = Length_mm)) +
  geom_jitter(width = 0.1, alpha = 0.15, size = 2) +
  geom_line(aes(y = fit)) +
  geom_ribbon(
    aes(x = FinalAge, ymin = lwr, ymax = upr, color = NULL), alpha = 0.3) +
  xlab("Age (years)") +
  ylab("Total length (mm)") +
  theme_bw()




female.vonbert <- function(x){
  617.2307131 * (1 - exp(-0.1582472 * (x + 4.2497746)))
}












### Male 

startsM <- vbStarts(formula = Length_mm~FinalAge,data=female)


# fit model function using least squares
mymodM <- nls(vbmod, data = male, start = startsM)

summary(mymodM)

pred <- predict(mymodM)

vbO <- vbFuns("typical")

vb_fitm <- nls(Length_mm~vbO(FinalAge,Linf,K, t0), data=male, start=startsM)

boot_fitm <- nlsBoot(vb_fitm)

boot_predsm <- data.frame(
  predict(boot_fitm, vbO, t = sort(unique(male$FinalAge))))


names(boot_predsm) <- c("FinalAge", "fit", "lwr", "upr")

male_preds <- merge(male, boot_preds, by = "FinalAge")




# plot
ggplot(male_preds, aes(x = FinalAge, y = Length_mm)) +
  geom_jitter(width = 0.1, alpha = 0.15, size = 2) +
  geom_line(aes(y = fit)) +
  geom_ribbon(
    aes(x = FinalAge, ymin = lwr, ymax = upr, color = NULL), alpha = 0.3) +
  xlab("Age (years)") +
  ylab("Total length (mm)") +
  theme_bw()
























##### Female curve with otolith/ scale aging #######

# otolith
female_otolith <- original_clean %>% 
  filter(Sex=="F",
         FinalAgeStructure == "Otolith") %>% 
  drop_na(Length_mm)
f.agesum <- female_otolith %>%
  summarize(minage=min(FinalAge),maxage=max(FinalAge))
f.agesum




vbmod <- Length_mm ~ Linf * (1 - exp(-K * (FinalAge - t0)))


starts <- vbStarts(formula = Length_mm~FinalAge,data=female_otolith)


# fit model function using least squares
mymod <- nls(vbmod, data = female_otolith, start = starts)

summary(mymod)

pred <- predict(mymod)

vbO <- vbFuns("typical")

vb_fit <- nls(Length_mm~vbO(FinalAge,Linf,K, t0), data=female_otolith, start=starts)

boot_fit <- nlsBoot(vb_fit)

boot_preds <- data.frame(
  predict(boot_fit, vbO, t = sort(unique(female_otolith$FinalAge))))


names(boot_preds) <- c("FinalAge", "fit", "lwr", "upr")

female_preds <- merge(female_otolith, boot_preds, by = "FinalAge")




# plot
ggplot(female_preds, aes(x = FinalAge, y = Length_mm)) +
  geom_jitter(width = 0.1, alpha = 0.15, size = 2) +
  geom_line(aes(y = fit)) +
  geom_ribbon(
    aes(x = FinalAge, ymin = lwr, ymax = upr, color = NULL), alpha = 0.3) +
  xlab("Age (years)") +
  ylab("Total length (mm)") +
  theme_bw()

## otolith VB curve is worse??






## Spine
female_spine <- original_clean %>% 
  filter(Sex=="F",
         FinalAgeStructure == "Spine") %>% 
  drop_na(Length_mm)
f.agesum <- female_spine %>%
  summarize(minage=min(FinalAge),maxage=max(FinalAge))
f.agesum




vbmod <- Length_mm ~ Linf * (1 - exp(-K * (FinalAge - t0)))


starts <- vbStarts(formula = Length_mm~FinalAge,data=female_spine)


# fit model function using least squares
mymod <- nls(vbmod, data = female_spine, start = starts)

summary(mymod)

pred <- predict(mymod)

vbO <- vbFuns("typical")

vb_fit <- nls(Length_mm~vbO(FinalAge,Linf,K, t0), data=female_spine, start=starts)

boot_fit <- nlsBoot(vb_fit)

boot_preds <- data.frame(
  predict(boot_fit, vbO, t = sort(unique(female_spine$FinalAge))))


names(boot_preds) <- c("FinalAge", "fit", "lwr", "upr")

female_preds <- merge(female_spine, boot_preds, by = "FinalAge")




# plot
ggplot(female_preds, aes(x = FinalAge, y = Length_mm)) +
  geom_jitter(width = 0.1, alpha = 0.15, size = 2) +
  geom_line(aes(y = fit)) +
  geom_ribbon(
    aes(x = FinalAge, ymin = lwr, ymax = upr, color = NULL), alpha = 0.3) +
  xlab("Age (years)") +
  ylab("Total length (mm)") +
  theme_bw()

## Spine aging is similarly ugly....








#### Non spine/otolith

female_scale <- original_clean %>% 
  filter(Sex=="F",
         FinalAgeStructure != "Spine",
         FinalAgeStructure !="Otolith") %>% 
  drop_na(Length_mm)
f.agesum <- female_scale %>%
  summarize(minage=min(FinalAge),maxage=max(FinalAge))
f.agesum




vbmod <- Length_mm ~ Linf * (1 - exp(-K * (FinalAge - t0)))


starts <- vbStarts(formula = Length_mm~FinalAge,data=female_scale)


# fit model function using least squares
mymod <- nls(vbmod, data = female_scale, start = starts)

summary(mymod)

pred <- predict(mymod)

vbO <- vbFuns("typical")

vb_fit <- nls(Length_mm~vbO(FinalAge,Linf,K, t0), data=female_scale, start=starts)

boot_fit <- nlsBoot(vb_fit)

boot_preds <- data.frame(
  predict(boot_fit, vbO, t = sort(unique(female_scale$FinalAge))))


names(boot_preds) <- c("FinalAge", "fit", "lwr", "upr")

female_preds <- merge(female_scale, boot_preds, by = "FinalAge")




# plot
ggplot(female_preds, aes(x = FinalAge, y = Length_mm)) +
  geom_jitter(width = 0.1, alpha = 0.15, size = 2) +
  geom_line(aes(y = fit)) +
  geom_ribbon(
    aes(x = FinalAge, ymin = lwr, ymax = upr, color = NULL), alpha = 0.3) +
  xlab("Age (years)") +
  ylab("Total length (mm)") +
  theme_bw()
















#### Decimal Ages
original_decimal <- original_clean %>%
  mutate(Cohort = SampleYear - FinalAge,
    SampleDate = mdy(SampleDate),
    Month = month(SampleDate)) %>% 
  mutate(decimal.length = case_when(
    Month %in% 9:12 ~ Length_mm,
    Month %in% 1:5 ~ Length_mm - 0.05*(Length_mm),
    Month %in% 6:8 ~ Length_mm - 0.025*(Length_mm)))



female_dec <- original_decimal %>% 
  filter(Sex=="F") %>% 
  drop_na(decimal.length)
f.agesum <- female_dec %>%
  summarize(minage=min(FinalAge),maxage=max(FinalAge))
f.agesum




vbmod <- decimal.length ~ Linf * (1 - exp(-K * (FinalAge - t0)))


starts <- vbStarts(formula = decimal.length~FinalAge,data=female_dec)


# fit model function using least squares
mymod <- nls(vbmod, data = female_dec, start = starts)

summary(mymod)

pred <- predict(mymod)

vbO <- vbFuns("typical")

vb_fit <- nls(decimal.length~vbO(FinalAge,Linf,K, t0), data=female_dec, start=starts)

boot_fit <- nlsBoot(vb_fit)

boot_preds <- data.frame(
  predict(boot_fit, vbO, t = sort(unique(female_dec$FinalAge))))


names(boot_preds) <- c("FinalAge", "fit", "lwr", "upr")

female_preds <- merge(female_dec, boot_preds, by = "FinalAge")




# plot
ggplot(female_preds, aes(x = FinalAge, y = decimal.length)) +
  geom_jitter(width = 0.1, alpha = 0.15, size = 2) +
  geom_line(aes(y = fit)) +
  geom_ribbon(
    aes(x = FinalAge, ymin = lwr, ymax = upr, color = NULL), alpha = 0.3) +
  xlab("Age (years)") +
  ylab("Total length (mm)") +
  theme_bw()




## only plot spawning season observations
original_spawners <- original_clean %>%
  mutate(Cohort = SampleYear - FinalAge,
         SampleDate = mdy(SampleDate),
         Month = month(SampleDate)) %>% 
  filter(Month %in% 9:12)



female_spawner <- original_spawners %>% 
  filter(Sex=="F") %>% 
  drop_na(Length_mm)
f.agesum <- female_spawner %>%
  summarize(minage=min(FinalAge),maxage=max(FinalAge))
f.agesum




vbmod <- Length_mm ~ Linf * (1 - exp(-K * (FinalAge - t0)))


starts <- vbStarts(formula = Length_mm~FinalAge,data=female_spawner)


# fit model function using least squares
mymod <- nls(vbmod, data = female_spawner, start = starts)

summary(mymod)

pred <- predict(mymod)

vbO <- vbFuns("typical")

vb_fit <- nls(Length_mm~vbO(FinalAge,Linf,K, t0), data=female_spawner, start=starts)

boot_fit <- nlsBoot(vb_fit)

boot_preds <- data.frame(
  predict(boot_fit, vbO, t = sort(unique(female_spawner$FinalAge))))


names(boot_preds) <- c("FinalAge", "fit", "lwr", "upr")

female_preds <- merge(female_spawner, boot_preds, by = "FinalAge")




# plot
ggplot(female_preds, aes(x = FinalAge, y = Length_mm)) +
  geom_jitter(width = 0.1, alpha = 0.15, size = 2) +
  geom_line(aes(y = fit)) +
  geom_ribbon(
    aes(x = FinalAge, ymin = lwr, ymax = upr, color = NULL), alpha = 0.3) +
  xlab("Age (years)") +
  ylab("Total length (mm)") +
  theme_bw()



# similarly messy when compared to the other non-/filtered curves....