library(FSA)     # for vbFuns(), vbStarts(), confint.bootCase()
library(nlstools)
library(plotrix)

female <- filter(original_clean,Sex=="F") %>% 
  drop_na(Length_mm)
f.agesum <- female %>%
  summarize(minage=min(FinalAge),maxage=max(FinalAge))
f.agesum


male <- filter(original_clean,Sex=="M") %>% 
  drop_na(Length_mm)
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
