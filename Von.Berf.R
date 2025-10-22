library(FSA)     # for vbFuns(), vbStarts(), confint.bootCase()
library(car)     # for Boot()
library(dplyr)   # for filter(), mutate()
library(ggplot2)

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


svtyp <- vbStarts(Length_mm~FinalAge,data=female)
unlist(svtyp)


?growthModelSsim
