##  chi-square tests and t-tests

load(file = '/export/home/xurren/WTCProject/Data/clinical533_reordered.RData')  ## load clinical data

race1 = rep(NA, dim(clinical)[1])
race1[clinical$race==1] = 1
race1[clinical$race>1] = 0
race1[is.na(clinical$race)] = 1  ## assume NA is white
race1 = factor(race1)
levels(race1) = c("non-white", "white")

smoke_status = factor(clinical$smoke_status)
levels(smoke_status) = c("non-smoker", "smoker")

age = clinical$age_RNA

PTSD_2levels = rep(NA, dim(clinical)[1])
PTSD_2levels[clinical$PTSD_SCID_3grp==1] = 1
PTSD_2levels[clinical$PTSD_SCID_3grp==3] = 0
PTSD_2levels = factor(PTSD_2levels)
levels(PTSD_2levels) = c("control", "current")

PTSD_3levels = clinical$PTSD_SCID_3grp
PTSD_3levels = factor(PTSD_3levels)
levels(PTSD_3levels) = c("current", "past", "control")

##  chi-square test on PTSD vs race 
PTSD2vsRace = table(PTSD_2levels, race1, useNA = "no") 
chisq.test(PTSD2vsRace)

PTSD3vsRace = table(PTSD_3levels, race1, useNA = "no")
chisq.test(PTSD3vsRace)

##  chi-square test on PTSD vs smoke_status
PTSD2vsSmoke = table(PTSD_2levels, smoke_status, useNA = "no") 
chisq.test(PTSD2vsSmoke)

PTSD3vsSmoke = table(PTSD_3levels, smoke_status, useNA = "no")
chisq.test(PTSD3vsSmoke)

##  t-test on PTSD(2 levels) vs age
x = age[which(PTSD_2levels == "control")]
y = age[which(PTSD_2levels == "current")]
t.test(x, y)

##  ANOVA on PTSD(3 levels) vs age
fit = lm(age~PTSD_3levels)
anova(fit)



