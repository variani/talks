### inc
library(plyr)

load_all("~/git/ugcd/solarius")
load_all("~/git/ugcd/gait")
load_all("~/git/hemostat/lme4qtl")

### model
K <- 0.05
dat <- mutate(subset(phen2, !is.na(ABOf3)), offset = log(K/(1 - K)))

m0 <- relmatGlmer(Throm ~ AGEf2 + SEXf + (1|ID), dat, 
  relmat = list(ID = dkin2), family = binomial)

# Sex-specificity only in the residual variance
m1 <- relmatGlmer(Throm ~ AGEsc + SEXf + (1|ID) + (0 + SEXf|RID), dat, 
  relmat = list(ID = dkin2), family = binomial,
  weights = rep(1e10, nrow(dat)), vcControl = list(rho0 = list(rid = 3)))

m1 <- relmatLmer(
  BMI ~ AGEsc + AGEsc2 + SEXf + (1|ID) + (0 + SEXf|RID) + (1|HHID), 
  phen, relmat = list(ID = dkin), 
  weights = rep(1e10, nrow(phen)), vcControl = list(rho0 = list(rid = 3)))





### BMI
dat <- phen2

m0 <- relmatLmer(
  BMI ~ AGEsc + AGEsc2 + SEXf + (1|HHID) + (1|ID), 
  dat, relmat = list(ID = dkin2))

# Sex-specificity only in the residual variance
m1 <- relmatLmer(
  BMI ~ AGEsc + AGEsc2 + SEXf + (1|ID) + (0 + AGEf2|RID) + (1|HHID), 
  dat, relmat = list(ID = dkin2), 
  weights = rep(1e10, nrow(dat)), vcControl = list(rho0 = list(rid = 3)))

# Sex-specificity in both polygenic and residual variances  
m2 <- relmatLmer(
  BMI ~ AGEsc + AGEsc2 + SEXf + (0 + AGEf2|ID) + (0 + AGEf2|RID) + (1|HHID), 
  dat, relmat = list(ID = dkin2), 
  weights = rep(1e10, nrow(dat)), vcControl = list(rho0 = list(rid = 5)))      


### APTT
dat <- phen2

m0 <- relmatLmer(
  APTT ~ AGEsc + AGEsc2 +(1|ID), 
  dat, relmat = list(ID = dkin2))

# Sex-specificity only in the residual variance
m1 <- relmatLmer(
  APTT ~ AGEsc + AGEsc2 + (1|ID) + (0 + AGEf2|RID), 
  dat, relmat = list(ID = dkin2), 
  weights = rep(1e10, nrow(dat)), vcControl = list(rho0 = list(rid = 3)))

# Sex-specificity in both polygenic and residual variances  
m2 <- relmatLmer(
  APTT ~ AGEsc + AGEsc2 + (0 + AGEf2|ID) + (0 + AGEf2|RID), 
  dat, relmat = list(ID = dkin2), 
  weights = rep(1e10, nrow(dat)), vcControl = list(rho0 = list(rid = 5)))
  
  
  
          
