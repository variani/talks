### inc
library(devtools)
load_all("~/git/hemostat/lme4qtl")

K <- 0.02
#dat <- mutate(phen, offset = log(K/(1 - K)))

#m <- relmatGlmer(aff ~ -1 + AGEsc + SEXfnum + ab0f2num + (1|ID), dat, offset = offset, relmat = list(ID = dkin), family = binomial)
#m0 <- relmatGlmer(aff ~ AGEsc + SEXfnum + ab0f2num + (1|ID), dat, relmat = list(ID = dkin), family = binomial) 
  
dat <- mutate(phen, offset = -qnorm(1 - K))

#m <- relmatGlmer(aff ~ -1 + AGEsc + AGEsc2 + SEXfnum + ab0f2num + (1|ID), dat, offset = offset, relmat = list(ID = dkin), family = binomial(probit))
#m0 <- relmatGlmer(aff ~ AGEsc + AGEsc2 + SEXfnum + ab0f2num + (1|ID), dat, relmat = list(ID = dkin), family = binomial(probit)) 

m <- relmatGlmer(aff ~ -1 + AGEf3num1 + AGEf3num2 + SEXfnum + ab0f2num + (1|ID), dat, offset = offset, relmat = list(ID = dkin), family = binomial(probit))
m0 <- relmatGlmer(aff ~ AGEf3 + SEXfnum + ab0f2num + (1|ID), dat, relmat = list(ID = dkin), family = binomial)

stop()

K <- 0.02
dat <- mutate(phen2, offset = -qnorm(1 - K))
m <- relmatGlmer(Throm ~ -1 + AGEsc + SEXfnum + ABOf3num + (1|ID), dat, offset = offset, relmat = list(ID = dkin2), family = binomial(probit))

K <- 0.05
dat <- mutate(phen2, offset = log(K/(1 - K)))
m <- relmatGlmer(Throm ~ -1 + AGEsc + SEXfnum + ABOf3num + (1|ID), dat, offset = offset, relmat = list(ID = dkin2), family = binomial)


stop()

### GAIT1 models 
m1 <- relmatGlmer(aff ~ (1|ID), phen, relmat = list(ID = dkin), family = binomial)

m2 <- relmatGlmer(aff ~ AGEsc + SEXf + (1|ID), phen, relmat = list(ID = dkin), family = binomial)

m3 <- relmatGlmer(aff ~ (1|ID), phen, relmat = list(ID = dkin), family = binomial(probit))

K <- 0.01
m4 <- relmatGlmer(aff ~ -1 + AGEsc + SEXf + (1|ID), mutate(phen, offset = log(K/(1 - K))), offset = offset, relmat = list(ID = dkin), family = binomial)

m40 <- relmatGlmer(aff ~ -1 + AGEsc + SEXf + (1|ID), phen, relmat = list(ID = dkin), family = binomial)
 
#m5 <- relmatGlmer(aff ~ -1 + (1|ID), phen, offset = rep(-qnorm(1 - 0.1), nrow(phen)), relmat = list(ID = dkin), family = binomial(probit))





