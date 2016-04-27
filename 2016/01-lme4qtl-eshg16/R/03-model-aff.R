### inc
library(devtools)
load_all("~/git/hemostat/lme4qtl")

### GAIT1 models 
m1 <- relmatGlmer(aff ~ (1|ID), phen, relmat = list(ID = dkin), family = binomial)

m2 <- relmatGlmer(aff ~ AGEsc + SEXf + (1|ID), phen, relmat = list(ID = dkin), family = binomial)

m3 <- relmatGlmer(aff ~ (1|ID), phen, relmat = list(ID = dkin), family = binomial(probit))

K <- 0.01
m4 <- relmatGlmer(aff ~ -1 + AGEsc + SEXf + (1|ID), mutate(phen, offset = log(K/(1 - K))), offset = offset, relmat = list(ID = dkin), family = binomial)

m40 <- relmatGlmer(aff ~ -1 + AGEsc + SEXf + (1|ID), phen, relmat = list(ID = dkin), family = binomial)
 
#m5 <- relmatGlmer(aff ~ -1 + (1|ID), phen, offset = rep(-qnorm(1 - 0.1), nrow(phen)), relmat = list(ID = dkin), family = binomial(probit))





