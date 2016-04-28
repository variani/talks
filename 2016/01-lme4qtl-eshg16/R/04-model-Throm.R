### inc
library(devtools)
load_all("~/git/hemostat/lme4qtl")

### GAIT2 models 
m1 <- relmatGlmer(Throm ~ (1|ID), phen2, relmat = list(ID = dkin2), family = binomial)

m2 <- relmatGlmer(Throm ~ AGEsc + SEXf + (1|ID), phen2, relmat = list(ID = dkin2), family = binomial)

m3 <- relmatGlmer(Throm ~ AGEsc + SEXf + (1|ID), phen2, relmat = list(ID = dkin2), family = binomial(probit))

K <- 0.01
m4 <- relmatGlmer(Throm ~ -1 + AGEsc + SEXfnum + (1|ID), mutate(phen2, offset = log(K/(1 - K))), offset = offset, relmat = list(ID = dkin2), family = binomial)

m40 <- relmatGlmer(Throm ~ -1 + AGEsc + SEXf + (1|ID), phen2, relmat = list(ID = dkin2), family = binomial)

m41 <- relmatGlmer(Throm ~ -1 + AGEsc + AGEsc2 + SEXfnum + (1|ID), mutate(phen2, offset = log(K/(1 - K))), offset = offset, relmat = list(ID = dkin2), family = binomial)

m5 <- relmatGlmer(Throm ~ -1 + AGEsc + SEXfnum + (1|ID), mutate(phen2, offset = -qnorm(1 - K)), offset = offset, relmat = list(ID = dkin2), family = binomial(probit))
 
m50 <- relmatGlmer(Throm ~ -1 + AGEsc + SEXf + (1|ID), phen2, relmat = list(ID = dkin2), family = binomial(probit)) 

m6 <- relmatGlmer(Throm ~ -1 + AGEsc + SEXfnum + ABOf2num + (1|ID), mutate(phen2, offset = log(K/(1 - K))), offset = offset, relmat = list(ID = dkin2), family = binomial)
#m6 <- relmatGlmer(Throm ~ -1 + AGEsc + AGEsc2 + SEXfnum + ABOf2num + (1|ID), mutate(phen2, offset = log(K/(1 - K))), offset = offset, relmat = list(ID = dkin2), family = binomial)

m60 <- relmatGlmer(Throm ~ AGEsc + SEXf + ABOf2 + (1|ID), phen2, relmat = list(ID = dkin2), family = binomial)

### print
printCoefmat(summary(m5)$coefficients, zap.ind = 3, digits = 4, signif.stars = TRUE)

m <- relmatGlmer(Throm ~ AGEf3 + SEXf + ABOf2 + (1|ID), phen2, relmat = list(ID = dkin2), family = binomial(probit))
m <- relmatGlmer(Throm ~ AGEf + SEXf + ABOf2 + (1|ID), phen2, relmat = list(ID = dkin2), family = binomial)

K <- 0.025
m <- relmatGlmer(Throm ~ -1 + AGEsc + SEXfnum + ABOf2num + (1|ID), mutate(phen2, offset = log(K/(1 - K))), offset = offset, relmat = list(ID = dkin2), family = binomial)

m <- relmatGlmer(Throm ~ -1 + AGEf3num1 + AGEf3num2 + SEXfnum + ABOf2num + (1|ID), mutate(phen2, offset = -qnorm(1 - K)), offset = offset, relmat = list(ID = dkin2), family = binomial(probit))

K <- 0.01
dat <- mutate(phen2, offset = -qnorm(1 - K))

m <- relmatGlmer(Throm ~ -1 + AGEsc + SEXfnum + ABOf2num + (1|ID), dat, offset = offset, relmat = list(ID = dkin2), family = binomial(probit)) 

